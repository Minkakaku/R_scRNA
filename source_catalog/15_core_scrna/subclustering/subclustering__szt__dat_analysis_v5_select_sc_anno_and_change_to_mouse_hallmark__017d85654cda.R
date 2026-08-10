library(tidyverse)
library(Seurat)


all_data <- readRDS("sp_fibro_rctd_anno.rds")
all_data <- all_data %>% JoinLayers()
DefaultAssay(all_data) <- "RNA"

all_data_sc2 <- readRDS('all_data_sc_subcluster2.rds')

out_dir_v5 <- "20251213 dat_anlysis"
dir.create(out_dir_v5)



path_ls <-
  readxl::excel_sheets("Supplementary Table 2_41586_2024_7806_MOESM4_ESM.xlsx") %>%
  setNames(., str_trim(str_extract(., "\\s[A-Za-z]*\\s"))) %>%
  setNames(., str_replace(names(.), "Ischemic", "IZ")) %>%
  map(function(x){
    dat <- readxl::read_excel("Supplementary Table 2_41586_2024_7806_MOESM4_ESM.xlsx", sheet = x)
    # if('avg_log2FC' %in% colnames(dat)){
    #   dat <-
    #     dat %>%
    #     dplyr::filter(avg_log2FC > .25, p_val_adj < 0.05, pct.1 > pct.2)
    # }
    dat$gene
  }) %>%
  map(str_to_title)

path_ls <- path_ls[c('BZ', 'IZ')]


# 0 findmarker of singlecell data and pathway score -----------
all_markers_sc <-
  FindAllMarkers(
    JoinLayers(all_data_sc2),
    group.by = 'annotation',
    return.thresh = 1.1
  )

de_gene_cm <- 
  all_markers_sc %>%
  dplyr::filter(cluster == 'CM') %>%
  dplyr::filter(
    avg_log2FC > 2,
    p_val_adj < 0.05,
    pct.1 - pct.2 > 0.3
  ) %>%
  dplyr::pull('gene')

# signature of papers
path_ls2 <-
  readxl::read_excel("44161_2022_160_MOESM2_ESM.xlsx", sheet = 1) %>%
  group_by(cluster) %>%
  nest() %>%
  deframe() %>%
  map('gene') %>%
  setNames(., fct_recode(names(.), RZ = 'Myh6', BZ1 = 'Ankrd1', BZ2 = 'Xirp2'))

# path_ls[['CM']] <- de_gene_cm
path_ls <- c(path_ls, path_ls2)

all_data_score_zone <- 
  all_data %>%
  AddModuleScore(path_ls)
colnames(all_data_score_zone@meta.data)[
    str_detect(colnames(all_data_score_zone@meta.data), 'Cluster')
  ] <- 
  names(path_ls)


# 1 nn analysis ----------
cell_coords_ls <- readRDS("cells_position_ls.rds")
cell_coords_ls2 <- cell_coords_ls %>% map(~  .[.$cell_id %in% colnames(all_data), ])

## find nn point -----------
k_nn <- 25
res_nn_ls <- 
  cell_coords_ls2 %>%
  map(function(dat){ # browser()
    cn_nn <- 
      RANN::nn2(
        dat[c('x', 'y')], 
        k = k_nn
      )
    
    list(
      nn = cn_nn,
      cell_id = dat$cell_id
    ) 
  })
  



## combined nn expression -----------
dat_zone_sum_ls <- 
  res_nn_ls %>%
  map(function(res_nn){ # browser()
    # dat
    print('subset data')
    dat_used <- FetchData(
      all_data_score_zone,
      vars = names(path_ls),
      cells = res_nn$cell_id
    )
    dat_used <- dat_used[res_nn$cell_id, ]
    
    stopifnot(all(rownames(dat_used) == res_nn$cell_id))

    
    # idx
    nn_idx <- res_nn$nn$nn.idx %>% t %>% c()
    group_idx <- rep(res_nn$cell_id, each = k_nn)
    
    # sum
    print('sum of nn cells')
    dat_merged <- dat_used %>% apply(2, function(x) tapply(x[nn_idx], group_idx, sum) )
    
    print('done.')
    
    # return
    dat_merged
  })


# # Find Cluster
# all_data_path_ls <-
#   dat_zone_sum_ls %>%
#   map(function(dat){
#     sc <- 
#       CreateSeuratObject(
#         CreateAssayObject(
#           # counts = t(dat),
#           data = t(dat)
#         ),
#         meta.data = all_data@meta.data
#       )
#     sc[['RNA']]$scale.data <- sc[['RNA']]$data
#     sc[['RNA']] <- as(sc[['RNA']], "Assay5")
#     sc
#   })
#   
# ll <- all_data_path_ls$A02994A6 %>% RunPCA(features = rownames(.)) %>% RunUMAP(dims = 1:3)


## kmeans cluster -------------

res_kmeans_used_ls <- 
  dat_zone_sum_ls %>%
  map(function(dat){
    set.seed(123)
    kmeans(dat, centers = 4) 
  })
  

p_heat_ls <-
  res_kmeans_used_ls %>%
  map('centers') %>%
  map(function(x){
    pheatmap::pheatmap(
      t(x),
      scale = "row",
      cluster_cols = F,
      # cluster_rows = F,
      border_color = NA,
      silent = T
    )
  })
  

pdf(
  file = file.path(out_dir_v5, "1 kmeans cluster IZ and BZ - heatmap.pdf"),
  height = 5,
  width  = 5
)
p_heat_ls %>% map('gtable') %>% map(cowplot::plot_grid)
dev.off()

df_kmeans_cluster <-
  res_kmeans_used_ls %>%
  map('cluster') %>%
  map(enframe, 'cell_id', 'cluster') %>%
  imap(function(dat, name){
    # cluster annotation by heatmap
    if(name == 'A02994A6'){
      dat <- dat %>%
        mutate(
          cluster = str_c('cluster', cluster),
          cluster = fct_recode(
            cluster, 
            "RZ_c1" = 'cluster1', 
            "RZ_c2" = 'cluster2',
            "BZ/IZ_c3" = 'cluster3',
            "BZ/IZ_c4" = 'cluster4'
          )
        )  
    }else{
      dat <- dat %>%
        mutate(
          cluster = str_c('cluster', cluster),
          cluster = fct_recode(
            cluster, 
            "BZ1_c1" = 'cluster1', 
            "IZ_c2" = 'cluster2',
            "RZ_c3" = 'cluster3',
            "BZ2_c4" = 'cluster4'
          )
        )  
    }
  }) %>%
  bind_rows() %>%
  column_to_rownames("cell_id")  



p_sp_iz_used_ls <-
  cell_coords_ls2 %>%
  imap(function(coord, name){
    coord %>%
      left_join(
        rownames_to_column(df_kmeans_cluster, "cell_id"),
        by = "cell_id"
      ) %>%
      
      slice_sample(n = 50000) %>%
      
      ggplot(aes(x = x, y = -y, color = as.character(cluster))) +
      geom_point(size = 0.3) +
      coord_equal() +
      theme_minimal() +
      # scale_color_manual(values = c(
      #   '1' = "red",
      #   '2' = 'blue'
      # )) +
      guides(color = guide_legend(override.aes = list(size = 2))) +
      labs(
        color = 'cluster',
        title = name
      )
  })
  

pdf(
  file = file.path(out_dir_v5, "1 kmeans cluster IZ and BZ - spatial plot.pdf"),
  height = 20,
  width  = 20
)
p_sp_iz_used_ls %>% walk(print)
dev.off()


## vlnplot of IZ and BZ ------------
 all_data_score_zone <-
  all_data_score_zone %>%
  AddMetaData(df_kmeans_cluster)

p_vln_zone_ls <-
  path_ls %>% 
  names() %>%
  map(function(x){
    all_data_score_zone %>%
      VlnPlot(
        features = x,
        group.by = "project",
        split.by = "cluster",
        pt.size = 0
      )  +
      geom_boxplot(width = 0.15, outliers = F, position = position_dodge(width = 0.9)) +
      ggtitle(x)
  
  })



pdf(
  file = file.path(out_dir_v5, "1 kmeans cluster IZ and BZ - vln plot.pdf"),
  height = 6,
  width  = 12
)
p_vln_zone_ls %>% walk(print)
dev.off()


# export
# saveRDS(all_data_score_zone, 'all_data_score_zone.rds')


# # test not used  --------------
# # combined data and cluster
# res_kmeans_used <- 
#   dat_bz_iz_sum_ls2$A02994A6 %>%
#   kmeans(centers = 2) 
# 
# 
# p_sp_iz_used <-
#   cell_coords_ls2$A02994A6 %>%
#   left_join(
#    enframe(res_kmeans_used$cluster, "cell_id", "cluster"),
#    by = "cell_id"
#   ) %>%
#   
#   slice_sample(n = 50000) %>%
#   
#   ggplot(aes(x = x, y = -y, color = as.character(cluster))) +
#   geom_point() +
#   theme_minimal() +
#   scale_color_manual(values = c(
#     '1' = "red",
#     '2' = 'blue'
#   ))
# 
# 
# ggsave(
#   p_sp_iz_used,
#   filename = file.path(out_dir_v5, "kmeans cluster IZ and BZ using marker gene - used.pdf"),
#   height = 20,
#   width  = 20
# )
#   
# 
# # not used combined data
# res_kmeans_not_used <- 
#   FetchData(
#     all_data,
#     vars = path_ls[1:2] %>% unlist,
#     cells = colnames(all_data)[all_data$project == 'A02994A6']
#   ) %>%
#   kmeans(centers = 2)
# 
# 
# 
# p_sp_iz_not_used <-
#   cell_coords_ls2$A02994A6 %>%
#   left_join(
#     enframe(res_kmeans_not_used$cluster, "cell_id", "cluster"),
#     by = "cell_id"
#   ) %>%
#   
#   slice_sample(n = 50000) %>%
#   
#   ggplot(aes(x = x, y = -y, color = as.character(cluster))) +
#   geom_point() +
#   # scale_color_continuous(trans = "reverse") +
#   theme_minimal() +
#   scale_color_manual(values = c(
#     '2' = "red",
#     '1' = 'blue'
#   ))
# 
# ggsave(
#   p_sp_iz_not_used,
#   filename = file.path(out_dir_v5, "kmeans cluster IZ and BZ using marker gene - not used combined.pdf"),
#   height = 20,
#   width  = 20
# )



# 2 gene expression in each cell -----------

gene_list <- c(
  'Cxcl10',
  'Hspa1b',
  'Hspa1a',
  'Slfn4',
  'Dyrk3',
  'Oasl1',
  'Il1rn'
)


all_dat_to_plot <- 
  all_data_score_zone %>%
  FetchData(
    vars =  c('project', 'group', 'roi', 'cluster', gene_list)
  )

all_dat_to_plot %>% head

all_data_score_zone[['SCT']] <- NULL
all_data_score_zone[['sketch']] <- NULL
all_data_score_zone <- JoinLayers(all_data_score_zone)
all_data_score_zone_slide1 <- all_data_score_zone[, all_data_score_zone$project == "A02994A6"]
all_data_score_zone_slide2 <- all_data_score_zone[, all_data_score_zone$project == "B04866A5"]

all_data_score_zone_slide1 <- CreateSeuratObject(
  CreateAssay5Object(
    counts = all_data_score_zone_slide1@assays$RNA$counts,
    data   = all_data_score_zone_slide1@assays$RNA$data
  ),
  meta.data = all_data_score_zone_slide1@meta.data
)

all_data_score_zone_slide2 <- CreateSeuratObject(
  CreateAssay5Object(
    counts = all_data_score_zone_slide2@assays$RNA$counts,
    data   = all_data_score_zone_slide2@assays$RNA$data
  ),
  meta.data = all_data_score_zone_slide2@meta.data
)



## groupby group ----------
p_dot_grp_slide1 <- 
  DotPlot(
    all_data_score_zone_slide1,
    features = gene_list,
    group.by = 'group'
  )


p_dot_grp_slide2 <- 
  DotPlot(
    all_data_score_zone_slide1,
    features = gene_list,
    group.by = 'group'
  )

# export 
ggsave(
  p_dot_grp_slide1, 
  filename = file.path(out_dir_v5, "2 selected gene dotplot of group in slide1.pdf"),
  height = 4,
  width = 8
)
ggsave(
  p_dot_grp_slide2, 
  filename = file.path(out_dir_v5, "2 selected gene dotplot of group in slide2.pdf"),
  height = 4,
  width = 8
)

## groupby roi -----------

p_dot_roi_slide1 <- 
  DotPlot(
    all_data_score_zone_slide1,
    features = gene_list,
    group.by = 'roi'
  )


p_dot_roi_slide2 <- 
  DotPlot(
    all_data_score_zone_slide2,
    features = gene_list,
    group.by = 'roi'
  )

# export
ggsave(
  p_dot_roi_slide1, 
  filename = file.path(out_dir_v5, "2 selected gene dotplot of roi in slide1.pdf"),
  height = 4,
  width = 8
)
ggsave(
  p_dot_roi_slide2, 
  filename = file.path(out_dir_v5, "2 selected gene dotplot of roi in slide2.pdf"),
  height = 4,
  width = 8
)




# groupby cluster

p_dot_cluster_slide1 <- 
  DotPlot(
    all_data_score_zone_slide1,
    features = gene_list,
    group.by = 'cluster'
  )


p_dot_cluster_slide2 <- 
  DotPlot(
    all_data_score_zone_slide2,
    features = gene_list,
    group.by = 'cluster'
  )

# export
ggsave(
  p_dot_cluster_slide1, 
  filename = file.path(out_dir_v5, "2 selected gene dotplot of cluster in slide1.pdf"),
  height = 4,
  width = 8
)
ggsave(
  p_dot_cluster_slide2, 
  filename = file.path(out_dir_v5, "2 selected gene dotplot of cluster in slide2.pdf"),
  height = 4,
  width = 8
)


# groupby annotation
p_dot_cell_slide1 <- 
  DotPlot(
    all_data_score_zone_slide1 %>% .[, !is.na(.$first_type)],
    features = gene_list,
    group.by = 'first_type'
  )


p_dot_cell_slide2 <- 
  DotPlot(
    all_data_score_zone_slide2 %>% .[, !is.na(.$first_type)],
    features = gene_list,
    group.by = 'first_type'
  )


# export
ggsave(
  p_dot_cell_slide1, 
  filename = file.path(out_dir_v5, "2 selected gene dotplot of annotation in slide1.pdf"),
  height = 8,
  width = 8
)
ggsave(
  p_dot_cell_slide2, 
  filename = file.path(out_dir_v5, "2 selected gene dotplot of annotation in slide2.pdf"),
  height = 8,
  width = 8
)


# 3 de analysis between CM -------

all_marker_ls_slide1 <-
  c( "CM") %>%
  setNames(., .) %>%
  map(function(cell){
    Idents( all_data_score_zone_slide1 ) <- 'first_type'
    
    all_data_score_zone_slide1 %>%
      FindMarkers(
        ident.1 = 'treat',
        ident.2 = 'Veh',
        group.by = 'group',
        subset.ident = cell
      )
  })

all_marker_ls_slide2 <-
  c( "CM") %>%
  setNames(., .) %>%
  map(function(cell){
    Idents( all_data_score_zone_slide2 ) <- 'first_type'
    
    tryCatch(
      all_data_score_zone_slide2 %>%
        FindMarkers(
          ident.1 = 'treat',
          ident.2 = 'Veh',
          group.by = 'group',
          subset.ident = cell
        ),
      error = function(e)return(NULL)
    )
  })

all_markers_ls_cm <- c(
    all_marker_ls_slide1 %>% setNames(., str_c(names(.), 'slide1', sep = "_")),
    all_marker_ls_slide2 %>% setNames(., str_c(names(.), 'slide2', sep = "_"))
  ) %>%
  purrr::discard(is.null)


# export de
all_markers_ls_cm %>%
  openxlsx::write.xlsx(
    file = file.path(out_dir_v5, "3 de analysis table.xlsx"),
    asTable = TRUE,
    withFilter = FALSE
  )


all_geneid_major_ls <-
  all_markers_ls_cm %>%
  map(function(dat){
    gene <-
      dat %>%
      rownames_to_column('gene') %>%
      dplyr::filter(avg_log2FC > 0.25, p_val_adj < 0.05, pct.1 > 0.1) %>%
      # group_by(cluster) %>%
      nest %>%
      deframe() %>%
      map(~{
        clusterProfiler::bitr(
          .x$gene,
          fromType = "SYMBOL",
          toType   = "ENTREZID",
          OrgDb    = "org.Mm.eg.db"
        ) %>%
          .$ENTREZID
      })
    
    gene  
  })

all_geneid_major_ls <- all_geneid_major_ls %>% map(1) %>% purrr::discard(is.null)

all_enrich_major_ls <- 
  all_geneid_major_ls %>%
  map(function(x){
    clusterProfiler::enrichGO(
      gene = x,
      ont = "BP",
      OrgDb = "org.Mm.eg.db",
      pvalueCutoff  = 1, # 0.05,
      pAdjustMethod = "BH",
      qvalueCutoff  = 1 # 0.2
    )
  })

# enrichment 
all_enrich_major_ls <-
  all_enrich_major_ls %>%
  map(function(x){ # browser()
    xx <- x@result$geneID %>% str_split("/")
    
    all_genes <- unlist(xx)
    gene_db <- 
      clusterProfiler::bitr(
        all_genes,
        fromType = "ENTREZID",
        toType   = "SYMBOL",
        OrgDb    = "org.Mm.eg.db"
      ) 
    
    xx2 <-
      xx %>%
      map_chr(function(g){
        str_c(gene_db$SYMBOL[match(g, gene_db$ENTREZID)], collapse = ",")
      })
    x@result$geneSymbol <- xx2
    x
  })


# enrichment plot
p_all_enrich_major_ls <-
  all_enrich_major_ls %>%
  imap(function(enrich, name){ # browser()
    enrichplot::dotplot(
      enrich,
      x = "GeneRatio",
      color = "p.adjust",
      showCategory = 10,
      label_format = 40
    ) + 
      # scale_fill_gradient2(low = "gray", high = "red") +
      ggtitle(name)
  })



pdf(
  file.path(out_dir_v5, "3 de enrichment analysis dotplot.pdf"), 
  height = 6.5, 
  width  = 8
)
p_all_enrich_major_ls %>%  walk(print)
dev.off()


# 4 hallmark ---------
DefaultAssay(all_data) <- "RNA"
gs_h <- msigdbr::msigdbr(
  db_species = "MM",
  species = "Mus musculus", # "Homo sapiens",
  collection = "MH"
) %>%
  with(split(gene_symbol, gs_name))

all_data_score_zone_hallmark <- all_data %>% AddModuleScore(gs_h)

p_hallmark_heat_ls <- 
  all_data_score_zone_hallmark[[]] %>%
  group_by(project) %>% 
  nest() %>%
  deframe() %>%
  imap(function(dat, name){
    dat %>%
      dplyr::select(
        first_type,
        starts_with("Cluster")
      ) %>%
      dplyr::filter(!is.na(first_type)) %>% 
      group_by(first_type) %>% 
      summarise(across(everything(), mean)) %>%
      column_to_rownames("first_type") %>%
      {pheatmap::pheatmap(
        t(.),
        # annotation_col = .['group'],
        scale = "row",
        border_color = NA,
        cluster_rows = TRUE,
        cluster_cols = TRUE,
        labels_row = names(gs_h),
        silent = TRUE
      )}
  })
  

pdf(
  file.path(out_dir_v5, "4 hallmark heatmap - cell.pdf"), 
  height = 10, 
  width  = 7
)
p_hallmark_heat_ls %>% map('gtable') %>%  imap(~ cowplot::plot_grid(.x, labels = .y, scale = 0.9))
dev.off()



# project  hallmark 
p_hallmark_heat_proj <-
  all_data_score_zone_hallmark[[]] %>%
  mutate(project = str_c(project, group, sep = "_")) %>%
  dplyr::select(
    project, 
    starts_with("Cluster")
  ) %>%
  # dplyr::filter(annotation_epi != "CEBPA Epithelial") %>% 
  group_by(project) %>%
  summarise(across(everything(), mean)) %>%
  column_to_rownames("project") %>%
  {pheatmap::pheatmap(
    t(.),
    # annotation_col = .['group'],
    scale = "row",
    border_color = NA,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    labels_row = names(gs_h),
    silent = TRUE
  )}

pdf(
  file.path(out_dir_v5, "4 hallmark heatmap - project.pdf"), 
  height = 10, 
  width  = 7
)
p_hallmark_heat_proj$gtable %>%  cowplot::plot_grid()
dev.off()



# group  hallmark 
p_hallmark_heat_grp_ls <-
  all_data_score_zone_hallmark[[]] %>%
  group_by(project) %>% 
  nest() %>%
  deframe() %>%
  imap(function(dat, name){
    dat %>%
      dplyr::select(
        group, 
        starts_with("Cluster")
      ) %>%
      # dplyr::filter(annotation_epi != "CEBPA Epithelial") %>% 
      group_by(group) %>%
      summarise(across(everything(), mean)) %>%
      column_to_rownames("group") %>%
      {pheatmap::pheatmap(
        t(.),
        # annotation_col = .['group'],
        scale = "row",
        border_color = NA,
        cluster_rows = TRUE,
        cluster_cols = FALSE,
        labels_row = names(gs_h),
        silent = TRUE
      )}
  })
  
pdf(
  file.path(out_dir_v5, "4 hallmark heatmap - group.pdf"), 
  height = 10, 
  width  = 7
)
p_hallmark_heat_grp_ls %>% map('gtable') %>%  imap(~ cowplot::plot_grid(.x, labels = .y, scale = 0.9))
dev.off()



 















