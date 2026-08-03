library(tidyverse)
library(Seurat)


all_data <- readRDS("all_data_projected.rds")
all_data <- all_data %>% JoinLayers()
DefaultAssay(all_data) <- "RNA"


out_dir_v4 <- "20251114 dat_anlysis"
dir.create(out_dir_v4)



path_ls <-
  readxl::excel_sheets("41586_2024_7806_MOESM5_ESM.xlsx") %>%
  setNames(., str_trim(str_extract(., "\\s[A-Z]*\\s"))) %>%
  map(function(x){
    dat <- readxl::read_excel("41586_2024_7806_MOESM5_ESM.xlsx", sheet = x) 
    if('avg_log2FC' %in% colnames(dat)){
      dat <-
        dat %>%
        dplyr::filter(avg_log2FC > .25, p_val_adj < 0.05, pct.1 > pct.2)
    }
    dat$gene
  }) %>%
  map(str_to_title)



# 1 nn analysis ----------
cell_coords_ls <- readRDS("cells_position_ls.rds")
cell_coords_ls2 <- cell_coords_ls %>% map(~  .[.$cell_id %in% colnames(all_data), ])

## find nn point -----------
k_nn <- 10
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
  

# dat_bz <- 
#   FetchData(
#     all_data,
#     vars = path_ls$BZ
#   )
# 
# dat_iz <- 
#   FetchData(
#     all_data,
#     vars = path_ls$IZ
#   )


## combined 10 nn expression -----------
dat_bz_iz_sum_ls <- 
  res_nn_ls %>%
  map(function(res_nn){ # browser()
    # dat
    print('BZ: subset data')
    dat_bz_used <- FetchData(
      all_data,
      vars = path_ls$BZ,
      cells = res_nn$cell_id
    )
    dat_bz_used <- dat_bz_used[res_nn$cell_id, ]
    
    print('IZ: subset data')
    dat_iz_used <- FetchData(
      all_data,
      vars = path_ls$IZ,
      cells = res_nn$cell_id
    )
    dat_iz_used <- dat_iz_used[res_nn$cell_id, ]
    
  
    stopifnot(all(rownames(dat_bz_used) == res_nn$cell_id))
    stopifnot(all(rownames(dat_iz_used) == res_nn$cell_id))
    # stopifnot(all(rownames(dat_bz_used) == rownames(dat_iz_used)))
    
    # cell_id_nn <-  res_nn$cell_id
    # cell_id_order <- match(rownames(dat_bz_used), res_nn$cell_id)
    
    
    # idx
    nn_idx <- res_nn$nn$nn.idx %>% t %>% c()
    group_idx <- rep(res_nn$cell_id, each = k_nn)
    
    # sum
    print('BZ: sum of nn cells')
    dat_bz <- dat_bz_used %>% apply(2, function(x) tapply(x[nn_idx], group_idx, sum) )
    
    print('IZ: sum of nn cells')
    dat_iz <- dat_iz_used %>% apply(2, function(x) tapply(x[nn_idx], group_idx, sum) )
    
    print('done.')
    
    # return
    list(
      dat_bz,
      dat_iz
    )
  })



dat_bz_iz_sum_ls2 <-
  dat_bz_iz_sum_ls %>%
  map(~ do.call(cbind, .))


## kmeans cluster -------------
res_kmeans_used_ls <- 
  dat_bz_iz_sum_ls2 %>%
  map(function(dat){
    set.seed(123)
    kmeans(dat, centers = 2) 
  })
  

p_heat_ls <-
  res_kmeans_used_ls %>%
  map('centers') %>%
  map(function(x){
    pheatmap::pheatmap(
      t(x),
      # scale = "row",
      cluster_cols = F,
      # cluster_rows = F,
      annotation_row = path_ls[1:2] %>% enframe("Group") %>% unnest() %>% column_to_rownames("value"),
      silent = T
    )
  })
  

df_kmeans_cluster <-
  res_kmeans_used_ls %>%
  map('cluster') %>%
  map(enframe, 'cell_id', 'cluster') %>%
  bind_rows() %>%
  mutate(
    cluster = str_c('cluster', cluster),
    cluster = fct_recode(cluster, "BZ" = 'cluster1', "IZ" = 'cluster2')
  ) %>%
  column_to_rownames("cell_id")  



p_sp_iz_used_ls <-
  cell_coords_ls2 %>%
  imap(function(coord, name){
    coord %>%
      left_join(
        enframe(res_kmeans_used_ls[[name]]$cluster, "cell_id", "cluster"),
        by = "cell_id"
      ) %>%
      
      slice_sample(n = 50000) %>%
      
      ggplot(aes(x = x, y = -y, color = as.character(cluster))) +
      geom_point() +
      theme_minimal() +
      scale_color_manual(values = c(
        '1' = "red",
        '2' = 'blue'
      )) +
      ggtitle(name)
  })
  

pdf(
  file = file.path(out_dir_v4, "kmeans cluster IZ and BZ using marker gene - used.pdf"),
  height = 20,
  width  = 20
)
p_sp_iz_used_ls %>% walk(print)
dev.off()


## vlnplot of IZ and BZ ------------
# all_data_score 来源于 `dat_analysis_v2.R`的`3 geneset score`
all_data_score <-
  all_data_score %>%
  AddMetaData(df_kmeans_cluster)


p_vln_bz <-
  all_data_score %>%
  VlnPlot(
    features = "BZ",
    group.by = "project",
    split.by = "cluster",
    pt.size = 0
  )  +
  geom_boxplot(width = 0.15, outliers = F, position = position_dodge(width = 0.9)) +
  facet_wrap(~ ident, scales = "free")

p_vln_iz <-
  all_data_score %>%
  VlnPlot(
    features = "IZ",
    group.by = "project",
    split.by = "cluster",
    pt.size = 0
  )  +
  geom_boxplot(width = 0.15, outliers = F, position = position_dodge(width = 0.9)) +
  facet_wrap(~ ident, scales = "free")

pdf(
  file = file.path(out_dir_v4, "Vlnplot IZ and BZ score in BZ and IZ area.pdf"),
  height = 6,
  width  = 12
)
p_vln_bz %>% print()
p_vln_iz %>% print()
dev.off()


# export
# saveRDS(all_data_score, 'all_data_score.rds')


# test not used  --------------
# combined data and cluster
res_kmeans_used <- 
  dat_bz_iz_sum_ls2$A02994A6 %>%
  kmeans(centers = 2) 


p_sp_iz_used <-
  cell_coords_ls2$A02994A6 %>%
  left_join(
   enframe(res_kmeans_used$cluster, "cell_id", "cluster"),
   by = "cell_id"
  ) %>%
  
  slice_sample(n = 50000) %>%
  
  ggplot(aes(x = x, y = -y, color = as.character(cluster))) +
  geom_point() +
  theme_minimal() +
  scale_color_manual(values = c(
    '1' = "red",
    '2' = 'blue'
  ))


ggsave(
  p_sp_iz_used,
  filename = file.path(out_dir_v4, "kmeans cluster IZ and BZ using marker gene - used.pdf"),
  height = 20,
  width  = 20
)
  

# not used combined data
res_kmeans_not_used <- 
  FetchData(
    all_data,
    vars = path_ls[1:2] %>% unlist,
    cells = colnames(all_data)[all_data$project == 'A02994A6']
  ) %>%
  kmeans(centers = 2)



p_sp_iz_not_used <-
  cell_coords_ls2$A02994A6 %>%
  left_join(
    enframe(res_kmeans_not_used$cluster, "cell_id", "cluster"),
    by = "cell_id"
  ) %>%
  
  slice_sample(n = 50000) %>%
  
  ggplot(aes(x = x, y = -y, color = as.character(cluster))) +
  geom_point() +
  # scale_color_continuous(trans = "reverse") +
  theme_minimal() +
  scale_color_manual(values = c(
    '2' = "red",
    '1' = 'blue'
  ))

ggsave(
  p_sp_iz_not_used,
  filename = file.path(out_dir_v4, "kmeans cluster IZ and BZ using marker gene - not used combined.pdf"),
  height = 20,
  width  = 20
)



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
  all_data_score %>%
  FetchData(
    vars =  c('project', 'group', 'roi', 'cluster', gene_list)
  )

all_dat_to_plot %>% head

all_data_score[['SCT']] <- NULL
all_data_score[['sketch']] <- NULL
all_data_score <- JoinLayers(all_data_score)
all_data_score_slide1 <- all_data_score[, all_data_score$project == "A02994A6"]
all_data_score_slide2 <- all_data_score[, all_data_score$project == "B04866A5"]

all_data_score_slide1 <- CreateSeuratObject(
  CreateAssay5Object(
    counts = all_data_score_slide1@assays$RNA$counts,
    data   = all_data_score_slide1@assays$RNA$data
  ),
  meta.data = all_data_score_slide1@meta.data
)

all_data_score_slide2 <- CreateSeuratObject(
  CreateAssay5Object(
    counts = all_data_score_slide2@assays$RNA$counts,
    data   = all_data_score_slide2@assays$RNA$data
  ),
  meta.data = all_data_score_slide2@meta.data
)



## groupby group ----------
p_dot_grp_slide1 <- 
  DotPlot(
    all_data_score_slide1,
    features = gene_list,
    group.by = 'group'
  )


p_dot_grp_slide2 <- 
  DotPlot(
    all_data_score_slide1,
    features = gene_list,
    group.by = 'group'
  )

# export 
ggsave(
  p_dot_grp_slide1, 
  filename = file.path(out_dir_v4, "1 selected gene dotplot of group in slide1.pdf"),
  height = 4,
  width = 8
)
ggsave(
  p_dot_grp_slide2, 
  filename = file.path(out_dir_v4, "1 selected gene dotplot of group in slide2.pdf"),
  height = 4,
  width = 8
)

## groupby roi -----------

p_dot_roi_slide1 <- 
  DotPlot(
    all_data_score_slide1,
    features = gene_list,
    group.by = 'roi'
  )


p_dot_roi_slide2 <- 
  DotPlot(
    all_data_score_slide2,
    features = gene_list,
    group.by = 'roi'
  )

# export
ggsave(
  p_dot_roi_slide1, 
  filename = file.path(out_dir_v4, "1 selected gene dotplot of roi in slide1.pdf"),
  height = 4,
  width = 8
)
ggsave(
  p_dot_roi_slide2, 
  filename = file.path(out_dir_v4, "1 selected gene dotplot of roi in slide2.pdf"),
  height = 4,
  width = 8
)




# groupby cluster

p_dot_cluster_slide1 <- 
  DotPlot(
    all_data_score_slide1,
    features = gene_list,
    group.by = 'cluster'
  )


p_dot_cluster_slide2 <- 
  DotPlot(
    all_data_score_slide2,
    features = gene_list,
    group.by = 'cluster'
  )

# export
ggsave(
  p_dot_cluster_slide1, 
  filename = file.path(out_dir_v4, "1 selected gene dotplot of cluster in slide1.pdf"),
  height = 4,
  width = 8
)
ggsave(
  p_dot_cluster_slide2, 
  filename = file.path(out_dir_v4, "1 selected gene dotplot of cluster in slide2.pdf"),
  height = 4,
  width = 8
)


# groupby annotation
p_dot_cell_slide1 <- 
  DotPlot(
    all_data_score_slide1,
    features = gene_list,
    group.by = 'annotation_full'
  )


p_dot_cell_slide2 <- 
  DotPlot(
    all_data_score_slide2,
    features = gene_list,
    group.by = 'annotation_full'
  )


# export
ggsave(
  p_dot_cell_slide1, 
  filename = file.path(out_dir_v4, "1 selected gene dotplot of annotation in slide1.pdf"),
  height = 8,
  width = 8
)
ggsave(
  p_dot_cell_slide2, 
  filename = file.path(out_dir_v4, "1 selected gene dotplot of annotation in slide2.pdf"),
  height = 8,
  width = 8
)


# 3 de analysis between CM -------

all_marker_ls_slide1 <-
  c( "CMs", "Ven-CMs", "Atrial-CMs") %>%
  setNames(., .) %>%
  map(function(cell){
    Idents( all_data_score_slide1 ) <- 'annotation_full'
    
    all_data_score_slide1 %>%
      FindMarkers(
        ident.1 = 'treat',
        ident.2 = 'Veh',
        group.by = 'group',
        subset.ident = cell
      )
  })

all_marker_ls_slide2 <-
  c( "CMs", "Ven-CMs", "Atrial-CMs") %>%
  setNames(., .) %>%
  map(function(cell){
    Idents( all_data_score_slide2 ) <- 'annotation_full'
    
    tryCatch(
      all_data_score_slide2 %>%
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
    file = file.path(out_dir_v4, "3 de enrichment table.xlsx"),
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
  file.path(out_dir_v4, "3 de enrichment analysis dotplot.pdf"), 
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

all_data_score_hallmark <- all_data %>% AddModuleScore(gs_h)

p_hallmark_heat_ls <- 
  all_data_score_hallmark[[]] %>%
  group_by(project) %>% 
  nest() %>%
  deframe() %>%
  imap(function(dat, name){
    dat %>%
      dplyr::select(
        annotation_full,
        starts_with("Cluster")
      ) %>%
      group_by(annotation_full) %>%
      summarise(across(everything(), mean)) %>%
      column_to_rownames("annotation_full") %>%
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
  file.path(out_dir_v4, "4 hallmark heatmap - cell.pdf"), 
  height = 10, 
  width  = 7
)
p_hallmark_heat_ls %>% map('gtable') %>%  imap(~ cowplot::plot_grid(.x, labels = .y, scale = 0.9))
dev.off()



# project  hallmark 
p_hallmark_heat_proj <-
  all_data_score_hallmark[[]] %>%
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
  file.path(out_dir_v4, "4 hallmark heatmap - project.pdf"), 
  height = 10, 
  width  = 7
)
p_hallmark_heat_proj$gtable %>%  cowplot::plot_grid()
dev.off()



# group  hallmark 
p_hallmark_heat_grp_ls <-
  all_data_score_hallmark[[]] %>%
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
  file.path(out_dir_v4, "4 hallmark heatmap - group.pdf"), 
  height = 10, 
  width  = 7
)
p_hallmark_heat_grp_ls %>% map('gtable') %>%  imap(~ cowplot::plot_grid(.x, labels = .y, scale = 0.9))
dev.off()



 















