library(tidyverse)
library(Seurat)


# flatDoubleList
flatDoubleList <- function(x, sep = "---"){
  stopifnot(is.vector(x))
  
  res_ls <- list()
  
  # layer 1
  name_1 <- names(x)
  id_1   <- seq_along(x)
  if(is.null(name_1)) name_1 <- id_1
  
  for(i in id_1){
    y <- x[[i]]
    stopifnot(is.vector(y))
    
    # layer 2
    name_2 <- names(y)
    id_2   <- seq_along(y)
    if(is.null(name_2)) name_2 <- id_2
    
    for(j in id_2){
      res  <- y[[j]]
      name <- paste0(name_1[i], sep, name_2[j])
      res_ls[[name]] <- res
    }
  }
  res_ls
}



all_data <- readRDS("all_data_projected.rds")

out_dir <- "20251017"
dir.create(out_dir)

# 1 cell porps ---------

p_col_grp <- 
  all_data[[]] %>%
  group_by(project, group, annotation_full) %>%
  summarise(count = n()) %>%
  group_by(project, group) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = group, y = freq, fill = annotation_full)) +
  geom_col() +
  facet_wrap(~ project) +
  ggsci::scale_fill_aaas() +
  theme_bw()


p_col_sample <- 
  all_data[[]] %>%
  mutate(roi = str_c(roi, group, sep = ":")) %>%
  group_by(project, roi, annotation_full) %>%
  summarise(count = n()) %>%
  group_by(project, roi) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = roi, y = freq, fill = annotation_full)) +
  geom_col() +
  facet_wrap(~ project, scales = "free_x") +
  ggsci::scale_fill_aaas() +
  theme_bw()

# export
ggsave(
  p_col_grp,
  filename = file.path(out_dir, "1. bar plot of cell props - group.pdf"),
  height = 6, 
  width  = 6
)

ggsave(
  p_col_sample,
  filename = file.path(out_dir, "1. bar plot of cell props - sample.pdf"),
  height = 6, 
  width  = 10
)


# 2 spatial plot -----------
cells_boundary_ls <- readRDS('cells_boundary_ls.rds')

sp_p_ls <-
  cells_boundary_ls %>%
  imap(function(cells_boundary, project){ # browser()
    
    all_meta <- all_data[[]] %>% rownames_to_column("cell_id") %>% dplyr::select(-x, -y)
    dat <- 
      cells_boundary %>%
      inner_join(all_meta, by = "cell_id")
    
    dat_sample <-
      dat %>% 
      mutate(roi = str_c(roi, group, sep = ":")) %>% 
      group_by(roi) %>% 
      summarise(x = min(x), y = min(y))
    
    p <- dat %>%
      ggplot(aes(x = x, y = y)) +
      geom_polygon(aes(fill = annotation_full, group = cell_id))  +
      geom_text(data = dat_sample, aes(x = x, y = y, label = roi)) +
      scale_y_reverse() +
      # coord_fixed(ratio = 1/image_wh_r) +
      coord_equal() +
      ggsci::scale_fill_igv() +
      theme_void()
    p
  })


# export image
sp_dir <- file.path(out_dir, "spatial plot") 
dir.create(sp_dir)

sp_p_ls %>% 
  iwalk(function(p, grp){
    pdf(
      file.path(sp_dir, str_glue("spatial plot {grp}.pdf")),
      height = 12,
      width  = 12
    )
    print(p)
    dev.off()
  })



# 3 geneset score ---------
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


all_data_score <- 
  all_data %>%
  AddModuleScore(path_ls)

colnames(all_data_score@meta.data)[str_detect(colnames(all_data_score@meta.data), "Cluster")] <- names(path_ls)

# path score
sp_p_ls_score <-
  cells_boundary_ls %>%
  imap(function(cells_boundary, project){ # browser()
    
    all_meta <- all_data_score[[]] %>% rownames_to_column("cell_id") %>% dplyr::select(-x, -y)
    dat <- 
      cells_boundary %>%
      inner_join(all_meta, by = "cell_id")
    
    dat_sample <-
      dat %>% 
      mutate(roi = str_c(roi, group, sep = ":")) %>% 
      group_by(roi) %>% 
      summarise(x = min(x), y = min(y))
    
    c('BZ', 'IZ', 'ISG') %>%
      setNames(., .) %>%
      map(function(path){
        val_rng <- quantile(dat[[path]], c(0.05, 0.95))
        print(val_rng)
        dat[[path]] <- pmin(dat[[path]], val_rng[2])
        dat[[path]] <- pmax(dat[[path]], val_rng[1])

        p <- dat %>%
          ggplot(aes(x = x, y = y)) +
          geom_polygon(aes(fill = .data[[path]], group = cell_id))  +
          geom_text(data = dat_sample, aes(x = x, y = y, label = roi)) +
          scale_y_reverse() +
          # coord_fixed(ratio = 1/image_wh_r) +
          coord_equal() +
          scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(Inf, "RdBu"))) + 
          theme_void() +
          ggtitle(path)
        p
      })
  })

# export geneset score
sp_p_ls_score %>% 
  iwalk(function(p_ls, grp){
    pdf(
      file.path(sp_dir, str_glue("spatial plot {grp} - score - v2.pdf")),
      height = 12,
      width  = 12
    )
    print(p_ls)
    dev.off()
  })


# 4 de analysis of celltype ------------
DefaultAssay(all_data) <- 'RNA'
all_markers_celltype <-
  FindAllMarkers(
    JoinLayers(all_data),
    group.by = "annotation_full", 
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    only.pos = TRUE
  )

all_geneid_ct_ls <-
  list(celltype = all_markers_celltype) %>%
  map(function(dat){
    gene <-
      dat %>%
      dplyr::filter(avg_log2FC > 0.25, p_val_adj < 0.05, pct.1 > 0.1) %>%
      group_by(cluster) %>%
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

all_geneid_ct_ls <- all_geneid_ct_ls$celltype

all_enrich_ct_ls <- 
  all_geneid_ct_ls %>%
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
all_enrich_ct_ls <-
  all_enrich_ct_ls %>%
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
p_all_enrich_ct_ls <-
  all_enrich_ct_ls %>%
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
  file.path(out_dir, "4 de enrichment analysis dotplot.pdf"), 
  height = 6.5, 
  width  = 8
)
p_all_enrich_ct_ls %>%  walk(print)
dev.off()

# 5 de analysis of group --------
dir_enrich <- file.path(out_dir, "enrichment analysis")
dir.create(dir_enrich)


DefaultAssay(all_data) <- 'RNA'
all_data <- JoinLayers(all_data)
all_marker_ls <-
  all_data$annotation_full %>%
  unique() %>%
  setNames(., .) %>%
  
  map(function(x){
    all_data %>% 
      {Idents(.) <- "annotation_full"; .} %>% 
      FindMarkers(
        ident.1 = 'treat',
        ident.2 = 'Veh',
        group.by = "group",
        subset.ident = x,
        logfc.threshold = 0
      )
  })


all_marker <- 
  all_marker_ls %>%
  map(rownames_to_column, "gene") %>% 
  enframe("cluster") %>%
  unnest("value")



all_de_genes_for_enrichmenet <-
  all_marker %>% 
  group_by(cluster) %>%
  nest() %>%
  deframe() %>%
  map(function(x){
    up_gene <- 
      x %>%
      dplyr::filter(
        p_val_adj < 0.05,
        # pct.1 - pct.2 > 0.1,
        avg_log2FC > 0.25
      ) %>%
      .$gene
    
    down_gene <- 
      x %>% 
      dplyr::filter(
        p_val_adj < 0.05,
        # pct.1 - pct.2 > 0.1,
        avg_log2FC < -0.25
      ) %>%
      .$gene
    
    # return
    list(
      treat = up_gene,
      Veh = down_gene
    )
  }) %>%
  flatDoubleList()


## GOBP -------
all_de_genes_id_for_enrichmenet <-
  all_de_genes_for_enrichmenet %>%
  map(function(x){
    clusterProfiler::bitr(
      x,
      OrgDb = "org.Mm.eg.db",
      fromType = "SYMBOL",
      toType = "ENTREZID"
    ) %>%
      .$ENTREZID
  }) %>%
  discard(function(x) length(x) < 1)

gobp_enrich_ls <-
  all_de_genes_id_for_enrichmenet %>%
  map(function(x){
    clusterProfiler::enrichGO(
      x,
      OrgDb = 'org.Mm.eg.db',
      ont = "BP",
      pvalueCutoff = 0.05, 
      pAdjustMethod = "BH", 
      minGSSize = 10,
      maxGSSize = 500, 
      qvalueCutoff = 0.2
    )
  })

# saveRDS(gobp_enrich_ls, file.path(dir_enrich, 'gobp_enrich_ls.rds'))

# expot enrichment table
gobp_enrich_ls %>%
  discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    
    write.csv(
      x@result,
      file = file.path(dir_enrich, str_glue("GOBP - enriment table - {name}.csv")),
      row.names = FALSE
    )
  })


# dotplot
gobp_enrich_dot_ls <-
  gobp_enrich_ls %>%
  imap(function(x, name){
    if(nrow(as.data.frame(x)) < 1){
      return(NULL)
    }
    enrichplot::dotplot(
      x,
      showCategory = 10,
      label_format = 50
    ) +
      ggtitle(name)
  })

# export dotplot
gobp_enrich_dot_ls %>%
  discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    ggsave(
      x,
      filename = file.path(dir_enrich, str_glue("BOBP - dotplot - {name}.pdf")),
      height = 6,
      width = 8
    )
  })



# 6 hallmark ---------
DefaultAssay(all_data) <- "RNA"
gs_h <- msigdbr::msigdbr(
  db_species = "HS",
  species = "Mus musculus", # "Homo sapiens",
  collection = "H"
) %>%
  with(split(gene_symbol, gs_name))

all_data_score_hallmark <- all_data %>% AddModuleScore(gs_h)

p_hallmark_heat <- 
  all_data_score_hallmark[[]] %>%
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
    cluster_cols = FALSE,
    labels_row = names(gs_h),
    silent = TRUE
  )}

pdf(
  file.path(out_dir, "6 hallmark heatmap - cell.pdf"), 
  height = 10, 
  width  = 7
)
p_hallmark_heat$gtable %>%  cowplot::plot_grid()
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
  file.path(out_dir, "6 hallmark heatmap - project.pdf"), 
  height = 10, 
  width  = 7
)
p_hallmark_heat_proj$gtable %>%  cowplot::plot_grid()
dev.off()



# group  hallmark 
p_hallmark_heat_grp <-
  all_data_score_hallmark[[]] %>%
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

pdf(
  file.path(out_dir, "6 hallmark heatmap - group.pdf"), 
  height = 10, 
  width  = 7
)
p_hallmark_heat_grp$gtable %>%  cowplot::plot_grid()
dev.off()



# group  hallmark 
p_hallmark_heat_time <-
  all_data_score_hallmark[[]] %>%
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
  file.path(out_dir, "6 hallmark heatmap - time.pdf"), 
  height = 10, 
  width  = 7
)
p_hallmark_heat_time$gtable %>%  cowplot::plot_grid()
dev.off()


# 7 spatial plot of gene ------------
cells_boundary_ls <- readRDS('cells_boundary_ls.rds')

gene_for_plot <- c(
  "Cxcr5", "Atf3", "Myc", "Notch2", "Cxcl10"
)

sp_p_ls_gene <-
  cells_boundary_ls %>%
  imap(function(cells_boundary, project){ # browser()
    
    all_meta <- FetchData(all_data, vars = c(gene_for_plot, "roi", "group", "project")) %>% rownames_to_column("cell_id")
    
    dat <- 
      cells_boundary %>%
      inner_join(all_meta, by = "cell_id")
    
    dat_sample <-
      dat %>% 
      mutate(roi = str_c(roi, group, sep = ":")) %>% 
      group_by(roi) %>% 
      summarise(x = min(x), y = min(y))
    
    gene_for_plot %>%
      setNames(., .) %>%
      map(function(path){
        p <- dat %>%
          ggplot(aes(x = x, y = y)) +
          geom_polygon(aes(fill = .data[[path]], group = cell_id))  +
          geom_text(data = dat_sample, aes(x = x, y = y, label = roi)) +
          scale_y_reverse() +
          # coord_fixed(ratio = 1/image_wh_r) +
          coord_equal() +
          scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(Inf, "RdBu"))) + 
          theme_void() +
          ggtitle(path)
        p
      })
  })

# export geneset score
sp_p_ls_gene %>% 
  iwalk(function(p_ls, grp){
    pdf(
      file.path(sp_dir, str_glue("spatial plot {grp} - gene.pdf")),
      height = 12,
      width  = 12
    )
    print(p_ls)
    dev.off()
  })







