library(Seurat)
library(tidyverse)

sp_epi2 <- readRDS("sp_subcluster/sp_epi2.rds")

epi_out_dir <- "epi - data analysis"
dir.create(epi_out_dir)

# sp_epi2 <- sp_epi[, sp_epi$annotation_epi != 'Others']

# 1 umap plot of epithelial --------------
p_umap_epi <-
  DimPlot(
    sp_epi2,
    group.by = 'annotation_epi',
    reduction = "harmonyumap"
  ) +
  ggsci::scale_color_aaas()

# export
ggsave(
  p_umap_epi,
  filename = file.path(epi_out_dir, "1 umap plot of epi cluster.pdf"),
  height = 6,
  width  = 8
)

p_umap_epi_sample <-
  DimPlot(
    sp_epi2,
    group.by = 'sample',
    reduction = "harmonyumap"
  ) +
  ggsci::scale_color_igv()

# export
ggsave(
  p_umap_epi_sample,
  filename = file.path(epi_out_dir, "1 umap plot of epi cluster - sample.pdf"),
  height = 6,
  width  = 8
)

# 1 cell props --------------------

# sp_epi$group %>% unique
# [1] PHT LRT LRN PHN
# 
# sp_epi$group2 %>% unique
# [1] "PH" "LR"
# 
# sp_epi$annotation_epi %>% unique
# [1] Tumor PH            Hepatocyte          Proliferating tumor Tumor               Tumor LR           
# [6] Cholangiocyte       Others             

p_col_group <-
  sp_epi2[[]] %>%
  group_by(group, annotation_epi) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  mutate(group = factor(group, levels = c('LRN', 'PHN', 'LRT', 'PHT'))) %>%
  
  ggplot(aes(x = group, y = freq, fill = annotation_epi)) +
  geom_col() +
  ggsci::scale_fill_aaas() +
  theme_classic()
  
p_col_group2 <-
  sp_epi2[[]] %>%
  group_by(group2, annotation_epi) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  mutate(group = factor(group2, levels = c('LR', 'PH'))) %>%
  
  ggplot(aes(x = group2, y = freq, fill = annotation_epi)) +
  geom_col() +
  ggsci::scale_fill_aaas() +
  theme_classic()


p_pie_group <-
  sp_epi2[[]] %>%
  group_by(group, annotation_epi) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  mutate(group = factor(group, levels = c('LRN', 'PHN', 'LRT', 'PHT'))) %>%
  
  ggplot(aes(x = 1, y = freq, fill = annotation_epi)) +
  geom_col() +
  ggsci::scale_fill_aaas() +
  theme_void() +
  coord_polar(theta = 'y') +
  facet_wrap(~ group)  


p_pie_group2 <-
  sp_epi2[[]] %>%
  group_by(group2, annotation_epi) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  mutate(group = factor(group2, levels = c('PH', 'LR'))) %>%
  
  ggplot(aes(x = 1, y = freq, fill = annotation_epi)) +
  geom_col() +
  ggsci::scale_fill_aaas() +
  theme_void() +
  coord_polar(theta = 'y') +
  facet_wrap(~ group2)

# export
ggsave(
  p_col_group, 
  filename = file.path(epi_out_dir, "2 cell props - barplot - group.pdf"),
  width = 6,
  height = 6
)

ggsave(
  p_col_group2, 
  filename = file.path(epi_out_dir, "2 cell props - barplot - group2.pdf"),
  width = 4,
  height = 6
)

ggsave(
  p_pie_group, 
  filename = file.path(epi_out_dir, "2 cell props - pie - group.pdf"),
  width = 6,
  height = 6
)

ggsave(
  p_pie_group2, 
  filename = file.path(epi_out_dir, "2 cell props - pie - group2.pdf"),
  width = 6,
  height = 4
)

# 2 cell de analysis --------------
# LR vs PH LR vs T PH vs T


all_markers_all <- 
  FindAllMarkers(
    sp_epi2,
    group.by = "annotation_epi",
    return.thresh = 1.1
  )

write.csv(
  all_markers_all,
  file = file.path(epi_out_dir, '2 de analysis - de table.csv')
)

## 2.1 volcano plot -----------
# scatter
p_volcanoplot_major_ls <- 
  list(
    Total = all_markers_all
  ) %>%
  imap(function(dat, name){
    
    set.seed(1234)
    rect_height <- 1
    dodge_width <- 0.1
    
    df_scatter <-
      dat %>%
      mutate(x_base = as.numeric(factor(cluster))) %>%
      slice_sample(prop = 1) %>%
      group_by(cluster) %>%
      
      mutate(
        x_delta = 1:n(),
        x_delta = scales::rescale(x_delta, to = c(dodge_width, 1 - dodge_width)),
        xx = x_base - 1 + x_delta,
        # y = ifelse(avg_log2FC > 0, avg_log2FC + rect_height, avg_log2FC - rect_height)
        y = avg_log2FC
      ) %>%
      mutate(
        group = ifelse(p_val_adj > 0.01, "p > 0.01", "p < 0.01"),
        group = factor(group, levels = c("p > 0.01", "p < 0.01"))
      ) %>%
      dplyr::arrange(group)
    
    # df_label 
    df_label_up <- df_scatter %>%
      dplyr::filter(as.character(group) == "p < 0.01") %>%
      group_by(cluster) %>% 
      slice_max(avg_log2FC, n = 10)
    
    df_label_down <- df_scatter %>%
      dplyr::filter(as.character(group) == "p < 0.01") %>%
      group_by(cluster) %>% 
      slice_min(avg_log2FC, n = 10)
    
    
    
    # grp rect  
    df_grp_rect <-
      dat %>%
      group_by(cluster) %>%
      summarise() %>%
      mutate(x_base = as.numeric(factor(cluster))) %>%
      mutate(xmin = x_base - 1, xmax = x_base) %>%
      mutate(ymin = -rect_height, ymax = rect_height)
    
    
    # background rect
    df_bg_rect <-
      dat %>%
      group_by(cluster) %>%
      summarise(ymax = max(avg_log2FC) + rect_height, ymin = min(avg_log2FC) - rect_height) %>%
      mutate(x_base = as.numeric(factor(cluster))) %>%
      mutate(xmin = x_base - 1 + dodge_width, xmax = x_base - dodge_width)
    
    
    
    df_scatter %>%
      ggplot(aes(x = xx, y = y)) +
      geom_rect(
        data = df_bg_rect, 
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
        fill = "gray90",
        inherit.aes = FALSE
      ) +
      geom_point(aes(color = group), size = 1) +
      geom_rect(
        data = df_grp_rect,
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = cluster),
        inherit.aes = FALSE
      ) +
      geom_text(
        data = df_grp_rect,
        aes(x = (xmin + xmax)/2, y = 0, label = cluster), 
        color = "white",
        show.legend = FALSE
      ) +
      ggrepel::geom_text_repel(data = df_label_up, aes(label = gene)) +
      ggrepel::geom_text_repel(data = df_label_down, aes(label = gene)) +
      scale_color_manual(values = c(
        'p < 0.01' = "red",
        'p > 0.01' = "gray70"
      )) +
      ggsci::scale_fill_aaas() +
      theme_minimal() +
      theme(
        panel.grid = element_blank(),
        axis.line.y = element_line(color = "black"),
        axis.text.x = element_blank()
      ) +
      labs(
        x = NULL,
        y = "Log2FC",
        title = name
      )
    
  })


pdf(
  file.path(epi_out_dir, "2 de analysis - volcano plot - all epi cluster.pdf"), 
  height = 5, 
  width  = 8
)
p_volcanoplot_major_ls %>%  walk(print)
dev.off()

## 2.2 enrichment analysis -------------

all_geneid_all_ls <-
  list(Total = all_markers_all) %>%
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
          OrgDb    = "org.Hs.eg.db"
        ) %>%
          .$ENTREZID
      })
    
    gene  
  })

all_geneid_all_ls <- all_geneid_all_ls$Total

all_enrich_all_ls <- 
  all_geneid_all_ls %>%
  map(function(x){
    clusterProfiler::enrichGO(
      gene = x,
      ont = "BP",
      OrgDb = "org.Hs.eg.db",
      pvalueCutoff  = 1, # 0.05,
      pAdjustMethod = "BH",
      qvalueCutoff  = 1 # 0.2
    )
  })


# enrichment 
all_enrich_all_ls <-
  all_enrich_all_ls %>%
  map(function(x){ # browser()
    xx <- x@result$geneID %>% str_split("/")
    
    all_genes <- unlist(xx)
    gene_db <- 
      clusterProfiler::bitr(
        all_genes,
        fromType = "ENTREZID",
        toType   = "SYMBOL",
        OrgDb    = "org.Hs.eg.db"
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
p_all_enrich_all_ls <-
  all_enrich_all_ls %>%
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


# export
all_enrich_all_ls %>%
  map('result') %>%
  openxlsx::write.xlsx(
    file = file.path(epi_out_dir, "2 de analysis - enrichment table.xlsx"),
    rowNames = FALSE
  )

pdf(
  file.path(epi_out_dir, "2 de analysis - dotplot - all epi cluster.pdf"), 
  height = 6.5, 
  width  = 8
)
p_all_enrich_all_ls %>%  walk(print)
dev.off()


## 2.3 PH vs Tumor and PH vs LR and LR vs Tumor -------------
ph_idx    <- colnames(sp_epi2)[sp_epi2$annotation_epi == "Tumor PH"]
lr_idx    <- colnames(sp_epi2)[sp_epi2$annotation_epi == "Tumor LR"]
tumor_idx <- colnames(sp_epi2)[sp_epi2$annotation_epi == "Tumor"]

de_ph_vs_tumor <-
  FindMarkers(
    sp_epi2[['RNA']],
    cells.1 = ph_idx,
    cells.2 = tumor_idx
  )

de_lr_vs_tumor <-
  FindMarkers(
    sp_epi2[['RNA']],
    cells.1 = lr_idx,
    cells.2 = tumor_idx
  )

de_lr_vs_ph <-
  FindMarkers(
    sp_epi2[['RNA']],
    cells.1 = lr_idx,
    cells.2 = ph_idx
  )

all_markers_tumor <- 
  list(
    'PH vs Tumor' = de_ph_vs_tumor,
    'LR vs Tumor' = de_lr_vs_tumor,
    'LR vs PH' = de_lr_vs_ph
  ) %>%
  map(rownames_to_column, 'gene') %>%
  enframe(name = "cluster") %>%
  unnest("value")

# export
write.csv(
  all_markers_tumor,
  file = file.path(epi_out_dir, '2 de analysis - de table - tumor.csv')
)

### 2.3.1 volcano plot ----------

# scatter
p_volcanoplot_tumor_ls <- 
  list(
    Tumor = all_markers_tumor
  ) %>%
  imap(function(dat, name){
    
    set.seed(1234)
    rect_height <- 1
    dodge_width <- 0.1
    
    df_scatter <-
      dat %>%
      mutate(x_base = as.numeric(factor(cluster))) %>%
      slice_sample(prop = 1) %>%
      group_by(cluster) %>%
      
      mutate(
        x_delta = 1:n(),
        x_delta = scales::rescale(x_delta, to = c(dodge_width, 1 - dodge_width)),
        xx = x_base - 1 + x_delta,
        # y = ifelse(avg_log2FC > 0, avg_log2FC + rect_height, avg_log2FC - rect_height)
        y = avg_log2FC
      ) %>%
      mutate(
        group = ifelse(p_val_adj > 0.01, "p > 0.01", "p < 0.01"),
        group = factor(group, levels = c("p > 0.01", "p < 0.01"))
      ) %>%
      dplyr::arrange(group)
    
    # df_label 
    df_label_up <- df_scatter %>%
      dplyr::filter(as.character(group) == "p < 0.01") %>%
      group_by(cluster) %>% 
      slice_max(avg_log2FC, n = 10)
    
    df_label_down <- df_scatter %>%
      dplyr::filter(as.character(group) == "p < 0.01") %>%
      group_by(cluster) %>% 
      slice_min(avg_log2FC, n = 10)
    
    
    
    # grp rect  
    df_grp_rect <-
      dat %>%
      group_by(cluster) %>%
      summarise() %>%
      mutate(x_base = as.numeric(factor(cluster))) %>%
      mutate(xmin = x_base - 1, xmax = x_base) %>%
      mutate(ymin = -rect_height, ymax = rect_height)
    
    
    # background rect
    df_bg_rect <-
      dat %>%
      group_by(cluster) %>%
      summarise(ymax = max(avg_log2FC) + rect_height, ymin = min(avg_log2FC) - rect_height) %>%
      mutate(x_base = as.numeric(factor(cluster))) %>%
      mutate(xmin = x_base - 1 + dodge_width, xmax = x_base - dodge_width)
    
    
    
    df_scatter %>%
      ggplot(aes(x = xx, y = y)) +
      geom_rect(
        data = df_bg_rect, 
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
        fill = "gray90",
        inherit.aes = FALSE
      ) +
      geom_point(aes(color = group), size = 1) +
      geom_rect(
        data = df_grp_rect,
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = cluster),
        inherit.aes = FALSE
      ) +
      geom_text(
        data = df_grp_rect,
        aes(x = (xmin + xmax)/2, y = 0, label = cluster), 
        color = "white",
        show.legend = FALSE
      ) +
      ggrepel::geom_text_repel(data = df_label_up, aes(label = gene)) +
      ggrepel::geom_text_repel(data = df_label_down, aes(label = gene)) +
      scale_color_manual(values = c(
        'p < 0.01' = "red",
        'p > 0.01' = "gray70"
      )) +
      ggsci::scale_fill_aaas() +
      theme_minimal() +
      theme(
        panel.grid = element_blank(),
        axis.line.y = element_line(color = "black"),
        axis.text.x = element_blank()
      ) +
      labs(
        x = NULL,
        y = "Log2FC",
        title = name
      )
    
  })


pdf(
  file.path(epi_out_dir, "2 de analysis - volcano plot - tumor cluster.pdf"), 
  height = 5, 
  width  = 10
)
p_volcanoplot_tumor_ls %>%  walk(print)
dev.off()

### 2.3.2 enrichment analysis -------------

all_markers_tumor2 <- 
  list(
    'PH vs Tumor' = de_ph_vs_tumor,
    'LR vs Tumor' = de_lr_vs_tumor,
    'LR vs PH' = de_lr_vs_ph,
    'PH vs LR' = de_lr_vs_ph %>% mutate(avg_log2FC = -avg_log2FC)
  ) %>%
  map(rownames_to_column, 'gene') %>%
  enframe(name = "cluster") %>%
  unnest("value")


all_geneid_tumor_ls <-
  list(Tumor = all_markers_tumor2) %>%
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
          OrgDb    = "org.Hs.eg.db"
        ) %>%
          .$ENTREZID
      })
    
    gene  
  })

all_geneid_tumor_ls <- all_geneid_tumor_ls$Tumor

all_enrich_tumor_ls <- 
  all_geneid_tumor_ls %>%
  map(function(x){
    clusterProfiler::enrichGO(
      gene = x,
      ont = "BP",
      OrgDb = "org.Hs.eg.db",
      pvalueCutoff  = 1, # 0.05,
      pAdjustMethod = "BH",
      qvalueCutoff  = 1 # 0.2
    )
  })


# enrichment 
all_enrich_tumor_ls <-
  all_enrich_tumor_ls %>%
  map(function(x){ # browser()
    xx <- x@result$geneID %>% str_split("/")
    
    all_genes <- unlist(xx)
    gene_db <- 
      clusterProfiler::bitr(
        all_genes,
        fromType = "ENTREZID",
        toType   = "SYMBOL",
        OrgDb    = "org.Hs.eg.db"
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
p_all_enrich_tumor_ls <-
  all_enrich_tumor_ls %>%
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


# export
all_enrich_tumor_ls %>%
  map('result') %>%
  openxlsx::write.xlsx(
    file = file.path(epi_out_dir, "2 de analysis - enrichment table - tumor.xlsx"),
    rowNames = FALSE
  )

pdf(
  file.path(epi_out_dir, "2 de analysis - dotplot - tumor cluster.pdf"), 
  height = 6.5, 
  width  = 8
)
p_all_enrich_tumor_ls %>%  walk(print)
dev.off()



# 3 infercnv score ------------------
infercnv_obj <- readRDS('infercnv - run3/infercnv_obj2.rds')

cnv_socre <- colMeans((infercnv_obj@expr.data - 1) ^ 2) %>% enframe('cell_id', "score")
cnv_socre_df <-
  c(
    infercnv_obj@observation_grouped_cell_indices,
    infercnv_obj@reference_grouped_cell_indices
  ) %>%
  map(function(x) colnames(infercnv_obj@expr.data)[x]) %>%
  enframe("cluster", "cell_id") %>%
  unnest(cell_id) %>%
  left_join(
    cnv_socre,
    by = c("cell_id")
  )

cnv_socre_df %>%
  ggplot(aes(x = cluster, y = log(score))) +
  geom_boxplot()




