library(Seurat)
library(tidyverse)

major_marker_ls <- c(
  'epithelials'  = 'EPCAM',
  'epithelials'  = 'KRT5',
  'epithelials'  = 'KRT14',
  
  "Foveolar cells" = "GKN1", #ACT
  "Foveolar cells" = "MUC5AC",
  "Foveolar cells" = "TFF1",
  
  'endothelials' = 'PECAM1', 
  'endothelials' = 'VWF', 
  
  "Fibroblasts" = "PDGFRA", 
  "Fibroblasts" = "CCDC80",
  "Fibroblasts" = "SFRP4",
  "Fibroblasts" = "COL1A1",
  "Fibroblasts" = "COL1A2",
  "Fibroblasts" = "DCN",
  "Fibroblasts" = "ACTA2",
  
  "B cells" = "CD19",
  "B cells" = "CD79A",
  "B cells" = "MS4A1",
  "Plasma" = "TNFRSF17",
  "Plasma" = "JCHAIN",
  "Plasma" = "MZB1",
  "Plasma" = "XBP1",
  
  'T cells'           = 'CD3D', 
  'T cells'           = 'CD3E', 
  'T cells'           = 'CD3G', 
  'T cells'           = 'TRAC', 
  
  "NK" = "GNLY",
  "NK" = "NKG7",
  "NK" = "NCAM1",
  
  'Proliferating' = 'MKI67',
  'Proliferating' = 'TOP2A',
  
  "Pericytes" = "PDGFRB",
  "Pericytes" = "RGS5",
  "Pericytes" = "MCAM",
  
  "Monocytes" = "CD14", 
  "Monocytes" = "FCGR3A",
  
  "Macrophage" = "CD68",
  "Macrophage" = "CD163",
  "Macrophage" = "C1QA",
  "Macrophage" = "C1QC", 
  'Macrophage'  = 'CSF1R',
  
  "Dendritic cell" = "CD83",
  'Dendritic cell' = "LILRA4", 
  'pDCs' = "IRF4", 
  'pDCs' = "IRF7", 
  'pDCs' = "IRF8", 
  'Dendritic cell' = "CLEC9A",
  'Dendritic cell' = "IDO1",
  'Dendritic cell' = "CD1C",
  'Dendritic cell' = "FCER1A",
  'Lamp3_DCs' = "LAMP3",
  
  'netrophil' = "CSF3R",
  'netrophil' = "G0S2",
  'netrophil' = "RETNLG",
  
  "Myocytes" = "TNNC1", 
  "Myocytes" = "MYH7", 
  
  "Mast cell" = "CPA3",
  "Mast cell" = "CTSG",
  "Mast cell" = "TPSAB1",
  "Mast cell" = "KIT"
)

theme_90angle <-
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
    strip.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )

Epithelial <- readRDS("Epithelial_subtype_remove_D2/Epithelial.rds")

Epithelial_out_dir <- 'Epithelial_subtype_remove_D2/reanno'
dir.create(Epithelial_out_dir)


## annotation -----------
Epithelial[[]] %>%
  group_by(RNA_snn_res.0.6, sample) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = RNA_snn_res.0.6, y = freq, fill = sample)) + 
  geom_col() +
  ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )
ggsave(file.path(Epithelial_out_dir, 'cluster prop by sample_scale.pdf'), width = 5, height = 6)

Epithelial[[]] %>%
  group_by(sample, RNA_snn_res.0.6) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = RNA_snn_res.0.6, y = freq, fill = sample)) + 
  geom_col() +
  ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )
ggsave(file.path(Epithelial_out_dir, 'cluster prop by sample.pdf'), width = 5, height = 6)

Epithelial[[]] %>%
  group_by(sample, group2, RNA_snn_res.0.6) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = RNA_snn_res.0.6, y = freq, fill = group2)) + 
  geom_boxplot() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )
ggsave(file.path(Epithelial_out_dir, 'boxplot_cluster prop by sample.pdf'), width = 8, height = 6)


Epithelial[[]] %>%
  group_by(RNA_snn_res.0.6, group2) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = RNA_snn_res.0.6, y = freq, fill = group2)) + 
  geom_col() +
  ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )
ggsave(file.path(Epithelial_out_dir, 'cluster prop by group2_scale.pdf'), width = 5, height = 6)

Epithelial[[]] %>%
  group_by(group2, RNA_snn_res.0.6) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = RNA_snn_res.0.6, y = freq, fill = group2)) + 
  geom_col() +
  ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )
ggsave(file.path(Epithelial_out_dir, 'cluster prop by group2.pdf'), width = 5, height = 6)

Epithelial$RNA_snn_res.0.6 %>% table()
table(Epithelial$RNA_snn_res.0.6, Epithelial$sample)

cell_cluster <- c(
  '0' = "Sensitive tumor cell", # ERBB2
  '1' = "Sensitive tumor cell", # ERBB2
  '2' = "Sensitive tumor cell", # EXO1, https://pubmed.ncbi.nlm.nih.gov/33569437/
  '3' = "Sensitive tumor cell", # AQP5, https://www.sciencedirect.com/science/article/abs/pii/S0027510724000186
  '4' = "low-nCounts", # BMP7, https://zhuanlan.zhihu.com/p/686585336
  '5' = "mix", # Fibro
  '6' = "Sensitive tumor cell", # CA9, https://rs.yiigle.com/CN371439201706/993396.htm
  '7' = "Sensitive tumor cell", # PTPRZ1, https://pmc.ncbi.nlm.nih.gov/articles/PMC6548002/
  '8' = "Sensitive tumor cell", # BMP7, https://www.sciencedirect.com/science/article/abs/pii/S0304383508007982
  '9' = "Sensitive tumor cell", # KIF20A, https://xb.xzhmu.edu.cn/cn/article/pdf/preview/10.3969/j.issn.2096-3882.2024.05.010.pdf
  '10' = "Metasis tumor cell", # NPY1R, https://www.nature.com/articles/s41598-022-05949-7
  '11' = "Resistant tumor cell", # ERBB2 
  '12' = "Sensitive tumor cell", # CLEC7A , 10.1007/s10238-025-01870-1
  '13' = "Sensitive tumor cell", # ERBB2 
  '14' = "low-nCounts", 
  '15' = "low-nCounts", 
  '16' = "Sensitive tumor cell", # PTPRZ1, https://pmc.ncbi.nlm.nih.gov/articles/PMC6548002/
  '17' = "noise",  # EPCAM不表达 
  '18' = "low-nCounts", 
  '19' = "low-nCounts", 
  '20' = "Sensitive tumor cell", # CA9, https://rs.yiigle.com/CN371439201706/993396.htm
  '21' = "noise"  # EPCAM不表达
)

Idents(Epithelial) <- "RNA_snn_res.0.6"

Epithelial <- RenameIdents(Epithelial, cell_cluster)
Epithelial$subtype <- Idents(Epithelial)

pdf(file.path(Epithelial_out_dir,str_glue("Dimplot dim{dims} annotation.pdf")), height = 6, width = 10)
DimPlot(Epithelial, group.by = "subtype", label = TRUE, reduction = 'harmonyumap',repel = TRUE)
dev.off()

pdf(file.path(Epithelial_out_dir, str_glue("Dimplot dim{dims} annotation split.pdf")), height = 7, width = 18)
DimPlot(Epithelial, group.by = "subtype",split.by = 'group', label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()


library(future)
plan(multisession, workers = 10)
options(future.globals.maxSize = 100 * 1024 ^ 3) # 80Gb

Epithelial_pure <- subset(JoinLayers(Epithelial), subtype %in% c("Sensitive tumor cell", "Resistant tumor cell", "Metasis tumor cell"))
Epithelial_pure <- FindVariableFeatures(Epithelial_pure)
Epithelial_pure <- ScaleData(Epithelial_pure, vars.to.regress = 'nCount_RNA')
Epithelial_pure <- RunPCA(Epithelial_pure)

pdf(file.path(Epithelial_out_dir, 'pca plot - pure.pdf'), height = 6, width = 8)
ElbowPlot(Epithelial_pure, ndims = 50)
dev.off()

dims <- 15
Epithelial_pure <- RunUMAP(Epithelial_pure, dims = 1:dims, return.model = TRUE)
DimPlot(Epithelial_pure, group.by = "subtype", label = TRUE, reduction = 'harmonyumap',repel = TRUE)

Epithelial_pure <- Epithelial_pure %>% 
  harmony::RunHarmony(
    dims.use = 1:dims,
    group.by.vars = "sample",
    reduction.save = "harmony"
    # theta = 3, # default 2, bigger value -- more diversity per cluster
    # lambda = 0.5, # default, smaller value -- more integrated
  )

Epithelial_pure <- RunUMAP(
  Epithelial_pure, 
  dims = 1:dims, 
  reduction = "harmony", 
  reduction.name = "harmonyumap", 
  reduction.key  = "harmonyUMAP_",
  return.model = TRUE
)

DimPlot(Epithelial_pure, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = 'harmonyumap',repel = TRUE)
DimPlot(Epithelial_pure, group.by = "subtype", label = TRUE, reduction = 'harmonyumap',repel = TRUE)

pdf(file.path(Epithelial_out_dir,str_glue("Dimplot dim{dims} annotation - pure.pdf")), height = 6, width = 10)
DimPlot(Epithelial_pure, group.by = "subtype", label = TRUE, reduction = 'harmonyumap',repel = TRUE)
dev.off()

pdf(file.path(Epithelial_out_dir, str_glue("Dimplot dim{dims} annotation split - pure.pdf")), height = 7, width = 16)
DimPlot(Epithelial_pure, group.by = "subtype",split.by = 'group', label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

saveRDS(Epithelial, file.path(Epithelial_out_dir, 'Epithelial_reanno.rds'))




## cell props barplot  --------
Epithelial_meta <- Epithelial_pure[[]]

p_col_proj <-
  Epithelial_meta %>%
  group_by(sample, subtype) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = sample, y = freq, fill = subtype)) + 
  geom_col() +
  ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )


# export
ggsave(
  p_col_proj,
  filename = file.path(Epithelial_out_dir, "cell_freq_sample.pdf"),
  height = 6,
  width  = 6
)

p_col_proj <-
  Epithelial_meta %>%
  group_by(group, subtype) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = group, y = freq, fill = subtype)) + 
  geom_col() +
  ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )

ggsave(
  p_col_proj,
  filename = file.path(Epithelial_out_dir, "cell_freq_group.pdf"),
  height = 6,
  width  = 4
)

p_col_proj <-
  Epithelial_meta %>%
  group_by(group2, subtype) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = group2, y = freq, fill = subtype)) + 
  geom_col() +
  ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )

ggsave(
  p_col_proj,
  filename = file.path(Epithelial_out_dir, "cell_freq_group2.pdf"),
  height = 6,
  width  = 4
)


cellType_marker_ls <- c(
  "Sensitive tumor cell" = 'ERBB2', 
  # 'Tumor cell' = 'XBP1', "Sensitive tumor cell", "Resistant tumor cell", "Metasis tumor cell"
  
  "Sensitive tumor cell" = 'EXO1', 
  "Sensitive tumor cell" = 'AQP5', 
  "Sensitive tumor cell" = 'BMP7', 
  "Sensitive tumor cell" = 'CA9', 
  "Sensitive tumor cell" = 'PTPRZ1', 
  "Sensitive tumor cell" = 'KIF20A', 
  "Sensitive tumor cell" = 'CLEC7A', 
  
  "Metasis tumor cell" = 'NPY1R', 
  "Resistant tumor cell" = 'MDM2'
)

pdf(file.path(Epithelial_out_dir, "dotplot celltype main_marker.pdf"), width = 10, height = 6)
DotPlot(Epithelial, group.by = "subtype", features = cellType_marker_ls, dot.scale = 10) + theme_90angle
dev.off()

pdf(file.path(Epithelial_out_dir, "dotplot celltype main_marker - pure.pdf"), width = 10, height = 5)
DotPlot(Epithelial_pure, group.by = "subtype", features = cellType_marker_ls, dot.scale = 10) + theme_90angle
dev.off()

saveRDS(Epithelial_pure, file.path(Epithelial_out_dir, 'Epithelial_pure_reanno.rds'))




# enrichment subtype -------------
Idents(Epithelial_pure) <- "subtype"
Epithelial_pure_markers <-
  FindAllMarkers(
    JoinLayers(Epithelial_pure),
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    only.pos = TRUE
  )
write.csv(Epithelial_pure_markers,file.path(Epithelial_out_dir, 'Epithelial_pure_DEGs_subtype.csv'))

library(EnhancedVolcano)
type<-unique(Epithelial_pure_markers$cluster)
volcano_ls <- type %>% 
  setNames(.,.) %>%
  imap(function(type,name){
    markers<- Epithelial_pure_markers[Epithelial_pure_markers$cluster == type, ]
    p<-EnhancedVolcano(markers,
                       lab = rownames(markers),
                       x = 'avg_log2FC',
                       y = 'p_val_adj',
                       title = str_glue('{name} volcano plot'),
                       pCutoff = 10e-2,
                       FCcutoff = 0,
                       pointSize = 3.0,
                       labSize = 6.0)
    p
  })
volcano_ls %>%
  iwalk(function(x, name){
    ggsave(
      x,
      filename = file.path(Epithelial_out_dir, str_glue("{name} diffgene valcano.pdf")),
      height = 6,
      width = 8
    )
  })


all_markers_major <- list(Epithelial_pure = Epithelial_pure_markers)
# scatter
p_volcanoplot_major_ls <- 
  all_markers_major %>%
  imap(function(dat, name){
    
    set.seed(1234)
    rect_height <- 1
    dodge_width <- 0.2
    
    df_scatter <-
      dat %>%
      mutate(x_base = as.numeric(factor(cluster))) %>%
      slice_sample(prop = 1) %>%
      group_by(cluster) %>%
      
      mutate(
        x_delta = 1:n(),
        x_delta = scales::rescale(x_delta, to = c(dodge_width, 1 - dodge_width)),
        xx = x_base - 1 + x_delta,
        y = ifelse(avg_log2FC > 0, avg_log2FC + rect_height, avg_log2FC - rect_height)
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
      geom_rect(
        data = df_grp_rect,
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = cluster),
        inherit.aes = FALSE
      ) +
      geom_point(aes(color = group), size = 1) +
      
      ggrepel::geom_text_repel(data = df_label_up, aes(label = gene)) +
      ggrepel::geom_text_repel(data = df_label_down, aes(label = gene)) +
      scale_color_manual(values = c(
        'p < 0.01' = "red",
        'p > 0.01' = "gray70"
      )) +
      # scale_fill_manual(values = grp_cols_sc) +
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

p_volcanoplot_major_ls %>%
  iwalk(function(x, name){
    ggsave(
      x,
      filename = file.path(Epithelial_out_dir, str_glue("{name} diffgene dot.pdf")),
      height = 8,
      width = 15
    )
  })


Epithelial_pure_de_genes_for_enrichmenet <- 
  Epithelial_pure_markers %>% 
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  dplyr::filter(pct.1-pct.2 > 0.075) %>%
  dplyr::filter(abs(avg_log2FC) > 0.5) %>% 
  group_by(cluster) %>%
  nest() %>%
  deframe() %>%
  map("gene")


enrich_ls <-
  Epithelial_pure_de_genes_for_enrichmenet %>%
  map(function(x){
    clusterProfiler::enrichGO(
      x,
      OrgDb = org.Hs.eg.db::org.Hs.eg.db, 
      ont='ALL',
      pAdjustMethod = 'BH',
      pvalueCutoff = 0.05, 
      qvalueCutoff = 0.05,
      keyType = 'SYMBOL'
    )
  })


Epithelial_pure_er_output <- file.path(Epithelial_out_dir, 'enrichment') 
dir.create(Epithelial_pure_er_output)


# expot enrichment table
enrich_ls %>% 
  # discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    dir.create(file.path(Epithelial_pure_er_output, name))
    
    write.csv(
      x@result,
      file = file.path(Epithelial_pure_er_output, name, "Go enriment table.csv"),
      row.names = FALSE
    )
  })

# dotplot
enrich_dot_ls <-
  enrich_ls %>%
  map(function(x){
    # if(nrow(as.data.frame(x)) < 1){
    #   return(NULL)
    # }
    enrichplot::dotplot(
      x,
      showCategory = 10,
      label_format = 50
    )
  })

enrich_bar_ls <-
  enrich_ls %>%
  map(function(x){
    # if(nrow(as.data.frame(x)) < 1){
    #   return(NULL)
    # }
    barplot(
      x,
      showCategory = 10,
      label_format = 50
    )
  })

# export dotplot
enrich_dot_ls %>%
  # discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    # cluster_name <- paste0("cluster_", name)
    # dir.create(file.path(er_output, cluster_name))
    
    ggsave(
      x,
      filename = file.path(Epithelial_pure_er_output, name, "go dotplot.pdf"),
      height = 6,
      width = 8
    )
  })

enrich_bar_ls %>%
  # discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    # cluster_name <- paste0("cluster_", name)
    # dir.create(file.path(er_output, cluster_name))
    
    ggsave(
      x,
      filename = file.path(Epithelial_pure_er_output, name, "go barplot.pdf"),
      height = 6,
      width = 8
    )
  })

Epithelial_pure_de_genes_id_for_enrichmenet <-
  Epithelial_pure_de_genes_for_enrichmenet %>%
  map(function(x){
    clusterProfiler::bitr(
      x,
      OrgDb = "org.Hs.eg.db",
      fromType = "SYMBOL",
      toType = "ENTREZID"
    ) %>%
      .$ENTREZID
  }) %>%
  purrr::discard(function(x) length(x) < 1)


enrich_ls <-
  Epithelial_pure_de_genes_id_for_enrichmenet %>%
  map(function(x){
    clusterProfiler::enrichKEGG(
      x,
      organism = "hsa", 
      keyType = "ncbi-geneid", 
      pvalueCutoff = 0.05, 
      pAdjustMethod = "BH", 
      minGSSize = 10,
      maxGSSize = 500, 
      qvalueCutoff = 0.2,
      use_internal_data = FALSE
    )
  })

# saveRDS(enrich_ls,file.path(Epithelial_pure_er_output, 'enrich_ls.rds'))

# expot enrichment table
enrich_ls %>%
  # discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    dir.create(file.path(Epithelial_pure_er_output, name))
    
    write.csv(
      x@result,
      file = file.path(Epithelial_pure_er_output, name, "KEGG enriment table.csv"),
      row.names = FALSE
    )
  })

# dotplot
enrich_dot_ls <-
  enrich_ls %>%
  map(function(x){
    if(nrow(as.data.frame(x)) < 1){
      return(NULL)
    }
    enrichplot::dotplot(
      x,
      showCategory = 10,
      label_format = 50
    )
  })

enrich_bar_ls <-
  enrich_ls %>%
  map(function(x){
    if(nrow(as.data.frame(x)) < 1){
      return(NULL)
    }
    barplot(
      x,
      showCategory = 10,
      label_format = 50
    )
  })

# export dotplot
enrich_dot_ls %>%
  purrr::discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    # cluster_name <- paste0("cluster_", name)
    # dir.create(file.path(er_output, cluster_name))
    
    ggsave(
      x,
      filename = file.path(Epithelial_pure_er_output, name, "kegg dotplot.pdf"),
      height = 6,
      width = 8
    )
  })

enrich_bar_ls %>%
  purrr::discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    # cluster_name <- paste0("cluster_", name)
    # dir.create(file.path(er_output, cluster_name))
    
    ggsave(
      x,
      filename = file.path(Epithelial_pure_er_output, name, "kegg barplot.pdf"),
      height = 6,
      width = 8
    )
  })




# enrichment RNA_snn_res.0.6 -------------
Idents(Epithelial_pure) <- "RNA_snn_res.0.6"
Epithelial_pure_markers <-
  FindAllMarkers(
    JoinLayers(Epithelial_pure),
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    only.pos = TRUE
  )
write.csv(Epithelial_pure_markers,file.path(Epithelial_out_dir, 'Epithelial_pure_DEGs_RNA_snn_res.0.6.csv'))

library(EnhancedVolcano)
type<-unique(Epithelial_pure_markers$cluster)
volcano_ls <- type %>% 
  setNames(.,.) %>%
  imap(function(type,name){
    markers<- Epithelial_pure_markers[Epithelial_pure_markers$cluster == type, ]
    p<-EnhancedVolcano(markers,
                       lab = rownames(markers),
                       x = 'avg_log2FC',
                       y = 'p_val_adj',
                       title = str_glue('{name} volcano plot'),
                       pCutoff = 10e-2,
                       FCcutoff = 0,
                       pointSize = 3.0,
                       labSize = 6.0)
    p
  })
volcano_ls %>%
  iwalk(function(x, name){
    ggsave(
      x,
      filename = file.path(Epithelial_out_dir, str_glue("{name} diffgene valcano.pdf")),
      height = 6,
      width = 8
    )
  })


all_markers_major <- list(Epithelial_pure = Epithelial_pure_markers)
# scatter
p_volcanoplot_major_ls <- 
  all_markers_major %>%
  imap(function(dat, name){
    
    set.seed(1234)
    rect_height <- 1
    dodge_width <- 0.2
    
    df_scatter <-
      dat %>%
      mutate(x_base = as.numeric(factor(cluster))) %>%
      slice_sample(prop = 1) %>%
      group_by(cluster) %>%
      
      mutate(
        x_delta = 1:n(),
        x_delta = scales::rescale(x_delta, to = c(dodge_width, 1 - dodge_width)),
        xx = x_base - 1 + x_delta,
        y = ifelse(avg_log2FC > 0, avg_log2FC + rect_height, avg_log2FC - rect_height)
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
      geom_rect(
        data = df_grp_rect,
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = cluster),
        inherit.aes = FALSE
      ) +
      geom_point(aes(color = group), size = 1) +
      
      ggrepel::geom_text_repel(data = df_label_up, aes(label = gene)) +
      ggrepel::geom_text_repel(data = df_label_down, aes(label = gene)) +
      scale_color_manual(values = c(
        'p < 0.01' = "red",
        'p > 0.01' = "gray70"
      )) +
      # scale_fill_manual(values = grp_cols_sc) +
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

p_volcanoplot_major_ls %>%
  iwalk(function(x, name){
    ggsave(
      x,
      filename = file.path(Epithelial_out_dir, str_glue("{name} diffgene dot.pdf")),
      height = 8,
      width = 15
    )
  })


Epithelial_pure_de_genes_for_enrichmenet <- 
  Epithelial_pure_markers %>% 
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  dplyr::filter(pct.1-pct.2 > 0.075) %>%
  dplyr::filter(abs(avg_log2FC) > 0.5) %>% 
  group_by(cluster) %>%
  nest() %>%
  deframe() %>%
  map("gene")


enrich_ls <-
  Epithelial_pure_de_genes_for_enrichmenet %>%
  map(function(x){
    clusterProfiler::enrichGO(
      x,
      OrgDb = org.Hs.eg.db::org.Hs.eg.db, 
      ont='ALL',
      pAdjustMethod = 'BH',
      pvalueCutoff = 0.05, 
      qvalueCutoff = 0.05,
      keyType = 'SYMBOL'
    )
  })


Epithelial_pure_er_output <- file.path(Epithelial_out_dir, 'enrichment_RNA_snn_res.0.6') 
dir.create(Epithelial_pure_er_output)


# expot enrichment table
enrich_ls %>% 
  # discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    dir.create(file.path(Epithelial_pure_er_output, name))
    
    write.csv(
      x@result,
      file = file.path(Epithelial_pure_er_output, name, "Go enriment table.csv"),
      row.names = FALSE
    )
  })

# dotplot
enrich_dot_ls <-
  enrich_ls %>%
  map(function(x){
    # if(nrow(as.data.frame(x)) < 1){
    #   return(NULL)
    # }
    enrichplot::dotplot(
      x,
      showCategory = 10,
      label_format = 50
    )
  })

enrich_bar_ls <-
  enrich_ls %>%
  map(function(x){
    # if(nrow(as.data.frame(x)) < 1){
    #   return(NULL)
    # }
    barplot(
      x,
      showCategory = 10,
      label_format = 50
    )
  })

# export dotplot
enrich_dot_ls %>%
  # discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    # cluster_name <- paste0("cluster_", name)
    # dir.create(file.path(er_output, cluster_name))
    
    ggsave(
      x,
      filename = file.path(Epithelial_pure_er_output, name, "go dotplot.pdf"),
      height = 6,
      width = 8
    )
  })

enrich_bar_ls %>%
  # discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    # cluster_name <- paste0("cluster_", name)
    # dir.create(file.path(er_output, cluster_name))
    
    ggsave(
      x,
      filename = file.path(Epithelial_pure_er_output, name, "go barplot.pdf"),
      height = 6,
      width = 8
    )
  })

Epithelial_pure_de_genes_id_for_enrichmenet <-
  Epithelial_pure_de_genes_for_enrichmenet %>%
  map(function(x){
    clusterProfiler::bitr(
      x,
      OrgDb = "org.Hs.eg.db",
      fromType = "SYMBOL",
      toType = "ENTREZID"
    ) %>%
      .$ENTREZID
  }) %>%
  purrr::discard(function(x) length(x) < 1)


enrich_ls <-
  Epithelial_pure_de_genes_id_for_enrichmenet %>%
  map(function(x){
    clusterProfiler::enrichKEGG(
      x,
      organism = "hsa", 
      keyType = "ncbi-geneid", 
      pvalueCutoff = 0.05, 
      pAdjustMethod = "BH", 
      minGSSize = 10,
      maxGSSize = 500, 
      qvalueCutoff = 0.2,
      use_internal_data = FALSE
    )
  })

# saveRDS(enrich_ls,file.path(Epithelial_pure_er_output, 'enrich_ls.rds'))

# expot enrichment table
enrich_ls %>%
  # discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    dir.create(file.path(Epithelial_pure_er_output, name))
    
    write.csv(
      x@result,
      file = file.path(Epithelial_pure_er_output, name, "KEGG enriment table.csv"),
      row.names = FALSE
    )
  })

# dotplot
enrich_dot_ls <-
  enrich_ls %>%
  map(function(x){
    if(nrow(as.data.frame(x)) < 1){
      return(NULL)
    }
    enrichplot::dotplot(
      x,
      showCategory = 10,
      label_format = 50
    )
  })

enrich_bar_ls <-
  enrich_ls %>%
  map(function(x){
    if(nrow(as.data.frame(x)) < 1){
      return(NULL)
    }
    barplot(
      x,
      showCategory = 10,
      label_format = 50
    )
  })

# export dotplot
enrich_dot_ls %>%
  purrr::discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    # cluster_name <- paste0("cluster_", name)
    # dir.create(file.path(er_output, cluster_name))
    
    ggsave(
      x,
      filename = file.path(Epithelial_pure_er_output, name, "kegg dotplot.pdf"),
      height = 6,
      width = 8
    )
  })

enrich_bar_ls %>%
  purrr::discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    # cluster_name <- paste0("cluster_", name)
    # dir.create(file.path(er_output, cluster_name))
    
    ggsave(
      x,
      filename = file.path(Epithelial_pure_er_output, name, "kegg barplot.pdf"),
      height = 6,
      width = 8
    )
  })





all_data <- readRDS(file.path("R_OUTPUT", "all_data.rds"))
all_data$sample %>% table()
all_data_remove_D2 <- subset(all_data, sample != 'D2')
all_data_remove_D2$sample %>% table()
saveRDS(all_data_remove_D2, 'all_data_remove_D2.rds')




# spatial plot -------------------
## 映射-------------------------------
library(tibble)
library(data.table)
all_data_remove_D2$annotation %>% table()
all_data_remove_D2$sample %>% table()
subtype <- 
  Epithelial[['RNA_snn_res.0.6']] %>% 
  rownames_to_column('cells')
subtype$RNA_snn_res.0.6 <- str_c('Epi_cluster_', subtype$RNA_snn_res.0.6)

meta.data.t <- 
  subtype %>%
  right_join(
    all_data_remove_D2[['annotation']] %>% rownames_to_column('cells'),
    by = "cells"
  ) %>% 
  mutate(
    annotation2 = ifelse((as.character(annotation) %in% c('Epithelial')),
                         as.character(RNA_snn_res.0.6), 
                         as.character(annotation)))%>%
  dplyr::select(cells, annotation2) %>%
  column_to_rownames('cells')

all_data_remove_D2 <- AddMetaData(all_data_remove_D2,meta.data.t['annotation2'])
all_data_remove_D2$annotation %>% table()
all_data_remove_D2$annotation2 %>% table()


## Spatial Dimplot -----------------
all_data_remove_D2_meta <- 
  all_data_remove_D2[[]] %>% 
  rownames_to_column('cell_id')
cells_boundary_ls <- readRDS('R_OUTPUT/data_for_analysis/cells_boundary_ls.rds')


ct_sp_out_dir <- file.path(Epithelial_out_dir, 'Spatial dimplot')
dir.create(ct_sp_out_dir)
ct_sp_p_ls <- 
  cells_boundary_ls %>% 
  imap(function(cells_boundary, sample){ # browser()
    if(sample %in% unique(all_data_remove_D2$sample)){
      anno_data   <- all_data_remove_D2_meta[c("sample", "group2")] %>% distinct()
      group_info  <- anno_data$group2[anno_data$sample == sample]
      
      dat <-
        cells_boundary %>%
        right_join(all_data_remove_D2_meta, by = "cell_id")
      
      p <- dat %>%
        ggplot(aes(x = vertex_x, y = vertex_y)) +
        geom_polygon(aes(fill = annotation2, group = cell_id))  +
        scale_y_reverse() +
        coord_equal() + 
        ggsci::scale_fill_igv() +
        theme_void() +
        labs(fill = NULL, title = paste(sample, group_info, sep = "_") )
    }
    else{
      return(NULL)
    }
    
    p
  }) %>%
  purrr::discard(function(x) length(x) < 1)

ct_sp_p_ls %>% 
  iwalk(function(p, name){
    if(all(is.na(p$data$vertex_x))) return(NULL)
    if(any(is.na(p$data$vertex_x))) warning("vertex_x has NAs ...", .immediate = T)
    pdf(
      file.path(ct_sp_out_dir, str_glue("{name} - spatial dimplot.pdf")),
      width  = 20,
      height = 20
    )
    print(p)
    dev.off()
  })


ct_sp_out_dir2 <- file.path(Epithelial_out_dir, 'Spatial dimplot_split')
dir.create(ct_sp_out_dir2)
ct_sp_p_ls2 <- 
  cells_boundary_ls %>% 
  imap(function(cells_boundary, sample){
    if(sample %in% unique(all_data_remove_D2$sample)){
      anno_data   <- all_data_remove_D2_meta[c("sample", "group2")] %>% distinct()
      group_info  <- anno_data$group2[anno_data$sample == sample]
      
      dat <-
        cells_boundary %>%
        right_join(all_data_remove_D2_meta, by = "cell_id")
      
      ct_color <- structure(
        names = dat$annotation2 %>% unique(),
        ggsci::pal_igv()(20)
      ) %>%
        c(
          Others = "gray90"
        )
      
      p_ls <-
        dat$annotation2 %>% unique() %>% # .[1] %>%
        set_names(., .) %>%
        map(function(ct){
          dat2 <- dat
          
          dat2$annotation2 <- as.character(dat2$annotation2)
          dat2$annotation2[dat2$annotation2 != ct] <- "Others"
          
          p <- dat2 %>%
            ggplot(aes(x = vertex_x, y = vertex_y)) +
            geom_polygon(aes(fill = annotation2, group = cell_id))  +
            scale_y_reverse() +
            coord_equal() +
            # ggsci::scale_fill_igv() +
            scale_fill_manual(values = ct_color) +
            theme_void() +
            labs(fill = NULL, title = paste(sample, group_info, ct, sep = "_") )
          p
        })
    }
    else{
      return(NULL)
    }
    
    p_ls
  }) %>% 
  purrr::discard(function(x) length(x) < 1)

ct_sp_p_ls2 %>% 
  iwalk(function(p, name){ # browser()
    if(all(is.na(p[[1]]$data$vertex_x))) return(NULL)
    if(any(is.na(p[[1]]$data$vertex_x))) warning("vertex_x has NAs ...", .immediate = T)
    pdf(
      file.path(ct_sp_out_dir2, str_glue("{name} - spatial dimplot - split by annotation2.pdf")),
      width  = 20,
      height = 20
    )
    print(p)
    dev.off()
  })

all_data_remove_D2_meta$annotation2 %>% 
  unique() %>% 
  setNames(., .) %>% 
  # .[1] %>% 
  iwalk(function(ct, ct_name){
    pdf(
      file.path(ct_sp_out_dir2, str_glue("{ct_name} - spatial dimplot - split by sample.pdf")),
      width  = 20,
      height = 20
    )
    ct_sp_p_ls2 %>% 
      walk(function(p){ # browser()
        if(all(is.na(p[[1]]$data$vertex_x))) return(NULL)
        if(any(is.na(p[[1]]$data$vertex_x))) warning("vertex_x has NAs ...", .immediate = T)
        print(p[[ct_name]])
      })
    dev.off()
  })

DotPlot(Epithelial_pure, features = 'TP53', group.by = 'subtype')
ggsave(file.path(Epithelial_out_dir, 'dotplot_TP53_subtype.pdf'), width = 6, height = 6)
