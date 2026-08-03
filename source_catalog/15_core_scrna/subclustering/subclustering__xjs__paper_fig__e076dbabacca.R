library(Seurat)
library(tidyverse)

allcolor <- c("#FFFF00","#DC143C","#20B2AA","#FFA500","#9370DB","#98FB98","#8B4513","#F08080","#1E90FF","#7CFC00", "#0000FF",
              "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB",
              "#40E0D0","#5F9EA0","#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4",
              "#32CD32","#F0E68C","#FFFFE0","#EE82EE","#FF6347","#6A5ACD","#9932CC","#8B008B","#DEB887")
out_dir <- 'paper_fig'
dir.create(out_dir)

# all_data<- readRDS('all_data_remove_Tumor3_annotation_full.rds')
# all_data$annotation %>% table()
# all_data_pure <- subset(all_data, annotation != 'Acinar-area' & annotation != 'Fibro-area')
# all_data_pure$annotation %>% table()
# saveRDS(all_data_pure, 'all_data_remove_Tumor3_annotation_full_pure.rds')


# UMAP annotation ---------------
all_data <- readRDS('all_data_remove_Tumor3_annotation_full_pure.rds')
all_data$annotation %>% table()

DimPlot(all_data, group.by = 'annotation', cols = allcolor, alpha = 0.4)
ggsave(file.path(out_dir, 'Dimplot_annotation.pdf'), width = 9, height = 6)

DimPlot(all_data, group.by = 'annotation', cols = allcolor, label = T, alpha = 0.4) + NoLegend()
ggsave(file.path(out_dir, 'Dimplot_annotation - no legnd.pdf'), width = 7, height = 6)


DimPlot(all_data, group.by = 'annotation', cols = allcolor, split.by = 'group', alpha = 0.4)
ggsave(file.path(out_dir, 'Dimplot_annotation - split.pdf'), width = 20, height = 6)

DimPlot(all_data, group.by = 'annotation', cols = allcolor, label = T, split.by = 'group', alpha = 0.4) + NoLegend()
ggsave(file.path(out_dir, 'Dimplot_annotation - no legnd - split.pdf'), width = 18, height = 6)




# dotplot annotation -------
theme_90angle <-
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
    strip.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )
cellTye_markers <- c(  
  'B cell' = "MS4A1", 
  'B cell' = "CD19", 
  
  'Epithelials' = "EPCAM", 
  'Ductal cell' = "CFTR", 
  
  'Endothelial cell' = "PECAM1", 
  
  'Tumor' = 'HMGA1', 
  'Tumor' = 'CEACAM1', 
  'Exocrine glandular cell' = "ERP27", 

  'Fibroblast' = "COL5A1", 
  'Fibroblast' = "CCDC80", 
  'Fibroblast' = "COL4A1", 
  
  'Mast cell' = "KIT", 
  'Mast cell' = "MS4A2", 
  'Mast cell' = "HDC", 
  
  'Myeloid cell' = "CD14", 
  'Myeloid cell' = "CD68", 
  'Myeloid cell' = "IDO1", 
  
  'Pancreatic β-cell' = "SLC30A8", 
  
  'Pericyte' = "MCAM", 
  'Pericyte' = "RGS5", 
  
  'Plasma cell' = "MZB1", 
  'Plasma cell' = "XBP1", 
  'Plasma cell' = "TNFRSF17", 
  
  'Schwann cell' = "MPZ", 
  'Schwann cell' = "SOX10", 
  
  'T cell' = "CD3E"
)

DotPlot(all_data, features = cellTye_markers, group.by = 'annotation') + theme_90angle
ggsave(file.path(out_dir, 'dotplot_annotation.pdf'), width = 13, height = 6)





# Myeloid pure ----------
# Myeloid <- readRDS('Myeloid_subtype/Myeloid.rds')
# Myeloid$subtype %>% table()
# Myeloid_pure <- subset(Myeloid, subtype %in% c('Neutropil', 'M2 macrophage', 'cDC2', 'M1 macrophage', 'Proliferating macrophage', 'IFN-responsive macrophage', 'LAMP3 DC'))
# Myeloid_pure$subtype %>% table()
# saveRDS(Myeloid_pure, 'Myeloid_subtype/Myeloid_pure.rds')

Myeloid <- readRDS('Myeloid_subtype/Myeloid_pure.rds')
Myeloid$subtype %>% table()

DimPlot(Myeloid, group.by = 'subtype', cols = allcolor, reduction = 'harmonyumap')

Myeloid <- FindVariableFeatures(Myeloid)
Myeloid <- ScaleData(Myeloid, vars.to.regress = 'nCount_RNA')
Myeloid <- RunPCA(Myeloid)
pdf(file.path(out_dir, 'Myeloid pca plot.pdf'), width = 8, height = 6)
ElbowPlot(Myeloid, ndims = 50)
dev.off()

dims <- 20

Myeloid <- RunUMAP(Myeloid, dims = 1:dims, return.model = T)
DimPlot(Myeloid, group.by = 'subtype', cols = allcolor)

Myeloid <- Myeloid %>%
  harmony::RunHarmony(
    dims.use = 1:dims, 
    group.by.vars = "sample", 
    reduction.save = "harmony"
  )

Myeloid <- RunUMAP(
  Myeloid, 
  dims = 1:dims, 
  reduction = 'harmony', 
  reduction.name = 'harmonyumap', 
  reduction.key = 'harmonyUMAP ', 
  return.model = T
  )

DimPlot(Myeloid, group.by = 'subtype', cols = allcolor, reduction = 'harmonyumap', alpha = 0.4)
ggsave(file.path(out_dir, 'Myeloid_Dimplot_subtype.pdf'), width = 9, height = 6)

DimPlot(Myeloid, group.by = 'subtype', cols = allcolor, label = T, reduction = 'harmonyumap', alpha = 0.4) + NoLegend()
ggsave(file.path(out_dir, 'Myeloid_Dimplot_subtype - no legnd.pdf'), width = 7, height = 6)


DimPlot(Myeloid, group.by = 'subtype', cols = allcolor, split.by = 'group', reduction = 'harmonyumap', alpha = 0.4)
ggsave(file.path(out_dir, 'Myeloid_Dimplot_subtype - split.pdf'), width = 20, height = 6)

DimPlot(Myeloid, group.by = 'subtype', cols = allcolor, label = T, split.by = 'group', reduction = 'harmonyumap', alpha = 0.4) + NoLegend()
ggsave(file.path(out_dir, 'Myeloid_Dimplot_subtype - no legnd - split.pdf'), width = 18, height = 6)


## dotplot --------
Myeloid_cellTye_markers <- c(  
  'M2 macrophage' = "CD163", 
  'M2 macrophage' = "MRC1", 
  'M2 macrophage' = "CTSL", 
  
  'Neutropil' = "CSF3R", 
  
  'DC' = "CD83", 
  
  'cDC2' = "FCER1A", 
  'cDC2' = "CD1C", 
  
  'M1 macrophage' = "CD300E", 
  'M1 macrophage' = "CCR2", 
  'M1 macrophage' = "CLEC5A", 
  
  'Proliferating macrophage' = "TOP2A", 
  'Proliferating macrophage' = "MKI67", 
  
  'IFN-responsive macrophage' = "IFIT3", 
  'IFN-responsive macrophage' = "MX1", 
  'IFN-responsive macrophage' = "IFIT1", 
  'IFN-responsive macrophage' = "IFIT2", 
  
  'LAMP3 DC' = "IDO1", 
  'LAMP3 DC' = "LAMP3"
)

DotPlot(Myeloid, features = Myeloid_cellTye_markers, group.by = 'subtype') + theme_90angle
ggsave(file.path(out_dir, 'Myeloid_dotplot_subtype.pdf'), width = 13, height = 6)

saveRDS(Myeloid, "~/Myeloid_subtype/Myeloid_pure.rds")




# Fibroblast -----------
Fibroblast <- readRDS('Fibroblast_subtype/Fibroblast.rds')

## volcano-------
Idents(Fibroblast) <- 'group'

groups <- Fibroblast$group %>% unique()
all_combns_ls <- groups %>% seq_along() %>% combn(2) %>% apply(2, c, simplify = FALSE)
all_combns_ls <- all_combns_ls[c(1, 3)]

Fibroblast_markers_ls <- 
  all_combns_ls %>% 
  map(function(x){ # browser()
    grp <- groups[x]
    # grp <- all_grp_ls[x]
    FindMarkers(
      # JoinLayers(Fibroblast),
      Fibroblast, 
      group.by = 'group',
      test.use = "wilcox", 
      ident.1 = grp[1],
      ident.2 = grp[2],
      # cells.1 = grp[[1]],
      # cells.2 = grp[[2]],
      logfc.threshold = 0.25,
      min.pct = 0.1,
      only.pos = TRUE
    ) %>% 
      rownames_to_column('gene')
  })
names(Fibroblast_markers_ls) <- 
  all_combns_ls %>% 
  map(function(x){ # browser()
    grp <- groups[x]
    str_glue('{grp[1]} vs {grp[2]}')
  }) %>% 
  unlist()

openxlsx::write.xlsx(Fibroblast_markers_ls, file.path(out_dir, 'Fibroblast_diff_gene.xlsx'))

library(EnhancedVolcano)
# type <- unique(Fibroblast_markers$cluster)

undebug(EnhancedVolcano)
volcano_ls <- 
  Fibroblast_markers_ls %>%
  imap(function(Fibroblast_markers, name){
    p <- EnhancedVolcano(Fibroblast_markers,
                         lab = Fibroblast_markers[['gene']], # Fibroblast_markers['gene'],
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
      filename = file.path(out_dir, str_glue("{name} fibroblast diffgene valcano.pdf")),
      height = 6,
      width = 8
    )
  })




## enrichment ----------------------------
Fibroblast_de_genes_for_enrichmenet_ls <- 
  Fibroblast_markers_ls %>% 
  map(function(Fibroblast_markers){ # browser()
    Fibroblast_markers %>% 
      dplyr::filter(p_val_adj < 0.05) %>%
      dplyr::filter(pct.1-pct.2 > 0.1) %>%
      dplyr::filter(abs(avg_log2FC) > 0.5) %>% 
      .[["gene"]]
  })
  

enrich_ls <- 
  Fibroblast_de_genes_for_enrichmenet_ls %>%
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


er_output <- file.path(out_dir, 'Fibroblast_enrichment') 
dir.create(er_output)

# expot enrichment table
enrich_ls %>%
  # discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    dir.create(file.path(er_output, name))
    
    write.csv(
      x@result,
      file = file.path(er_output, name, "Go enrichment table.csv"),
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
      filename = file.path(er_output, name, "go dotplot.pdf"),
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
      filename = file.path(er_output, name, "go barplot.pdf"),
      height = 6,
      width = 8
    )
  })

Fibroblast_de_genes_id_for_enrichmenet_ls <-
  Fibroblast_de_genes_for_enrichmenet_ls %>%
  map(function(x){ # browser()
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
  Fibroblast_de_genes_id_for_enrichmenet_ls %>%
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

# saveRDS(enrich_ls, file.path(er_output, 'enrich_ls.rds'))

# expot enrichment table
enrich_ls %>%
  # discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    dir.create(file.path(er_output, name))
    
    write.csv(
      x@result,
      file = file.path(er_output, name, "KEGG enriment table.csv"),
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
      filename = file.path(er_output, name, "kegg dotplot.pdf"),
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
      filename = file.path(er_output, name, "kegg barplot.pdf"),
      height = 6,
      width = 8
    )
  })

# saveRDS(Fibroblast, file.path(out_dir, 'Fibroblast.rds'))




# T cell ---------- 
All_T_anno <- readRDS('T_subtype/All_T_anno.rds')
All_T_anno$subtype %>% table()
DimPlot(All_T_anno, group.by = 'subtype', reduction = 'harmonyumap')


## dotplot --------
All_T_anno_cellTye_markers <- c(  
  'T' = 'CD3E', 
  
  'CD8' = 'CD8A',
  
  'CD4' = 'CD4',
  
  'CD8 Tem' = 'GZMK', 
  'CD8 Tem' = 'GZMA',
  'CD8 Tem' = 'KLRG1',
  
  'CD4 Trm'='CXCR6',
  'CD4 Trm'='ITGAE',
  
  "CD8 Tnaive" = "SELL",
  "CD8 Tnaive" = "LEF1",
  "CD8 Tnaive" = "TCF7",
  "CD8 Tnaive" = "CCR7",
  
  "gdT" = "TRDC", # https://apps.allenimmunology.org/aifi/resources/imm-health-atlas/cell-type-descriptions/cd8-t-cells-gdt-cells-and-mait-cells/
  
  'CD4 Treg'='FOXP3',
  'CD4 Treg'='CTLA4',
  'CD4 Treg'='TIGIT',
  
  'Cycling T'='MKI67',
  'Cycling T'='TOP2A'
)

DimPlot(All_T_anno, group.by = 'subtype', cols = allcolor, reduction = 'harmonyumap', alpha = 0.4)
ggsave(file.path(out_dir, 'All_T_anno_Dimplot_subtype.pdf'), width = 8, height = 6)

DimPlot(All_T_anno, group.by = 'subtype', cols = allcolor, label = T, reduction = 'harmonyumap', alpha = 0.4) + NoLegend()
ggsave(file.path(out_dir, 'All_T_anno_Dimplot_subtype - no legnd.pdf'), width = 7, height = 6)

All_T_anno$group %>% table()
All_T_anno$group <- factor(All_T_anno$group, levels = c('N', 'CP', 'T'))
DimPlot(All_T_anno, group.by = 'subtype', cols = allcolor, split.by = 'group', reduction = 'harmonyumap', alpha = 0.4)
ggsave(file.path(out_dir, 'All_T_anno_Dimplot_subtype - split.pdf'), width = 18, height = 6)

DimPlot(All_T_anno, group.by = 'subtype', cols = allcolor, label = T, split.by = 'group', reduction = 'harmonyumap', alpha = 0.4) + NoLegend()
ggsave(file.path(out_dir, 'All_T_anno_Dimplot_subtype - no legnd - split.pdf'), width = 16, height = 6)

DotPlot(All_T_anno, features = All_T_anno_cellTye_markers, group.by = 'subtype') + theme_90angle
ggsave(file.path(out_dir, 'All_T_anno_dotplot_subtype.pdf'), width = 10, height = 4)




# B subtype ----------- 
b_out_dir <- file.path(out_dir, "B_subtype") 
dir.create(b_out_dir)
all_data$annotation %>% table()
B <- subset(JoinLayers(all_data), annotation %in% c('B cell', 'Plasma cell'))

## PCA ----------------------------------
DefaultAssay(B) <- "RNA"
B <- FindVariableFeatures(B)
B <- ScaleData(B, vars.to.regress = 'nCount_RNA')
B <- RunPCA(B)

### pca inspector ------
pdf(file.path(b_out_dir, 'pca plot.pdf'), height = 6, width = 8)
ElbowPlot(B, ndims = 50)
dev.off()


## dimension reduction --------
dims <- 10 

library(future)
plan(multisession, workers = 10)
options(future.globals.maxSize = 10 * 1024 ^ 3) # 10Gb

B <- RunUMAP(B, dims = 1:dims,return.model = TRUE)

pdf(file.path(b_out_dir, str_glue("Dimplot dim{dims} sample.pdf")), height = 6, width = 8)
DimPlot(B, group.by = "sample", label = TRUE) + NoLegend()
dev.off()

B <- B %>% 
  harmony::RunHarmony(
    dims.use = 1:dims,
    group.by.vars = "sample",
    reduction.save = "harmony"
    # theta = 3, # default 2, bigger value -- more diversity per cluster
    # lambda = 0.5, # default, smaller value -- more integrated
  )

B <- RunUMAP(
  B, 
  dims = 1:dims, 
  reduction = "harmony", 
  reduction.name = "harmonyumap", 
  reduction.key  = "harmonyUMAP_",
  return.model = TRUE
)

pdf(file.path(b_out_dir, str_glue("Dimplot dim{dims} sample - harmony.pdf")), height = 6, width = 8)
DimPlot(B, group.by = "sample", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()


reduction.use <- 'harmony' # 'harmony' 'pca'
### cluster -------
B <-
  FindNeighbors(
    B ,
    reduction = reduction.use,
    dims = 1:dims,
    k.param = 20
  )

B <-
  FindClusters(
    B,
    resolution = c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.5)
  )

# dimplot
pdf(file.path(b_out_dir, str_glue("Dimplot dim{dims} cluster.pdf")), height = 6, width = 8)
DimPlot(B, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(B, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(B, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(B, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(B, group.by = "RNA_snn_res.1", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(B, group.by = "RNA_snn_res.1.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(B, group.by = "RNA_snn_res.1.5", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()

# Idents(B)<-'RNA_snn_res.0.8'
# tmp <- FindSubCluster(B, cluster = 0, graph.name = "RNA_snn", algorithm = 1,resolution = 0.3)
# Idents(tmp) <- "sub.cluster"
# DimPlot(tmp, group.by = "sub.cluster", label = TRUE, reduction = 'harmonyumap')
# DotPlot(tmp, group.by = "sub.cluster", features = B_major_marker_ls) + RotatedAxis()
# B <- tmp

### dotplot -------
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

#https://zhuanlan.zhihu.com/p/25031149689
B_major_marker_ls <- c(
  "B cells" = "CD19",
  "B cells" = "CD79A",
  "B cells" = "MS4A1",
  'naive_B' = 'CD20',
  'naive_B' = 'CD24',
  'naive_B' = 'IGHM',
  'naive_B' = 'BANK1',
  'naive_B' = 'FCRL1',
  'naive_B' = 'HLA-DQB1',
  'naive_B' = 'HLA-DQA1',
  'naive_B' = 'HVCN1',
  'naive_B' = 'IGHD',
  'naive_B' = 'FCER2',
  'naive_B' = 'IL4R',
  'Follicular B cells' = "CD40", 
  'Follicular B cells' = "BCL6", 
  'memory_B'= 'CD27',
  'memory_B'= 'FCRL4',
  'memory_B'= 'TACI',
  'memory_B'= 'AIM2',
  'memory_B'= 'TXNIP',
  'memory_B'= 'GPR183',
  'memory_B'= 'TNFRSF13B',
  'memory_B'= 'CD22',
  'memory_B'= 'FCRL2',
  'memory_B'= 'ADAM28',
  'memory_B'= 'IGLL5',
  'memory_B'= 'P2RX5',
  'memory_B'= 'BLK',
  # 'Marginal zone B cells' = "CD27", 
  'Marginal zone B cells' = "CR2",
  'Marginal zone B cells' = "CD1C",
  'cyclingB'='MKI67',
  'cyclingB'='STMN1',
  'cyclingB'='HMGB2',
  'cyclingB'='BIRC5',
  "Plasma" = "CD138",
  "Plasma" = "PRDM1",
  "Plasma" = "BLIMP1",
  "Plasma" = "IGHA1",
  "Plasma" = "IGHA2",
  "Plasma" = "IGHG1",
  "Plasma" = "IGHG2",
  "Plasma" = "IGHG3",
  "Plasma" = "IGHG4",
  "Plasma" = "IGKC",
  "Plasma" = "TNFRSF17",
  "Plasma" = "JCHAIN",
  "Plasma" = "MZB1",
  "Plasma" = "XBP1",
  "Plasma" = "SDC1",
  #'Plasmablasts' = "MKI67", 
  'Plasmablasts' = "CD38", 
  'Plasmablasts' = "IRF4", 
  'Bregs' = "IL10", 
  'Bregs' = "CD1D", 
  'Bregs' = "TGFB1", 
  'GC B'= 'S1PI2',
  # 'GC B'='BCL6',
  'GC B'='RGS13',
  'GC B'='TCL1A',
  'GC B'= 'LRMP',
  'GC B'= 'SUGCT',
  'GC B'= 'MME',
  'GC B'='AICDA'
)

Idents(B)<-'RNA_snn_res.0.2'
B_markers <-
  FindAllMarkers(
    JoinLayers(B),
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    only.pos = TRUE
  )

B_top5_markers <- 
  B_markers %>% 
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 30) %>%
  arrange(desc(pct.1)) %>%
  slice_head(n = 5) %>%
  nest() %>%
  deframe() %>% 
  map(pull, 'gene') %>%
  imap(function(genes, name){
    re <- genes %>%
      map(function(gene){
        # browser()
        names(gene) = str_glue("cluster_{name}")
        gene
      })
    re
  }) %>%
  unlist() %>% 
  .[!duplicated(.)]

names(B_top5_markers) <- sub('.*\\.','',names(B_top5_markers))


# dotplot
pdf(file.path(b_out_dir, str_glue("dotplot dims{dims} major_marker.pdf")), width = 25, height = 6)
DotPlot(B, group.by = "RNA_snn_res.0.2", features = major_marker_ls) + theme_90angle
DotPlot(B, group.by = "RNA_snn_res.0.2", features = B_major_marker_ls) + RotatedAxis()
DotPlot(B, group.by = "RNA_snn_res.0.2", features = B_top5_markers) + theme_90angle
# DotPlot(B, group.by = "RNA_snn_res.0.2", features = Fibro_major_marker_ls2) + RotatedAxis()
# DotPlot(B, group.by = "RNA_snn_res.0.2", features = CD4_major_marker_ls) + RotatedAxis()
# DotPlot(B, group.by = "RNA_snn_res.0.2", features = CD8_major_marker_ls) + RotatedAxis()
dev.off()


## annotation -------------------------------------------------------------------------
# cell_cluster <- c(
#   '0_0'  = "Memory B cells", 
#   '0_1'  = "Naive B cells",
#   '1'  = "Plasma cells", # MZB1, XBP1
#   '2'  = "B cells", # no other markers
#   '3'  = "Plasma cells", # MZB1, XBP1
#   '4'  = "Plasma cells", # MZB1, XBP1
#   '5'  = "Epithelials", # Epcam
#   '6'  = "Marginal zone B cells", # CR2 https://www.science.org/doi/10.1126/science.abf7470
#   '7' = "rCAFs", # CCDC80, CCL19
#   '8' = "Plasma cells", 
#   '9' = "Plasma cells",  
#   '10' = "Macrophage-like B cells", # CD68, CD19
#   '11' = "Cycling B cells"
# )

#Idents(B) <- "RNA_snn_res.0.8"

cell_cluster <- c(
  '0' = "Naive B cell", #
  '1' = "mix", # Fibro
  '2' = "Plasma cell",
  '3' = "Plasmablast",
  '4' = "mix", # Myeloid
  '5' = "mix" # Epi
)

Idents(B) <- "RNA_snn_res.0.2"

B <- RenameIdents(B, cell_cluster)
B$subtype <- Idents(B)


pdf(file.path(b_out_dir, str_glue("Dimplot dim{dims} annotation.pdf")), height = 6, width = 8)
DimPlot(B, group.by = "subtype", label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

pdf(file.path(b_out_dir, str_glue("Dimplot dim{dims} annotation split.pdf")), height = 5, width = 16)
DimPlot(B, group.by = "subtype",split.by = 'group', label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

saveRDS(B, file.path(b_out_dir, 'B.rds'))


## B - pure --------------
B_pure <- subset(JoinLayers(B), subtype %in% c('Plasma cell', 'Plasmablast', 'Naive B cell'))

DimPlot(B_pure, group.by = 'subtype', cols = allcolor, label = T, alpha = 0.4) + NoLegend()
ggsave(file.path(b_out_dir, 'Bcell_umap_subtype_noLegend_orig.pdf'), width = 6, height = 5)
DimPlot(B_pure, group.by = 'subtype', cols = allcolor, label = F, alpha = 0.4)
ggsave(file.path(b_out_dir, 'Bcell_umap_subtype_orig.pdf'), width = 7.5, height = 5)

### PCA ----------------------------------
DefaultAssay(B_pure) <- "RNA"
B_pure <- FindVariableFeatures(B_pure)
B_pure <- ScaleData(B_pure, vars.to.regress = 'nCount_RNA')
B_pure <- RunPCA(B_pure)
pdf(file.path(b_out_dir, 'Bcell_pca plot.pdf'), height = 6, width = 8)
ElbowPlot(B_pure, ndims = 50)
dev.off()

dims <- 10

library(future)
plan(multisession, workers = 10)
options(future.globals.maxSize = 10 * 1024 ^ 3) # 10Gb

B_pure <- RunUMAP(B_pure, dims = 1:dims,return.model = TRUE)
DimPlot(B_pure, group.by = 'subtype', cols = allcolor, label = T, alpha = 0.4) + NoLegend()

B_pure <- B_pure %>% 
  harmony::RunHarmony(
    dims.use = 1:dims,
    group.by.vars = "sample",
    reduction.save = "harmony"
    # theta = 3, # default 2, bigger value -- more diversity per cluster
    # lambda = 0.5, # default, smaller value -- more integrated
  )

B_pure <- RunUMAP(
  B_pure, 
  dims = 1:dims, 
  reduction = "harmony", 
  reduction.name = "harmonyumap", 
  reduction.key  = "harmonyUMAP_",
  return.model = TRUE
)

DimPlot(B_pure, group.by = 'subtype', cols = allcolor, label = T, alpha = 0.4) + NoLegend()
ggsave(file.path(b_out_dir, 'Bcell_umap_subtype_noLegend.pdf'), width = 6, height = 5)
DimPlot(B_pure, group.by = 'subtype', cols = allcolor, label = F, alpha = 0.4)
ggsave(file.path(b_out_dir, 'Bcell_umap_subtype.pdf'), width = 7.5, height = 5)




B <- B_pure
b_out_dir <- file.path(b_out_dir, 'B_pure')
dir.create(b_out_dir)
### cluster -------
B <- 
  FindNeighbors(
    B ,
    reduction = reduction.use,
    dims = 1:dims,
    k.param = 10
  )

B <-
  FindClusters(
    B,
    resolution = c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.5)
  )

# dimplot
pdf(file.path(b_out_dir, str_glue("Dimplot dim{dims} cluster.pdf")), height = 6, width = 8)
DimPlot(B, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(B, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(B, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(B, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(B, group.by = "RNA_snn_res.1", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(B, group.by = "RNA_snn_res.1.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(B, group.by = "RNA_snn_res.1.5", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()


Idents(B)<-'RNA_snn_res.0.2'
B_markers <-
  FindAllMarkers(
    JoinLayers(B),
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    only.pos = TRUE
  )

B_top5_markers <- 
  B_markers %>% 
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 30) %>%
  arrange(desc(pct.1)) %>%
  slice_head(n = 5) %>%
  nest() %>%
  deframe() %>% 
  map(pull, 'gene') %>%
  imap(function(genes, name){
    re <- genes %>%
      map(function(gene){
        # browser()
        names(gene) = str_glue("cluster_{name}")
        gene
      })
    re
  }) %>%
  unlist() %>% 
  .[!duplicated(.)]

names(B_top5_markers) <- sub('.*\\.','',names(B_top5_markers))


# dotplot
pdf(file.path(b_out_dir, str_glue("dotplot dims{dims} major_marker.pdf")), width = 25, height = 6)
DotPlot(B, group.by = "RNA_snn_res.0.2", features = major_marker_ls) + theme_90angle
DotPlot(B, group.by = "RNA_snn_res.0.2", features = B_major_marker_ls) + RotatedAxis()
DotPlot(B, group.by = "RNA_snn_res.0.2", features = B_top5_markers) + theme_90angle
# DotPlot(B, group.by = "RNA_snn_res.0.2", features = Fibro_major_marker_ls2) + RotatedAxis()
# DotPlot(B, group.by = "RNA_snn_res.0.2", features = CD4_major_marker_ls) + RotatedAxis()
# DotPlot(B, group.by = "RNA_snn_res.0.2", features = CD8_major_marker_ls) + RotatedAxis()
dev.off()


## annotation -------------------------------------------------------------------------
# cell_cluster <- c(
#   '0_0'  = "Memory B cells", 
#   '0_1'  = "Naive B cells",
#   '1'  = "Plasma cells", # MZB1, XBP1
#   '2'  = "B cells", # no other markers
#   '3'  = "Plasma cells", # MZB1, XBP1
#   '4'  = "Plasma cells", # MZB1, XBP1
#   '5'  = "Epithelials", # Epcam
#   '6'  = "Marginal zone B cells", # CR2 https://www.science.org/doi/10.1126/science.abf7470
#   '7' = "rCAFs", # CCDC80, CCL19
#   '8' = "Plasma cells", 
#   '9' = "Plasma cells",  
#   '10' = "Macrophage-like B cells", # CD68, CD19
#   '11' = "Cycling B cells"
# )

#Idents(B) <- "RNA_snn_res.0.8"

cell_cluster <- c(
  '0' = "Naive B cell", #
  '1' = "Plasmablast", 
  '2' = "Plasma cell",
  '3' = "Plasmablast",
  '4' = "Naive B cell"
)

Idents(B) <- "RNA_snn_res.0.2"

B <- RenameIdents(B, cell_cluster)
B$subtype <- Idents(B)


pdf(file.path(b_out_dir, str_glue("Dimplot dim{dims} annotation.pdf")), height = 6, width = 8)
DimPlot(B, group.by = "subtype", label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

pdf(file.path(b_out_dir, str_glue("Dimplot dim{dims} annotation split.pdf")), height = 5, width = 16)
DimPlot(B, group.by = "subtype",split.by = 'group', label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

saveRDS(B, file.path(b_out_dir, 'B_pure.rds'))


cellTye_markers <- c(
  "B cells" = "CD19",
  "B cells" = "CD79A",
  "B cells" = "MS4A1",
  
  'naive_B' = 'CD20',
  'naive_B' = 'CD24',
  'naive_B' = 'IGHM',
  'naive_B' = 'BANK1',
  'naive_B' = 'FCRL1',
  'naive_B' = 'HLA-DQB1',
  'naive_B' = 'HLA-DQA1',
  'naive_B' = 'HVCN1',
  'naive_B' = 'IGHD',
  'naive_B' = 'FCER2',
  'naive_B' = 'IL4R', 
  
  "Plasma" = "MZB1",
  "Plasma" = "XBP1", 
  # "Plasma" = "PRDM1", 
  # "Plasma" = "TNFRSF17",
  
  'Plasmablasts' = "CD38", 
  'Plasmablasts' = "IRF4"
)

pdf(file.path(b_out_dir, "dotplot celltype main_marker.pdf"), width = 10, height = 4)
DotPlot(B, group.by = "subtype", features = cellTye_markers, dot.scale = 10) + theme_90angle
dev.off()


## cell props barplot  --------
B_meta <- B[[]]

B$sample <- factor(B$sample, levels = c(B_meta$sample %>% unique()))

p_col_proj <-
  B_meta %>%
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
  filename = file.path(b_out_dir, "cell_freq_sample.pdf"),
  height = 6,
  width  = 6
)

p_col_proj <-
  B_meta %>%
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
  filename = file.path(b_out_dir, "cell_freq_group.pdf"),
  height = 6,
  width  = 4
)
