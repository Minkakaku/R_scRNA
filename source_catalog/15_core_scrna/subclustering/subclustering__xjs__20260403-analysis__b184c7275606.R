# 20260403 Normal Epi subcluster ------------

result_dir <- "20260403-NormalEpi_subcluster" 
dir.create(result_dir)

## read RDS----
library(Seurat)
library(tidyverse)
library(BPCells)
library(clusterProfiler)
library(org.Hs.eg.db)

all_data <- readRDS(file.path("all_data_RCTD_subtype_1023.rds"))

# markers: https://cloud.tencent.com/developer/article/2115987
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

valid_genes <- major_marker_ls[major_marker_ls %in% rownames(spatial)]
cat("有效基因:", length(valid_genes), "/", length(major_marker_ls))

# 提取差集：virus_genes中存在但valid_Vriusgenes中缺失的基因
missing_gene <- setdiff(major_marker_ls, valid_genes)
missing_gene

theme_90angle <-
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
    strip.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )


## Epithelials -------------------------------------------------------------------
Normal <- all_data[, all_data@meta.data$RCTD == "Normal"]

### PCA --------
DefaultAssay(Normal) <- "RNA"
Normal <- FindVariableFeatures(Normal)
Normal <- ScaleData(Normal, vars.to.regress = 'nCount_RNA')
Normal <- RunPCA(Normal)


### pca inspector ------
pdf(file.path(result_dir, 'pca plot.pdf'), height = 6, width = 8)
ElbowPlot(Normal, ndims = 50)
dev.off()

pdf(file.path(result_dir, 'dimheatmap 1-5.pdf'), height = 6, width = 10)
DimHeatmap(Normal, dims = 1:5)
dev.off()
pdf(file.path(result_dir, 'dimheatmap 6-10.pdf'), height = 6, width = 10)
DimHeatmap(Normal, dims = 6:10)
dev.off()
pdf(file.path(result_dir, 'dimheatmap 11-15.pdf'), height = 6, width = 10)
DimHeatmap(Normal, dims = 11:15)
dev.off()
pdf(file.path(result_dir, 'dimheatmap 16-20.pdf'), height = 6, width = 10)
DimHeatmap(Normal, dims = 16:20)
dev.off()


### dimension reduction --------
dims <- 17

library(future)
plan(multisession, workers = 10)
options(future.globals.maxSize = 10 * 1024 ^ 3) # 10Gb


Normal <- RunUMAP(Normal, dims = 1:dims,return.model = TRUE)



pdf(file.path(result_dir, str_glue("Dimplot dim{dims} sample.pdf")), height = 6, width = 8)
DimPlot(Normal, group.by = "sample", label = TRUE) + NoLegend()
dev.off()

pdf(file.path(result_dir, str_glue("Dimplot dim{dims} sample - split.pdf")), height = 5, width = 12)
DimPlot(Normal, group.by = "sample", reduction = "umap", split.by = "sample")
dev.off()



Normal <- Normal %>% 
  harmony::RunHarmony(
    dims.use = 1:dims,
    group.by.vars = "sample",
    reduction.save = "harmony"
    # theta = 3, # default 2, bigger value -- more diversity per cluster
    # lambda = 0.5, # default, smaller value -- more integrated
  )

Normal <- RunUMAP(
  Normal, 
  dims = 1:dims, 
  reduction = "harmony", 
  reduction.name = "harmonyumap", 
  reduction.key  = "harmonyUMAP_",
  return.model = TRUE
)


pdf(file.path(result_dir, str_glue("Dimplot dim{dims} sample - harmony.pdf")), height = 6, width = 8)
DimPlot(Normal, group.by = "sample", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()

pdf(file.path(result_dir, str_glue("Dimplot dim{dims} sample - harmony - split.pdf")), height = 5, width = 12)
DimPlot(Normal, group.by = "sample", reduction = "harmonyumap", split.by = "sample")
dev.off()



reduction.use <- 'harmony' # 'harmony' 'pca'


### cluster -------
Normal  <-
  FindNeighbors(
    Normal ,
    reduction = reduction.use,
    dims = 1:dims,
    k.param = 20
  )


Normal <-
  FindClusters(
    Normal,
    resolution = c(0.2, 0.4, 0.6, 0.8)
  )



# dimplot
pdf(file.path(result_dir, str_glue("Dimplot dim{dims} cluster.pdf")), height = 6, width = 8)
DimPlot(Normal, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Normal, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Normal, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Normal, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()




### spatial plot------------
spatial_outdir <- file.path(result_dir, "spatial plot")
dir.create(spatial_outdir)

cells_boundary_ls <- readRDS("R_OUTPUT/data_for_analysis/cells_boundary_ls.rds")

Normal_meta <-
  Normal[[]] %>%
  rownames_to_column("cell_id")

all_data_meta <- 
  all_data[[]] %>% 
  rownames_to_column("cell_id")

sp_p_ls <- 
  cells_boundary_ls %>% 
  imap(function(cells_boundary, project){ # browser()
    # image_wh_r <- image_wh_r_ls[[project]]
    
    dat <-
      cells_boundary %>%
      right_join(all_data_meta, by = "cell_id")
    
    dat2 <-
      cells_boundary %>%
      right_join(Normal_meta, by = "cell_id")
    
    p <- dat[which(dat$RCTD != 'Normal'),] %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) +
      geom_polygon(aes(group = cell_id), fill = "gray90")  +
      scale_y_reverse() +
      coord_equal() + 
      # coord_fixed(ratio = 1/image_wh_r) +
      theme_void() +
      ggnewscale::new_scale_fill() +
      geom_polygon(data = dat2, aes(x = vertex_x, y = vertex_y, fill = RNA_snn_res.0.6, group = cell_id)) +
      ggsci::scale_fill_igv() +
      theme_void()
    p
  })

sp_p_ls %>% 
  iwalk(function(p, grp){
    pdf(file.path(spatial_outdir, str_glue("Normal RNA_snn_res.0.6 spatial plot {grp}.pdf")), height = 18, width = 18)
    print(p)
    dev.off()
  })




### dotplot -------
# 1. Epi2: Xiao S, Zhou L. Gastric Stem Cells: Physiological and Pathological Perspectives. Front Cell Dev Biol. 2020
# 2. Epi3: Kim TH et al., Stomach development, stem cells and disease. Development. 2016
# 3. Yaser Greish et al., Biodegradable Scaffolds for Gastric Tissue Regeneration. 2017
# 4. Epi4: Busslinger GA et al., Human gastrointestinal epithelial of the esophagus, stomach, and duodenum resolved at single-cell resolution. Cell Rep. 2021
# 5. Epi5: Engevik AC et al., The Physiology of the Gastric Parietal Cell. Physiol Rev. 2020
# 6. Epi6: Zhang P et al.,Dissecting the Single-Cell Transcriptome Network Underlying Gastric Premalignant Lesions and Early Gastric Cancer. Cell Rep. 2019
# 7. Epi7: Sathe A et al., Single-Cell Genomic Characterization Reveals the Cellular Reprogramming of the Gastric Tumor Microenvironment. Clin Cancer Res. 2020
# 8. Epi8: Zhang M et al., Dissecting transcriptional heterogeneity in primary gastric adenocarcinoma by single cell RNA sequencing. Gut. 2021
Epithelial_major_marker_ls <- c(
  'Surface mucous cells' = "MUC5AC", # 1
  'Surface mucous cells' = "TFF1", # 1 
  'Surface mucous cells' = "PAS", # 1 
  'Surface mucous cells' = "GKN1", # 9 
  'Surface mucous cells' = "GKN2", # 9 
  
  'Mucous gland cells' = "MUC6", # 1 # 5:TFF2 
  
  'Mucous neck cells' = "TFF2", # 1 
  
  'Parietal cells' = "GIF", # 1 
  'Parietal cells' = "ESRRG", # 2
  'Parietal cells' = "ATP4B", # 2
  
  'Enteroendocrine cells' = "CHGA", # 3,6
  'Enteroendocrine cells' = "CHGB", # 6
  'Enteroendocrine cells' = "GAST", # 9
  'Enteroendocrine cells' = "PROX1", # 9
  
  'Acid-secreting parietal cells' = "ATP4A", # 4 
  #'Acid-secreting parietal cells' = "ATP4B", # 4 
  
  # 'Proliferating isthmus cells' = "MKI67", # 4 
  'Proliferating isthmus cells' = "STMN1", # 4
  
  'Proliferative cells' = "MKI67", # 6
  'Proliferative cells' = "BRIC5", # 6
  'Proliferative cells' = "CDKN2A", # 9
  'Proliferative cells' = "RBP4", # 9
  
  'Chief cells' = "LIPF", # 4
  'Chief cells' = "PGA5", # 4
  'Chief cells' = "MIST1", # 5
  'Chief cells' = "PGC", # 5
  'Chief cells' = "PGA3", # 6
  'Chief cells' = "PGA4", # 6
  
  # 'Endocrine cells' = "GAST", # 5
  'Endocrine cells' = "GHRL", # 5
  'Endocrine cells' = "SST", # 5
  'Endocrine cells' = "HDC", # 5
  
  'Tuft cells' = "TRPM5", # 5
  'Tuft cells' = "SH2D6", # 5
  
  'Tumor cells' = "CEACAM5", # 6
  'Tumor cells' = "CEACAM6", # 6
  'Tumor cells' = "EPCAM", # 9
  'Tumor cells' = "CDH17", # 9
  'Tumor cells' = "MDM2", # P53 pathway
  
  'Enterocytes' = "FABP1", # 6
  'Enterocytes' = "APOA1", # 6
  'Enterocytes' = "FABP2", # 8
  
  'Primary cancer cells' = "NOTCH2", # 7
  'Primary cancer cells' = "NOTCH2NL", # 7
  'Primary cancer cells' = "KIF5B", # 7
  'Primary cancer cells' = "ERBB4", # 7
  
  'Metastatic cancer cells' = "CDK12", # 7
  'Metastatic cancer cells' = "ERBB2", # 7
  'Metastatic cancer cells' = "CLDN11", # 7 
  
  'Driver gene' = "FOS", # 7
  'Driver gene' = "JUN", # 7
  
  'Metaplasia genes' = "CD44", # 8 癌前
  # 'Metaplasia genes' = "TFF2", # 8
  # 'Metaplasia genes' = "SOX9", # 8
  'Metaplasia genes' = "HLA-DPA1", # 8
  'Metaplasia genes' = "HLA-DRA", # 8
  
  'Intestinal metaplasia cells' = "TFF3", # 9
  'Intestinal metaplasia cells' = "CDX1", # 9
  'Intestinal metaplasia cells' = "CDX2", # 9
  
  'Goblet cells' = "MUC2", # 6, 9
  'Goblet cells' = "ITLN1", # 6, 9
  'Goblet cells' = "HES6", # 9
  
  'Metaplastic-stem cells' = "EPHB2", # 6
  'Metaplastic-stem cells' = "SOX9", # 6  
  'Metaplastic-stem cells' = "OLFM4", # 6, 9
  'Metaplastic-stem cells' = "REG1A", # 9
  'Metaplastic-stem cells' = "CLDN3" # 9
)

# find markers
Idents(Normal) <- "RNA_snn_res.0.6"
Normal_markers <-
  FindAllMarkers(
    JoinLayers(Normal),
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    only.pos = TRUE
  )

Normal_top5_markers <- 
  Normal_markers %>% 
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

names(Normal_top5_markers) <- sub('.*\\.','',names(Normal_top5_markers))



# dotplot
pdf(file.path(result_dir, str_glue("dotplot dims{dims} major_marker.pdf")), width = 35, height = 6)
#DotPlot(Normal, group.by = "RNA_snn_res.0.6", features = major_marker_ls)+RotatedAxis()
#DotPlot(Normal, group.by = "RNA_snn_res.0.6", features = Epithelial_major_marker_ls)+RotatedAxis()
DotPlot(Normal, group.by = "RNA_snn_res.0.6", features = Normal_top5_markers)+RotatedAxis()
# DotPlot(Normal, group.by = "RNA_snn_res.0.6", features = Fibro_major_marker_ls2)+RotatedAxis()
# DotPlot(Normal, group.by = "RNA_snn_res.0.6", features = B_major_marker_ls)+RotatedAxis()
# DotPlot(Normal, group.by = "RNA_snn_res.0.6", features = macrophage_major_marker_ls)+RotatedAxis()
# DotPlot(Normal, group.by = "RNA_snn_res.0.6", features = T_major_marker_ls)+RotatedAxis()
# DotPlot(Normal, group.by = "RNA_snn_res.0.6", features = CD8_major_marker_ls)+RotatedAxis()
# DotPlot(Normal, group.by = "RNA_snn_res.0.6", features = CD4_major_marker_ls)+RotatedAxis()
dev.off()



pdf(file.path(result_dir, str_glue("dotplot dims{dims} major_marker_ls.pdf")), width = 25, height = 6)
DotPlot(Normal, group.by = "RNA_snn_res.0.6", features = major_marker_ls) + RotatedAxis()
dev.off()

pdf(file.path(result_dir, str_glue("dotplot dims{dims} Epithelial_major_marker_ls.pdf")), width = 25, height = 6)
DotPlot(Normal, group.by = "RNA_snn_res.0.6", features = Epithelial_major_marker_ls) + RotatedAxis()
dev.off()



# DotPlot(Normal, group.by = "RNA_snn_res.0.8", features = Fibro_major_marker_ls2)+RotatedAxis()

VlnPlot(Normal, features = 'nCount_RNA', pt.size = 0)




### annotation -----------
cell_cluster <- c(
  '0' = "Normal_C0", 
  '1' = "Normal_C1", # MUC5AC, GKN1
  '2' = "Normal_C2",
  '3' = "Normal_C3",
  '4' = "Normal_C4",
  '5' = "Normal_C5", 
  '6' = "Normal_C6", 
  '7' = "Normal_C7",
  '8' = "Normal_C8", # 
  '9' = "Normal_C9", # 
  '10' = "Normal_C10", # 
  '11' = "Normal_C11"
)


Idents(Normal) <- "RNA_snn_res.0.6"

Normal <- RenameIdents(Normal, cell_cluster)
Normal$subtype <- Idents(Normal)

pdf(file.path(result_dir,str_glue("Dimplot dim{dims} annotation.pdf")), height = 6, width = 8)
DimPlot(Normal, group.by = "subtype", label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

pdf(file.path(result_dir, str_glue("Dimplot dim{dims} annotation split.pdf")), height = 7, width = 15)
DimPlot(Normal, group.by = "subtype",split.by = 'group', label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

saveRDS(Normal, file.path(result_dir, 'Normal.rds'))






### spatial plot------------
spatial_outdir <- file.path(result_dir, "spatial plot")
dir.create(spatial_outdir)

cells_boundary_ls <- readRDS("R_OUTPUT/data_for_analysis/cells_boundary_ls.rds")

Normal_meta <-
  Normal[[]] %>%
  rownames_to_column("cell_id")

all_data_meta <- 
  all_data[[]] %>% 
  rownames_to_column("cell_id")

sp_p_ls <- 
  cells_boundary_ls %>% 
  imap(function(cells_boundary, project){ # browser()
    # image_wh_r <- image_wh_r_ls[[project]]
    
    dat <-
      cells_boundary %>%
      right_join(all_data_meta, by = "cell_id")
    
    dat2 <-
      cells_boundary %>%
      right_join(Normal_meta, by = "cell_id")
    
    p <- dat[which(dat$RCTD != 'Normal'),] %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) +
      geom_polygon(aes(group = cell_id), fill = "gray90")  +
      scale_y_reverse() +
      coord_equal() + 
      # coord_fixed(ratio = 1/image_wh_r) +
      theme_void() +
      ggnewscale::new_scale_fill() +
      geom_polygon(data = dat2, aes(x = vertex_x, y = vertex_y, fill = subtype, group = cell_id)) +
      ggsci::scale_fill_igv() +
      theme_void()
    p
  })

sp_p_ls %>% 
  iwalk(function(p, grp){
    pdf(file.path(spatial_outdir, str_glue("Normal subtype spatial plot {grp}.pdf")), height = 18, width = 18)
    print(p)
    dev.off()
  })







## cluster7 是混杂细胞，剔除后重新保存
Normal_filt <- Normal[, Normal$subtype != "Normal_C7"]
unique(Normal_filt$subtype)


saveRDS(Normal_filt, file.path(result_dir, "Normal_filt.rds"))


pdf(file.path(result_dir,str_glue("Dimplot dim{dims} annotation - trimmed.pdf")), height = 6, width = 8)
DimPlot(Normal_filt, group.by = "subtype", label = TRUE, reduction = 'harmonyumap', repel = TRUE) + NoLegend()
dev.off()

pdf(file.path(result_dir, str_glue("Dimplot dim{dims} annotation split - trimmed.pdf")), height = 7, width = 15)
DimPlot(Normal_filt, group.by = "subtype", split.by = 'group', label = TRUE, reduction = 'harmonyumap', repel = TRUE) + NoLegend()
dev.off()



### cell props barplot  --------
Normal_meta <- Normal_filt[[]]

Normal_meta$sample <- factor(Normal_meta$sample, 
                             levels = c("PA_1", "PA_2", "PC"))

p_col_proj <-
  Normal_meta %>%
  group_by(sample, subtype) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = sample, y = freq, fill = subtype)) + 
  geom_col() +
  # ggsci::scale_fill_igv() + 
  scale_fill_manual(values = plan) + 
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )


# export
ggsave(
  p_col_proj,
  filename = file.path(result_dir, "cell_freq_sample.pdf"),
  height = 6,
  width  = 5
)

p_col_proj <- 
  Normal_meta %>%
  group_by(group, subtype) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = group, y = freq, fill = subtype)) + 
  geom_col() +
  # ggsci::scale_fill_igv() + 
  scale_fill_manual(values = plan) + 
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )

ggsave(
  p_col_proj,
  filename = file.path(result_dir, "cell_freq_group.pdf"),
  height = 6,
  width  = 4
)

# cellType_marker_ls <- c(
#   'Intestinal metaplasia cells' = "CDX1", # 9
#   'Intestinal metaplasia cells' = "CDX2", # 9
#   
#   'Surface mucous cells' = "MUC5AC", # 1
#   'Surface mucous cells' = "GKN1", # 9 
#   
#   'Proliferating genes' = "STMN1", # 4
#   'Proliferating genes' = "MKI67", # 6
#   
#   'Metastatic cancer cells' = "ERBB2", # 7
#   
#   'Tumor cells' = "EPCAM", # 9
#   'Tumor cells' = "MDM2", # P53 pathway
#   
#   "Fibroblasts" = "PDGFRA", 
#   "Fibroblasts" = "CCDC80",
#   "Fibroblasts" = "SFRP4",
#   
#   'myCAF' = 'MMP11',
#   'myCAF' = 'CTHRC1',
#   'myCAF' = 'COL11A1',
#   'myCAF' = 'POSTN',
#   'myCAF' = 'COL10A1',
#   'myCAF' = 'TNC',
#   
#   'Chief cells' = "LIPF", # 4
#   'Chief cells' = "PGA3", # 6
#   
#   'T cells' = 'CD3E', 
#   'T cells' = 'CD3G', 
#   'CD8_Tex'='HAVCR2',
#   'CD8_Tex'='CTLA4',
#   'CD8_Tex'='LAG3',
#   'CD8_Tex'='TIGIT',
#   
#   'Goblet cells' = "MUC2", # 6, 9
#   'Goblet cells' = "ITLN1", # 6, 9
#   'Goblet cells' = "HES6", # 9
#   
#   'Enteroendocrine cells' = "PROX1", # 9
#   
#   "B cells" = "CD19",
#   "B cells" = "CD79A",
#   "Plasma" = "MZB1",
#   "Plasma" = "XBP1",
#   
#   'Primary cancer cells' = "NOTCH2", # 7
#   'Primary cancer cells' = "KIF5B", # 7
#   
#   'Metaplasia genes' = "CD44", # 8 癌前
#   
#   'Metaplastic-stem cells' = "EPHB2", # 6
#   'Metaplastic-stem cells' = "SOX9" # 6  
# )
# 
# pdf(file.path(Epithelial_out_dir, "dotplot celltype main_marker.pdf"), width = 20, height = 8)
# DotPlot(Epithelial, group.by = "subtype", features = cellType_marker_ls, dot.scale = 10) + theme_90angle
# dev.off()




### volcano-------
Epithelial_markers <-
  FindAllMarkers(
    JoinLayers(Normal_filt),
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    only.pos = TRUE
  )
write.csv(Epithelial_markers, file.path(result_dir, 'Epithelial_diff_gene.csv'))




### Normal_filt 映射------
## cluster7 是混杂细胞，剔除后重新保存
library(tibble)
library(dplyr)
meta.data.t <- 
  Normal[['subtype']] %>% 
  rownames_to_column('cells') %>% 
  right_join(
    all_data[['RCTD_subtype']] %>% rownames_to_column('cells'), 
    by = "cells"
  ) %>% 
  mutate(
    RCTD_subtype2 = ifelse((as.character(RCTD_subtype) == 'Normal'), as.character(subtype), as.character(RCTD_subtype))) %>% 
  dplyr::select(cells, RCTD_subtype2) %>% 
  column_to_rownames('cells')


all_data <- AddMetaData(all_data, metadata = meta.data.t['RCTD_subtype2'])


saveRDS(all_data, "all_data_RCTD_subtype2.rds")









library(EnhancedVolcano)
type <- unique(Epithelial_markers$cluster)
volcano_ls <- type %>% 
  setNames(.,.) %>%
  imap(function(type,name){
    markers<- Epithelial_markers[Epithelial_markers$cluster == type, ]
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
      filename = file.path(result_dir, str_glue("{name} diffgene valcano.pdf")),
      height = 6,
      width = 8
    )
  })

Epithelial_markers_ls <- list(Epithelial_markers)
names(Epithelial_markers_ls) <- 'Normal'

# scatter
p_volcanoplot_major_ls <- 
  Epithelial_markers_ls %>%
  imap(function(dat, name){
    
    set.seed(1234)
    rect_height <- 0.5
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
      filename = file.path(result_dir, str_glue("{name} diffgene dot.pdf")),
      height = 10,
      width = 20
    )
  })




### enrichment ----------------------------
Epithelial_de_genes_for_enrichmenet <- 
  Epithelial_markers %>% 
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  dplyr::filter(pct.1 > 0.15) %>%
  dplyr::filter(abs(avg_log2FC) > 0.5) %>% 
  group_by(cluster) %>%
  nest() %>%
  deframe() %>%
  map("gene")


enrich_ls <- 
  Epithelial_de_genes_for_enrichmenet %>%
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


er_output<-'Epithelial_subtype/enrichment' # dir.create(er_output)

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

Epithelial_de_genes_id_for_enrichmenet <-
  Epithelial_de_genes_for_enrichmenet %>%
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
  Epithelial_de_genes_id_for_enrichmenet %>%
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







# 20260427 RCTD 导出 ----------
annotation_outdir <- "20260427-annotation"
dir.create(annotation_outdir)

all_data <- readRDS("all_data_RCTD_subtype2.rds")

# all_cols <- c(
#   "B cell" = "#5050FFFF",
#   "Endothelium" = "#CE3D32FF",
#   "Fibroblast" = "#749B58FF",
#   "Macrophage" = "#F0E685FF",
#   "Mast cell" = "#466983FF",
#   "Monocyte" = "#BA6338FF",
#   "Neutrophil" = "#5DB1DDFF",
#   "Normal" = "#802268FF",
#   "T cell" = "#6BD76BFF",
#   "Tumor" = "#D595A7FF",
#   "Unknown" = "#924822FF"
# )
all_cols <- c(
  "B cell" = "#5050FF",
  "Endothelium" = "#CE3D32",
  "Fibroblast" = "#749B58",
  "Macrophage" = "#F0E685",
  "Mast cell" = "#466983",
  "Monocyte" = "#BA6338",
  "Neutrophil" = "#5DB1DD",
  "Normal" = "#802268", 
  "T cell" = "#6BD76B",
  "Tumor" = "#D595A7",
  "Unknown" = "#924822"
)
rctd_metadata <- all_data[[]] %>% 
  rownames_to_column("cell_id") %>% 
  mutate(cell_id = str_remove(cell_id, paste0(project, "_"))) %>% 
  dplyr::select(cell_id, group = RCTD)

rctd_metadata$color <- all_cols[rctd_metadata$group]

write.csv(rctd_metadata, file = file.path(annotation_outdir, "RCTD_annotation.csv"), row.names = FALSE)





rctd_metadata <- all_data[[]] %>% 
  rownames_to_column("cell_id") %>% 
  mutate(cell_id = str_remove(cell_id, paste0(project, "_"))) %>% 
  dplyr::select(cell_id, group = RCTD_subtype2)


write.csv(rctd_metadata, file = file.path(annotation_outdir, "Normal_Tumor_sub_annotation.csv"), row.names = FALSE)


