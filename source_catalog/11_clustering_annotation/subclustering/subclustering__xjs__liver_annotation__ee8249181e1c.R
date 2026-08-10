library(tidyverse)
library(Seurat)

major_marker_ls <- c(
  'Epithelial'  = 'Epcam',
  'Epithelial'  = 'Krt5',
  'Epithelial'  = 'Krt14',
  'Epithelial'  = 'Krt20',
  
  'Endothelial' = 'Pecam1', 
  'Endothelial' = 'Vwf', 
  'Endothelial' = 'Aqp1', 
  'Endothelial' = 'Cdh5',
  'Endothelial' = 'Plvap',
  
  "Fibroblast" = "Pdgfra", 
  "Fibroblast" = "Ccdc80",
  "Fibroblast" = "Col1a1",
  "Fibroblast" = "Col1a2",
  "Fibroblast" = "Dcn",
  "Fibroblast" = "Acta2",
  
  "B cell" = "Cd19",
  "B cell" = "Cd79a",
  "B cell" = "Ms4a1",
  "Plasma" = "Tnfrsf17",
  "Plasma" = "Jchain", 
  "Plasma" = "Igj", 
  "Plasma" = "Mzb1",
  "Plasma" = "Xbp1",
  
  'T cell' = 'Cd3d', 
  'T cell' = 'Cd3e', 
  'T cell' = 'Cd3g', 
  'T cell' = 'Trac', 
  'T cell' = 'Bcl11b', 
  
  "NK" = "Gnly",
  "NK" = "Nkg7",
  "NK" = "Ncam1", 
  "NK" = "Unc13d", 
  
  'Proliferating' = 'Mki67',
  'Proliferating' = 'Top2a',
  
  "Pericyte" = "Pdgfrb",
  "Pericyte" = "Rgs5",
  "Pericyte" = "Mcam",
  
  "Monocyte" = "Cd14", 
  "Monocyte" = "Fcgr3a", 
  "Monocyte" = "Itgam",
  
  "Macrophage" = "Cd68",
  "Macrophage" = "Cd163",
  "Macrophage" = "C1qa",
  "Macrophage" = "C1qc", 
  'Macrophage'  = 'Csf1r',
  
  "Dendritic cell" = "Cd83",
  'Dendritic cell' = "Lilra4", 
  'Dendritic cell' = "Clec9a",
  'Dendritic cell' = "Ido1",
  'Dendritic cell' = "Cd1c",
  'Dendritic cell' = "Fcer1a", 
  'Dendritic cell' = "Cst3", 
  'pDC' = "Irf4", 
  'pDC' = "Irf7", 
  'pDC' = "Irf8", 
  'Lamp3_DC' = "Lamp3",
  
  'Netrophil' = "Csf3r", 
  'Netrophil' = "Fcgr3b",
  'Netrophil' = "S100a8",
  'Netrophil' = "S100a9",
  'Netrophil' = "G0s2",
  'Netrophil' = "Retnlg",
  
  "Mast cell" = "Ms4a2",
  "Mast cell" = "Kit",
  "Mast cell" = "Hdc",
  "Mast cell" = "Tpsab1",
  "Mast cell" = "Ctsg",
  "Mast cell" = "Cpa3", 
  
  # 'Smooth muscle cell' = 'Pdgfrb', 
  'Smooth muscle cell' = 'Pln', 
  'Smooth muscle cell' = 'Mylk', 
  "Smooth muscle cell" = 'Myh11', 
  "Smooth muscle cell" = 'Actg2', 
  
  "Schwann cell" = 'Mpz', 
  "Schwann cell" = 'Sox10', 
  
  'Neuron' = 'Tubb3', 
  'Neuron' = 'Snap25', 
  
  "Mesothelial" = 'Msln', 
  
  'Erythrocytes' = "Hba-a1",
  'Erythrocytes' = "Hba-a2"
)

theme_90angle <- 
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
    strip.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )

all_data <- readRDS(file = file.path("R_OUTPUT", "all_data.rds"))
all_data$group %>% table()

liver <- subset(JoinLayers(all_data), group %in% c('HS_liver', 'NC_liver'))

## PCA --------
library(future)
plan(multisession, workers = 10)
options(future.globals.maxSize = 100 * 1024 ^ 3) # 100Gb

DefaultAssay(liver) <- "RNA"
liver <- FindVariableFeatures(liver)
liver <- ScaleData(liver, vars.to.regress = 'nCount_RNA')
liver <- RunPCA(liver)

liver_out_dir <-'liver_annotation' 
dir.create(liver_out_dir)


### pca inspector ------
pdf(file.path(liver_out_dir, 'pca plot.pdf'), height = 6, width = 8)
ElbowPlot(liver, ndims = 50)
dev.off()


## dimension reduction --------
dims <- 20

liver <- RunUMAP(liver, dims = 1:dims,return.model = TRUE)

pdf(file.path(liver_out_dir, str_glue("Dimplot dim{dims} sample.pdf")), height = 6, width = 8)
DimPlot(liver, group.by = "sample", label = TRUE) + NoLegend()
dev.off()

pdf(file.path(liver_out_dir, str_glue("Dimplot dim{dims} group.pdf")), height = 6, width = 8)
DimPlot(liver, group.by = "group", label = TRUE) + NoLegend()
dev.off()


liver <- liver %>% 
  harmony::RunHarmony(
    dims.use = 1:dims,
    group.by.vars = "sample",
    reduction.save = "harmony"
    # theta = 3, # default 2, bigger value -- more diversity per cluster
    # lambda = 0.5, # default, smaller value -- more integrated
  )

liver <- RunUMAP(
  liver, 
  dims = 1:dims, 
  reduction = "harmony", 
  reduction.name = "harmonyumap", 
  reduction.key  = "harmonyUMAP_",
  return.model = TRUE
)


pdf(file.path(liver_out_dir, str_glue("Dimplot dim{dims} sample - harmony.pdf")), height = 6, width = 8)
DimPlot(liver, group.by = "sample", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()

pdf(file.path(liver_out_dir, str_glue("Dimplot dim{dims} group - harmony.pdf")), height = 6, width = 8)
DimPlot(liver, group.by = "group", label = TRUE, reduction = "harmonyumap")
dev.off() 



reduction.use <- 'harmony' # 'harmony' 'pca'
### cluster -------
liver <- 
  FindNeighbors(
    liver ,
    reduction = reduction.use,
    dims = 1:dims,
    k.param = 20
  )


liver <-
  FindClusters(
    liver,
    resolution = c(0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.5)
  )

# dimplot
pdf(file.path(liver_out_dir, str_glue("Dimplot dim{dims} cluster.pdf")), height = 6, width = 8)
DimPlot(liver, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(liver, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(liver, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(liver, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(liver, group.by = "RNA_snn_res.1", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(liver, group.by = "RNA_snn_res.1.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(liver, group.by = "RNA_snn_res.1.5", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()

DimPlot(liver, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = 'harmonyumap')
DotPlot(liver, group.by = "RNA_snn_res.0.8", features = major_marker_ls)+RotatedAxis()




liver_major_marker_ls <- c(
  'Kupffer' = 'Cd5l', 
  'Kupffer' = 'Vsig4', 
  
  'LAM' = 'Gpnmb', 
  'LAM' = 'Fabp5', 
  
  'Cholangiocyte' = 'Sox9', 
  'Cholangiocyte' = 'Spp1', 
  'Cholangiocyte' = 'Cd24a', 
  
  'LSEC1' = 'Kdr', 
  'LSEC1' = 'Lyve1', 
  'LSEC2' = 'Fbxl7', 
  'LSEC2' = 'Lrmda', 
  
  'Hepatocyte' = 'Cyp2e1', 
  'Hepatocyte' = 'Rora', 
  'Hepatocyte' = 'Gpr183', 
  'Hepatocyte' = 'Socs3', 
  'Hepatocyte' = 'Maff', 
  'Hepatocyte' = 'Dusp5', 
  'Contral_hep' = 'Oat', 
  'Contral_hep' = 'Glul', 
  'Mid_hep' = 'Igfbp2', 
  'Mid_hep' = 'Hamp2', 
  'Potral_hep' = 'Hal', 
  'Potral_hep' = 'Sds'
)

# find markers
Idents(liver) <- "RNA_snn_res.0.8"
liver_markers <-
  FindAllMarkers(
    JoinLayers(liver),
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    only.pos = TRUE
  )

liver_top5_markers <- 
  liver_markers %>% 
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>% 
  dplyr::filter(avg_log2FC > 0.75) %>% 
  dplyr::filter(pct.1 - pct.2 > 0.1) %>% 
  arrange(desc(avg_log2FC)) %>%
  # slice_head(n = 30) %>%
  # arrange(desc(pct.1)) %>%
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

names(liver_top5_markers) <- sub('.*\\.','',names(liver_top5_markers))



# dotplot
pdf(file.path(liver_out_dir, str_glue("dotplot dims{dims} major_marker.pdf")), width = 35, height = 6)
DotPlot(liver, group.by = "RNA_snn_res.0.8", features = major_marker_ls)+RotatedAxis()
DotPlot(liver, group.by = "RNA_snn_res.0.8", features = liver_major_marker_ls)+RotatedAxis()
DotPlot(liver, group.by = "RNA_snn_res.0.8", features = liver_top5_markers)+RotatedAxis()
# DotPlot(liver, group.by = "RNA_snn_res.0.6", features = Fibro_major_marker_ls2)+RotatedAxis()
# DotPlot(liver, group.by = "RNA_snn_res.0.6", features = B_major_marker_ls)+RotatedAxis()
# DotPlot(liver, group.by = "RNA_snn_res.0.6", features = macrophage_major_marker_ls)+RotatedAxis()
# DotPlot(liver, group.by = "RNA_snn_res.0.6", features = T_major_marker_ls)+RotatedAxis()
# DotPlot(liver, group.by = "RNA_snn_res.0.6", features = CD8_major_marker_ls)+RotatedAxis()
# DotPlot(liver, group.by = "RNA_snn_res.0.6", features = CD4_major_marker_ls)+RotatedAxis()
dev.off()

# DotPlot(liver, group.by = "RNA_snn_res.0.8", features = Fibro_major_marker_ls2)+RotatedAxis()
VlnPlot(liver, features = 'nCount_RNA', pt.size = 0)




## annotation -----------
liver$RNA_snn_res.0.8 %>% table()
# cell_cluster <- c(
#   '0' = "", # 
#   '1' = "", # 
#   '2' = "", # 
#   '3' = "", # 
#   '4' = "", # 
#   '5' = "", # 
#   '6' = "", # 
#   '7' = "", # 
#   '8' = "", # 
#   '9' = "", # 
#   '10' = "", # 
#   '11' = "", # 
#   '12' = "", # 
#   '13' = "", # 
#   '14' = "", # 
#   '15' = "", # 
#   '16' = "", # 
#   '17' = "", # 
#   '18' = "", # 
#   '19' = "", # 
#   '20' = "", # 
#   '21' = "", # 
#   '22' = "" # 
# )

cell_cluster <- c(
  '0' = "Hepatocyte", # Hal
  '1' = "Hepatocyte", # Cyp2e1
  '2' = "Hepatocyte", # Hal
  '3' = "Hepatocyte", # Cyp2e1
  '4' = "Endothelial", # Pecam1
  '5' = "Pericyte", # Plvap
  '6' = "Hepatocyte", # Hal
  '7' = "Macrophage", # Cd68
  '8' = "Hepatocyte", # 
  '9' = "Hepatocyte", # Cyp2e1
  '10' = "Endothelial", # Pecam1
  '11' = "Macrophage", # Cd68
  '12' = "Fibroblast", # Col1a1
  '13' = "Cholangiocyte", # Epcam, Sox9, Cd24a
  '14' = "Pericyte", # Plvap
  '15' = "Netrophil", # 
  '16' = "Endothelial", # Pecam1
  '17' = "T cell", # 
  '18' = "Macrophage", # Cd68
  '19' = "B cell"
)


Idents(liver) <- "RNA_snn_res.0.8"

liver <- RenameIdents(liver, cell_cluster)
liver$subtype <- Idents(liver)

pdf(file.path(liver_out_dir,str_glue("Dimplot dim{dims} annotation.pdf")), height = 6, width = 10)
DimPlot(liver, group.by = "subtype", label = TRUE, reduction = 'harmonyumap',repel = TRUE)
dev.off()

pdf(file.path(liver_out_dir, str_glue("Dimplot dim{dims} annotation split.pdf")), height = 7, width = 18)
DimPlot(liver, group.by = "subtype",split.by = 'group', label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()




## cell props barplot  --------
liver_meta <- liver[[]]


p_col_proj <-
  liver_meta %>%
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
  filename = file.path(liver_out_dir, "cell_freq_sample.pdf"),
  height = 6,
  width  = 8
)

p_col_proj <- 
  liver_meta %>%
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
  filename = file.path(liver_out_dir, "cell_freq_group.pdf"),
  height = 6,
  width  = 4
)

cellType_marker_ls <- c(
  'Hepatocyte' = 'Hal', 
  'Hepatocyte' = 'Cyp2e1', 
  
  'Endothelial' = 'Pecam1', 
  
  "Pericyte" = 'Pdgfrb', 
  "Pericyte" = 'Plvap', 
  
  'Macrophage' = 'Cd68', 
  'Macrophage' = 'Cd163', 
  'Macrophage' = 'Csf1r', 
  
  'Fibroblast' = 'Pdgfra',
  'Fibroblast' = 'Col1a1',
  'Fibroblast' = 'Col1a2',

  "Cholangiocyte" = 'Epcam', 
  "Cholangiocyte" = 'Sox9',
  "Cholangiocyte" = 'Cd24a', 
  
  "Netrophil" = 'Cd14', 
  "Netrophil" = 'Csf3r',
  
  "T cell" = 'Cd3d', 
  "T cell" = 'Cd3e', 
  "T cell" = 'Cd3g', 
  
  'B cell' = 'Cd19', 
  'B cell' = 'Cd79a', 
  'B cell' = 'Ms4a1'
)

pdf(file.path(liver_out_dir, "dotplot celltype main_marker.pdf"), width = 12, height = 6)
DotPlot(liver, group.by = "subtype", features = cellType_marker_ls, dot.scale = 10) + theme_90angle
dev.off()

saveRDS(liver, file.path(liver_out_dir, 'liver.rds'))
