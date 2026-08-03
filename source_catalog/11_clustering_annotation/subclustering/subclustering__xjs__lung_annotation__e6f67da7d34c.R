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

lung <- subset(JoinLayers(all_data), group %in% c('HS_lung', 'NC_lung'))

## PCA --------
library(future)
plan(multisession, workers = 10)
options(future.globals.maxSize = 100 * 1024 ^ 3) # 100Gb

DefaultAssay(lung) <- "RNA"
lung <- FindVariableFeatures(lung)
lung <- ScaleData(lung, vars.to.regress = 'nCount_RNA')
lung <- RunPCA(lung)

lung_out_dir <-'lung_annotation' 
dir.create(lung_out_dir)


### pca inspector ------
pdf(file.path(lung_out_dir, 'pca plot.pdf'), height = 6, width = 8)
ElbowPlot(lung, ndims = 50)
dev.off()


## dimension reduction --------
dims <- 20

lung <- RunUMAP(lung, dims = 1:dims, return.model = TRUE)

pdf(file.path(lung_out_dir, str_glue("Dimplot dim{dims} sample.pdf")), height = 6, width = 8)
DimPlot(lung, group.by = "sample", label = TRUE) + NoLegend()
dev.off()

pdf(file.path(lung_out_dir, str_glue("Dimplot dim{dims} group.pdf")), height = 6, width = 8)
DimPlot(lung, group.by = "group", label = TRUE) + NoLegend()
dev.off()


lung <- lung %>% 
  harmony::RunHarmony(
    dims.use = 1:dims,
    group.by.vars = "sample",
    reduction.save = "harmony"
    # theta = 3, # default 2, bigger value -- more diversity per cluster
    # lambda = 0.5, # default, smaller value -- more integrated
  )

lung <- RunUMAP(
  lung, 
  dims = 1:dims, 
  reduction = "harmony", 
  # n.neighbors = 50,
  reduction.name = "harmonyumap", 
  reduction.key  = "harmonyUMAP_",
  return.model = TRUE
)


pdf(file.path(lung_out_dir, str_glue("Dimplot dim{dims} sample - harmony.pdf")), height = 6, width = 8)
DimPlot(lung, group.by = "sample", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()

pdf(file.path(lung_out_dir, str_glue("Dimplot dim{dims} group - harmony.pdf")), height = 6, width = 8)
DimPlot(lung, group.by = "group", label = TRUE, reduction = "harmonyumap")
dev.off() 



reduction.use <- 'harmony' # 'harmony' 'pca'
### cluster -------
lung <- 
  FindNeighbors(
    lung ,
    reduction = reduction.use,
    k.param = 50, 
    dims = 1:dims
  )


lung <-
  FindClusters(
    lung,
    resolution = c(0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.5)
  )

# dimplot
pdf(file.path(lung_out_dir, str_glue("Dimplot dim{dims} cluster.pdf")), height = 6, width = 8)
DimPlot(lung, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(lung, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(lung, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(lung, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(lung, group.by = "RNA_snn_res.1", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(lung, group.by = "RNA_snn_res.1.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(lung, group.by = "RNA_snn_res.1.5", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()

DimPlot(lung, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = 'harmonyumap')
DotPlot(lung, group.by = "RNA_snn_res.0.6", features = major_marker_ls)+RotatedAxis()




lung_major_marker_ls <- c(
  'AT1 cells' = "Ager", # 
  'AT1 cells' = "Sema3b", # 
  
  'AT2 cells' = "Abca3", #  
  'AT2 cells' = "Sftpa1", #  
  'AT2 cells' = "Sftpb", #  
  'AT2 cells' = "Sftpc", #  
  'AT2 cells' = "Hhip", #  
  
  'Basal cells' = "Krt5", # 
  'Basal cells' = "Ngfr", # 
  'Basal cells' = "Tp63", # 
  'Basal cells' = "Krt15", # 
  
  'Deuterosomal cells' = 'Cdc208', 
  
  'Goblet cells' = 'Muc5ac', 
  'Goblet cells' = 'Ceacam5', 
  
  'Pulmonary neuroendocrine cells' = 'Grp', 
  'Pulmonary neuroendocrine cells' = 'Calca', 
  
  'Tuft cells' = 'Pou2f3', 
  'Tuft cells' = 'Ascl2', 
  
  'Lonocyes' = 'Ascl3', 
  'Lonocyes' = 'Bsnd', 
  'Lonocyes' = 'Foxi1', 
  
  'Hillock cells' = 'Krt14', 
  'Hillock cells' = 'Krt6a', 
  
  'Clumnar cells' = 'Foxj1', 
  'Clumnar cells' = 'Rsph1', 
  'Clumnar cells' = 'Cdhr3', 
  'Clumnar cells' = 'Cdhr4', 
  
  'Mucous cells' = 'Muc5b', 
  'Mucous cells' = 'Spdef', 
  'Mucous cells' = 'Bpifb2', 
  
  'Serous cells' = 'Lyz', 
  'Serous cells' = 'Ltf', 
  
  'Secretory cells' = 'Scgb1a1', 
  'Secretory cells' = 'Scgb3a2', 
  
  'Ciliated cells' = "Akap14", # 
  'Ciliated cells' = "Aldh3b1", # 
  'Ciliated cells' = "Ankrd66", # 
  'Ciliated cells' = "C11orf97", # 
  'Ciliated cells' = "Dnai1", # 
  
  # 'Club cells' = "Scgb1a1", # 
  'Club cells' = "Scgb3a1", # 
  'Club cells' = "Pigr"
)

# find markers
Idents(lung) <- "RNA_snn_res.0.6"
lung_markers <-
  FindAllMarkers(
    JoinLayers(lung),
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    only.pos = TRUE
  )

lung_top5_markers <- 
  lung_markers %>% 
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

names(lung_top5_markers) <- sub('.*\\.','',names(lung_top5_markers))



# dotplot
pdf(file.path(lung_out_dir, str_glue("dotplot dims{dims} major_marker.pdf")), width = 35, height = 6)
DotPlot(lung, group.by = "RNA_snn_res.0.6", features = major_marker_ls)+RotatedAxis()
DotPlot(lung, group.by = "RNA_snn_res.0.6", features = lung_major_marker_ls)+RotatedAxis()
DotPlot(lung, group.by = "RNA_snn_res.0.6", features = lung_top5_markers)+RotatedAxis()
# DotPlot(lung, group.by = "RNA_snn_res.0.6", features = Fibro_major_marker_ls2)+RotatedAxis()
# DotPlot(lung, group.by = "RNA_snn_res.0.6", features = B_major_marker_ls)+RotatedAxis()
# DotPlot(lung, group.by = "RNA_snn_res.0.6", features = macrophage_major_marker_ls)+RotatedAxis()
# DotPlot(lung, group.by = "RNA_snn_res.0.6", features = T_major_marker_ls)+RotatedAxis()
# DotPlot(lung, group.by = "RNA_snn_res.0.6", features = CD8_major_marker_ls)+RotatedAxis()
# DotPlot(lung, group.by = "RNA_snn_res.0.6", features = CD4_major_marker_ls)+RotatedAxis()
dev.off()

# DotPlot(lung, group.by = "RNA_snn_res.0.8", features = Fibro_major_marker_ls2)+RotatedAxis()
VlnPlot(lung, features = 'nCount_RNA', group.by = "RNA_snn_res.0.6", pt.size = 0)




## annotation -----------
lung$RNA_snn_res.0.8 %>% table()
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
  '0' = "Endothelial", # Pecam1
  '1' = "Endothelial", # Cdh5
  '2' = "AT2 cell", # Abca3
  '3' = "Fibroblast", # Pdgfra, Col*
  '4' = "Myeloid", # Cd14, Cd68, Cd83
  '5' = "Clumnar cell", # Foxj1, Rsph1
  '6' = "Pericyte", # Pdgfrb, Mcam
  '7' = "Endothelial", # Pecam1
  '8' = "AT1 cell", # Ager
  '9' = "Fibroblast", # Pdgfra, Col*
  '10' = "Endothelial", # Pecam1
  '11' = "Smooth muscle cell", # Myh11, Actg2
  '12' = "Clumnar cell", # Foxj1, Rsph1
  '13' = "B cell", # Cd19, Cd79a, Ms4a1
  '14' = "Neutrophil", # Csf3r
  '15' = "Macrophage", # Cd68
  '16' = "T cell", # Cd3d, Cd3e, Cd3g, Bcl11b
  '17' = "Lymphatic endothelial cell", # Lyve1, Ccl21a, Ccl21b
  '18' = "Smooth muscle cell", # Pln
  '19' = "Myeloid" # Cd14, Cd68, Cd83
)

Idents(lung) <- "RNA_snn_res.0.6"

lung <- RenameIdents(lung, cell_cluster)
lung$subtype <- Idents(lung)

pdf(file.path(lung_out_dir,str_glue("Dimplot dim{dims} annotation.pdf")), height = 6, width = 10)
DimPlot(lung, group.by = "subtype", label = TRUE, reduction = 'harmonyumap',repel = TRUE)
dev.off()

pdf(file.path(lung_out_dir, str_glue("Dimplot dim{dims} annotation split.pdf")), height = 6, width = 13)
DimPlot(lung, group.by = "subtype",split.by = 'group', label = TRUE, reduction = 'harmonyumap', repel = TRUE) + NoLegend()
dev.off()




## cell props barplot  --------
lung_meta <- lung[[]]


p_col_proj <-
  lung_meta %>%
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
  filename = file.path(lung_out_dir, "cell_freq_sample.pdf"),
  height = 6,
  width  = 8
)

p_col_proj <- 
  lung_meta %>%
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
  filename = file.path(lung_out_dir, "cell_freq_group.pdf"),
  height = 6,
  width  = 4
)

cellType_marker_ls <- c(
  'Endothelial' = 'Pecam1', 
  'Endothelial' = 'Cdh5', 
  
  'AT2 cell' = 'Abca3', 
  
  'Fibroblast' = 'Col1a1', 
  'Fibroblast' = 'Col1a2', 
  'Fibroblast' = 'Pdgfra', 
  
  'Myeloid' = 'Cd14', 
  
  'Clumnar cell' = 'Foxj1', 
  'Clumnar cell' = 'Rsph1', 
  
  'Pericyte' = 'Pdgfrb', 
  'Pericyte' = 'Mcam', 
  
  'AT1 cell' = 'Ager', 
  
  'Smooth muscle cell' = 'Myh11', 
  'Smooth muscle cell' = 'Actg2', 
  'Smooth muscle cell' = 'Pln', 
  
  'B cell' = 'Cd19', 
  'B cell' = 'Cd79a', 
  'B cell' = 'Ms4a1', 
  
  'Neutrophil' = 'Csf3r', 
  
  'Macrophage' = 'Cd68', 
  
  'T cell' = 'Cd3d', 
  'T cell' = 'Cd3e', 
  'T cell' = 'Cd3g', 
  'T cell' = 'Bcl11b', 
  
  'Lymphatic endothelial cell' = 'Lyve1', 
  'Lymphatic endothelial cell' = 'Ccl21a', 
  'Lymphatic endothelial cell' = 'Ccl21b'
)

pdf(file.path(lung_out_dir, "dotplot celltype main_marker.pdf"), width = 15, height = 7)
DotPlot(lung, group.by = "subtype", features = cellType_marker_ls, dot.scale = 10) + theme_90angle
dev.off()

saveRDS(lung, file.path(lung_out_dir, 'lung.rds'))
