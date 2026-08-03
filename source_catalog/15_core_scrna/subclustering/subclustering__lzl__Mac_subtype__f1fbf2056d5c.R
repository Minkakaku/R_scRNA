library(Seurat)
library(tidyverse)
library(future)
# set default ---------
out_dir <- 'Mac_subtype'
dir.create(out_dir)
theme_90angle <-
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
    strip.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )

all_data <- readRDS('/mnt/NexData_file/home/lzl/case/XE2602120662-01/output/data/adata_filtered_add_celltype.rds')

# change assay
all_data[['RNA']] <- as(all_data[['RNA']], "Assay5")
all_data@assays$RNA$data <- NULL
all_data[['RNA']] <- split(all_data[['RNA']], all_data$project, layers = "counts")
all_data$annotation <- gsub("Netrophil", "Neutrophil", all_data$annotation)
plan(sequential) 
all_data <- all_data %>% NormalizeData()
#####
Mac_data <- subset(all_data,annotation %in% c("Macrophage"))
Mac_data <- JoinLayers(Mac_data)
## PCA --------
Mac_data <- FindVariableFeatures(Mac_data)
length(VariableFeatures(Mac_data)) #2000
Mac_data <- ScaleData(Mac_data)
Mac_data <- RunPCA(Mac_data)

### pca inspector ------
pdf(file.path(out_dir, 'pca plot.pdf'), height = 6, width = 8)
ElbowPlot(Mac_data, ndims = 50)
dev.off()

pdf(file.path(out_dir, 'dimheatmap 1-5.pdf'), height = 6, width = 10)
DimHeatmap(Mac_data, dims = 1:5)
dev.off()
pdf(file.path(out_dir, 'dimheatmap 6-10.pdf'), height = 6, width = 10)
DimHeatmap(Mac_data, dims = 6:10)
dev.off()
pdf(file.path(out_dir, 'dimheatmap 11-15.pdf'), height = 6, width = 10)
DimHeatmap(Mac_data, dims = 11:15)
dev.off()
pdf(file.path(out_dir, 'dimheatmap 16-20.pdf'), height = 6, width = 10)
DimHeatmap(Mac_data, dims = 16:20)
dev.off()

## dimension reduction --------
dims <- 15 # 18

library(future)
options(future.globals.maxSize = 10 * 1024 ^ 3) # 10Gb
plan(multisession, workers = 10)

Mac_data <- RunUMAP(Mac_data, dims = 1:dims,return.model = TRUE)

pdf(file.path(out_dir, str_glue("Dimplot dim{dims} sample.pdf")), height = 8, width = 8)
DimPlot(Mac_data, group.by = "sample", label = TRUE) + NoLegend()
dev.off()

Mac_data <- Mac_data %>% 
  harmony::RunHarmony(
    dims.use = 1:dims,
    group.by.vars = "project",
    # theta = 4, # default 2, bigger value -- more diversity per cluster
    # sigma = 0.5,
    # lambda = 0.5, # default, smaller value -- more integrated
    reduction.save = "harmony"
  )

Mac_data <- RunUMAP(
  Mac_data, 
  dims = 1:dims, 
  reduction = "harmony", 
  reduction.name = "harmonyumap", 
  reduction.key  = "harmonyUMAP_",
  return.model = TRUE
)

pdf(file.path(out_dir, str_glue("Dimplot dim{dims} sample - harmony.pdf")), height = 8, width = 8)
DimPlot(Mac_data, group.by = "sample", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()

reduction.use <- 'harmony' # 'harmony' 'pca'


### cluster -------
Mac_data  <-
  FindNeighbors(
    Mac_data ,
    reduction = reduction.use,
    dims = 1:dims,
    k.param = 20
  )


Mac_data <-
  FindClusters(
    Mac_data,
    resolution = c(0.1,0.2, 0.4, 0.6, 0.8)
  )



# dimplot
pdf(file.path(out_dir, str_glue("Dimplot dim{dims} cluster.pdf")), height = 6, width = 6)
DimPlot(Mac_data, group.by = "RNA_snn_res.0.1", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Mac_data, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Mac_data, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Mac_data, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Mac_data, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(Macrophage, group.by = "RNA_snn_res.1", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(Macrophage, group.by = "RNA_snn_res.1.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(Macrophage, group.by = "RNA_snn_res.1.5", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()

Idents(Mac_data) <- "RNA_snn_res.0.4"
###marker_ls
major_marker_ls <- c(
  # 第0行：B cell
  "B cell" = "Cd19", 
  "B cell" = "Cd79a", 
  "B cell" = "Ms4a1", 
  "Plasma" = "Tnfrsf17", 
  "Plasma" = "Jchain", 
  "Plasma" = "Igj", 
  "Plasma" = "Mzb1", 
  "Plasma" = "Xbp1", 
  
  "T cell" = "Cd3d", 
  "T cell" = "Cd3e", 
  "T cell" = "Cd3g", 
  "T cell" = "Trac", 
  "T cell" = "Bcl11b",
  "T cell" = "Cd4", 
  "T cell" = "Foxp3", 
  "T cell" = "Cd8a", 
  "T cell" = "Gzmb",
  
  "NK" = "Gnly", 
  "NK" = "Nkg7", 
  "NK" = "Ncam1", 
  "NK" = "Unc13d", 
  
  # 第11行：Plasma
  
  "Adipocyte" = "Plin1",
  "Adipocyte" = "Plin4",
  "Adipocyte" = "Adipoq",
  "Adipocyte" = "Fabp4",
  "Adipocyte" = "Lpl",
  
  
  "Endothelial" = "Pecam1",
  "Endothelial" = "Vwf", 
  "Endothelial" = "Aqp1", 
  "Endothelial" = "Cdh5", 
  "Endothelial" = "Plvap", 
  
  # 第4行：Fibroblast
  "Fibroblast" = "Pdgfra", 
  "Fibroblast" = "Ccdc80", 
  "Fibroblast" = "Col1a1", 
  "Fibroblast" = "Col1a2", 
  "Fibroblast" = "Dcn", 
  "Fibroblast" = "Acta2",
  
  "Pericyte" = "Pdgfrb", 
  "Pericyte" = "Rgs5", 
  "Pericyte" = "Mcam", 
  
  # 第3行：Macrophage
  "Macrophage" = "Cd68", 
  "Macrophage" = "Cd163", 
  "Macrophage" = "C1qa", 
  "Macrophage" = "C1qc", 
  "Macrophage" = "Csf1r", 
  #"Macrophage" = "Il1b", 
  
  # 第10行：Lamp3_DC
  "Lamp3_DC" = "Lamp3",
  
  # 第9行：Smooth muscle cell
  "Smooth muscle cell" = "Pln", 
  "Smooth muscle cell" = "Mylk", 
  "Smooth muscle cell" = "Myh11", 
  "Smooth muscle cell" = "Actg2", 
  "Smooth muscle cell" = "Tagln", 
  
  # 其他已在图中出现的细胞类型
  "Dendritic cell" = "Cd83", 
  "Dendritic cell" = "Lilra4", 
  #"Dendritic cell" = "Clec9a", 
  "Dendritic cell" = "Ido1", 
  "Dendritic cell" = "Cd1c", 
  "Dendritic cell" = "Fcer1a", 
  "Dendritic cell" = "Cst3", 
  
  "Epithelial" = "Epcam",
  "Epithelial" = "Krt20",
  
  "Monocyte" = "Cd14",
  "Monocyte" = "Fcgr3a",
  "Monocyte" = "Itgam",
  
  # 修正拼写错误 Netrophil → Neutrophil
  "Neutrophil" = "Csf3r", 
  "Neutrophil" = "Fcgr3b", 
  "Neutrophil" = "S100a8", 
  "Neutrophil" = "S100a9", 
  "Neutrophil" = "G0s2", 
  "Neutrophil" = "Retnlg", 
  
  
  "Neuron" = "Tubb3",
  "Neuron" = "Snap25",
  
  
  "pDC" = "Irf4",
  "pDC" = "Irf7",
  "pDC" = "Irf8",
  
  "cDC1" = "Xcr1",
  "cDC1" = "Clec9a",
  
  "Proliferating" = "Mki67",
  "Proliferating" = "Top2a",
  
  "Schwann cell" = "Mpz",
  "Schwann cell" = "Sox10",
  "Schwann cell" = "Ngfr",
  
  "M1" = "Il1b",
  "M2" = "Mrc1",
  "M2" = "Arg1",
  "M2" = "Clec7a",
  
  "Mesothelial" = "Msln",
  
  "Mast cell" = "Ms4a2", 
  "Mast cell" = "Kit", 
  "Mast cell" = "Hdc", 
  "Mast cell" = "Tpsab1", 
  "Mast cell" = "Ctsg", 
  "Mast cell" = "Cpa3",
  
  # 小鼠红细胞血红蛋白规范命名
  "Erythrocytes" = "Hba-a1", 
  "Erythrocytes" = "Hba-a2"
)
Macrophage_major_marker_ls <- c(
  "M1"= "Nos2",
  "M1"= "Il12",
  "M1"= "Fcgr1a",
  "M1"= "Fcgr1b",
  "M1"= "Fcgr1c",
  "M1"= "Cd80",
  "M1"= "Cxcr10",
  "M1"= "Il23",
  "M1"= "Cxcl9",
  "M1"= "Cxcl10",
  "M1"= "Cxcl11",
  "M1"= "Cd86",
  "M1"= "Il1a",
  "M1"= "Il1b",
  "M1"= "Il6",
  "M1"= "Tnf",
  "M1"= "Mhcii",
  "M1"= "Ccl5",
  "M1"= "Irf5",
  "M1"= "Irf1",
  "M1"= "Cd40",
  "M1"= "Kynu",
  "M1"= "Ccr7",
  "M1"= "Vcan",
  "M1"= "Cd300e",
  "M1"= "Inhba",
  "M1"= "Tnip3",
  
  "M2"= "Arg1",
  "M2"= "Il10",
  "M2"= "Cd32",
  "M2"= "Cd163",
  "M2"= "Fcer2",
  "M2"= "Cd200r1",
  "M2"= "Pdcd1lg2",
  "M2"= "Cd274",
  "M2"= "Marco",
  "M2"= "Csf1r",
  "M2"= "Mrc1",
  "M2"= "Il1rn",
  "M2"= "Il1r2",
  "M2"= "Il4r",
  "M2"= "Ccl4",
  "M2"= "Ccl13",
  "M2"= "Ccl20",
  "M2"= "Ccl17",
  "M2"= "Ccl18",
  "M2"= "Ccl22",
  "M2"= "Ccl24",
  "M2"= "Lyve1",
  "M2"= "Vegfa",
  "M2"= "Vegfb",
  "M2"= "Vegfc",
  "M2"= "Vegfd",
  "M2"= "Egf",
  "M2"= "Ctsa",
  "M2"= "Ctsb",
  "M2"= "Cstc",
  "M2"= "Ctsd",
  "M2"= "Tgfb1",
  "M2"= "Cd209",
  "M2"= "Cd163l1",
  "M2"= "Folr2",
  
  "M0"= "Cd68",
  "M0"= "C1qa",
  "M0"= "C1qb",
  "M0"= "Selenop",
  "M0"= "Lyz",
  "M0"= "Hla-dpb1",
  "M0"= "Hla-dpa1",
  "M0"= "Aif1", 
  
  "TAM" = "Spp1", 
  "TAM" = "Cqc"
)

existing_genes <- rownames(all_data)

# 2. 过滤 marker 列表，只保留存在的基因
Macrophage_major_marker_ls <- Macrophage_major_marker_ls[
  Macrophage_major_marker_ls %in% existing_genes
]

Macrophage_markers <-
  FindAllMarkers(
    JoinLayers(Mac_data),
    test.use = "wilcox",
    logfc.threshold = 0,
    min.pct = 0,
    only.pos = TRUE
  )
write.csv(Macrophage_markers, file = "Mac_subtype/Macrophage_all_markers.csv", row.names = TRUE)


Macrophage_top5_markers <- 
  Macrophage_markers %>% 
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

names(Macrophage_top5_markers) <- sub('.*\\.','',names(Macrophage_top5_markers))


# dotplot
# 如果已做其他细胞类型的亚群分析，这里就可以顺便画一下其他类型的细胞亚群基因，比如Fibro_major_marker_ls2，否则只用画T细胞
pdf(file.path(out_dir, str_glue("dotplot dims{dims} major_marker.pdf")), width = 25, height = 6)
DotPlot(Mac_data, group.by = "RNA_snn_res.0.4", features = major_marker_ls)+RotatedAxis()+theme_90angle
DotPlot(Mac_data, group.by = "RNA_snn_res.0.4", features = Macrophage_major_marker_ls)+RotatedAxis()+theme_90angle
DotPlot(Mac_data, group.by = "RNA_snn_res.0.4", features = Macrophage_top5_markers)+RotatedAxis()+theme_90angle
dev.off()

Mac_data$RNA_snn_res.0.1 %>% table()

saveRDS(Mac_data ,file = "Mac_subtype/Macrophage_cluster.rds")
##
## annotation ----------------------------------------------------------------
cell_cluster <- c(
  '0'  = "Unknown", # 
  '1'  = "Fibro-like macrophage", #
  '2'  = "M2-like macrophage",
  '3' ="Proliferating cell",
  '4' = "Unknown",
  "5" = "M2-like macrophage",
  "6" = "Mrc1+ CD163+ M2-like macrophage",
  "7" = "Irf5+ Irf1+ M1-like macrophages",
  "8" = "M1-like macrophage",
  "9" = "Proliferating cell"
)

Idents(Mac_data) <- "RNA_snn_res.0.4"



Mac_data <- RenameIdents(Mac_data, cell_cluster)
Mac_data$subtype <- Idents(Mac_data)


pdf(file.path(out_dir, str_glue("Dimplot dim{dims} annotation.pdf")), height = 6, width = 8)
DimPlot(Mac_data, group.by = "subtype", label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

pdf(file.path(out_dir, str_glue("Dimplot dim{dims} annotation split.pdf")), height = 6, width = 15)
DimPlot(Mac_data, group.by = "subtype",split.by = 'group', label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()




keep <- c(
  "M2-like macrophage",
  "Mrc1+ CD163+ M2-like macrophage",
  "M1-like macrophage",
  "Irf5+ Irf1+ M1-like macrophages"
)

Mac_clean <- subset(Mac_data, subtype %in% keep)

Mac_clean$subtype <- droplevels(Mac_clean$subtype)
Idents(Mac_clean) <- Mac_clean$subtype

table(Mac_clean$subtype)
saveRDS(Mac_clean, "Mac_subtype/Mac_clean.rds")
#
Mac_clean <- readRDS("Mac_subtype/Mac_clean.rds")
Mac_clean <- RunUMAP(Mac_clean, dims = 1:dims,n.neighbors = 50,return.model = TRUE)

pdf(file.path(out_dir, str_glue("Dimplot dim{dims} Mac annotation.pdf")), height = 6, width = 8)
DimPlot(Mac_clean, group.by = "subtype", label = TRUE, reduction = 'umap',repel = TRUE) + NoLegend()
dev.off()

pdf(file.path(out_dir, str_glue("Dimplot dim{dims} Mac annotation split.pdf")), height = 6, width = 15)
DimPlot(Mac_clean, group.by = "subtype",split.by = 'group', label = TRUE, reduction = 'umap',repel = TRUE) + NoLegend()
dev.off()

###############
#######
Mac_data_filter <- subset(Mac_clean,subtype%in%c(
  "M2-like macrophage",
  "Mrc1+ CD163+ M2-like macrophage",
  "M1-like macrophage",
  "Irf5+ Irf1+ M1-like macrophages"
))
## PCA --------
Mac_data_filter <- FindVariableFeatures(Mac_data_filter)
length(VariableFeatures(Mac_data_filter)) #2000
Mac_data_filter <- ScaleData(Mac_data_filter)
Mac_data_filter <- RunPCA(Mac_data_filter)

out_dir <- "Mac_subtype_filter"
dir.create(out_dir)
### pca inspector ------
pdf(file.path(out_dir, 'pca plot.pdf'), height = 6, width = 8)
ElbowPlot(Mac_data_filter, ndims = 50)
dev.off()

pdf(file.path(out_dir, 'dimheatmap 1-5.pdf'), height = 6, width = 10)
DimHeatmap(Mac_data_filter, dims = 1:5)
dev.off()
pdf(file.path(out_dir, 'dimheatmap 6-10.pdf'), height = 6, width = 10)
DimHeatmap(Mac_data_filter, dims = 6:10)
dev.off()
pdf(file.path(out_dir, 'dimheatmap 11-15.pdf'), height = 6, width = 10)
DimHeatmap(Mac_data_filter, dims = 11:15)
dev.off()
pdf(file.path(out_dir, 'dimheatmap 16-20.pdf'), height = 6, width = 10)
DimHeatmap(Mac_data_filter, dims = 16:20)
dev.off()

## dimension reduction --------
dims <- 15# 18

library(future)
options(future.globals.maxSize = 10 * 1024 ^ 3) # 10Gb
plan(multisession, workers = 10)

Mac_data_filter <- RunUMAP(Mac_data_filter, dims = 1:dims,return.model = TRUE)

pdf(file.path(out_dir, str_glue("Dimplot dim{dims} sample.pdf")), height = 8, width = 8)
DimPlot(Mac_data_filter, group.by = "sample", label = TRUE) + NoLegend()
dev.off()

Mac_data_filter <- Mac_data_filter %>% 
  harmony::RunHarmony(
    dims.use = 1:dims,
    group.by.vars = "project",
    # theta = 4, # default 2, bigger value -- more diversity per cluster
    # sigma = 0.5,
    # lambda = 0.5, # default, smaller value -- more integrated
    reduction.save = "harmony"
  )

Mac_data_filter <- RunUMAP(
  Mac_data_filter, 
  dims = 1:dims, 
  reduction = "harmony", 
  reduction.name = "harmonyumap", 
  reduction.key  = "harmonyUMAP_",
  return.model = TRUE
)

pdf(file.path(out_dir, str_glue("Dimplot dim{dims} sample - harmony.pdf")), height = 8, width = 8)
DimPlot(Mac_data_filter, group.by = "sample", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()

reduction.use <- 'harmony' # 'harmony' 'pca'


### cluster -------
Mac_data_filter  <-
  FindNeighbors(
    Mac_data_filter ,
    reduction = reduction.use,
    dims = 1:dims,
    k.param = 15  #20改15
  )


Mac_data_filter <-
  FindClusters(
    Mac_data_filter,
    resolution = c(0.1,0.2, 0.4, 0.6, 0.8)
  )



# dimplot
pdf(file.path(out_dir, str_glue("Dimplot dim{dims} cluster.pdf")), height = 6, width = 6)
DimPlot(Mac_data_filter, group.by = "RNA_snn_res.0.1", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Mac_data_filter, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Mac_data_filter, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Mac_data_filter, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Mac_data_filter, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(Macrophage, group.by = "RNA_snn_res.1", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(Macrophage, group.by = "RNA_snn_res.1.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(Macrophage, group.by = "RNA_snn_res.1.5", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()
Idents(Mac_data_filter) <- "RNA_snn_res.0.2"

####

Mac_f_markers <-
  FindAllMarkers(
    JoinLayers(Mac_data_filter),
    test.use = "wilcox",
    logfc.threshold = 0,
    min.pct = 0,
    only.pos = TRUE
  )
write.csv(Mac_f_markers, file = "Mac_subtype_filter/Mac_f_all_markers.csv", row.names = TRUE)


Mac_f_top5_markers <- 
  Mac_f_markers %>% 
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

names(Mac_f_top5_markers) <- sub('.*\\.','',names(Mac_f_top5_markers))


# dotplot
# 如果已做其他细胞类型的亚群分析，这里就可以顺便画一下其他类型的细胞亚群基因，比如Fibro_major_marker_ls2，否则只用画T细胞
pdf(file.path(out_dir, str_glue("dotplot dims{dims} major_marker.pdf")), width = 25, height = 6)
DotPlot(Mac_data_filter, group.by = "RNA_snn_res.0.2", features = major_marker_ls)+RotatedAxis()+theme_90angle
DotPlot(Mac_data_filter, group.by = "RNA_snn_res.0.2", features = Macrophage_major_marker_ls)+RotatedAxis()+theme_90angle
DotPlot(Mac_data_filter, group.by = "RNA_snn_res.0.2", features = Mac_f_top5_markers)+RotatedAxis()+theme_90angle
dev.off()

cell_cluster <- c(
  '0'  = "M2-like macrophage", # 
  '1'  = "Mrc1+ CD163+ M2-like macrophage", #
  '2'  = "M0-like macrophage",
  '3'  = "M1/M2-like macrophage",
  '4'  = "pDC-like macrophage"
)



Mac_data_filter <- RenameIdents(Mac_data_filter, cell_cluster)
Mac_data_filter$subtype <- Idents(Mac_data_filter)


pdf(file.path(out_dir, str_glue("Dimplot dim{dims} annotation.pdf")), height = 6, width = 8)
DimPlot(Mac_data_filter, group.by = "subtype", label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

pdf(file.path(out_dir, str_glue("Dimplot dim{dims} annotation split.pdf")), height = 6, width = 15)
DimPlot(Mac_data_filter, group.by = "subtype",split.by = 'group', label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()
saveRDS(Mac_data_filter, "Mac_subtype_filter/Mac_data_filter.rds")


Mac_plot <- subset(Mac_data_filter, subtype != "pDC-like macrophage")

saveRDS(Mac_plot, "Mac_subtype_filter/Mac_plot.rds")
pdf(file.path(out_dir, str_glue("Dimplot mac annotation_rmDC.pdf")), height = 6, width = 8)
DimPlot(
  Mac_plot,
  group.by = "subtype",
  label = TRUE,
  reduction = 'harmonyumap',
  repel = TRUE
) + NoLegend()
dev.off()

# 第二张分组拆分UMAP
pdf(file.path(out_dir, str_glue("Dimplot mac annotation split_rmDC.pdf")), height = 6, width = 15)
DimPlot(
  Mac_plot,
  group.by = "subtype",
  split.by = 'group',
  label = TRUE,
  reduction = 'harmonyumap',
  repel = TRUE
) + NoLegend()
dev.off()
#####prop
df = Mac_plot@meta.data
#df = read.csv("/lvdata/wzb/scRNA/FWSC20240734/ReAll/celltype/meta.xls",sep = "\t",row.names=1)
summary <- table(df[,c('subtype','group')])
summary_t <- as.data.frame(summary)
summary_prop <- summary_t %>%
  group_by(group) %>%             # 按分组(C/T)分组计算
  mutate(total = sum(Freq), prop = Freq/total) %>%  # 计算组内总数、占比
  ungroup()


cell_colors <- c(
  "M2-like macrophage" = "#E64B35",           # 珊瑚红
  "Mrc1+ CD163+ M2-like macrophage" = "#4DBBD5",  # 青蓝
  "M0-like macrophage" = "#00A087",           #  teal
  "M1/M2-like macrophage" = "#3C5488"        # 藏青
)

p2 <- ggplot(summary_prop, aes(x = group, y = prop, fill = subtype)) +
  geom_col(position = "stack", width = 0.7) +  # 堆叠柱状图
  scale_y_continuous(labels = scales::percent_format(), expand = c(0,0)) + # y轴转为百分比
  # 配色
  scale_fill_manual(values = cell_colors) +
  # 坐标轴、标签主题设置
  labs(
    x = "Group",
    y = "Cell proportion",
    fill = "Cell type"
  ) +
  # 论文级主题美化
  theme_classic() +
  theme(
    axis.title = element_text(size = 12, color = "black"),
    axis.text = element_text(size = 11, color = "black"),
    legend.title = element_text(size = 13, color = "black"),
    legend.text = element_text(size = 12, color = "black"),
    legend.position = "right",
    plot.margin = margin(10,20,10,10)
  )

pdf(file.path(out_dir, str_glue("mac celltype_prop.pdf")),width = 8,height = 8)
p2
dev.off()
####
df = Mac_plot@meta.data

# ========== 第一步：在循环外，汇总所有样本的数据 ==========

# 计算所有样本 × 所有亚型的频数
summary <- table(df[, c('subtype', 'project')])
summary_t <- as.data.frame(summary)

# 按样本(project)分组，计算每个样本内的比例
summary_prop <- summary_t %>%
  group_by(project) %>%
  mutate(total = sum(Freq), prop = Freq/total) %>%
  ungroup()

pdf(
  file.path(out_dir, "mac celltype_prop - all samples.pdf"),
  width = 10,   # 加宽以容纳多个样本
  height = 8
)

p <- ggplot(summary_prop, aes(x = project, y = prop, fill = subtype)) +
  geom_col(position = "stack", width = 0.7) +
  scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0)) +
  scale_fill_manual(values = cell_colors) +
  labs(
    x = "",
    y = "Cell proportion",
    fill = "Cell type",
    title = "Macrophage cell subtype proportion across samples"
  ) +
  theme_classic() +
  theme(
    axis.title = element_text(size = 12, color = "black"),
    axis.text.x = element_text(size = 11, color = "black", angle = 45, hjust = 1),  # x轴标签倾斜防重叠
    axis.text.y = element_text(size = 11, color = "black"),
    legend.title = element_text(size = 13, color = "black"),
    legend.text = element_text(size = 12, color = "black"),
    legend.position = "right",
    plot.margin = margin(10, 20, 10, 10),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
  )

print(p)
dev.off()
#############
##################
all_data$subtype <- as.character(all_data$annotation)  # 或 NA

# 提取 KC_data 的细胞名和对应注释
Mac_cells <- colnames(Mac_plot)
Mac_subtype <- Mac_plot$subtype  # 如 KC_data$cell_type 或 KC_data$RNA_snn_res.0.4

unique(Mac_plot$subtype)
# 创建命名向量
subtype_map <- setNames(as.character(Mac_subtype), Mac_cells)

# 赋值到 all_data
all_data$subtype[match(Mac_cells, colnames(all_data))] <- subtype_map[Mac_cells]

# 非 KC 细胞保留原注释
non_kc <- setdiff(colnames(all_data), Mac_cells)
all_data$subtype[non_kc] <- as.character(all_data$annotation[non_kc])
unique(all_data$subtype)


######
Mac_celltype <- c(
  "M2-like macrophage",
  "Mrc1+ CD163+ M2-like macrophage",
  "M0-like macrophage",
  "M1/M2-like macrophage"
)

plot_dat <-
  FetchData(
    all_data, 
    vars = c("subtype")  # 获取原始注释，不是二分的 immune
  ) %>%
  rownames_to_column('cell_id')

p_sp_immnue_ls <- cells_position_ls %>%
  imap(function(coords, sam){
    dat2 <- coords %>% left_join(plot_dat, by = "cell_id")
    
    # 创建颜色映射：免疫细胞用原色，其他用 gray60
    dat2$color_group <- ifelse(
      dat2$subtype %in% Mac_celltype,
      dat2$subtype,           # 保留原始细胞类型名
      "Others"                   # 其他统一
    )
    
    # 构建颜色向量
    plot_colors <- cell_colors[Mac_celltype]  # 提取免疫细胞颜色
    plot_colors["Others"] <- "gray60"             # 其他颜色
    
    # 只保留 dat2 中实际存在的类型
    present_types <- unique(dat2$color_group)
    plot_colors <- plot_colors[names(plot_colors) %in% present_types]
    
    # 透明度：免疫细胞不透明，其他半透明
    alpha_values <- setNames(
      ifelse(names(plot_colors) == "Others", 0.3, 1),
      names(plot_colors)
    )
    
    # 确保因子顺序：免疫细胞在前，Others 最后
    type_levels <- c(setdiff(present_types, "Others"), "Others")
    dat2$color_group <- factor(dat2$color_group, levels = type_levels)
    
    p <- ggplot(dat2, aes(x = SPATIAL_1, y = SPATIAL_2)) +
      geom_tile(aes(fill = color_group, alpha = color_group), height = 8, width = 8) +
      theme_void() +
      coord_equal() +
      scale_y_reverse() +
      
      scale_alpha_manual(values = alpha_values) +
      
      scale_fill_manual(
        values = plot_colors,
        guide = guide_legend(override.aes = list(alpha = 1))  # 图例统一不透明
      ) +
      
      theme(
        panel.background = element_rect(fill = NA, color = NA),
        plot.title = element_text(hjust = 0.5, size = 16)
      ) +
      ggtitle(sam)
    
    return(p)
  })

pdf(
  file.path(out_dir, "spatial plot of Mac subtype.pdf"),
  height = 8, width = 8
)
walk(p_sp_immnue_ls, print)
dev.off()
