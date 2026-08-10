library(Seurat)
library(tidyverse)
library(BPCells)
library(clusterProfiler)
library(org.Hs.eg.db)
library(scDblFinder)


all_data <- readRDS('/mnt/NexData_file/home/lzl/case/XE2510316140-03/output/data/adata_filtered_add_celltype.rds')

# change assay
all_data[['RNA']] <- as(all_data[['RNA']], "Assay5")
all_data@assays$RNA$data <- NULL
all_data[['RNA']] <- split(all_data[['RNA']], all_data$project, layers = "counts")

# normalize
library(future)
plan("sequential") 
all_data <- all_data %>% NormalizeData()

dir.create("data")
saveRDS(all_data ,file = "/mnt/NexData_file/home/lzl/20260515_XE2510316140-03/data/all_data.rds")

all_data$celltype %>% table()


major_marker_ls <- c(
  # 第0行：B cell
  "B cell" = "CD19", 
  "B cell" = "CD79A", 
  "B cell" = "MS4A1", 
  "Plasma" = "TNFRSF17", 
  "Plasma" = "JCHAIN", 
  "Plasma" = "IGJ", 
  "Plasma" = "MZB1", 
  "Plasma" = "XBP1", 
  
  "T cell" = "CD3D", 
  "T cell" = "CD3E", 
  "T cell" = "CD3G", 
  "T cell" = "TRAC", 
  "T cell" = "BCL11B",
  "T cell" = "CD4", 
  "T cell" = "FOXP3", 
  "T cell" = "CD8A", 
  "T cell" = "GZMB",
  
  "NK" = "GNLY", 
  "NK" = "NKG7", 
  "NK" = "NCAM1", 
  "NK" = "UNC13D", 
  
  # 第11行：Plasma
  
  "Adipocyte" = "PLIN1",
  "Adipocyte" = "PLIN4",
  "Adipocyte" = "ADIPOQ",
  "Adipocyte" = "FABP4",
  "Adipocyte" = "LPL",
  
  
  "Endothelial" = "PECAM1",
  "Endothelial" = "VWF", 
  "Endothelial" = "AQP1", 
  "Endothelial" = "CDH5", 
  "Endothelial" = "PLVAP", 
  
  # 第4行：Fibroblast
  "Fibroblast" = "PDGFRA", 
  "Fibroblast" = "CCDC80", 
  "Fibroblast" = "COL1A1", 
  "Fibroblast" = "COL1A2", 
  "Fibroblast" = "DCN", 
  "Fibroblast" = "ACTA2",
  
  "Pericyte" = "PDGFRB", 
  "Pericyte" = "RGS5", 
  "Pericyte" = "MCAM", 
  
  # 第3行：Macrophage
  "Macrophage" = "CD68", 
  "Macrophage" = "CD163", 
  "Macrophage" = "C1QA", 
  "Macrophage" = "C1QC", 
  "Macrophage" = "CSF1R", 
  #"Macrophage" = "IL1B", 
  
  # 第10行：Lamp3_DC
  "Lamp3_DC" = "LAMP3",
  
  # 第9行：Smooth muscle cell
  "Smooth muscle cell" = "PLN", 
  "Smooth muscle cell" = "MYLK", 
  "Smooth muscle cell" = "MYH11", 
  "Smooth muscle cell" = "ACTG2", 
  "Smooth muscle cell" = "TAGLN", 
  
  # 其他已在图中出现的细胞类型
  "Dendritic cell" = "CD83", 
  "Dendritic cell" = "LILRA4", 
  #"Dendritic cell" = "CLEC9A", 
  "Dendritic cell" = "IDO1", 
  "Dendritic cell" = "CD1C", 
  "Dendritic cell" = "FCER1A", 
  "Dendritic cell" = "CST3", 
  
  "Epithelial" = "EPCAM",
  "Epithelial" = "KRT20",
  
  "Monocyte" = "CD14",
  "Monocyte" = "FCGR3A",
  "Monocyte" = "ITGAM",
  
  "Neutrophil" = "CSF3R", 
  "Neutrophil" = "FCGR3B", 
  "Neutrophil" = "S100A8", 
  "Neutrophil" = "S100A9", 
  "Neutrophil" = "G0S2", 
  "Neutrophil" = "RETNLG", 
  
  
  "Neuron" = "TUBB3",
  "Neuron" = "SNAP25",
  
  
  "pDC" = "IRF4",
  "pDC" = "IRF7",
  "pDC" = "IRF8",
  
  "cDC1" = "XCR1",
  "cDC1" = "CLEC9A",
  
  "Proliferating" = "MKI67",
  "Proliferating" = "TOP2A",
  
  "Schwann cell" = "MPZ",
  "Schwann cell" = "SOX10",
  "Schwann cell" = "NGFR",
  
  "M1" = "IL1B",
  "M2" = "MRC1",
  "M2" = "ARG1",
  "M2" = "CLEC7A",
  
  "Mesothelial" = "MSLN",
  
  "Mast cell" = "MS4A2", 
  "Mast cell" = "KIT", 
  "Mast cell" = "HDC", 
  "Mast cell" = "TPSAB1", 
  "Mast cell" = "CTSG", 
  "Mast cell" = "CPA3",
  
  "Erythrocytes" = "HBA-A1", 
  "Erythrocytes" = "HBA-A2",
  
  'AM'="ITGAX",
  'AM'="FCGR1A",
  'AM'='SIGLEC1'
  
)

# Macrophage -------------------------------------------------------------------

Macrophage_out_dir <- "Macrophage_subtype" 
dir.create(Macrophage_out_dir)
Macrophage <- subset(all_data, annotation == 'Macrophage')

## PCA --------
Macrophage <- FindVariableFeatures(Macrophage)

options(future.globals.maxSize = 50 * 1024^3) # 把并行全局变量上限改成 50GB（足够用）
Macrophage <- ScaleData(Macrophage, vars.to.regress = 'nCount_RNA')
Macrophage <- RunPCA(Macrophage)

### pca inspector ------
pdf(file.path(Macrophage_out_dir, 'pca plot.pdf'), height = 6, width = 8)
ElbowPlot(Macrophage, ndims = 50)
dev.off()


Macrophage$annotation %>% table()
#Macrophage$celltype %>% table()
## dimension reduction --------
dims <- 20 # 18

library(future)
#options(future.globals.maxSize = 10 * 1024 ^ 3) # 10Gb
plan(multisession, workers = 10)

Macrophage <- RunUMAP(Macrophage, dims = 1:dims,return.model = TRUE)

pdf(file.path(Macrophage_out_dir, str_glue("Dimplot dim{dims} sample.pdf")), height = 6, width = 8)
DimPlot(Macrophage, group.by = "sample", label = TRUE) + NoLegend()
dev.off()

Macrophage <- Macrophage %>% 
  harmony::RunHarmony(
    dims.use = 1:dims,
    group.by.vars = "project",
    # theta = 4, # default 2, bigger value -- more diversity per cluster
    # sigma = 0.5,
    # lambda = 0.5, # default, smaller value -- more integrated
    reduction.save = "harmony"
  )

Macrophage <- RunUMAP(
  Macrophage, 
  dims = 1:dims, 
  reduction = "harmony", 
  reduction.name = "harmonyumap", 
  reduction.key  = "harmonyUMAP_",
  return.model = TRUE
)

pdf(file.path(Macrophage_out_dir, str_glue("Dimplot dim{dims} sample - harmony.pdf")), height = 6, width = 8)
DimPlot(Macrophage, group.by = "sample", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()

reduction.use <- 'harmony' # 'harmony' 'pca'


### cluster -------
Macrophage  <-
  FindNeighbors(
    Macrophage ,
    reduction = reduction.use,
    dims = 1:dims,
    k.param = 20
  )


Macrophage <-
  FindClusters(
    Macrophage,
    resolution = c(0.2, 0.4, 0.6, 0.8)
  )



# dimplot
pdf(file.path(Macrophage_out_dir, str_glue("Dimplot dim{dims} cluster.pdf")), height = 6, width = 8)
DimPlot(Macrophage, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Macrophage, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Macrophage, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Macrophage, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(Macrophage, group.by = "RNA_snn_res.1", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(Macrophage, group.by = "RNA_snn_res.1.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(Macrophage, group.by = "RNA_snn_res.1.5", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()

Idents(Macrophage) <- "RNA_snn_res.0.2"
###marker_ls
Macrophage_major_marker_ls <- c(
  "M1"= "NOS2",
  "M1"= "IL12",
  "M1"= "FCGR1A",
  "M1"= "FCGR1B",
  "M1"= "FCGR1C",
  "M1"= "CD80",
  "M1"= "CXCR10",
  "M1"= "IL23",
  "M1"= "CXCL9",
  "M1"= "CXCL10",
  "M1"= "CXCL11",
  "M1"= "CD86",
  "M1"= "IL1A",
  "M1"= "IL1B",
  "M1"= "IL6",
  "M1"= "TNF",
  "M1"= "MHCII",
  "M1"= "CCL5",
  "M1"= "IRF5",
  "M1"= "IRF1",
  "M1"= "CD40",
  "M1"= "KYNU",
  "M1"= "CCR7",
  "M1"= "VCAN",
  "M1"= "CD300E",
  "M1"= "INHBA",
  "M1"= "TNIP3",　
  
  "M2"= "ARG1",
  "M2"= "IL10",
  "M2"= "CD32",
  "M2"= "CD163",
  "M2"= "FCER2",
  "M2"= "CD200R1",
  "M2"= "PDCD1LG2",
  "M2"= "CD274",
  "M2"= "MARCO",
  "M2"= "CSF1R",
  "M2"= "MRC1",
  "M2"= "IL1RN",
  "M2"= "Il1R2",
  "M2"= "IL4R",
  "M2"= "CCL4",
  "M2"= "CCL13",
  "M2"= "CCL20",
  "M2"= "CCL17",
  "M2"= "CCL18",
  "M2"= "CCL22",
  "M2"= "CCL24",
  "M2"= "LYVE1",
  "M2"= "VEGFA",
  "M2"= "VEGFB",
  "M2"= "VEGFC",
  "M2"= "VEGFD",
  "M2"= "EGF",
  "M2"= "CTSA",
  "M2"= "CTSB",
  "M2"= "CSTC",
  "M2"= "CTSD",
  "M2"= "TGFB1",
  "M2"= "CD209",
  "M2"= "CD163L1",
  "M2"= "FOLR2",　
  
  "M0"= "CD68",
  "M0"= "C1QA",
  "M0"= "C1QB",
  "M0"= "SELENOP",
  "M0"= "LYZ",
  "M0"= "HLA-DPB1",
  "M0"= "HLA-DPA1",
  "M0"= "AIF1", 
  
  "TAM" = "SPP1", 
  "TAM" = "CQC"
)
existing_genes <- rownames(all_data)

# 2. 过滤 marker 列表，只保留存在的基因
Macrophage_major_marker <- Macrophage_major_marker_ls[
  Macrophage_major_marker_ls %in% existing_genes
]


Macrophage_markers <-
  FindAllMarkers(
    JoinLayers(Macrophage),
    test.use = "wilcox",
    logfc.threshold = 0,
    min.pct = 0,
    only.pos = TRUE
  )
write.csv(Macrophage_markers, file = "Macrophage_subtype/Macrophage_all_markers.csv", row.names = TRUE)


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
pdf(file.path(Macrophage_out_dir, str_glue("dotplot dims{dims} major_marker.pdf")), width = 25, height = 6)
DotPlot(Macrophage, group.by = "RNA_snn_res.0.2", features = major_marker_ls)+RotatedAxis()+theme_90angle
DotPlot(Macrophage, group.by = "RNA_snn_res.0.2", features = Macrophage_major_marker_ls)+RotatedAxis()+theme_90angle
#DotPlot(Macrophage, group.by = "RNA_snn_res.0.4", features = Macrophage_major_marker_ls_1)+RotatedAxis()+theme_90angle
#DotPlot(Macrophage, group.by = "RNA_snn_res.0.6", features = Macrophage_markers)+RotatedAxis()
DotPlot(Macrophage, group.by = "RNA_snn_res.0.2", features = Macrophage_top5_markers)+RotatedAxis()+theme_90angle
#DotPlot(Macrophage, group.by = "RNA_snn_res.0.6", features = Fibro_major_marker_ls2)+RotatedAxis()
dev.off()

Macrophage$RNA_snn_res.0.2 %>% table()

saveRDS(Macrophage ,file = "Macrophage_subtype/Macrophage_cluster.rds")
Macrophage <- readRDS("Macrophage_subtype/Macrophage_cluster.rds")

Macrophage$RNA_snn_res.0.2 %>% table()
## annotation ----------------------------------------------------------------
cell_cluster <- c(
  '0'  = "M0-like_macrophage", # 
  '1'  = "M2-like_macrophage", #
  '2'  = "Unkown", #
  '3'  = "CD1C+DC", # 
  '4'  = "Cytotoxic T-like_macrophage", # 
  '5'  = "IL1B+M1-like_macrophage", # 
  '6'  = "M2-like_macrophage", # 
  '7'  = "M1-like_macrophage", # 
  '8'  = "MMP9+macrophage",# 偏向CD8+TEM
  '9'  = "Epi-like_macrophage", # 
  '10' = "Monocyte", # 
  '11' = "Lamp3+DC", # 
  '12' = "CLEC9A+DC"
)

Idents(Macrophage) <- "RNA_snn_res.0.2"



Macrophage <- RenameIdents(Macrophage, cell_cluster)
Macrophage$subtype <- Idents(Macrophage)


pdf(file.path(Macrophage_out_dir, str_glue("Dimplot dim{dims} annotation.pdf")), height = 6, width = 8)
DimPlot(Macrophage, group.by = "subtype", label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

pdf(file.path(Macrophage_out_dir, str_glue("Dimplot dim{dims} annotation split.pdf")), height = 6, width = 15)
DimPlot(Macrophage, group.by = "subtype",split.by = 'group', label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()


all_data$Macrophage_annotation <- as.character(all_data$annotation)

# 2. 创建映射：细胞名 → 精细mac亚型
cell_subtype <- Macrophage$subtype
names(cell_subtype) <- colnames(Macrophage)

# 3. 只把 "mac" 替换成精细亚型，其他细胞保持原样
all_data$Macrophage_annotation[colnames(Macrophage)] <- as.character(cell_subtype)

# 查看结果（正确）
table(all_data$Macrophage_annotation)
table(all_data$annotation)
table(Macrophage$subtype)
saveRDS(all_data, "Macrophage_subtype/all_data_Macrophage_annotation.rds")

pdf(file.path(Macrophage_out_dir, str_glue("Macrophage ncount.pdf")), height = 6, width = 8)
VlnPlot(object = Macrophage,
        features = "nCount_RNA",  # 你要画的指标
        group.by = "RNA_snn_res.0.2",  # 按cluster分组
        pt.size = 0,  # 去掉散点，更干净
        fill = "seurat_clusters"  # 每个cluster自动配色
)
dev.off()

keep <- c(
  "M0-like_macrophage",
  "M2-like_macrophage",
  "CD1C+DC",
  "IL1B+M1-like_macrophage",
  "M1-like_macrophage",
  "MMP9+macrophage",
  "Monocyte",
  "Lamp3+DC",
  "CLEC9A+DC"
)

Macrophage_clean <- subset(Macrophage, subtype %in% keep)

Macrophage_clean$subtype <- droplevels(Macrophage_clean$subtype)
Idents(Macrophage_clean) <- Macrophage_clean$subtype

table(Macrophage_clean$subtype)
saveRDS(Macrophage_clean, "Macrophage_subtype/Macrophage_clean.rds")

Macrophage_clean

Macrophage_clean <- RunUMAP(Macrophage_clean, dims = 1:dims,n.neighbors = 50,return.model = TRUE)

pdf(file.path(Macrophage_out_dir, str_glue("Dimplot dim{dims} Macrophage annotation.pdf")), height = 6, width = 8)
DimPlot(Macrophage_clean, group.by = "subtype", label = TRUE, reduction = 'umap',repel = TRUE) + NoLegend()
dev.off()

pdf(file.path(Macrophage_out_dir, str_glue("Dimplot dim{dims} Macrophage annotation split.pdf")), height = 6, width = 15)
DimPlot(Macrophage_clean, group.by = "subtype",split.by = 'group', label = TRUE, reduction = 'umap',repel = TRUE) + NoLegend()
dev.off()


#画一个注释完的点图
Macrophage_major_marker_ls_2 <- c(
  "M0-like_macrophage"="CD68",
  
  "M2-like_macrophage"="MRC1",
  "M2-like_macrophage"="CCL18",
  "M2-like_macrophage"="CD163",
  "M2-like_macrophage"="VEGFB",
  
  "CD1C+DC"="CD1C",
  "CD1C+DC"="CD1D",
  "CD1C+DC"="CLEC10A",
  
  
  "IL1B+M1-like_macrophage"="IL1B", 
  
  "M1-like_macrophag"="CXCL9", 
  "M1-like_macrophag"="CXCL10", 
  "M1-like_macrophag"="IRF1", 
  
  
  "MMP9+macrophage"="MMP9", 
  

  "Monocyte"="CD14",
  "Monocyte"="ITGAM",
  "Monocyte"="FCGR3A",

  
  "Lamp3+DC"="LAMP3",
  "Lamp3+DC"="CD83",
  
  
  "CLEC9A+DC" = "CLEC9A",
  "CLEC9A+DC" = "XCR1"
)

pdf(file.path(Macrophage_out_dir, str_glue("dotplot Macrophage subtype major_marker.pdf")), width = 25, height = 6)
DotPlot(Macrophage_clean, group.by = "subtype", features = Macrophage_major_marker_ls_2)+RotatedAxis()+theme_90angle
dev.off()
###########提取mac单独做
Macrophage_clean <- readRDS(file.path(Macrophage_out_dir,"Macrophage_clean.rds"))

keep2 <- c(
  "M0-like_macrophage",
  "M2-like_macrophage",
  "IL1B+M1-like_macrophage",
  "M1-like_macrophage",
  "MMP9+macrophage"
)
Macrophage_clean_2 <-subset(Macrophage_clean, subtype %in% keep2)

Macrophage_clean_2 <- RunUMAP(Macrophage_clean_2, dims = 1:dims,n.neighbors = 50,return.model = TRUE)

## PCA --------
Macrophage_clean_2 <- FindVariableFeatures(Macrophage_clean_2)
length(VariableFeatures(Macrophage_clean_2)) #2000
Macrophage_clean_2 <- ScaleData(Macrophage_clean_2)
Macrophage_clean_2 <- RunPCA(Macrophage_clean_2)

out_dir <- "Macrophage_subtype_filter"
dir.create(out_dir)
### pca inspector ------
pdf(file.path(out_dir, 'pca plot.pdf'), height = 6, width = 8)
ElbowPlot(Macrophage_clean_2, ndims = 50)
dev.off()

pdf(file.path(out_dir, 'dimheatmap 1-5.pdf'), height = 6, width = 10)
DimHeatmap(Macrophage_clean_2, dims = 1:5)
dev.off()
pdf(file.path(out_dir, 'dimheatmap 6-10.pdf'), height = 6, width = 10)
DimHeatmap(Macrophage_clean_2, dims = 6:10)
dev.off()
pdf(file.path(out_dir, 'dimheatmap 11-15.pdf'), height = 6, width = 10)
DimHeatmap(Macrophage_clean_2, dims = 11:15)
dev.off()
pdf(file.path(out_dir, 'dimheatmap 16-20.pdf'), height = 6, width = 10)
DimHeatmap(Macrophage_clean_2, dims = 16:20)
dev.off()

## dimension reduction --------
dims <- 15# 18

library(future)
options(future.globals.maxSize = 10 * 1024 ^ 3) # 10Gb
plan(multisession, workers = 10)

Macrophage_clean_2 <- RunUMAP(Macrophage_clean_2, dims = 1:dims,return.model = TRUE)

pdf(file.path(out_dir, str_glue("Dimplot dim{dims} sample.pdf")), height = 8, width = 8)
DimPlot(Macrophage_clean_2, group.by = "sample", label = TRUE) + NoLegend()
dev.off()

Macrophage_clean_2 <- Macrophage_clean_2 %>% 
  harmony::RunHarmony(
    dims.use = 1:dims,
    group.by.vars = "project",
    # theta = 4, # default 2, bigger value -- more diversity per cluster
    # sigma = 0.5,
    # lambda = 0.5, # default, smaller value -- more integrated
    reduction.save = "harmony"
  )

Macrophage_clean_2 <- RunUMAP(
  Macrophage_clean_2, 
  dims = 1:dims, 
  reduction = "harmony", 
  reduction.name = "harmonyumap", 
  reduction.key  = "harmonyUMAP_",
  return.model = TRUE
)

pdf(file.path(out_dir, str_glue("Dimplot dim{dims} sample - harmony.pdf")), height = 8, width = 8)
DimPlot(Macrophage_clean_2, group.by = "sample", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()

reduction.use <- 'harmony' # 'harmony' 'pca'


### cluster -------
Macrophage_clean_2  <-
  FindNeighbors(
    Macrophage_clean_2 ,
    reduction = reduction.use,
    dims = 1:dims,
    k.param = 15  #20改15
  )


Macrophage_clean_2 <-
  FindClusters(
    Macrophage_clean_2,
    resolution = c(0.1,0.2, 0.4, 0.6, 0.8)
  )



# dimplot
pdf(file.path(out_dir, str_glue("Dimplot dim{dims} cluster.pdf")), height = 6, width = 6)
DimPlot(Macrophage_clean_2, group.by = "RNA_snn_res.0.1", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Macrophage_clean_2, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Macrophage_clean_2, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Macrophage_clean_2, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Macrophage_clean_2, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(Macrophage, group.by = "RNA_snn_res.1", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(Macrophage, group.by = "RNA_snn_res.1.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(Macrophage, group.by = "RNA_snn_res.1.5", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()
Idents(Macrophage_clean_2) <- "RNA_snn_res.0.4"


target_cells <- WhichCells(Macrophage_clean_2, idents = "1")
pdf(file.path(out_dir, str_glue("Dimplot dim{dims} cluster1.pdf")), height = 6, width = 6)
# 高亮显示
DimPlot(Macrophage_clean_2, 
        reduction = 'harmonyumap',
        cells.highlight = target_cells,
        cols.highlight = "red",
        cols = "grey90") + 
  NoLegend()
dev.off()
###

####

Mac_f_markers <-
  FindAllMarkers(
    JoinLayers(Macrophage_clean_2),
    test.use = "wilcox",
    logfc.threshold = 0,
    min.pct = 0,
    only.pos = TRUE
  )
write.csv(Mac_f_markers, file = "Macrophage_subtype_filter/Mac_f_all_markers.csv", row.names = TRUE)


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
DotPlot(Macrophage_clean_2, group.by = "RNA_snn_res.0.8", features = major_marker_ls)+RotatedAxis()+theme_90angle
DotPlot(Macrophage_clean_2, group.by = "RNA_snn_res.0.8", features = Macrophage_major_marker_ls)+RotatedAxis()+theme_90angle
DotPlot(Macrophage_clean_2, group.by = "RNA_snn_res.0.8", features = Mac_f_top5_markers)+RotatedAxis()+theme_90angle
dev.off()

pdf(file.path(out_dir, str_glue("vln_ncount.pdf")), width = 12, height = 6)
VlnPlot(Macrophage_clean_2, 
        features = c("nCount_RNA", "nFeature_RNA"), 
        group.by = "RNA_snn_res.0.4",
        ncol = 2,
        pt.size = 0) + 
  NoLegend()
DotPlot(Macrophage_clean_2, features = c("CCL18", "VEGFB", "CTSA")) + 
  RotatedAxis()
dev.off()
#####
cell_cluster <- c(
  '0'  = "Macrophage", # 
  '1'  = "Fibro-like_macrophage", #
  '2'  = "CYP27A1+ Macro", #
  '3'  = "STAB1+ Macro", # 
  '4'  = "MRC1+ Macro", # 
  '5'  = "Mix", # 
  '6'  = "CXCL10+ Macro", # 
  '7'  = "Neutrophil", # 
  '8'  = "Endothelial",# 
  '9'  = "MMP9+ Macro", # 
  '10' = "Mix", # 
  '11' = "Immune mixed" # 
)

Idents(Macrophage_clean_2) <- "RNA_snn_res.0.4"



Macrophage_clean_2 <- RenameIdents(Macrophage_clean_2, cell_cluster)
Macrophage_clean_2$subtype <- Idents(Macrophage_clean_2)


pdf(file.path(out_dir, str_glue("Dimplot dim{dims} annotation.pdf")), height = 6, width = 8)
DimPlot(Macrophage_clean_2, group.by = "subtype", label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

pdf(file.path(out_dir, str_glue("Dimplot dim{dims} annotation split.pdf")), height = 6, width = 15)
DimPlot(Macrophage_clean_2, group.by = "subtype",split.by = 'group', label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()
##
keep <- c(
  "CYP27A1+ Macro", #
  "STAB1+ Macro", # 
  "MRC1+ Macro", # 
  "CXCL10+ Macro", # 
  "MMP9+ Macro",
  "Macrophage"
)

Macrophage_clean_3 <- subset(Macrophage_clean_2, subtype %in% keep)

Macrophage_clean_3$subtype <- droplevels(Macrophage_clean_3$subtype)
Idents(Macrophage_clean_3) <- Macrophage_clean_3$subtype

table(Macrophage_clean_3$subtype)



Macrophage_clean_3 <- RunUMAP(Macrophage_clean_3, dims = 1:dims,n.neighbors = 50,return.model = TRUE)

pdf(file.path(out_dir, str_glue("Dimplot dim{dims} Macrophage annotation.pdf")), height = 6, width = 8)
DimPlot(Macrophage_clean_3, group.by = "subtype", label = TRUE, reduction = 'umap',repel = TRUE) + NoLegend()
dev.off()

pdf(file.path(out_dir, str_glue("Dimplot dim{dims} Macrophage annotation split.pdf")), height = 6, width = 15)
DimPlot(Macrophage_clean_3, group.by = "subtype",split.by = 'group', label = TRUE, reduction = 'umap',repel = TRUE) + NoLegend()
dev.off()
saveRDS(Macrophage_clean_3, file.path(out_dir,"Macrophage_clean_3"))
#####再删
Macrophage_clean_3 <- readRDS(file.path(out_dir,"Macrophage_clean_3"))
Macrophage_clean_3 <- subset(Macrophage_clean_3, RNA_snn_res.0.4 %in% c("0","2","3","4","5"))


Macrophage_clean_3 <- RunUMAP(Macrophage_clean_3, dims = 1:dims,n.neighbors = 50,return.model = TRUE)


Macrophage_clean_3$subtype <- gsub("M0", "Macrophage", Macrophage_clean_3$subtype)
# 修复 STAB1+ M2
Macrophage_clean_3$subtype <- gsub("STAB1+ M2", "STAB1+ Macro", Macrophage_clean_3$subtype, fixed = TRUE)

# 修复 MRC1+ M2
Macrophage_clean_3$subtype <- gsub("MRC1+ M2", "CXCL10+ Macro", Macrophage_clean_3$subtype, fixed = TRUE)

table(Macrophage_clean_3$subtype)
pdf(file.path(out_dir, str_glue("Dimplot dim{dims} Macrophage annotation.pdf")), height = 6, width = 8)
DimPlot(Macrophage_clean_3, group.by = "subtype", label = TRUE, reduction = 'umap',repel = TRUE) + NoLegend()
dev.off()

pdf(file.path(out_dir, str_glue("Dimplot dim{dims} Macrophage annotation split.pdf")), height = 6, width = 15)
DimPlot(Macrophage_clean_3, group.by = "subtype",split.by = 'group', label = TRUE, reduction = 'umap',repel = TRUE) + NoLegend()
dev.off()

#####重聚类
Macrophage_clean_4 <- Macrophage_clean_3
#Macrophage_clean_4 <- RunUMAP(Macrophage_clean_4, dims = 1:dims,n.neighbors = 50,return.model = TRUE)

## PCA --------
Macrophage_clean_4 <- FindVariableFeatures(Macrophage_clean_4)
length(VariableFeatures(Macrophage_clean_4)) #2000
Macrophage_clean_4 <- ScaleData(Macrophage_clean_4)
Macrophage_clean_4 <- RunPCA(Macrophage_clean_4)

out_dir <- "Macrophage_subtype_filter2"
dir.create(out_dir)
### pca inspector ------
pdf(file.path(out_dir, 'pca plot.pdf'), height = 6, width = 8)
ElbowPlot(Macrophage_clean_4, ndims = 50)
dev.off()

pdf(file.path(out_dir, 'dimheatmap 1-5.pdf'), height = 6, width = 10)
DimHeatmap(Macrophage_clean_4, dims = 1:5)
dev.off()
pdf(file.path(out_dir, 'dimheatmap 6-10.pdf'), height = 6, width = 10)
DimHeatmap(Macrophage_clean_4, dims = 6:10)
dev.off()
pdf(file.path(out_dir, 'dimheatmap 11-15.pdf'), height = 6, width = 10)
DimHeatmap(Macrophage_clean_4, dims = 11:15)
dev.off()
pdf(file.path(out_dir, 'dimheatmap 16-20.pdf'), height = 6, width = 10)
DimHeatmap(Macrophage_clean_4, dims = 16:20)
dev.off()

## dimension reduction --------
dims <- 10# 18

library(future)
options(future.globals.maxSize = 10 * 1024 ^ 3) # 10Gb
plan(multisession, workers = 10)

Macrophage_clean_4 <- RunUMAP(Macrophage_clean_4, dims = 1:dims,return.model = TRUE)

pdf(file.path(out_dir, str_glue("Dimplot dim{dims} sample.pdf")), height = 8, width = 8)
DimPlot(Macrophage_clean_4, group.by = "sample", label = TRUE) + NoLegend()
dev.off()

Macrophage_clean_4 <- Macrophage_clean_4 %>% 
  harmony::RunHarmony(
    dims.use = 1:dims,
    group.by.vars = "project",
    # theta = 4, # default 2, bigger value -- more diversity per cluster
    # sigma = 0.5,
    # lambda = 0.5, # default, smaller value -- more integrated
    reduction.save = "harmony"
  )

Macrophage_clean_4 <- RunUMAP(
  Macrophage_clean_4, 
  dims = 1:dims, 
  reduction = "harmony", 
  reduction.name = "harmonyumap", 
  reduction.key  = "harmonyUMAP_",
  return.model = TRUE
)

pdf(file.path(out_dir, str_glue("Dimplot dim{dims} sample - harmony.pdf")), height = 8, width = 8)
DimPlot(Macrophage_clean_4, group.by = "sample", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()

reduction.use <- 'harmony' # 'harmony' 'pca'


### cluster -------
Macrophage_clean_4  <-
  FindNeighbors(
    Macrophage_clean_4 ,
    reduction = reduction.use,
    dims = 1:dims,
    k.param = 15  #20改15
  )


Macrophage_clean_4 <-
  FindClusters(
    Macrophage_clean_4,
    resolution = c(0.1,0.2, 0.4, 0.6, 0.8)
  )



# dimplot
pdf(file.path(out_dir, str_glue("Dimplot dim{dims} cluster.pdf")), height = 6, width = 6)
DimPlot(Macrophage_clean_4, group.by = "RNA_snn_res.0.1", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Macrophage_clean_4, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Macrophage_clean_4, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Macrophage_clean_4, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Macrophage_clean_4, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(Macrophage, group.by = "RNA_snn_res.1", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(Macrophage, group.by = "RNA_snn_res.1.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(Macrophage, group.by = "RNA_snn_res.1.5", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()


Idents(Macrophage_clean_4) <- "RNA_snn_res.0.4"


####

Mac_f_markers <-
  FindAllMarkers(
    JoinLayers(Macrophage_clean_4),
    test.use = "wilcox",
    logfc.threshold = 0,
    min.pct = 0,
    only.pos = TRUE
  )
write.csv(Mac_f_markers, file = "Macrophage_subtype_filter2/Mac_f_all_markers.csv", row.names = TRUE)


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
DotPlot(Macrophage_clean_4, group.by = "RNA_snn_res.0.4", features = major_marker_ls)+RotatedAxis()+theme_90angle
DotPlot(Macrophage_clean_4, group.by = "RNA_snn_res.0.4", features = Macrophage_major_marker_ls)+RotatedAxis()+theme_90angle
DotPlot(Macrophage_clean_4, group.by = "RNA_snn_res.0.4", features = Mac_f_top5_markers)+RotatedAxis()+theme_90angle
dev.off()


cell_cluster <- c(
  '0'  = "Macrophage", # 
  '1'  = "Macrophage", #
  '2'  = "CHIT1+ Macro", #
  '3'  = "STAB1+ Macro", # 
  '4'  = "HSPA6+ Macro", # 
  '5'  = "CYP27A1+ Macro", # 
  '6'  = "MMP9+ Macro", # 
  '7'  = "MARCO+ Macro", # 
  '8'  = "Epi-like"
)

Idents(Macrophage_clean_4) <- "RNA_snn_res.0.4"



Macrophage_clean_4 <- RenameIdents(Macrophage_clean_4, cell_cluster)
Macrophage_clean_4$subtype <- Idents(Macrophage_clean_4)


pdf(file.path(out_dir, str_glue("Dimplot dim{dims} annotation.pdf")), height = 6, width = 8)
DimPlot(Macrophage_clean_4, group.by = "subtype", label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

pdf(file.path(out_dir, str_glue("Dimplot dim{dims} annotation split.pdf")), height = 6, width = 15)
DimPlot(Macrophage_clean_4, group.by = "subtype",split.by = 'group', label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()
saveRDS(Macrophage_clean_4, file.path(out_dir,"Macrophage_clean_4"))



Macrophage_clean_4_pure <- subset(Macrophage_clean_4,subtype!="Epi-like")
pdf(file.path(out_dir, str_glue("Dimplot dim{dims} annotation pure.pdf")), height = 6, width = 8)
DimPlot(Macrophage_clean_4_pure, group.by = "subtype", label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

pdf(file.path(out_dir, str_glue("Dimplot dim{dims} annotation split pure.pdf")), height = 6, width = 15)
DimPlot(Macrophage_clean_4_pure, group.by = "subtype",split.by = 'group', label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

saveRDS(Macrophage_clean_4_pure, file.path(out_dir,"Macrophage_clean_4_pure.rds"))
## cell props barplot  ----------------------------------------------------------
Macrophage_meta <- Macrophage_clean_4_pure[[]]
Macrophage_meta$subtype <- as.character(Macrophage_meta$subtype)
Macrophage_meta$subtype <- as.factor(Macrophage_meta$subtype)
Macrophage_meta <- Macrophage_meta %>%
  mutate(
    group_2 = case_when(
      grepl("S", sample) ~ "S",
      grepl("P", sample) ~ "P",
      grepl("R", sample) ~ "R",
      TRUE ~ NA_character_  # 都不匹配时设为 NA
    )
  )
Macrophage_meta <- Macrophage_meta  %>%
  mutate(
    group_3 = case_when(
      grepl("A", sample)&grepl("S", sample) ~ "A_S",
      grepl("A", sample)&grepl("P", sample) ~ "A_P",
      grepl("A", sample)&grepl("R", sample) ~ "A_R",
      grepl("C", sample)&grepl("S", sample) ~ "C_S",
      grepl("C", sample)&grepl("P", sample) ~ "C_P",
      TRUE ~ NA_character_  # 都不匹配时设为 NA
    )
  )

Macrophage_meta$sample <- factor(Macrophage_meta$sample,levels=c(Macrophage_meta$sample %>% unique()))

p_col_proj <-
  Macrophage_meta %>%
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
  filename = file.path(out_dir, "cell_freq_sample.pdf"),
  height = 6,
  width  = 8
)

p_col_proj <-
  Macrophage_meta %>%
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
  filename = file.path(out_dir, "cell_freq_group.pdf"),
  height = 6,
  width  = 4
)


p_col_proj <-
  Macrophage_meta %>%
  group_by(group_3, subtype) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = group_3, y = freq, fill = subtype)) + 
  geom_col() +
  ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )


ggsave(
  p_col_proj,
  filename = file.path(out_dir, "cell_freq_group3.pdf"),
  height = 6,
  width  = 4
)
###
######################
df = Macrophage_meta
#df = read.csv("/lvdata/wzb/scRNA/FWSC20240734/ReAll/celltype/meta.xls",sep = "\t",row.names=1)
summary <- table(df[,c('subtype','group')])
roe <- as.data.frame(ROIE(summary))
long_df <- tidyr::gather(roe, key = "group", value = "Ro/e", colnames(roe)[1:ncol(roe)-1])
write.csv(roe, file = "/mnt/NexData_file/home/lzl/20260515_XE2510316140-03/Macrophage_subtype_filter2/Roe.csv",quote=F)


p = pheatmap::pheatmap(roe, 
                       color = colorRampPalette(brewer.pal(n = 7, name = "Reds"))(10),
                       cluster_rows = F, 
                       cluster_cols = F, 
                       display_numbers = TRUE,
                       number_color = "black",
                       cellwidth = 30, cellheight = 16, fontsize = 10,                    # fontsize字号大小
                       border_color = '#ffffff',
                       angle_col='45')
dev.off()
p$gtable$heights[2] = unit(15,"bigpts")

p$gtable[[2]][p$gtable[[2]]$name =="legend",1:4] = c(4,5,4,6)
p$gtable <- gtable_add_grob(p$gtable, textGrob("Ro/e",x = unit(0.32,"npc"),y = unit(0.7,"npc")), t = 2, l = 5, r= 5, b=2)
#pdf("heatmap.pdf",width = ncol(roe)*1.5,height = nrow(roe)*0.3)
pdf("/mnt/NexData_file/home/lzl/20260515_XE2510316140-03/Macrophage_subtype_filter2/heatmap.pdf",width = ncol(roe)*2,height = nrow(roe)*0.3)
p
dev.off()
png("/mnt/NexData_file/home/lzl/20260515_XE2510316140-03/Macrophage_subtype_filter2/heatmap.png",width = ncol(roe)*200,height = nrow(roe)*30)
p
dev.off()


df_p = df[grepl("P",df$project),]
summary <- table(df_p[,c('subtype','group')])
roe <- as.data.frame(ROIE(summary))
long_df <- tidyr::gather(roe, key = "group", value = "Ro/e", colnames(roe)[1:ncol(roe)-1])
write.csv(roe, file = "/mnt/NexData_file/home/lzl/20260515_XE2510316140-03/Macrophage_subtype_filter2/Roe_P.csv",quote=F)


p = pheatmap::pheatmap(roe, 
                       color = colorRampPalette(brewer.pal(n = 7, name = "Reds"))(10),
                       cluster_rows = F, 
                       cluster_cols = F, 
                       display_numbers = TRUE,
                       number_color = "black",
                       cellwidth = 30, cellheight = 16, fontsize = 10,                    # fontsize字号大小
                       border_color = '#ffffff',
                       angle_col='45')
dev.off()
p$gtable$heights[2] = unit(15,"bigpts")

p$gtable[[2]][p$gtable[[2]]$name =="legend",1:4] = c(4,5,4,6)
p$gtable <- gtable_add_grob(p$gtable, textGrob("Ro/e",x = unit(0.32,"npc"),y = unit(0.7,"npc")), t = 2, l = 5, r= 5, b=2)
#pdf("heatmap.pdf",width = ncol(roe)*1.5,height = nrow(roe)*0.3)
pdf("/mnt/NexData_file/home/lzl/20260515_XE2510316140-03/Macrophage_subtype_filter2/heatmap_P.pdf",width = ncol(roe)*2,height = nrow(roe)*0.3)
p
dev.off()
png("/mnt/NexData_file/home/lzl/20260515_XE2510316140-03/Macrophage_subtype_filter2/heatmap_P.png",width = ncol(roe)*200,height = nrow(roe)*30)
p
dev.off()

###S,P,R
df = Macrophage_meta
summary <- table(df[,c('subtype','group_2')])
roe <- as.data.frame(ROIE(summary))
long_df <- tidyr::gather(roe, key = "group_2", value = "Ro/e", colnames(roe)[1:ncol(roe)-1])
write.csv(roe, file = "/mnt/NexData_file/home/lzl/20260515_XE2510316140-03/Macrophage_subtype_filter2/Roe_S_P_R.csv",quote=F)


p = pheatmap::pheatmap(roe, 
                       color = colorRampPalette(brewer.pal(n = 7, name = "Reds"))(10),
                       cluster_rows = F, 
                       cluster_cols = F, 
                       display_numbers = TRUE,
                       number_color = "black",
                       cellwidth = 30, cellheight = 16, fontsize = 10,                    # fontsize字号大小
                       border_color = '#ffffff',
                       angle_col='45')
p$gtable$heights[2] = unit(15,"bigpts")

p$gtable[[2]][p$gtable[[2]]$name =="legend",1:4] = c(4,5,4,6)
p$gtable <- gtable_add_grob(p$gtable, textGrob("Ro/e",x = unit(0.32,"npc"),y = unit(0.7,"npc")), t = 2, l = 5, r= 5, b=2)
#pdf("heatmap.pdf",width = ncol(roe)*1.5,height = nrow(roe)*0.3)
pdf("/mnt/NexData_file/home/lzl/20260515_XE2510316140-03/Macrophage_subtype_filter2/heatmap_S_P_R.pdf",width = ncol(roe)*2,height = nrow(roe)*0.4)
p
dev.off()

###S,P,R,A,C

df = Macrophage_meta
summary <- table(df[,c('subtype','group_3')])
roe <- as.data.frame(ROIE(summary))
long_df <- tidyr::gather(roe, key = "group_3", value = "Ro/e", colnames(roe)[1:ncol(roe)-1])
write.csv(roe, file = "/mnt/NexData_file/home/lzl/20260515_XE2510316140-03/Macrophage_subtype_filter2/Roe_group_3.csv",quote=F)


p = pheatmap::pheatmap(roe, 
                       color = colorRampPalette(brewer.pal(n = 7, name = "Reds"))(10),
                       cluster_rows = F, 
                       cluster_cols = F, 
                       display_numbers = TRUE,
                       number_color = "black",
                       cellwidth = 30, cellheight = 16, fontsize = 10,                    # fontsize字号大小
                       border_color = '#ffffff',
                       angle_col='45')
p$gtable$heights[2] = unit(15,"bigpts")

p$gtable[[2]][p$gtable[[2]]$name =="legend",1:4] = c(4,5,4,6)
p$gtable <- gtable_add_grob(p$gtable, textGrob("Ro/e",x = unit(0.32,"npc"),y = unit(0.7,"npc")), t = 2, l = 5, r= 5, b=2)
#pdf("heatmap.pdf",width = ncol(roe)*1.5,height = nrow(roe)*0.3)
pdf("/mnt/NexData_file/home/lzl/20260515_XE2510316140-03/Macrophage_subtype_filter2/heatmap_group_3.pdf",width = ncol(roe)*2,height = nrow(roe)*0.4)
p
dev.off()

###############
Idents(Macrophage_clean_4_pure) <- "subtype"

####

Mac_f_markers <-
  FindAllMarkers(
    JoinLayers(Macrophage_clean_4_pure),
    test.use = "wilcox",
    logfc.threshold = 0,
    min.pct = 0,
    only.pos = TRUE
  )
write.csv(Mac_f_markers, file = "Macrophage_subtype_filter2/Mac_f_all_markers subtype.csv", row.names = TRUE)


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
pdf(file.path(out_dir, str_glue("dotplot subtype major_marker.pdf")), width = 25, height = 6)
DotPlot(Macrophage_clean_4_pure, group.by = "subtype", features = major_marker_ls)+RotatedAxis()+theme_90angle
DotPlot(Macrophage_clean_4_pure, group.by = "subtype", features = Macrophage_major_marker_ls)+RotatedAxis()+theme_90angle
DotPlot(Macrophage_clean_4_pure, group.by = "subtype", features = Mac_f_top5_markers)+RotatedAxis()+theme_90angle
dev.off()





#####功能

library(clusterProfiler)
library(org.Hs.eg.db)

# 按 cluster 分组做富集
marker_list <- Mac_f_markers %>%
  filter(p_val_adj < 0.05, avg_log2FC > 0.5) %>%
  split(.$cluster)

enrichGO_results <- imap(marker_list, function(df, cluster_name) {
  message("正在富集：", cluster_name, " (", nrow(df), " genes)")
  
  entrez <- bitr(df$gene, "SYMBOL", "ENTREZID", org.Hs.eg.db)$ENTREZID
  
  if (length(entrez) < 10) return(NULL)
  
  enrichGO(
    gene = entrez,
    OrgDb = org.Hs.eg.db,
    ont = "BP",
    readable = TRUE
  )
}) %>% compact()

names(enrichGO_results)
enrichGO_results
# 可视化
iwalk(enrichGO_results, function(ego, cluster_name) {
  
  # 检查是否有显著结果
  if (is.null(ego) || nrow(as.data.frame(ego)) == 0) {
    message(cluster_name, " 无显著富集结果")
    return(NULL)
  }
  
  # 点图
  p <- dotplot(ego, showCategory = 15) +
    ggtitle(paste0(cluster_name, " - GO BP Enrichment")) +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.text.y = element_text(size = 10)
    )
  
  # 保存
  ggsave(
    filename = file.path(out_dir, paste0("GO_dotplot_", cluster_name, ".pdf")),
    plot = p,
    width = 10,
    height = 8
  )
  
  message("已保存: ", cluster_name)
})
