library(Seurat)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(BPCells)
library(clusterProfiler)
library(org.Hs.eg.db)


all_data <- readRDS('/mnt/NexData_file/home/lxy/case/XE2603262474-01/stRNA_2/output/data/all_data.rds')

# Macrophages -------------------------------------------------------------------
#all_data$annotation <- str_replace(all_data$annotation, 'Proliferating T cell', 'T cell')
all_data$annotation %>% table()
Macrophage_out_dir <- "Macrophage_subtype" 
dir.create(Macrophage_out_dir)
Macrophage <- subset(all_data, annotation == 'Macrophage')

all_data$annotation %>% table()
all_data$celltype %>% table()
## PCA --------
Macrophage <- FindVariableFeatures(Macrophage)

Macrophage <- ScaleData(Macrophage, vars.to.regress = 'nCount_RNA')
Macrophage <- RunPCA(Macrophage)

### pca inspector ------
pdf(file.path(Macrophage_out_dir, 'pca plot.pdf'), height = 6, width = 8)
ElbowPlot(Macrophage, ndims = 50)
dev.off()


## dimension reduction --------
dims <- 18 # 18

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

Idents(Macrophage) <- "RNA_snn_res.0.4"

### dotplot -------
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


Idents(Macrophage)<-'RNA_snn_res.0.4'
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
DotPlot(Macrophage, group.by = "RNA_snn_res.0.4", features = major_marker_ls)+RotatedAxis()+theme_90angle
DotPlot(Macrophage, group.by = "RNA_snn_res.0.4", features = Macrophage_major_marker_ls)+RotatedAxis()+theme_90angle
#DotPlot(Macrophage, group.by = "RNA_snn_res.0.6", features = Fibro_major_marker_ls_1)+RotatedAxis()+theme_90angle
#DotPlot(Macrophage, group.by = "RNA_snn_res.0.6", features = Fibro_major_marker_ls_2)+RotatedAxis()+theme_90angle
#DotPlot(Macrophage, group.by = "RNA_snn_res.0.6", features = Macrophage_markers)+RotatedAxis()
DotPlot(Macrophage, group.by = "RNA_snn_res.0.4", features = Macrophage_top5_markers)+RotatedAxis()+theme_90angle
#DotPlot(Macrophage, group.by = "RNA_snn_res.0.6", features = Fibro_major_marker_ls2)+RotatedAxis()
dev.off()


saveRDS(Macrophage ,file = "Macrophage_subtype/Macrophage_cluster.rds")
Macrophage <- readRDS("Macrophage_subtype/Macrophage_cluster.rds")
## annotation ----------------------------------------------------------------
cell_cluster <- c(
  '0'  = "M2", # 
  '1'  = "Fibroblsat-like_Macrophage", #
  '2'  = "cDC2", # 
  '3'  = "LOH-like_Macrophage", # 
  '4'  = "PTS2-like_Macrophage", # 
  '5'  = "T_B_cell-like_Macrophage",#  
  '6'  = "Endothelial-like_Macrophage", #  
  '7'  = "M1", #  
  '8'  = "Endo_Neutrophil-like_Macrophage", #  
  '9'  = "cDC1", #  
  '10'  = "Mast cell", #  
  '11'  = "Podocyte", #  
  '12'  = "Neutrophil", #  
  '13'  = "Proliferating T cell", #  
  '14'  = "LAMP3_DC"
)

Idents(Macrophage) <- "RNA_snn_res.0.4"




Macrophage <- RenameIdents(Macrophage, cell_cluster)
Macrophage$subtype <- Idents(Macrophage)


pdf(file.path(Macrophage_out_dir, str_glue("Dimplot dim{dims} annotation.pdf")), height = 6, width = 8)
DimPlot(Macrophage, group.by = "subtype", label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

pdf(file.path(Macrophage_out_dir, str_glue("Dimplot dim{dims} annotation split.pdf")), height = 6, width = 15)
DimPlot(Macrophage, group.by = "subtype",split.by = 'group', label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()


keep <- c(
  "M1",
  "M2",
  "cDC1",
  "cDC2",
  "LAMP3_DC"
)

Macrophage_clean <- subset(Macrophage, subtype %in% keep)

Macrophage_clean$subtype <- droplevels(Macrophage_clean$subtype)
Idents(Macrophage_clean) <- Macrophage_clean$subtype

table(Macrophage_clean$subtype)
saveRDS(Macrophage_clean, "Macrophage_subtype/Macrophage_clean.rds")

Macrophage <- Macrophage_clean


#Macrophage <- readRDS("Macrophage_subtype/Macrophage_clean.rds")
Macrophage <- RunUMAP(Macrophage, dims = 1:18,n.neighbors = 100,return.model = TRUE)

pdf(file.path(Macrophage_out_dir, str_glue("Dimplot dim{dims} Macrophage annotation.pdf")), height = 6, width = 8)
DimPlot(Macrophage, group.by = "subtype", label = TRUE, reduction = 'umap',repel = TRUE) + NoLegend()
dev.off()

pdf(file.path(Macrophage_out_dir, str_glue("Dimplot dim{dims} Macrophage annotation split.pdf")), height = 6, width = 15)
DimPlot(Macrophage, group.by = "subtype",split.by = 'group', label = TRUE, reduction = 'umap',repel = TRUE) + NoLegend()
dev.off()




# 保留：0,1,7
Macrophage <- subset(Macrophage, RNA_snn_res.0.4 %in% c(0,2,7,9,14))
Macrophage$RNA_snn_res.0.4 %>% table()


## PCA --------
Macrophage <- FindVariableFeatures(Macrophage)

Macrophage <- ScaleData(Macrophage, vars.to.regress = 'nCount_RNA')
Macrophage <- RunPCA(Macrophage)

### pca inspector ------
pdf(file.path(Macrophage_out_dir, 'pca plot_1.pdf'), height = 6, width = 8)
ElbowPlot(Macrophage, ndims = 50)
dev.off()

## dimension reduction --------
dims <- 15 # 18

Macrophage <- RunUMAP(Macrophage, dims = 1:dims,return.model = TRUE)


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
pdf(file.path(Macrophage_out_dir, str_glue("Dimplot dim{dims} cluster_1.pdf")), height = 6, width = 8)
DimPlot(Macrophage, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Macrophage, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Macrophage, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Macrophage, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(Macrophage, group.by = "RNA_snn_res.1", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(Macrophage, group.by = "RNA_snn_res.1.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(Macrophage, group.by = "RNA_snn_res.1.5", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()


Idents(Macrophage)<-'RNA_snn_res.0.4'
Macrophage_markers <-
  FindAllMarkers(
    JoinLayers(Macrophage),
    test.use = "wilcox",
    logfc.threshold = 0,
    min.pct = 0,
    only.pos = TRUE
  )
write.csv(Macrophage_markers, file = "Macrophage_subtype/Macrophage_all_markers_1.csv", row.names = TRUE)

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

#Idents(Macrophage) <- "RNA_snn_res.0.4"
pdf(file.path(Macrophage_out_dir, str_glue("dotplot dims{dims} major_marker_1.pdf")), width = 25, height = 6)
DotPlot(Macrophage, group.by = "RNA_snn_res.0.4", features = major_marker_ls)+RotatedAxis()+theme_90angle
DotPlot(Macrophage, group.by = "RNA_snn_res.0.4", features = Macrophage_major_marker_ls)+RotatedAxis()+theme_90angle
DotPlot(Macrophage, group.by = "RNA_snn_res.0.4", features = Macrophage_top5_markers)+RotatedAxis()+theme_90angle
#DotPlot(Macrophage, group.by = "RNA_snn_res.0.6", features = Fibro_major_marker_ls2)+RotatedAxis()
dev.off()


saveRDS(Macrophage ,file = "Macrophage_subtype/Macrophage_cluster_1.rds")
Macrophage <- readRDS("Macrophage_subtype/Macrophage_cluster_1.rds")
## annotation ----------------------------------------------------------------
cell_cluster <- c(
  '0'  = "M2", # 
  '1'  = "cDC2", #
  '2'  = "M2", #
  '3'  = "M1", #
  '4'  = "M2", #
  '5'  = "cDC1", #
  '6'  = "LAMP3_DC", #
  '7'  = "pDC"
)


Idents(Macrophage) <- "RNA_snn_res.0.4"




Macrophage <- RenameIdents(Macrophage, cell_cluster)
Macrophage$subtype <- Idents(Macrophage)


pdf(file.path(Macrophage_out_dir, str_glue("Dimplot dim{dims} annotation_1.pdf")), height = 6, width = 8)
DimPlot(Macrophage, group.by = "subtype", label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

pdf(file.path(Macrophage_out_dir, str_glue("Dimplot dim{dims} annotation split_1.pdf")), height = 6, width = 15)
DimPlot(Macrophage, group.by = "subtype",split.by = 'group', label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

table(Macrophage$subtype)
saveRDS(Macrophage, "Macrophage_subtype/Macrophage_clean_1.rds")

#keep <- c(
#  "CD4+",
#  "CD8+TEM",
#  "CD8+Effector",
#  "Treg",
#  "Proliferating T cell"
#)

#Macrophage_clean <- subset(Macrophage, subtype %in% keep)

#Macrophage_clean$subtype <- droplevels(Macrophage_clean$subtype)
#Idents(Macrophage_clean) <- Macrophage_clean$subtype

#table(Macrophage_clean$subtype)
#saveRDS(Macrophage_clean, "Macrophage_subtype/Macrophage_clean.rds")

#Macrophage <- Macrophage_clean



#Macrophage <- RunUMAP(Macrophage, dims = 1:16,n.neighbors = 50,return.model = TRUE)

#pdf(file.path(Macrophage_out_dir, str_glue("Dimplot dim{dims} T annotation.pdf")), height = 6, width = 8)
#DimPlot(Macrophage, group.by = "subtype", label = TRUE, reduction = 'umap',repel = TRUE) + NoLegend()
#dev.off()

#pdf(file.path(Macrophage_out_dir, str_glue("Dimplot dim{dims} T annotation split.pdf")), height = 6, width = 15)
#DimPlot(Macrophage, group.by = "subtype",split.by = 'group', label = TRUE, reduction = 'umap',repel = TRUE) + NoLegend()
#dev.off()



Macrophage <- readRDS("Macrophage_subtype/Macrophage_clean_1.rds")
Macrophage_1 <- readRDS("Macrophage_subtype/Macrophage_cluster.rds")
cell_cluster <- c(
  '0'  = "M2", # 
  '1'  = "Fibroblsat-like_Macrophage", #
  '2'  = "cDC2", # 
  '3'  = "LOH-like_Macrophage", # 
  '4'  = "PTS2-like_Macrophage", # 
  '5'  = "T_B_cell-like_Macrophage",#  
  '6'  = "Endothelial-like_Macrophage", #  
  '7'  = "M1", #  
  '8'  = "Endo_Neutrophil-like_Macrophage", #  
  '9'  = "cDC1", #  
  '10'  = "Mast cell", #  
  '11'  = "Podocyte", #  
  '12'  = "Neutrophil", #  
  '13'  = "Proliferating T cell", #  
  '14'  = "LAMP3_DC"
)


Idents(Macrophage_1) <- "RNA_snn_res.0.4"
Macrophage_1 <- RenameIdents(Macrophage_1, cell_cluster)
Macrophage_1$subtype <- Idents(Macrophage_1)

Macrophage_1$subtype2 <- as.character(Macrophage_1$subtype)
cell_subtype <- Macrophage$subtype
names(cell_subtype) <- colnames(Macrophage)
Macrophage_1$subtype2[colnames(Macrophage)] <- as.character(cell_subtype)
table(Macrophage_1$subtype2)
table(Macrophage_1$subtype)
table(Macrophage$subtype)
saveRDS(Macrophage_1, "Macrophage_subtype/Macrophage_allannotation.rds")
######下面合并到alldata的注释先不运行
all_data$annotation %>% table()
#all_data <- readRDS("output/data/all_data.rds")
# 1. 先把 all_data$annotation 复制到新列 T_annotation
all_data$annotation1 <- as.character(all_data$annotation)

# 2. 创建映射：细胞名 → 精细T亚型
cell_subtype <- Macrophage_1$subtype2
names(cell_subtype) <- colnames(Macrophage_1)

# 3. 只把 "T cell" 替换成精细亚型，其他细胞保持原样
all_data$annotation1[colnames(Macrophage_1)] <- as.character(cell_subtype)

# 查看结果（正确）
table(all_data$annotation1)
table(all_data$annotation)
table(Macrophage_1$subtype2)
saveRDS(all_data, "Macrophage_subtype/all_data_T_annotation.rds")
saveRDS(Macrophage_1 ,file = "Macrophage_subtype/Macrophage_annotation.rds")



Macrophage <- readRDS("Macrophage_subtype/Macrophage_clean_1.rds")
#画一个注释完的点图
Macrophage_major_marker_ls_3 <- c(
  "Macrophage"= "CD68",
  
  #"M1"= "FCGR1A",
  "M1"= "CXCL9",
  "M1"= "CXCL10",
  "M1"= "CXCL11",
  "M1"= "IL1B",
  "M1"= "IRF1",
  "M1"= "CD40",
  
  "M2"= "CD163",
  "M2"= "CSF1R",
  "M2"= "MRC1",
  "M2"= "CD209",
  "M2"= "FOLR2",　
  
  "pDC" = "IRF4", 
  "pDC" = "IRF7", 
  "pDC" = "IRF8", 
  
  "cDC1" = "XCR1",
  "cDC1" = "CLEC9A",
  
  "cDC2" = "CD1C", 
  "cDC2" = "FCER1A", 
  
  "LAMP3_DC" = "LAMP3"
)
pdf(file.path(Macrophage_out_dir, str_glue("dotplot Macrophage subtype major_marker.pdf")), width = 13, height = 6)
DotPlot(Macrophage, group.by = "subtype", features = Macrophage_major_marker_ls_3)+RotatedAxis()+theme_90angle
dev.off()


Macrophage$subtype %>% table()
Macrophage$sample %>% table()
Macrophage$group %>% table()
## cell props barplot  ----------------------------------------------------------
Macrophage_meta <- Macrophage[[]]

Macrophage_meta$sample <- factor(Macrophage_meta$sample,levels=c(Macrophage_meta$sample %>% unique()))

p_col_proj <-
  Macrophage_meta %>%
  group_by(sample, subtype) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = sample, y = freq, fill = subtype)) + 
  geom_col() +
  scale_fill_brewer(palette = "Set2", name = "Macrophage Subtype") +
  #ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )


# export
ggsave(
  p_col_proj,
  filename = file.path("Macrophage_subtype", "cell_freq_sample.pdf"),
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
  filename = file.path(Macrophage_out_dir, "cell_freq_group.pdf"),
  height = 6,
  width  = 4
)

#改图
# 提取meta
Macrophage_meta <- Macrophage[[]]

group_map <- c(
  "N"       = "Normal",
  "M_igAN"  = "Mild",
  "S_igAN"  = "Severe"
)
# 新分组显示顺序（左→右：Normal, Mild, Severe）
group_levels <- c("Normal", "Mild", "Severe")
# 对应配色（保持你指定色值不变）
group_col <- c(
  "Normal"  = "#2E86AB",
  "Mild"    = "#A23B72",
  "Severe"  = "#F18F01"
)


# 确保meta里的group是因子，顺序和颜色对应
Macrophage_meta$Group <- recode(Macrophage_meta$group, !!!group_map)
Macrophage_meta$Group <- factor(Macrophage_meta$Group, levels = group_levels)


# ========== 2. 计算样本水平百分比（用于箱线图 & 均值堆叠图） ==========
# 每个样本内各细胞亚群百分比
data_long <- Macrophage_meta %>%
  group_by(sample, Group, subtype) %>%
  summarise(Count = n(), .groups = "drop_last") %>%
  group_by(sample) %>%
  mutate(Percentage = Count / sum(Count) * 100) %>%
  ungroup()

# 保存原始百分比数据
write.csv(data_long, 
          file.path("Macrophage_subtype", "Macrophage_sample_percentage.csv"), 
          row.names = FALSE)

# 按Group + CellType 计算均值（用于堆叠柱状图）
category_avg <- data_long %>%
  group_by(Group, subtype) %>%
  summarise(MeanPercentage = mean(Percentage), .groups = "drop")

# 保存均值统计数据
write.csv(category_avg, 
          file.path("Macrophage_subtype", "Macrophage_group_mean_percentage.csv"), 
          row.names = FALSE)

# ========== 堆叠柱状图 ==========
p_category_stack <- ggplot(category_avg, 
                           aes(x = Group, y = MeanPercentage, fill = subtype)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  scale_fill_brewer(palette = "Set2", name = "Macrophage Subtype") +
  labs(title = "Macrophage Subtype Composition by Group",
       x = "Group", y = "Percentage (%)") +
  coord_cartesian(ylim = c(0, 100)) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 10, face = "bold"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
  )

ggsave(file.path("Macrophage_subtype", "Macrophage_StackedBar.pdf"), 
       p_category_stack, width = 8, height = 6, dpi = 300)
print(p_category_stack)


# ==================== 按细胞类型分组的箱线图 ====================
comp_list <- list(c("Normal","Mild"), c("Normal","Severe"), c("Mild","Severe"))
p <- ggplot(data_long, aes(x = Group, y = Percentage, fill = Group)) +
  geom_boxplot(width = 0.7, alpha = 0.8, outlier.shape = 19, outlier.size = 1.5) +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.7, aes(color = Group)) +
  facet_wrap(~subtype, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = group_col) +
  scale_color_manual(values = group_col) +
  labs(title = "Macrophage - Cell Type Percentages",
       x = "Group", y = "Percentage (%)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 10, face = "bold"),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5))+
  stat_compare_means(
    comparisons = comp_list,
    method = "wilcox.test",
    label = "p.format",
    size = 3
  )

# 保存+展示
ggsave(file.path("Macrophage_subtype","Macrophage_Boxplot.pdf"), p, width = 12, height = 8, dpi = 300)
print(p)