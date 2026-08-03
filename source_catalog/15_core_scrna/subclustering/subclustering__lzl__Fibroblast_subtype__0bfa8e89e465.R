all_data$celltype %>% table()


library(Seurat)
library(tidyverse)
library(BPCells)
library(clusterProfiler)
library(org.Hs.eg.db)



all_data <- adata
all_data$celltype %>% table()



#all_data <- readRDS(file.path("R_OUTPUT", "all_data.rds"))

theme_90angle <-
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
    strip.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )


# Fibroblasts -------------------------------------------------------------------

Fibroblast_out_dir <- "Fibroblast_subtype" 
dir.create(Fibroblast_out_dir)
Fibroblast <- subset(all_data, annotation == 'Fibroblast')

## PCA --------
Fibroblast <- FindVariableFeatures(Fibroblast)
future::plan("sequential")
options(future.globals.maxSize = 50 * 1024^3) # 把并行全局变量上限改成 50GB（足够用）
Fibroblast <- ScaleData(Fibroblast, vars.to.regress = 'nCount_RNA')
Fibroblast <- RunPCA(Fibroblast)

### pca inspector ------
pdf(file.path(Fibroblast_out_dir, 'pca plot.pdf'), height = 6, width = 8)
ElbowPlot(Fibroblast, ndims = 50)
dev.off()


Fibroblast$annotation %>% table()
Fibroblast$celltype %>% table()
## dimension reduction --------
dims <- 17 # 18

library(future)
#options(future.globals.maxSize = 10 * 1024 ^ 3) # 10Gb
plan(multisession, workers = 10)

Fibroblast <- RunUMAP(Fibroblast, dims = 1:dims,return.model = TRUE)

pdf(file.path(Fibroblast_out_dir, str_glue("Dimplot dim{dims} sample.pdf")), height = 6, width = 8)
DimPlot(Fibroblast, group.by = "sample", label = TRUE) + NoLegend()
dev.off()

Fibroblast <- Fibroblast %>% 
  harmony::RunHarmony(
    dims.use = 1:dims,
    group.by.vars = "project",
    # theta = 4, # default 2, bigger value -- more diversity per cluster
    # sigma = 0.5,
    # lambda = 0.5, # default, smaller value -- more integrated
    reduction.save = "harmony"
  )

Fibroblast <- RunUMAP(
  Fibroblast, 
  dims = 1:dims, 
  reduction = "harmony", 
  reduction.name = "harmonyumap", 
  reduction.key  = "harmonyUMAP_",
  return.model = TRUE
)

pdf(file.path(Fibroblast_out_dir, str_glue("Dimplot dim{dims} sample - harmony.pdf")), height = 6, width = 8)
DimPlot(Fibroblast, group.by = "sample", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()

reduction.use <- 'harmony' # 'harmony' 'pca'


### cluster -------
Fibroblast  <-
  FindNeighbors(
    Fibroblast ,
    reduction = reduction.use,
    dims = 1:dims,
    k.param = 20
  )


Fibroblast <-
  FindClusters(
    Fibroblast,
    resolution = c(0.2, 0.4, 0.6, 0.8)
  )



# dimplot
pdf(file.path(Fibroblast_out_dir, str_glue("Dimplot dim{dims} cluster.pdf")), height = 6, width = 8)
DimPlot(Fibroblast, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Fibroblast, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Fibroblast, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Fibroblast, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(Fibroblast, group.by = "RNA_snn_res.1", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(Fibroblast, group.by = "RNA_snn_res.1.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(Fibroblast, group.by = "RNA_snn_res.1.5", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()

### dotplot -------
Fibro_major_marker_ls <- c(
  "NF" = "IGFBP6", 
  #"NF" = "PI16", 
  # "NF" = "DPT", 
  # "NF" = "CXCL12", 
  "NF" = "GPX3", 
  "NF" = "TNXB", 
  
  "CAF" = "FAP",
  "CAF" = "PDGFRA",
  "CAF" = "PDGFRB",
  "CAF" = "HSF1",
  "CAF" = "STAT3",
  "CAF" = "MYC",
  "CAF" = "CD10",
  "CAF" = "GPR77",
  
  "mCAF" = "MMP11",
  "mCAF" = "MYL9",
  "mCAF" = "POSTN",
  "mCAF" = "TPM2",
  "mCAF" = "TAGLN",
  "mCAF" = "HOPX",
  "mCAF" = "COL1A2",
  
  "iCAF" = "PLA2G2A",
  "iCAF" = "CFD",
  "iCAF" = "CD34",
  "iCAF" = "CXCL14",
  "iCAF" = "PI16",
  "iCAF" = "PRG4",
  "iCAF" = "APOD",
  "iCAF" = "DPT",
  "iCAF" = "CXCL12",
  "iCAF" = "CXCL1",
  "iCAF" = "IL6",
  "iCAF" = "CCL2", 
  
  "vCAF" = "NOTCH3",
  "vCAF" = "COL18A1",
  "vCAF" = "ACTA2",
  "vCAF" = "MCAM",
  "vCAF" = "RGS5",
  "vCAF" = "MYH11",
  "vCAF" = "ADIRF", 
  "vCAF" = "RERGL", 
  
  "tCAF" = "PDPN",
  "tCAF" = "MME",
  "tCAF" = "NDRG1",
  "tCAF" = "TMEM158",
  "tCAF" = "VEGFA",
  "tCAF" = "ENO1",
  "tCAF" = "GAPDH",
  "tCAF" = "CAIX",
  "tCAF" = "HSPH1", 
  "tCAF" = "HSP90AA1",
  
  "ifnCAF" = "IL32",
  "ifnCAF" = "IDO1",
  "ifnCAF" = "CXCL9",
  "ifnCAF" = "CXCL10",
  "ifnCAF" = "CXCL11",
  
  "apCAF" = "HLA-DRA",
  "apCAF" = "HLA-DRB1",
  "apCAF" = "CD74",
  "apCAF" = "IGLC2",
  "apCAF" = "IGKC",
  "apCAF" = "CXCR4",
  
  "rCAF" = "CCL19",
  "rCAF" = "CCL21", 
  
  "dCAF" = "MKI67",
  "dCAF" = "TUBA1B", 
  
  "meCAF" = "SLC2A1",
  "meCAF" = "BNIP3",
  "meCAF" = "ADM",
  "meCAF" = "ANGPTL4"
)

Idents(Fibroblast)<-'RNA_snn_res.0.4'
Fibroblast_markers <-
  FindAllMarkers(
    JoinLayers(Fibroblast),
    test.use = "wilcox",
    logfc.threshold = 0,
    min.pct = 0,
    only.pos = TRUE
  )
write.csv(Fibroblast_markers, file = "Fibroblast_subtype/Fibroblast_all_markers.csv", row.names = TRUE)

Fibroblast_top5_markers <- 
  Fibroblast_markers %>% 
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

names(Fibroblast_top5_markers) <- sub('.*\\.','',names(Fibroblast_top5_markers))


# dotplot
# 如果已做其他细胞类型的亚群分析，这里就可以顺便画一下其他类型的细胞亚群基因，比如Fibro_major_marker_ls2，否则只用画T细胞
pdf(file.path(Fibroblast_out_dir, str_glue("dotplot dims{dims} major_marker.pdf")), width = 25, height = 6)
DotPlot(Fibroblast, group.by = "RNA_snn_res.0.4", features = major_marker_ls)+RotatedAxis()+theme_90angle
#DotPlot(Fibroblast, group.by = "RNA_snn_res.0.4", features = T_major_marker_ls)+RotatedAxis()+theme_90angle
#DotPlot(Fibroblast, group.by = "RNA_snn_res.0.4", features = T_major_marker_ls_1)+RotatedAxis()+theme_90angle
#DotPlot(Fibroblast, group.by = "RNA_snn_res.0.6", features = Fibroblast_markers)+RotatedAxis()
DotPlot(Fibroblast, group.by = "RNA_snn_res.0.4", features = Fibro_major_marker_ls)+RotatedAxis()+theme_90angle
DotPlot(Fibroblast, group.by = "RNA_snn_res.0.4", features = Fibroblast_top5_markers)+RotatedAxis()+theme_90angle

dev.off()




saveRDS(Fibroblast ,file = "Fibroblast_subtype/Fibroblast_cluster.rds")
Fibroblast <- readRDS("Fibroblast_subtype/Fibroblast_cluster.rds")

Fibroblast$RNA_snn_res.0.4 %>% table()
## annotation ----------------------------------------------------------------
cell_cluster <- c(
  '0'  = "NF", # 
  '1'  = "Macrophage-like_Fib", #
  '2'  = "Fib_others", #
  '3'  = "iCAF", # 
  '4'  = "mCAF", # 
  '5'  = "T_cell-like_Fib", # 
  '6'  = "Plasma-like_Fib", # 
  '7'  = "rCAF", # 
  '8'  = "Endothelial",# 
  '9'  = "ifnCAF", # 
  '10' = "Epithelial-like_Fib", # 
  '11' = "Immunity_cell-like_Fib", # 
  '12' = "Immunity_cell-like_Fib"
)

Idents(Fibroblast) <- "RNA_snn_res.0.4"

# 基础版小提琴图
png(file.path(Fibroblast_out_dir, "Fibroblast_nCount_RNA_vln.png"), width = 1200, height = 800, res = 150)
VlnPlot(object = Fibroblast,
        features = c("nCount_RNA", "nFeature_RNA"),  # 你要画的指标
        group.by = "RNA_snn_res.0.4",  # 按cluster分组
        pt.size = 0,  # 去掉散点，更干净
        fill = "RNA_snn_res.0.4",
        raster = FALSE  # 每个cluster自动配色
)
dev.off()

getwd()
Fibroblast <- RenameIdents(Fibroblast, cell_cluster)
Fibroblast$subtype <- Idents(Fibroblast)


pdf(file.path(Fibroblast_out_dir, str_glue("Dimplot dim{dims} annotation.pdf")), height = 6, width = 8)
DimPlot(Fibroblast, group.by = "subtype", label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

pdf(file.path(Fibroblast_out_dir, str_glue("Dimplot dim{dims} annotation split.pdf")), height = 6, width = 15)
DimPlot(Fibroblast, group.by = "subtype",split.by = 'group', label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()


#all_data <- readRDS("output/data/all_data.rds")
# 1. 先把 all_data$annotation 复制到新列 T_annotation
all_data$annotation2 <- as.character(all_data$T_annotation)

# 2. 创建映射：细胞名 → 精细T亚型
cell_subtype <- Fibroblast$subtype
names(cell_subtype) <- colnames(Fibroblast)

# 3. 只把 "T cell" 替换成精细亚型，其他细胞保持原样
all_data$annotation2[colnames(Fibroblast)] <- as.character(cell_subtype)

# 查看结果（正确）
table(all_data$T_annotation)
table(all_data$annotation2)
table(Fibroblast$subtype)
saveRDS(all_data, "Fibroblast_subtype/all_data_Fibroblast_annotation.rds")

saveRDS(Fibroblast, "Fibroblast_subtype/Fibroblast_annotation.rds")

keep <- c(
  "NF",
  "iCAF",
  "mCAF",
  "rCAF",
  "ifnCAF"
)

Fibroblast_clean <- subset(Fibroblast, subtype %in% keep)

Fibroblast_clean$subtype <- droplevels(Fibroblast_clean$subtype)
Idents(Fibroblast_clean) <- Fibroblast_clean$subtype

table(Fibroblast_clean$subtype)
saveRDS(Fibroblast_clean, "Fibroblast_subtype/Fibroblast_clean.rds")
#saveRDS(Fibroblast, "Fibroblast_subtype/Fibroblast_clean.rds")
Fibroblast <- Fibroblast_clean


Fibroblast <- RunUMAP(Fibroblast, dims = 1:dims,n.neighbors = 30,return.model = TRUE)

pdf(file.path(Fibroblast_out_dir, str_glue("Dimplot dim{dims} Fibroblast annotation.pdf")), height = 6, width = 8)
DimPlot(Fibroblast, group.by = "subtype", label = TRUE, reduction = 'umap',repel = TRUE) + NoLegend()
dev.off()

pdf(file.path(Fibroblast_out_dir, str_glue("Dimplot dim{dims} Fibroblast annotation split.pdf")), height = 6, width = 15)
DimPlot(Fibroblast, group.by = "subtype",split.by = 'group', label = TRUE, reduction = 'umap',repel = TRUE) + NoLegend()
dev.off()



#图形效果不好，需要重新过滤一下
Fibroblast <- Fibroblast_clean
Fibroblast <- readRDS("Fibroblast_subtype/Fibroblast_clean.rds")

## PCA --------
Fibroblast <- FindVariableFeatures(Fibroblast)

Fibroblast <- ScaleData(Fibroblast, vars.to.regress = 'nCount_RNA')
Fibroblast <- RunPCA(Fibroblast)

### pca inspector ------
pdf(file.path(Fibroblast_out_dir, 'pca plot_1.pdf'), height = 6, width = 8)
ElbowPlot(Fibroblast, ndims = 50)
dev.off()

## dimension reduction --------
dims <- 16 # 18

library(future)
#options(future.globals.maxSize = 10 * 1024 ^ 3) # 10Gb
plan(multisession, workers = 10)

Fibroblast <- RunUMAP(Fibroblast, dims = 1:dims,return.model = TRUE)


Fibroblast <- Fibroblast %>% 
  harmony::RunHarmony(
    dims.use = 1:dims,
    group.by.vars = "project",
    # theta = 4, # default 2, bigger value -- more diversity per cluster
    # sigma = 0.5,
    # lambda = 0.5, # default, smaller value -- more integrated
    reduction.save = "harmony"
  )

Fibroblast <- RunUMAP(
  Fibroblast, 
  dims = 1:dims, 
  reduction = "harmony", 
  reduction.name = "harmonyumap", 
  reduction.key  = "harmonyUMAP_",
  return.model = TRUE
)



reduction.use <- 'harmony' # 'harmony' 'pca'


### cluster -------
Fibroblast  <-
  FindNeighbors(
    Fibroblast ,
    reduction = reduction.use,
    dims = 1:dims,
    k.param = 20
  )


Fibroblast <-
  FindClusters(
    Fibroblast,
    resolution = c(0.2, 0.4, 0.6, 0.8)
  )



# dimplot
pdf(file.path(Fibroblast_out_dir, str_glue("Dimplot dim{dims} cluster_1.pdf")), height = 6, width = 8)
DimPlot(Fibroblast, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Fibroblast, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Fibroblast, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Fibroblast, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(Fibroblast, group.by = "RNA_snn_res.1", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(Fibroblast, group.by = "RNA_snn_res.1.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(Fibroblast, group.by = "RNA_snn_res.1.5", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()


Idents(Fibroblast)<-'RNA_snn_res.0.8'
Fibroblast_markers <-
  FindAllMarkers(
    JoinLayers(Fibroblast),
    test.use = "wilcox",
    logfc.threshold = 0,
    min.pct = 0,
    only.pos = TRUE
  )
write.csv(Fibroblast_markers, file = "Fibroblast_subtype/Fibroblast_all_markers_1.csv", row.names = TRUE)

Fibroblast_top5_markers <- 
  Fibroblast_markers %>% 
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

names(Fibroblast_top5_markers) <- sub('.*\\.','',names(Fibroblast_top5_markers))


pdf(file.path(Fibroblast_out_dir, str_glue("dotplot dims{dims} major_marker_2.pdf")), width = 25, height = 6)
DotPlot(Fibroblast, group.by = "RNA_snn_res.0.8", features = major_marker_ls)+RotatedAxis()+theme_90angle
#DotPlot(Fibroblast, group.by = "RNA_snn_res.0.8", features = T_major_marker_ls)+RotatedAxis()+theme_90angle
#DotPlot(Fibroblast, group.by = "RNA_snn_res.0.8", features = T_major_marker_ls_1)+RotatedAxis()+theme_90angle
DotPlot(Fibroblast, group.by = "RNA_snn_res.0.8", features = Fibro_major_marker_ls)+RotatedAxis()
DotPlot(Fibroblast, group.by = "RNA_snn_res.0.8", features = Fibroblast_top5_markers)+RotatedAxis()+theme_90angle
#DotPlot(Fibroblast, group.by = "RNA_snn_res.0.6", features = Fibro_major_marker_ls2)+RotatedAxis()
dev.off()



saveRDS(Fibroblast ,file = "Fibroblast_subtype/Fibroblast_cluster_1.rds")
Fibroblast <- readRDS("Fibroblast_subtype/Fibroblast_cluster_1.rds")
## annotation ----------------------------------------------------------------
cell_cluster <- c(
  '0'  = "NF", # 
  '1'  = "iCAF", #
  '2'  = "mCAF", # 
  '3'  = "tCAF", # CD8+TEM
  '4'  = "Fib_others", # 
  '5'  = "Fib_others",#  
  '6'  = "mCAF", #  
  '7'  = "Immunity_cell-like_Fib", #  
  '8'  = "ifnCAF",# 
  '9'  = "Fib_others", #  
  '10' = "Immunity_cell-like_Fib", # 
  '11' = "Immunity_cell-like_Fib"
)

Idents(Fibroblast) <- "RNA_snn_res.0.8"




Fibroblast <- RenameIdents(Fibroblast, cell_cluster)
Fibroblast$subtype <- Idents(Fibroblast)


pdf(file.path(Fibroblast_out_dir, str_glue("Dimplot dim{dims} annotation_1.pdf")), height = 6, width = 8)
DimPlot(Fibroblast, group.by = "subtype", label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

pdf(file.path(Fibroblast_out_dir, str_glue("Dimplot dim{dims} annotation split_1.pdf")), height = 6, width = 15)
DimPlot(Fibroblast, group.by = "subtype",split.by = 'group', label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

keep <- c(
  "NF",
  "iCAF",
  "mCAF",
  "tCAF",
  "ifnCAF"
)

Fibroblast_clean <- subset(Fibroblast, subtype %in% keep)

Fibroblast_clean$subtype <- droplevels(Fibroblast_clean$subtype)
Idents(Fibroblast_clean) <- Fibroblast_clean$subtype

table(Fibroblast_clean$subtype)
saveRDS(Fibroblast_clean, "Fibroblast_subtype/Fibroblast_clean_1.rds")

Fibroblast <- Fibroblast_clean

Fibroblast <- RunUMAP(Fibroblast, dims = 1:16,n.neighbors = 50,return.model = TRUE)

pdf(file.path(Fibroblast_out_dir, str_glue("Dimplot dim{dims} Fibroblast annotation_1.pdf")), height = 6, width = 8)
DimPlot(Fibroblast, group.by = "subtype", label = TRUE, reduction = 'umap',repel = TRUE) + NoLegend()
dev.off()

pdf(file.path(Fibroblast_out_dir, str_glue("Dimplot dim{dims} Fibroblast annotation split_1.pdf")), height = 6, width = 15)
DimPlot(Fibroblast, group.by = "subtype",split.by = 'group', label = TRUE, reduction = 'umap',repel = TRUE) + NoLegend()
dev.off()
#运行到这里




#画一个注释完的点图
Fibro_major_marker_ls_2 <- c(
  "NF" = "TNXB", 
  
  "CAF" = "FAP",
  "CAF" = "PDGFRA",
  "CAF" = "PDGFRB",
  
  "mCAF" = "MMP11",
  "mCAF" = "POSTN",
  "mCAF" = "HOPX",
  
  #"iCAF" = "PLA2G2A",
  #"iCAF" = "CFD",
  "iCAF" = "CD34",
  #"iCAF" = "CXCL14",
  "iCAF" = "PI16",
  #"iCAF" = "PRG4",
  #"iCAF" = "APOD",
  "iCAF" = "DPT",
  "iCAF" = "CXCL12",
  #"iCAF" = "CXCL1",
  #"iCAF" = "IL6",
  #"iCAF" = "CCL2", 
  "tCAF" = "PDPN",
  "tCAF" = "VEGFA",
  

  
  #"ifnCAF" = "IL32",
  #"ifnCAF" = "IDO1",
  "ifnCAF" = "CXCL9",
  "ifnCAF" = "CXCL10"
)

pdf(file.path(Fibroblast_out_dir, str_glue("dotplot Fibroblast subtype major_marker.pdf")), width = 25, height = 6)
DotPlot(Fibroblast, group.by = "subtype", features = Fibro_major_marker_ls_2)+RotatedAxis()+theme_90angle
dev.off()



## cell props barplot  ----------------------------------------------------------
Fibroblast_meta <- Fibroblast_clean[[]]

Fibroblast_meta$sample <- factor(Fibroblast_meta$sample,levels=c(Fibroblast_meta$sample %>% unique()))

p_col_proj <-
  Fibroblast_meta %>%
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
  filename = file.path("Fibroblast_subtype", "cell_freq_sample.pdf"),
  height = 6,
  width  = 8
)

p_col_proj <-
  Fibroblast_meta %>%
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
  filename = file.path(Fibroblast_out_dir, "cell_freq_group.pdf"),
  height = 6,
  width  = 4
)


######################
df = Fibroblast_meta
#df = read.csv("/lvdata/wzb/scRNA/FWSC20240734/ReAll/celltype/meta.xls",sep = "\t",row.names=1)
summary <- table(df[,c('subtype','group')])
roe <- as.data.frame(ROIE(summary))
long_df <- tidyr::gather(roe, key = "group", value = "Ro/e", colnames(roe)[1:ncol(roe)-1])
write.csv(roe, file = "/mnt/NexData_file/home/lzl/20260515_XE2510316140-03/Fibroblast_subtype/Roe.csv",quote=F)


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
pdf("/mnt/NexData_file/home/lzl/20260515_XE2510316140-03/Fibroblast_subtype/heatmap.pdf",width = 2.5,height = 2)
p
dev.off()
png("/mnt/NexData_file/home/lzl/20260515_XE2510316140-03/Fibroblast_subtype/heatmap.png",width = 250,height = 200)
p
dev.off()


df_p = df[grepl("P",df$project),]
summary <- table(df_p[,c('subtype','group')])
roe <- as.data.frame(ROIE(summary))
long_df <- tidyr::gather(roe, key = "group", value = "Ro/e", colnames(roe)[1:ncol(roe)-1])
write.csv(roe, file = "/mnt/NexData_file/home/lzl/20260515_XE2510316140-03/Fibroblast_subtype/Roe_P.csv",quote=F)


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
pdf("/mnt/NexData_file/home/lzl/20260515_XE2510316140-03/Fibroblast_subtype/heatmap_P.pdf",width = 2.5,height = 2)
p
dev.off()
png("/mnt/NexData_file/home/lzl/20260515_XE2510316140-03/Fibroblast_subtype/heatmap_P.png",width = 250,height = 200)
p
dev.off()
