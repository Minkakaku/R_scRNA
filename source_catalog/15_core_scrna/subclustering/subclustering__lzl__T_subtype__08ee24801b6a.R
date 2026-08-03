library(Seurat)
library(tidyverse)
library(future)
# set default ---------
out_dir <- 'T_subtype'
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
table(all_data$annotation)
#####
Tcell_data <- subset(all_data,annotation %in% c("T cell"))
Tcell_data <- JoinLayers(Tcell_data)
## PCA --------
Tcell_data <- FindVariableFeatures(Tcell_data)
length(VariableFeatures(Tcell_data)) #2000
Tcell_data <- ScaleData(Tcell_data)
Tcell_data <- RunPCA(Tcell_data)

### pca inspector ------
pdf(file.path(out_dir, 'pca plot.pdf'), height = 6, width = 8)
ElbowPlot(Tcell_data, ndims = 50)
dev.off()

pdf(file.path(out_dir, 'dimheatmap 1-5.pdf'), height = 6, width = 10)
DimHeatmap(Tcell_data, dims = 1:5)
dev.off()
pdf(file.path(out_dir, 'dimheatmap 6-10.pdf'), height = 6, width = 10)
DimHeatmap(Tcell_data, dims = 6:10)
dev.off()
pdf(file.path(out_dir, 'dimheatmap 11-15.pdf'), height = 6, width = 10)
DimHeatmap(Tcell_data, dims = 11:15)
dev.off()
pdf(file.path(out_dir, 'dimheatmap 16-20.pdf'), height = 6, width = 10)
DimHeatmap(Tcell_data, dims = 16:20)
dev.off()

## dimension reduction --------
dims <- 15# 18

library(future)
options(future.globals.maxSize = 10 * 1024 ^ 3) # 10Gb
plan(multisession, workers = 10)

Tcell_data <- RunUMAP(Tcell_data, dims = 1:dims,return.model = TRUE)

pdf(file.path(out_dir, str_glue("Dimplot dim{dims} sample.pdf")), height = 8, width = 8)
DimPlot(Tcell_data, group.by = "sample", label = TRUE) + NoLegend()
dev.off()

Tcell_data <- Tcell_data %>% 
  harmony::RunHarmony(
    dims.use = 1:dims,
    group.by.vars = "project",
    # theta = 4, # default 2, bigger value -- more diversity per cluster
    # sigma = 0.5,
    # lambda = 0.5, # default, smaller value -- more integrated
    reduction.save = "harmony"
  )

Tcell_data <- RunUMAP(
  Tcell_data, 
  dims = 1:dims, 
  reduction = "harmony", 
  reduction.name = "harmonyumap", 
  reduction.key  = "harmonyUMAP_",
  return.model = TRUE
)

pdf(file.path(out_dir, str_glue("Dimplot dim{dims} sample - harmony.pdf")), height = 8, width = 8)
DimPlot(Tcell_data, group.by = "sample", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()

reduction.use <- 'harmony' # 'harmony' 'pca'


### cluster -------
Tcell_data  <-
  FindNeighbors(
    Tcell_data ,
    reduction = reduction.use,
    dims = 1:dims,
    k.param = 15  #20改15
  )


Tcell_data <-
  FindClusters(
    Tcell_data,
    resolution = c(0.1,0.2, 0.4, 0.6, 0.8)
  )



# dimplot
pdf(file.path(out_dir, str_glue("Dimplot dim{dims} cluster.pdf")), height = 6, width = 6)
DimPlot(Tcell_data, group.by = "RNA_snn_res.0.1", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Tcell_data, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Tcell_data, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Tcell_data, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Tcell_data, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(Macrophage, group.by = "RNA_snn_res.1", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(Macrophage, group.by = "RNA_snn_res.1.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(Macrophage, group.by = "RNA_snn_res.1.5", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()

Idents(Tcell_data) <- "RNA_snn_res.0.2"

#####
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
###marker_ls
Tcell_major_marker_ls <- c(
  # ==================== 基础大类 (Base Lineages) ====================
  # T cells (Pan-T markers)
  "T_cell" = "Cd3d",
  "T_cell" = "Cd3e",
  "T_cell" = "Cd3g",
  "T_cell" = "Trac",
  "T_cell" = "Trgc",
  "T_cell" = "Trbc2",
  
  # Exhausted T cells
  "Exhausted_T_cell" = "Lag3",
  
  # CD4+ T cell (Pan-CD4 markers)
  "Cd4_T_cell" = "Cd4",
  "Cd4_T_cell" = "Cxcl13",
  "Cd4_T_cell" = "Il32",
  "Cd4_T_cell" = "Ltb",
  
  # Cytotoxic / CD8+ T cell (Pan-CD8 markers)
  "Cytotoxic_T_cell" = "Cd8",
  "Cd8_T_cell" = "Cd8a",
  "Cd8_T_cell" = "Ccl4l2",
  "Cd8_T_cell" = "Cd8b",
  "Cd8_T_cell" = "Gzmk",
  "Cd8_T_cell" = "Hla-dra",
  
  # NK cell (Pan-NK markers)
  "Nk_cell" = "Nkg7",
  "Nk_cell" = "Klrd1",
  "Nk_cell" = "Klrf1",
  "Nk_cell" = "Cd56",
  "Nk_cell" = "Fgfbp2",
  "Nk_cell" = "Gnly",
  "Nk_cell" = "Klrb1",
  "Nk_cell" = "Ncam1",
  
  # Treg (Pan-Treg markers)
  "Treg" = "Foxp3",
  "Treg" = "Il2ra",
  "Treg" = "Il32",
  "Treg" = "Tsys",
  
  # ==================== CD4+ 精细亚型 ====================
  # Gzma+ effector Cd4+ T cell
  "Gzma_effector_Cd4_T_cell" = "Cd40lg",
  "Gzma_effector_Cd4_T_cell" = "Gzma",
  "Gzma_effector_Cd4_T_cell" = "Gzmk",
  
  # Naive Cd4+ T cell
  "Naive_Cd4_T_cell" = "Lef1",
  "Naive_Cd4_T_cell" = "Ccr7",
  "Naive_Cd4_T_cell" = "Sell",
  "Naive_Cd4_T_cell" = "Mcr7",
  "Naive_Cd4_T_cell" = "Tcf7",
  
  # Effector memory Cd4+ T cell
  "Effector_memory_Cd4_T_cell" = "Abxa1",
  "Effector_memory_Cd4_T_cell" = "Ccl20",
  "Effector_memory_Cd4_T_cell" = "Ccr6",
  "Effector_memory_Cd4_T_cell" = "Gpr183",
  "Effector_memory_Cd4_T_cell" = "Il2",
  
  # Memory Cd4+ T cell
  "Memory_Cd4_T_cell" = "Abxa2",
  "Memory_Cd4_T_cell" = "S100a4",
  
  # Exhausted Cd4+ T cell
  "Exhausted_Cd4_T_cell" = "ICA1",
  "Exhausted_Cd4_T_cell" = "Tnfrsf4",
  
  # ==================== CD8+ 精细亚型 ====================
  # Gzmb+ effector Cd8+ T cell
  "Gzmk_effector_Cd8_T_cell" = "Gzma",
  "Gzmk_effector_Cd8_T_cell" = "Gzmk",
  "Gzmk_effector_Cd8_T_cell" = "Eomes",
  "Gzmk_effector_Cd8_T_cell" = "Klrg1",
  
  # Naive Cd8+ T cell
  "Naive_Cd8_T_cell" = "Ccr7",
  "Naive_Cd8_T_cell" = "Lef1",
  "Naive_Cd8_T_cell" = "Sell",
  "Naive_Cd8_T_cell" = "Tcf7",
  
  # Pre-effector Cd8+ T cell
  "Pre_effector_Cd8_T_cell" = "Tcf7",
  "Pre_effector_Cd8_T_cell" = "Gzmk",
  
  # Effector memory Cd8+ T cell
  "Effector_memory_Cd8_T_cell" = "Ccl5",
  "Effector_memory_Cd8_T_cell" = "Cxcr3",
  "Effector_memory_Cd8_T_cell" = "Il2",
  "Effector_memory_Cd8_T_cell" = "Il7r",
  "Effector_memory_Cd8_T_cell" = "Zfp683",
  
  # Memory Cd8+ T cell
  "Memory_Cd8_T_cell" = "Zfp683",
  
  # Exhausted Cd8+ T cell
  "Exhausted_Cd8_T_cell" = "Cxcl13",
  "Exhausted_Cd8_T_cell" = "Havcr2",
  "Exhausted_Cd8_T_cell" = "Tigit",
  
  # ==================== NK 精细亚型 ====================
  # Fcgr3ahi_Gzmbhi1_Nk_cell
  "Fcgr3ahi_Gzmbhi1_Nk_cell" = "Fcgr3ahi",
  "Fcgr3ahi_Gzmbhi1_Nk_cell" = "Gzmbhi",
  "Fcgr3ahi_Gzmbhi1_Nk_cell" = "Prf1h1",
  
  # Fcgr3alo_Gzmhhi1_Nk_cell
  "Fcgr3alo_Gzmhhi1_Nk_cell" = "Fcgr3alo",
  "Fcgr3alo_Gzmhhi1_Nk_cell" = "Gzmhh1",
  
  # Fcgr3alo_Xck1hi_Nk_cell
  "Fcgr3alo_Xck1hi_Nk_cell" = "Fcgr3alo",
  "Fcgr3alo_Xck1hi_Nk_cell" = "Xcl1",
  "Fcgr3alo_Xck1hi_Nk_cell" = "Xcl2",
  
  # ==================== Treg 精细亚型 ====================
  # Resolve Treg
  "Suppressive_Treg" = "Ctla4",
  "Suppressive_Treg" = "Foxp3hi",
  "Suppressive_Treg" = "Ikzf2hi",
  "Suppressive_Treg" = "Il2rahi",
  "Suppressive_Treg" = "Tigit",
  
  # Resting Treg
  "Resting_Treg" = "Ctla4",
  "Resting_Treg" = "Foxp3lo",
  "Resting_Treg" = "Ikzf2hi",
  "Resting_Treg" = "Il2ralo",
  "Resting_Treg" = "Tigit",
  
  # ==================== 特殊/应激状态 ====================
  # Hspa1b+ stress-response T cell
  "Hspa1b_stress_response_T_cell" = "Hspa1b",
  "Hspa1b_stress_response_T_cell" = "Cd40lg",
  "Hspa1b_stress_response_T_cell" = "Hspa1a",
  "Hspa1b_stress_response_T_cell" = "Il2",
  "Hspa1b_stress_response_T_cell" = "Tnf",
  
  # NKT cell
  "MAIT_cell" = "Slc4a10"
)
existing_genes <- rownames(all_data)

# 2. 过滤 marker 列表，只保留存在的基因
Tcell_major_marker_ls <- Tcell_major_marker_ls[
  Tcell_major_marker_ls %in% existing_genes
]

keep_idx <- !duplicated(Tcell_major_marker_ls)
Tcell_major_marker_ls <- Tcell_major_marker_ls[keep_idx]
####

Tcell_markers <-
  FindAllMarkers(
    JoinLayers(Tcell_data),
    test.use = "wilcox",
    logfc.threshold = 0,
    min.pct = 0,
    only.pos = TRUE
  )
write.csv(Tcell_markers, file = "T_subtype/Tcell_all_markers.csv", row.names = TRUE)


Tcell_top5_markers <- 
  Tcell_markers %>% 
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

names(Tcell_top5_markers) <- sub('.*\\.','',names(Tcell_top5_markers))


# dotplot
# 如果已做其他细胞类型的亚群分析，这里就可以顺便画一下其他类型的细胞亚群基因，比如Fibro_major_marker_ls2，否则只用画T细胞
pdf(file.path(out_dir, str_glue("dotplot dims{dims} major_marker.pdf")), width = 25, height = 6)
DotPlot(Tcell_data, group.by = "RNA_snn_res.0.2", features = major_marker_ls)+RotatedAxis()+theme_90angle
DotPlot(Tcell_data, group.by = "RNA_snn_res.0.2", features = Tcell_major_marker_ls)+RotatedAxis()+theme_90angle
DotPlot(Tcell_data, group.by = "RNA_snn_res.0.2", features = Tcell_top5_markers)+RotatedAxis()+theme_90angle
dev.off()

Tcell_data$RNA_snn_res.0.1 %>% table()

saveRDS(Tcell_data ,file = "T_subtype/T_cluster.rds")
##
Tcell_data <- readRDS("T_subtype/T_cluster.rds")
## annotation ----------------------------------------------------------------
cell_cluster <- c(
  '0'  = "Cd8+ T cell", # 
  '1'  = "Fibroblast-like T cell", #
  '2'  = "Exhausted Cd8+ T cell",
  '3'  = "Macrophage-like T cell",
  '4'  = "Cd4+ T cell",
  '5'  = "Hspa1b stress response T cell",
  '6'  = "Endothelial-like T cell"
)

Idents(Tcell_data) <- "RNA_snn_res.0.2"



Tcell_data <- RenameIdents(Tcell_data, cell_cluster)
Tcell_data$subtype <- Idents(Tcell_data)


pdf(file.path(out_dir, str_glue("Dimplot dim{dims} annotation.pdf")), height = 6, width = 8)
DimPlot(Tcell_data, group.by = "subtype", label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

pdf(file.path(out_dir, str_glue("Dimplot dim{dims} annotation split.pdf")), height = 6, width = 15)
DimPlot(Tcell_data, group.by = "subtype",split.by = 'group', label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

#######
Tcell_data_filter <- subset(Tcell_data,subtype%in%c("Cd8+ T cell","Exhausted Cd8+ T cell","Cd4+ T cell","Hspa1b stress response T cell"))
## PCA --------
Tcell_data_filter <- FindVariableFeatures(Tcell_data_filter)
length(VariableFeatures(Tcell_data_filter)) #2000
Tcell_data_filter <- ScaleData(Tcell_data_filter)
Tcell_data_filter <- RunPCA(Tcell_data_filter)

out_dir <- "T_subtype_filter"
dir.create(out_dir)
### pca inspector ------
pdf(file.path(out_dir, 'pca plot.pdf'), height = 6, width = 8)
ElbowPlot(Tcell_data_filter, ndims = 50)
dev.off()

pdf(file.path(out_dir, 'dimheatmap 1-5.pdf'), height = 6, width = 10)
DimHeatmap(Tcell_data_filter, dims = 1:5)
dev.off()
pdf(file.path(out_dir, 'dimheatmap 6-10.pdf'), height = 6, width = 10)
DimHeatmap(Tcell_data_filter, dims = 6:10)
dev.off()
pdf(file.path(out_dir, 'dimheatmap 11-15.pdf'), height = 6, width = 10)
DimHeatmap(Tcell_data_filter, dims = 11:15)
dev.off()
pdf(file.path(out_dir, 'dimheatmap 16-20.pdf'), height = 6, width = 10)
DimHeatmap(Tcell_data_filter, dims = 16:20)
dev.off()

## dimension reduction --------
dims <- 10# 18

library(future)
options(future.globals.maxSize = 10 * 1024 ^ 3) # 10Gb
plan(multisession, workers = 10)

Tcell_data_filter <- RunUMAP(Tcell_data_filter, dims = 1:dims,return.model = TRUE)

pdf(file.path(out_dir, str_glue("Dimplot dim{dims} sample.pdf")), height = 8, width = 8)
DimPlot(Tcell_data_filter, group.by = "sample", label = TRUE) + NoLegend()
dev.off()

Tcell_data_filter <- Tcell_data_filter %>% 
  harmony::RunHarmony(
    dims.use = 1:dims,
    group.by.vars = "project",
    # theta = 4, # default 2, bigger value -- more diversity per cluster
    # sigma = 0.5,
    # lambda = 0.5, # default, smaller value -- more integrated
    reduction.save = "harmony"
  )

Tcell_data_filter <- RunUMAP(
  Tcell_data_filter, 
  dims = 1:dims, 
  reduction = "harmony", 
  reduction.name = "harmonyumap", 
  reduction.key  = "harmonyUMAP_",
  return.model = TRUE
)

pdf(file.path(out_dir, str_glue("Dimplot dim{dims} sample - harmony.pdf")), height = 8, width = 8)
DimPlot(Tcell_data_filter, group.by = "sample", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()

reduction.use <- 'harmony' # 'harmony' 'pca'


### cluster -------
Tcell_data_filter  <-
  FindNeighbors(
    Tcell_data_filter ,
    reduction = reduction.use,
    dims = 1:dims,
    k.param = 15  #20改15
  )


Tcell_data_filter <-
  FindClusters(
    Tcell_data_filter,
    resolution = c(0.1,0.2, 0.4, 0.6, 0.8)
  )



# dimplot
pdf(file.path(out_dir, str_glue("Dimplot dim{dims} cluster.pdf")), height = 6, width = 6)
DimPlot(Tcell_data_filter, group.by = "RNA_snn_res.0.1", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Tcell_data_filter, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Tcell_data_filter, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Tcell_data_filter, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Tcell_data_filter, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(Macrophage, group.by = "RNA_snn_res.1", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(Macrophage, group.by = "RNA_snn_res.1.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(Macrophage, group.by = "RNA_snn_res.1.5", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()
Idents(Tcell_data_filter) <- "RNA_snn_res.0.4"

####

Tcell_f_markers <-
  FindAllMarkers(
    JoinLayers(Tcell_data_filter),
    test.use = "wilcox",
    logfc.threshold = 0,
    min.pct = 0,
    only.pos = TRUE
  )
write.csv(Tcell_f_markers, file = "T_subtype_filter/Tcell_f_all_markers.csv", row.names = TRUE)


Tcell_f_top5_markers <- 
  Tcell_f_markers %>% 
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

names(Tcell_f_top5_markers) <- sub('.*\\.','',names(Tcell_f_top5_markers))


# dotplot
# 如果已做其他细胞类型的亚群分析，这里就可以顺便画一下其他类型的细胞亚群基因，比如Fibro_major_marker_ls2，否则只用画T细胞
pdf(file.path(out_dir, str_glue("dotplot dims{dims} major_marker.pdf")), width = 25, height = 6)
DotPlot(Tcell_data_filter, group.by = "RNA_snn_res.0.4", features = major_marker_ls)+RotatedAxis()+theme_90angle
DotPlot(Tcell_data_filter, group.by = "RNA_snn_res.0.4", features = Tcell_major_marker_ls)+RotatedAxis()+theme_90angle
DotPlot(Tcell_data_filter, group.by = "RNA_snn_res.0.4", features = Tcell_f_top5_markers)+RotatedAxis()+theme_90angle
dev.off()

## annotation ----------------------------------------------------------------
cell_cluster <- c(
  '0'  = "Effector memory Cd8+ T cell", # 
  '1'  = "NK cell", #
  '2'  = "Exhausted Cd8+ T cell",
  '3'  = "Fibroblast-like T cell",
  '4'  = "Cd4+ T cell",
  '5'  = "Exhausted Cd8+ T cell",
  '6'  = "Hspa1b stress response T cell",
  '7' = "Cd8+ T cell"
  
)

Idents(Tcell_data) <- "RNA_snn_res.0.4"



Tcell_data_filter <- RenameIdents(Tcell_data_filter, cell_cluster)
Tcell_data_filter$subtype <- Idents(Tcell_data_filter)


pdf(file.path(out_dir, str_glue("Dimplot dim{dims} annotation.pdf")), height = 6, width = 8)
DimPlot(Tcell_data_filter, group.by = "subtype", label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

pdf(file.path(out_dir, str_glue("Dimplot dim{dims} annotation split.pdf")), height = 6, width = 15)
DimPlot(Tcell_data_filter, group.by = "subtype",split.by = 'group', label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

pdf(file.path(out_dir, str_glue("ncount.pdf")), height = 6, width = 15)

saveRDS(Tcell_data_filter ,file = "T_subtype/Tcell_data_filter.rds")


Tcell_plot <- subset(Tcell_data_filter, subtype != "Fibroblast-like T cell")


saveRDS(Tcell_plot ,file = "T_subtype/Tcell_plot.rds")
pdf(file.path(out_dir, str_glue("Dimplot dim{dims} annotation_rmFibro.pdf")), height = 6, width = 8)
DimPlot(
  Tcell_plot,
  group.by = "subtype",
  label = TRUE,
  reduction = 'harmonyumap',
  repel = TRUE
) + NoLegend()
dev.off()

# 第二张分组拆分UMAP
pdf(file.path(out_dir, str_glue("Dimplot dim{dims} annotation split_rmFibro.pdf")), height = 6, width = 15)
DimPlot(
  Tcell_plot,
  group.by = "subtype",
  split.by = 'group',
  label = TRUE,
  reduction = 'harmonyumap',
  repel = TRUE
) + NoLegend()
dev.off()

table(Tcell_plot$subtype)


df = Tcell_plot@meta.data
#df = read.csv("/lvdata/wzb/scRNA/FWSC20240734/ReAll/celltype/meta.xls",sep = "\t",row.names=1)
summary <- table(df[,c('subtype','group')])
summary_t <- as.data.frame(summary)
summary_prop <- summary_t %>%
  group_by(group) %>%             # 按分组(C/T)分组计算
  mutate(total = sum(Freq), prop = Freq/total) %>%  # 计算组内总数、占比
  ungroup()


cell_colors <- c(
  "NK cell" = "#E7A428",                      # 金黄
  "Cd8+ T cell" = "#7B68EE",                  # 中紫蓝
  "Effector memory Cd8+ T cell" = "#2E8B57",  # 海绿
  "Exhausted Cd8+ T cell" = "#C71585",        # 深洋红
  "Cd4+ T cell" = "#87CEEB",                  # 天蓝
  "Hspa1b stress response T cell" = "#8B4513" #  saddle棕
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

pdf(file.path(out_dir, str_glue("T celltype_prop.pdf")),width = 8,height = 8)
p2
dev.off()



df = Tcell_plot@meta.data

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
  file.path(out_dir, "T celltype_prop - all samples.pdf"),
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
    title = "T cell subtype proportion across samples"
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
##################
all_data$subtype <- as.character(all_data$annotation)  # 或 NA

# 提取 KC_data 的细胞名和对应注释
T_cells <- colnames(Tcell_plot)
T_subtype <- Tcell_plot$subtype  # 如 KC_data$cell_type 或 KC_data$RNA_snn_res.0.4

unique(Tcell_plot$subtype)
# 创建命名向量
subtype_map <- setNames(as.character(T_subtype), T_cells)

# 赋值到 all_data
all_data$subtype[match(T_cells, colnames(all_data))] <- subtype_map[T_cells]

# 非 KC 细胞保留原注释
non_kc <- setdiff(colnames(all_data), T_cells)
all_data$subtype[non_kc] <- as.character(all_data$annotation[non_kc])
unique(all_data$subtype)


######
T_celltype <- c(
  "NK cell",
  "Cd8+ T cell",
  "Effector memory Cd8+ T cell",
  "Exhausted Cd8+ T cell",
  "Cd4+ T cell",
  "Hspa1b stress response T cell"
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
      dat2$subtype %in% T_celltype,
      dat2$subtype,           # 保留原始细胞类型名
      "Others"                   # 其他统一
    )
    
    # 构建颜色向量
    plot_colors <- cell_colors[T_celltype]  # 提取免疫细胞颜色
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
  file.path(out_dir, "spatial plot of T subtype.pdf"),
  height = 8, width = 8
)
walk(p_sp_immnue_ls, print)
dev.off()
