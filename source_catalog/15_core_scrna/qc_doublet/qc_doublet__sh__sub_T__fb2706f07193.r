.libPaths(c('/home/szt/miniforge3/envs/rpy/lib/R/library', .libPaths()))
library('harmony',lib.loc = '/home/szt/miniforge3/envs/rpy/lib/R/library' )
library('ggplot2')
library('Seurat')
library('tidyverse')
# library(qs,lib.loc = '/mnt/NexData_file/home/wmx/envs/scRNA/lib/R/library')

set.seed(123)
outdir = "/mnt/NexData_file/home/wmx/Project/AE2507033447-08/20260413/subcluster/T_cells_R/"

dir.create(outdir)

# meta_file = "/mnt/NexData_file/home/sh/Projects/AE2507033447-08/06-Cell_type_annotation/all_cells_features.xls"
# t = read.delim2(meta_file,sep = '\t',header = T)
# rownames(t) = t[,1]
# t = t[,-1]

# # 读入数据，并进行标准化降维聚类处理，使用harmony去批次 -------
# rdsfile = '/mnt/NexData_file/home/wmx/Project/AE2507033447-08/20260413/result.rds'
# sce = readRDS(rdsfile)
# > identical(rownames(t),rownames(sce@meta.data))
# [1] TRUE
# sce@meta.data = t
# T_sub <- subset(sce, subset = cell_type == "T cells")

# if ("RNA" %in% Assays(T_sub)) {
#   try({
#     T_sub[["RNA"]] <- JoinLayers(T_sub[["RNA"]])
#   }, silent = TRUE)
# }

# T_sub <- NormalizeData(T_sub)
# T_sub <- FindVariableFeatures(T_sub, selection.method = "vst", nfeatures = 2000)
# T_sub <- ScaleData(T_sub)

# T_sub <- RunPCA(T_sub, npcs = 50, seed.use = 123)
# # 肘击图
# p_elbow <- ElbowPlot(T_sub, ndims = 50)
# ggsave(
#   filename = file.path(outdir, "T_sub_pca_elbowplot.pdf"),
#   plot = p_elbow,
#   width = 6,
#   height = 4
# )
# ggsave(
#   filename = file.path(outdir, "T_sub_pca_elbowplot.png"),
#   plot = p_elbow,
#   width = 6,
#   height = 4,
#   dpi = 300
# )

# pca_sample <- DimPlot(T_sub,reduction='pca',group.by = 'sample')
# ggsave(
#   filename = file.path(outdir, "T_sub_pca_sample.pdf"),
#   plot = pca_sample,
#   width = 6,
#   height = 4
# )
# ggsave(
#   filename = file.path(outdir, "T_sub_pca_sample.png"),
#   plot = pca_sample,
#   width = 6,
#   height = 4,
#   dpi = 300
# )
# # Harmony去批次
# T_sub <- RunHarmony(T_sub, group.by.vars = "sample")
# p_elbow <- ElbowPlot(T_sub, ndims = 50,reduction = 'harmony')
# ggsave(
#   filename = file.path(outdir, "T_sub_harmony_elbowplot.pdf"),
#   plot = p_elbow,
#   width = 6,
#   height = 4
# )
# ggsave(
#   filename = file.path(outdir, "T_sub_harmony_elbowplot.png"),
#   plot = p_elbow,
#   width = 6,
#   height = 4,
#   dpi = 300
# )
# # 聚类
# T_sub <- FindNeighbors(T_sub, reduction = "harmony", dims = 1:20)
# T_sub <- RunUMAP(T_sub, reduction = "harmony", dims = 1:20)
# T_sub <- RunTSNE(T_sub, reduction = "harmony", dims = 1:20)
# T_sub <- FindClusters(T_sub, resolution = 0.5)

# # 保存数据
# saveRDS(T_sub ,paste0(outdir, "T_sub_harmony_umap_tsne.rds"))
T_sub <- readRDS(paste0(outdir, "T_sub_harmony_umap_tsne.rds"))

# 进行non-T，CD4，CD8 划分 ----------
# 1 簇这里 CD4 CD8 都表达，需要分辨率提升一下，把它分开
# 调一下降维聚类的参数
T_sub <- FindNeighbors(T_sub, reduction = "harmony", dims = 1:16, verbose = TRUE, k.param = 20)
T_sub <- RunUMAP(T_sub, reduction = "harmony", dims = 1:16 , verbose = TRUE,n.neighbors = 20,min.dist = 0.5)
T_sub <- RunTSNE(T_sub, reduction = "harmony", dims = 1:16, verbose = TRUE)
T_sub <- FindClusters(T_sub, resolution = 1.0)
T_sub <- FindClusters(T_sub, resolution = 0.8)

# DimPlot(T_sub,
#                reduction = 'umap',
#                group.by = "RNA_snn_res.0.8",
#                label = T
#                )
# 1簇有点奇怪，看一下差异基因 ，定为ILC
# tets <- FindMarkers(T_sub,ident.1 = "1")
# tets <- tets %>%
#   dplyr::filter(avg_log2FC > 0.25 & p_val_adj < 0.05) %>%
#   dplyr::arrange(desc(avg_log2FC))
# head(tets,n =100)
# T_sub$pct_counts_mt <- as.numeric(T_sub$pct_counts_mt)
# T_sub$pct_counts_ribo <- as.numeric(T_sub$pct_counts_ribo)
# FeaturePlot(
#     T_sub,
#     features = c(
#         'pct_counts_mt','pct_counts_ribo'
#     ),
#     ncol = 2
# )


p_group <- DimPlot(T_sub , 
                    reduction = 'umap' , 
                    split.by = 'group',
                    group.by = 'RNA_snn_res.0.8'
                    )
p_sample <- DimPlot(T_sub , 
                    reduction = 'umap' , 
                    split.by = 'sample',
                    ncol = 6,
                    group.by = 'RNA_snn_res.0.8'
                    )
ggsave(
  filename = file.path(outdir, "T_sub_cluster_group.pdf"),
  plot = p_group,
  width = 16,
  height = 6,
  dpi = 300
)

ggsave(
  filename = file.path(outdir, "T_sub_cluster_sample.pdf"),
  plot = p_sample,
  width = 24,
  height = 24,
  dpi = 300
)

count_data_sample <- T_sub@meta.data %>%
  group_by(sample, RNA_snn_res.0.8) %>%
  summarise(cell_count = n(), .groups = "drop")

p <- ggplot(count_data_sample, aes(x = RNA_snn_res.0.8, y = cell_count, color = RNA_snn_res.0.8)) +
  geom_point(size = 2) +
  facet_wrap(~ sample, scales = "free_y") +
  labs(x = "Cluster", y = "Cell count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(p ,filename = paste0(outdir,"test.pdf"),width = 25,height = 24)

count_data_group <- T_sub@meta.data %>%
  group_by(group, RNA_snn_res.0.8) %>%
  summarise(cell_count = n(), .groups = "drop")

p1 <- ggplot(count_data_group, aes(x = RNA_snn_res.0.8, y = cell_count, color = RNA_snn_res.0.8)) +
  geom_point(size = 2) +
  facet_wrap(~ group, scales = "free_y") +
  labs(x = "Cluster", y = "Cell count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(p1 ,filename = paste0(outdir,"test1.pdf"),width = 16,height = 6)


P <- FeaturePlot(
    T_sub,
    features = c(
        "PTPRC","EPCAM",
        "MS4A1","C1QA","PECAM1","COL1A1","CD14","S100A9","S100A8",
        "MKI67","TOP2A","PCNA",
        # T 
        "CD3D","CD3E","CD3G",
        "TRAC","TRDC","TRGC1",
        # NK
        "GNLY","NKG7","NCR1",
        # ILC
        'IL2RA', 'ID2', 'IL7R', 'GATA3',
        # Others 
        "CD4","CD8A","CD8B","doublet_score"
    ),
    ncol = 6
)

# DotPlot(T_sub,features = c('IL2RA', 'ID2', 'RORA', 'GATA3',"IL7R"),group.by = 'RNA_snn_res.0.8')
# FeaturePlot(T_sub,features = c('IL2RA', 'ID2', 'RORA', 'GATA3',"IL7R"))
ggsave(
  filename = file.path(outdir, "T_sub_features.pdf"),
  plot = P,
  width = 30,
  height = 20,
  dpi = 300
)

ggsave(
  filename = file.path(outdir, "T_sub_features.png"),
  plot = P,
  width = 30,
  height = 20,
  dpi = 300
)

T_sub$doublet_score <- as.numeric(T_sub$doublet_score)
P1 <- DotPlot(T_sub,   
        features = c(
        # T 
        "CD3D","CD3E","CD3G",
        "TRDC",'TRGC1', 'TRGC2',
        "S100A9","S100A8",
        # NK
        "GNLY","NKG7","NCR1",
        # Cycling
        "MKI67","TOP2A","PCNA",
        # ILC
        'IL2RA', 'ID2', 'IL7R', 'GATA3',
        # Others 
        "CD4","CD8A","CD8B","doublet_score"
    ),
    group.by ='RNA_snn_res.0.8'
    )
ggsave(
  filename = file.path(outdir, "T_sub_features_dot.pdf"),
  plot = P1,
  width = 20,
  height = 14,
  dpi = 300
)
ggsave(
  filename = file.path(outdir, "T_sub_features_dot.png"),
  plot = P1,
  width = 20,
  height = 14,
  dpi = 300
)

# 绘制降维聚类的图 --------
p_d <- DimPlot(T_sub,
               reduction = 'umap',
               group.by = 'group'
               )

p_d1 <- DimPlot(T_sub,
               reduction = 'umap',
               group.by = "RNA_snn_res.0.8",
               label = T
               )

p_d2 <- DimPlot(T_sub,
               reduction = 'tsne',
               group.by = "RNA_snn_res.0.8",
               label = T
               )

p_d3 <- DimPlot(T_sub,
               reduction = 'umap',
               group.by = "RNA_snn_res.0.8",
               )

p_d4 <- DimPlot(T_sub,
               reduction = 'tsne',
               group.by = "RNA_snn_res.0.8",
               )
ggsave(
  filename = file.path(outdir, "T_sub_umap_group.pdf"),
  plot = p_d,
  width = 10,
  height = 7,
)
ggsave(
  filename = file.path(outdir, "T_sub_umap_group.png"),
  plot = p_d,
  width = 10,
  height = 7,
  dpi = 300
)

ggsave(
  filename = file.path(outdir, "T_sub_umap.pdf"),
  plot = p_d3,
  width = 10,
  height = 7,
)
ggsave(
  filename = file.path(outdir, "T_sub_umap.png"),
  plot = p_d3,
  width = 10,
  height = 7,
  dpi = 300
)

ggsave(
  filename = file.path(outdir, "T_sub_umap_label.pdf"),
  plot = p_d1,
  width = 10,
  height = 7,
)
ggsave(
  filename = file.path(outdir, "T_sub_umap_label.png"),
  plot = p_d1,
  width = 10,
  height = 7,
  dpi = 300
)

ggsave(
  filename = file.path(outdir, "T_sub_tsne_label.pdf"),
  plot = p_d4,
  width = 10,
  height = 7,
)
ggsave(
  filename = file.path(outdir, "T_sub_tsne_label.png"),
  plot = p_d4,
  width = 10,
  height = 7,
  dpi = 300
)

ggsave(
  filename = file.path(outdir, "T_sub_tsne.pdf"),
  plot = p_d4,
  width = 10,
  height = 7,
)
ggsave(
  filename = file.path(outdir, "T_sub_tsne.png"),
  plot = p_d4,
  width = 10,
  height = 7,
  dpi = 300
)

T_sub$total_counts <- as.numeric(T_sub$total_counts)
T_sub$n_genes_by_counts <- as.numeric(T_sub$n_genes_by_counts)
pp <- VlnPlot(T_sub, features="total_counts", group.by="RNA_snn_res.0.8", pt.size=0) + theme(legend.position = 'none')
pp1 <- VlnPlot(T_sub, features="n_genes_by_counts", group.by="RNA_snn_res.0.8", pt.size=0) + theme(legend.position = 'none')
library(patchwork)
pp2 <- pp + pp1
ggsave(
  filename = file.path(outdir, "T_sub_qc.pdf"),
  plot = pp2,
  width = 10,
  height = 10,
)
ggsave(
  filename = file.path(outdir, "T_sub_qc.png"),
  plot = pp2,
  width = 10,
  height = 10,
  dpi = 300
)

# 根据结果注释细胞
result <- cbind(seq(0,19),c("CD8_T",
                            "CD4_T",
                            "CD8_T",
                            "ILC",
                            "CD4_T",
                            "CD8_T",
                            "CD4_T",
                            "NK",
                            "NKT",
                            "CD4_T",
                            "Cycling",
                            "ILC",
                            "NK",
                            "Doublets",
                            "ILC",
                            "Cycling",
                            "CD8_T",
                            "NK",
                            "Cycling",
                            "Doublets"
)) %>% as.data.frame()

# 选择结果写入seurat对象
new.cluster.ids <- result$V2
T_sub$RNA_snn_res.0.8 <-  factor(T_sub$RNA_snn_res.0.8,levels = as.character(c(0:19)))
names(new.cluster.ids) <- levels(T_sub$RNA_snn_res.0.8)
Idents(T_sub) <- "RNA_snn_res.0.8"
sce2 <- RenameIdents(T_sub, new.cluster.ids)
sce2@meta.data[["CellType"]] <- sce2@active.ident

p1 <- DimPlot(sce2, reduction = "umap",
              label = TRUE, pt.size = 0.5) + NoLegend()
p2 <- DimPlot(sce2, reduction = "umap",group.by = "RNA_snn_res.0.8",
              label = TRUE, pt.size = 0.5) + NoLegend()

p1+p2

saveRDS(sce2 ,paste0(outdir, "T_sub_harmony_umap_tsne_cell_type.rds"))


p_d5 <- DimPlot(sce2,
               reduction = 'umap',
               group.by = "CellType",
               label = T
               )

p_d6 <- DimPlot(sce2,
               reduction = 'umap',
               group.by = "CellType",
               )


ggsave(
  filename = file.path(outdir, "T_sub_umap_celltype_label.pdf"),
  plot = p_d5,
  width = 10,
  height = 7,
)
ggsave(
  filename = file.path(outdir, "T_sub_umap_celltype_label.png"),
  plot = p_d5,
  width = 10,
  height = 7,
  dpi = 300
)

ggsave(
  filename = file.path(outdir, "T_sub_umap_celltype.pdf"),
  plot = p_d6,
  width = 10,
  height = 7,
)
ggsave(
  filename = file.path(outdir, "T_sub_umap_celltype.png"),
  plot = p_d6,
  width = 10,
  height = 7,
  dpi = 300
)


# 根据CellType 将 各个分开后再做降维聚类，并且保存
cellnames <- c('CD4_T','CD8_T')
for(cells in cellnames){
    outdir1 = paste0(outdir,cells)
    dir.create(outdir1)
    message("Processing: ", cells)
    sub <- subset(sce2,subset = CellType == cells)
    DefaultAssay(sub) <- "RNA"
    sub <- NormalizeData(sub, verbose = FALSE)
    sub <- FindVariableFeatures(sub, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
    sub <- ScaleData(sub, verbose = FALSE)
    sub <- RunPCA(sub, npcs = 50, verbose = FALSE)
    p_elbow <- ElbowPlot(sub, ndims = 50)
    ggsave(
        filename = file.path(outdir1, "T_sub_pca_elbowplot.pdf"),
        plot = p_elbow,
        width = 6,
        height = 4
    )
    ggsave(
        filename = file.path(outdir1, "T_sub_pca_elbowplot.png"),
        plot = p_elbow,
        width = 6,
        height = 4,
        dpi = 300
    )

    pca_sample <- DimPlot(sub,reduction='pca',group.by = 'sample')
    ggsave(
        filename = file.path(outdir1, "T_sub_pca_sample.pdf"),
        plot = pca_sample,
        width = 6,
        height = 4
    )
    ggsave(
        filename = file.path(outdir1, "T_sub_pca_sample.png"),
        plot = pca_sample,
        width = 6,
        height = 4,
        dpi = 300
    )
    # Harmony去批次
    sub <- RunHarmony(sub, group.by.vars = "sample")
    p_elbow <- ElbowPlot(sub, ndims = 50,reduction = 'harmony')
    ggsave(
        filename = file.path(outdir1, "T_sub_harmony_elbowplot.pdf"),
        plot = p_elbow,
        width = 6,
        height = 4
    )
    ggsave(
        filename = file.path(outdir1, "T_sub_harmony_elbowplot.png"),
        plot = p_elbow,
        width = 6,
        height = 4,
        dpi = 300
    )
    # sub <- FindNeighbors(sub, reduction = "harmony", dims = 1:21)
    # sub <- RunUMAP(sub, reduction = "harmony", dims = 1:21, verbose = FALSE)
    # sub <- RunTSNE(sub, reduction = "harmony", dims = 1:21, verbose = FALSE)
    # sub <- FindClusters(sub, resolution = 0.5, verbose = FALSE)


    saveRDS(sub, file = file.path(outdir1, paste0(cells, "_sub.rds")))
  
    # 提取 metadata 和降维坐标
    meta_df <- sub@meta.data
    meta_df$cell_barcode <- rownames(meta_df)
  
    pca_df <- as.data.frame(Embeddings(sub, "pca"))
    pca_df$cell_barcode <- rownames(pca_df)
  
    umap_df <- as.data.frame(Embeddings(sub, "umap"))
    umap_df$cell_barcode <- rownames(umap_df)
  
    tsne_df <- as.data.frame(Embeddings(sub, "tsne"))
    tsne_df$cell_barcode <- rownames(tsne_df)
  
    # 保存表格
    write.csv(meta_df, file = file.path(outdir1, paste0(cells, "_metadata.csv")), row.names = FALSE)
    write.csv(pca_df,  file = file.path(outdir1, paste0(cells, "_pca.csv")), row.names = FALSE)
    write.csv(umap_df, file = file.path(outdir1, paste0(cells, "_umap.csv")), row.names = FALSE)
    write.csv(tsne_df, file = file.path(outdir1, paste0(cells, "_tsne.csv")), row.names = FALSE)
}
