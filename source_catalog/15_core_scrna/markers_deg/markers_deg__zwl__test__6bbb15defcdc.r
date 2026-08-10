.libPaths()
library(Seurat)
library(harmony)
library(ggplot2)
library(magrittr)
library(dplyr)

setwd('/mnt/NexData_file/home/zwl/Project/SC2509044733/require2/R/')
getwd()
S<-readRDS('../result.rds')
DimPlot(S,reduction='umap',raster=FALSE)

# S <- NormalizeData(S, normalization.method = "LogNormalize", scale.factor = 10000)
S <- FindVariableFeatures(S, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(S)
S <- ScaleData(S, features = all.genes)
S <- RunPCA(S, features = VariableFeatures(object = S))
set.seed(666)
S <- RunHarmony(S, group.by.vars = "sample")

ElbowPlot(S,ndims = 50,reduction = "harmony")
S <- FindNeighbors(S,reduction = "harmony",dims = 1:30)
S <- FindClusters(S,resolution = 0.5)
S <- RunUMAP(S,reduction = "harmony",dims = 1:30)
# S <- RunTSNE(S,reduction = "harmony",dims = 1:30)
# bdata_cell_names <- rownames(S@meta.data)
# adata_cell_names <- all_cells_features$cell_id
# adata_cell_types <- as.character(all_cells_features$cell_type)
# bdata_cell_types <- as.character(S@meta.data$cell_type)
# mapped_indices <- match(bdata_cell_names, adata_cell_names)
# matched_cell_types <- adata_cell_types[mapped_indices]
# bdata_cell_types[!is.na(mapped_indices)] <- matched_cell_types[!is.na(mapped_indices)]
# S@meta.data$cell_type <- as.factor(bdata_cell_types)

DimPlot(S,reduction='umap',raster=FALSE,group.by = c('RNA_snn_res.0.8','group'),label = T)
gene = c('MSX1','SHH','MSX2','PAX7','SOX2','SOX21','MYBL2','ZFP36L2','NKX6-2','S100B','MAP2','ZFHX3',
                  'ISL1','LHX1','PITX2','NKX2-2','SIM1')

#FeaturePlot(object = S, features = gene,ncol = 6)
p = DotPlot(S,features = gene,group.by = 'RNA_snn_res.0.8')+theme(axis.text.x = element_text(angle = 45,hjust = 1))
ggsave(filename = '1.marker_dotplot1.png',p,w=10,h=7)
ggsave(filename = '1.marker_dotplot1.pdf',p,w=10,h=7)
gene = c(# 通用神经前体/祖细胞标记
  "SOX2", "PAX6", "NES", "SOX1",
  # 脊髓背侧祖细胞域标记
  "PAX3", "PAX7",
  # 脊髓腹侧祖细胞域标记
  # "NKX6-1", "OLIG2", "NKX2-2", "FOXA2",
  # 运动神经元 MN 标记
  "MNX1", "ISL1", "LHX3", "SLC5A7",
  # V2 中间神经元标记
  "GATA3", "VSX2",
  # 神经递质功能亚型标记
  # GABA能（抑制性）
  "GAD1", "GAD2", "SLC32A1",
  # 谷氨酸能（兴奋性）
  "SLC17A6", "SLC17A7",
  # 甘氨酸能
  "SLC6A5",
  # 胆碱能
  "CHAT", "SLC18A3",
  # # HOX 家族（脊髓头尾轴节段身份标记）
  # "HOXA2", "HOXB2",  # 上段颈段
  # "HOXA3", "HOXB3", "HOXA4", "HOXB4",  # 中下段颈段
  "HOXA5", "HOXB5",# "HOXC4", "HOXC5", "HOXD3", "HOXD4",  # 颈胸交界 / 上段胸段
  # "HOXB6", "HOXB7", "HOXB8",  # 中下段胸段
  # 皮层类器官分层标记（文中皮层部分）
  "SATB2","RELN")
#FeaturePlot(object = S, features = gene,cols = 6)
p = DotPlot(S,features = gene,group.by = 'RNA_snn_res.0.8')+theme(axis.text.x = element_text(angle = 45,hjust = 1))
ggsave(filename = '1.marker_dotplot2.png',p,w=10,h=7)
ggsave(filename = '1.marker_dotplot2.pdf',p,w=10,h=7)
saveRDS(S,'D:/项目/FWAE20254246/result_harmony_pc30.rds')

cluster_map <- c(
  "0" = "Glycinergic",
  "1" = "Cholinergic",
  "2" = "HOXB5+ ImmN",
  "3" = "HOXB5+ ImmN",
  "4" = "Glycinergic",
  "5" = "Cholinergic",
  "6" = "MN",
  "7" = "Cholinergic",
  "8" = "pMN",
  "9" = "HOXB5+ ImmN",
  "10" = "MN",
  "11" = "NPC",
  "12" = "Cholinergic",
  "13" = "HOXB5+ ImmN",
  "14" = "Glycinergic",
  "15" = "dI4-6 InN"
)

S$cell_type <- unname(cluster_map[as.character(S$RNA_snn_res.0.8)])
p = DimPlot(S,reduction='umap',raster=FALSE,group.by = c('cell_type','RNA_snn_res.0.8'),label = T)
ggsave(filename = '1.cell_type_umap.png',p,w=10,h=5)
ggsave(filename = '1.cell_type_umap.pdf',p,w=10,h=5)

markers <- FindAllMarkers(
  object = S,
  group.by = 'RNA_snn_res.0.8',
  only.pos = TRUE, # 仅保留阳性表达基因
  min.pct = 0.25, # 基因在细胞中表达的最小比例
  logfc.threshold = 0.25 # 最小 logFC 阈值
)
write.csv(file = '3.cluster_markers.csv',markers)
markers2 <- FindAllMarkers(
  object = S,
  group.by = 'cell_type',
  only.pos = TRUE, # 仅保留阳性表达基因
  min.pct = 0.25, # 基因在细胞中表达的最小比例
  logfc.threshold = 0.25 # 最小 logFC 阈值
)
write.csv(file = '3.celltype_markers.csv',markers2)
