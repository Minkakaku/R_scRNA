setwd('/mnt/NexData_file/home/zwl/Project/FWSC20241362/')
rm(list = ls())
gc()
library(Seurat)
library(data.table)
Sys.setenv("VROOM_CONNECTION_SIZE" = 10 * 1024 * 1024)

setwd('./sample12/cellchat/')
sceasy::convertFormat('./adata_raw.h5ad', from="anndata", to="seurat",outFile='adata_raw.rds')

# 读取 UMI 计数文件
umi_counts <- fread("/mnt/NexData_file/home/zwl/Project/FWSC20241362/normal_data/GSE181919/GSE181919_UMI_counts.txt.gz")
umi_counts[1:4,1:4]
# 读取条形码元数据文件
barcode_metadata <- fread("/mnt/NexData_file/home/zwl/Project/FWSC20241362/normal_data/GSE181919/GSE181919_Barcode_metadata.txt.gz")

# 将 UMI 计数数据转换为矩阵，假设第一列是基因名，其余列是细胞条形码
count_matrix <- umi_counts[,-1]
rownames(count_matrix) <- umi_counts$V1
count_matrix[1:4,1:4]
# 创建 Seurat 对象
seurat_obj <- CreateSeuratObject(counts = count_matrix, 
                                 meta.data = barcode_metadata,
                                 min.cells = 3, #在不少于三个细胞中表达
                                 min.features = 200)
saveRDS(seurat_obj,file = './ALL_data.rds',compress = F)
seurat_obj <- readRDS('./ALL_data.rds')
# 查看 Seurat 对象
sce <- subset(seurat_obj,sample.id %in% c('N06','N15','N21'))
write.csv(file = 'sample3_celltype.csv',table(sce$sample.id,sce$cell.type))
# 
# test = Read10X('./normal_data/LB/')
# obj <- CreateSeuratObject(counts = test, min.cells = 3,min.features = 200)

#test:group DEG muscat##################################################################
library(dplyr)
library(ggplot2)
library(limma)
library(muscat)
library(purrr)
library(SingleCellExperiment)
setwd('/mnt/NexData_file/home/zwl/Project/FWSC20241362/')
seurat_obj <- readRDS('./ALL_data.rds')
obj <- subset(seurat_obj,sample.id %in% c('N06','N15','N21', #normal
                                          'C04','C15','C21'))#cancer 
obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
## Filter #####################################################
# remove undetected or lowly expressed genes
obj <- obj[rowSums(GetAssayData(obj,layer = 'counts') > 0) > 10, ]
sce <- as.SingleCellExperiment(obj)
# remove cells

# compute sum-factors & normalize
library(scater)
sce <- scater::computeLibraryFactors(sce)
sce <- logNormCounts(sce)

## Data preparation ###########################################
(sce <- prepSCE(sce, 
                kid = "cell.type", # subpopulation assignments
                gid = "tissue.type",  # group IDs (ctrl/stim)
                sid = "sample.id",   # sample IDs (ctrl/stim.1234)
                drop = TRUE))  # drop all other colData columns

nk <- length(kids <- levels(sce$cluster_id))
ns <- length(sids <- levels(sce$sample_id))
names(kids) <- kids; names(sids) <- sids

## Data overview ##############################################
t(table(sce$cluster_id, sce$sample_id))
sce <- runUMAP(sce, pca = 20)

# wrapper to prettify reduced dimension plots
.plot_dr <- function(sce, dr, col)
  plotReducedDim(sce, dimred = dr, colour_by = col) +
  guides(fill = guide_legend(override.aes = list(alpha = 1, size = 3))) +
  theme_minimal() + theme(aspect.ratio = 1)

# downsample to max. 100 cells per cluster
cs_by_k <- split(colnames(sce), sce$cluster_id)
cs100 <- unlist(sapply(cs_by_k, function(u) 
  sample(u, min(length(u), 100))))

# plot t-SNE & UMAP colored by cluster & group ID
for (dr in c( "UMAP"))
  for (col in c("cluster_id", "group_id"))
    .plot_dr(sce[, cs100], dr, col)

## DS analysis ##############################################
pb <- aggregateData(sce,
                    assay = "counts", fun = "sum",
                    by = c("cluster_id", "sample_id"))
# one sheet per subpopulation
assayNames(pb)
t(head(assay(pb)))   #统计每个样本基因的bin
(pb_mds <- pbMDS(pb))

# run DS analysis
res <- pbDS(pb, verbose = FALSE)
# access results table for 1st comparison
tbl <- res$table[[1]]
# one data.frame per cluster
names(tbl)
# view results for 1st cluster
k1 <- tbl[[1]]
head(format(k1[, -ncol(k1)], digits = 2))
# construct design & contrast matrix
ei <- metadata(sce)$experiment_info
mm <- model.matrix(~ 0 + ei$group_id)
dimnames(mm) <- list(ei$sample_id, levels(ei$group_id))
contrast <- makeContrasts("CA-NL", levels = mm)

# run DS analysis
pbDS(pb, design = mm, contrast = contrast)

# filter FDR < 5%, abs(logFC) > 1 & sort by adj. p-value
tbl_fil <- lapply(tbl, function(u) {
  u <- dplyr::filter(u, p_adj.loc < 0.05, abs(logFC) > 1)
  dplyr::arrange(u, p_adj.loc)
})

# nb. of DS genes & % of total by cluster
n_de <- vapply(tbl_fil, nrow, numeric(1))
p_de <- format(n_de / nrow(sce) * 100, digits = 3)
data.frame("#DS" = n_de, "%DS" = p_de, check.names = FALSE)

# ALL DEG in each cluster bind to dataframe
ALL_DEG <- bind_rows(tbl_fil)
write.csv(file = './normal_data/muscat_allcluster_DEG.csv',ALL_DEG)

## Visualizing results ##############################################
library(UpSetR)
de_gs_by_k <- map(tbl_fil, "gene")
upset(fromList(de_gs_by_k))

library(dplyr)
length(intersect(de_gs_by_k[["Fibroblasts"]],de_gs_by_k[["T.cells"]]))
length(intersect(de_gs_by_k[["Fibroblasts"]],de_gs_by_k[["Endothelial.cells"]]))







