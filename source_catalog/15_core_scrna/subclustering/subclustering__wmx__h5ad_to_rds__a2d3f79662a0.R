# install.packages('reticulate')
# BiocManager::install(c(
#   "LoomExperiment",
#   "SingleCellExperiment"
# ), ask = FALSE, update = FALSE)
# remotes::install_github("cellgeni/sceasy")

library(reticulate)
library(sceasy,lib.loc = "C:/Users/wmx/AppData/Local/R/win-library/4.5")
library(Seurat)
library(Matrix)
use_condaenv('D:/Softwares/Anaconda/envs/scRNA')

# # 简单转换，这里的rds文件中已经有obs的信息以及umap等的信息
h5ad <- "D:/Project/SC2510235895/shouhou7/Fibroblasts_raw.h5ad"
rds <- "D:/Project/SC2510235895/shouhou7/Fibroblasts.rds"
sceasy::convertFormat(h5ad,from='anndata',
                      to='seurat',outFile=rds)

h5ad <- "D:/Project/SC2510235895/shouhou7/Macrophages_raw.h5ad"
rds <- "D:/Project/SC2510235895/shouhou7/Macrophages.rds"
sceasy::convertFormat(h5ad,from='anndata',
                      to='seurat',outFile=rds)


h5ad <- "D:/Project/SC2510235895/shouhou7/Neutrophils_raw.h5ad"
rds <- "D:/Project/SC2510235895/shouhou7/Neutrophils.rds"
sceasy::convertFormat(h5ad,from='anndata',
                      to='seurat',outFile=rds)


h5ad <- "D:/Project/SC2510235895/shouhou7/Tcells_raw.h5ad"
rds <- "D:/Project/SC2510235895/shouhou7/Tcells.rds"
sceasy::convertFormat(h5ad,from='anndata',
                      to='seurat',outFile=rds)
# ----------------
# 读入转换完成的rds,接下来将信息加入进去(最后生成的结果的scale.data,count矩阵)
# ----------------
# 读 h5ad
ad <- import("anndata", convert = FALSE)
adata <- ad$read_h5ad("D:/Project/SC2510235895/shouhou7/adata_counts.h5ad")

# 取 adata.X
x <- py_to_r(adata$X)

# 如果不是稀疏矩阵，可转为稀疏矩阵
if (!inherits(x, "sparseMatrix")) {
  x <- as(x, "dgCMatrix")
}
# AnnData 通常是 cell x gene
# Seurat 需要 gene x cell，所以要转置
x <- t(x)
x <- as(x, "TsparseMatrix")
x <- as(x, "dgCMatrix")

# 补上行列名
rownames(x) <- py_to_r(adata$var_names$to_list())
colnames(x) <- py_to_r(adata$obs_names$to_list())


rdsfile = list.files("D:/Project/SC2510235895/shouhou7/",pattern = "s.rds",recursive = T,full.names = T)
for(files in rdsfile){
  test <- readRDS(files) 
  dim(test@assays$RNA@data)

  # 和当前 Seurat 对象对齐
  common_cells <- intersect(colnames(test), colnames(x))
  common_genes <- intersect(rownames(test),rownames(x))
  obj <- subset(test,features=common_genes, cells = common_cells)
  x1 <- x[common_genes, common_cells]

  identical(colnames(x1),colnames(obj))
  identical(rownames(x1),rownames(obj))

  # 放到 count 层
  obj[["RNA"]]$counts <- x1
  finalname =  gsub( x = files,pattern  = ".rds",replacement = "_final.rds")
  saveRDS(object = obj,finalname)
}
