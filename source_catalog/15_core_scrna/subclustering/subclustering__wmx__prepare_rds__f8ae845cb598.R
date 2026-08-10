# 按照要求提取对应的rds
# 因为要做monocle2，需要提取对应的count数据
# 最终注释的数据来自 "/mnt/NexData_file/home/sh/Projects/AE2511046228-02/注释数据/Main_mapping.rds"
# 原始counts矩阵来自"/mnt/NexData_file/scRNA/scVision_output/AE2511046228-02/03-Cell_filter/adata_counts.h5ad"
rm(list = ls())
library(reticulate)
library(sceasy)
library(Seurat)
library(Matrix)
library(tidyverse)
setwd("/mnt/NexData_file/home/wmx/Project/AE2511046228-02/shouhou12/")

# 读入rds数据
rds <- readRDS("/mnt/NexData_file/home/sh/Projects/AE2511046228-02/注释数据/Main_mapping.rds")

# 读入对应的count矩阵
library(hdf5r)
library(Matrix)

h5file <-"/mnt/NexData_file/scRNA/scVision_output/AE2511046228-02/03-Cell_filter/adata_counts.h5ad"
# 读取并整理的函数
library(hdf5r)
library(Matrix)

read_h5ad_X_counts <- function(h5file) {
  f <- H5File$new(h5file, mode = "r")
  on.exit(f$close_all())
  
  cells <- as.character(f[["obs/_index"]]$read())
  
  if ("var" %in% names(f) && "_index" %in% names(f[["var"]])) {
    genes <- as.character(f[["var/_index"]]$read())
  } else if ("var" %in% names(f) && "features" %in% names(f[["var"]])) {
    genes <- as.character(f[["var/features"]]$read())
  } else if ("var" %in% names(f) && "gene_ids" %in% names(f[["var"]])) {
    genes <- as.character(f[["var/gene_ids"]]$read())
  } else {
    stop("没有找到基因名，请先运行 f[['var']]$ls() 查看 var 结构")
  }
  
  message("n cells: ", length(cells))
  message("n genes: ", length(genes))
  
  g <- f[["X"]]
  
  data    <- as.numeric(g[["data"]]$read())
  indices <- as.integer(g[["indices"]]$read())
  indptr  <- as.integer(g[["indptr"]]$read())
  
  message("length(data): ", length(data))
  message("length(indices): ", length(indices))
  message("length(indptr): ", length(indptr))
  
  # AnnData CSR: cells x genes
  # 直接构造成转置后的 Seurat counts: genes x cells, dgCMatrix
  counts <- new(
    "dgCMatrix",
    x = data,
    i = indices,
    p = indptr,
    Dim = c(length(genes), length(cells))
  )
  
  rownames(counts) <- make.unique(genes)
  colnames(counts) <- make.unique(cells)
  
  return(counts)
}

counts <- read_h5ad_X_counts(h5file)

# 和当前 Seurat 对象对齐
common_cells <- intersect(colnames(rds), colnames(counts))
common_genes <- intersect(rownames(rds),rownames(counts))
obj <- subset(rds,features=common_genes, cells = common_cells)
count_data <- counts[common_genes, common_cells]
# 验证
identical(colnames(count_data),colnames(obj))
identical(rownames(count_data),rownames(obj))

# 放到 count 层
obj[["RNA"]]$counts <- count_data

# 保存
saveRDS(object = obj,"/mnt/NexData_file/home/wmx/Project/AE2511046228-02/shouhou12/Main_mapping_count.rds")


# 创建新列 cell_type1，基于原有列 cell_type 和 sub_cell_type
obj$cell_type1 <- case_when(
  # 规则1：B细胞或浆细胞 -> BandPlasma
  obj$cell_type %in% c("B cells", "Plasma cells") ~ "BandPlasma",
  
  # 规则2：指定 CD8T 亚群 -> CD8T
  obj$sub_cell_type %in% c("CD8T naive", "CD8T effector", 
                           "terminally exhausted CD8T cells", 
                           "transitory exhausted CD8T cells ", 
                           "precursor exhausted CD8T cells", 
                           "ProliferativeT") ~ "CD8T",
  
  # 规则3：指定 CD4T 亚群 -> CD4T
  obj$sub_cell_type %in% c("CD4T memory", "CD4T naive", "Tfh", "Th17","Treg") ~ "CD4T",
  
  # 其余情况：保留原有的 cell_type（或你可以设为 NA 或其他）
  TRUE ~ as.character(obj$cell_type)
)

# 检查结果
table(obj$cell_type1)

# 依次提取各个细胞类型的数据，保存rds以便后续的拟时序分析
celltypelist = c("BandPlasma","CD4T","CD8T","Endothelial cells","Fibroblasts")
grouplist = unique(obj$group)
outdir = "/mnt/NexData_file/home/wmx/Project/AE2511046228-02/shouhou12/Monocle_rds/"
dir.create(outdir,showWarnings = F)
safe_name <- function(x) {
  gsub("[/\\\\:;,*?\"<>| ]+", "_", x)
}

for (cellname in celltypelist) {
  
  message("========================================")
  message("[Celltype] Processing: ", cellname)
  
  subobj <- subset(obj, subset = cell_type1 == cellname)
  message("[Celltype] ", cellname, " cells: ", ncol(subobj))
  
  if (ncol(subobj) == 0) {
    message("[Skip] No cells for cell_type = ", cellname)
    next
  }
  
  for (groups in grouplist) {
    
    message("----------------------------------------")
    message("[Group] Processing: ", groups, " | Celltype: ", cellname)
    
    subobj_f <- subset(subobj, subset = group == groups)
    message("[Subset] ", groups, " / ", cellname, " cells: ", ncol(subobj_f))
    
    if (ncol(subobj_f) == 0) {
      message("[Skip] No cells for group = ", groups, ", cell_type = ", cellname)
      next
    }
    
    safe_group <- safe_name(groups)
    safe_cellname <- safe_name(cellname)
    
    outfile <- file.path(outdir, paste0(safe_group, "_", safe_cellname, ".rds"))
    
    message("[Save] ", outfile)
    saveRDS(object = subobj_f, file = outfile)
    message("[Done] Saved: ", groups, " / ", cellname)
  }
}

message("========================================")
message("[All Done] Finished saving subgroup RDS files.")


