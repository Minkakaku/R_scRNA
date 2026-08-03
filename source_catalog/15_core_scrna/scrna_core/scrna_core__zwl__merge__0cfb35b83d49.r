setwd('/mnt/NexData_file/home/zwl/Project/FWAE20252150')
library(fs)
library(Seurat)

# 1. 获取所有样本目录
s_files <- list.files(
  path = "/mnt/NexData_file/nas/outputdata/cellranger_out/FWAE20252150-02",
  pattern = "^S",  # 正则表达式：^表示开头，S匹配大写S
  full.names = F)

# 2. 创建空列表存储Seurat对象
seurat_list <- list()

# 3. 循环读取每个样本
for (i in s_files) {
  sample_name <- i
  matrix_dir <- paste0("/mnt/NexData_file/nas/outputdata/cellranger_out/FWAE20252150-02/",i,
                       "/outs/", "filtered_feature_bc_matrix/")
  
  # 读取数据并创建Seurat对象
  data <- Read10X(data.dir = matrix_dir)
  seurat_obj <- CreateSeuratObject(
    counts = data,
    project = sample_name,
    min.cells = 3,  # 过滤低表达基因
    min.features = 200  # 过滤低质量细胞
  )
  
  # 添加样本标识
  seurat_obj[["sample"]] <- sample_name
  
  # 添加到列表
  seurat_list[[sample_name]] <- seurat_obj
}

# 4. 合并所有Seurat对象
if (length(seurat_list) > 0) {
  combined_seurat <- merge(
    x = seurat_list[[1]],
    y = seurat_list[-1],
    add.cell.ids = names(seurat_list)
  )
  
  # 查看合并后对象信息
  print(combined_seurat)
  
  # 可选：保存合并后的对象
  saveRDS(combined_seurat, file = "combined_seurat02.rds",compress = F)
} else {
  warning("No valid samples found to merge")
}
