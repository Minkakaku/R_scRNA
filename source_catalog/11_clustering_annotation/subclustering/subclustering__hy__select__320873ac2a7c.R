# ==================== R 脚本 ====================
# 功能：从三个 rds 对象中提取 "matri-prod CAF" 细胞并保存
# 输入文件：
#   All_filter_state_samp_new_with_CNs10-15_C.rds
#   All_filter_state_samp_new_with_CNs10-15_A.rds
#   All_filter_state_samp_new_with_CNs10-15_N.rds
# 输出文件：
#   All_filter_state_samp_new_with_CNs10-15_C_mcaf.rds
#   All_filter_state_samp_new_with_CNs10-15_A_mcaf.rds
#   All_filter_state_samp_new_with_CNs10-15_N_mcaf.rds
# ================================================

# 设置工作目录（若文件不在当前目录，请修改为实际路径）
# setwd("/path/to/your/files")

# 定义样本标识
sample_ids <- c("C", "A", "N")

# 循环处理每个样本
for (id in sample_ids) {
  # 构造输入输出文件名
  infile  <- paste0("All_filter_state_samp_new_with_CNs10-15_", id, ".rds")
  outfile <- paste0("All_filter_state_samp_new_with_CNs10-15_", id, "_mcaf.rds")
  
  # 检查输入文件是否存在
  if (!file.exists(infile)) {
    warning("文件不存在: ", infile, "，跳过")
    next
  }
  
  cat("\n正在处理: ", infile, "\n")
  
  # 读取 RDS 对象
  obj <- readRDS(infile)
  
  # ----- 提取元数据（meta.data / colData）-----
  if (inherits(obj, "Seurat")) {
    meta <- obj@meta.data
  } else if (inherits(obj, "SingleCellExperiment")) {
    meta <- as.data.frame(SummarizedExperiment::colData(obj))
  } else if (is.data.frame(obj)) {
    # 若对象本身是数据框，假定行名为细胞名
    meta <- obj
  } else {
    stop("未知的对象类型，请检查 RDS 文件内容")
  }
  
  # 检查 'celltype' 列是否存在
  if (!"celltype" %in% colnames(meta)) {
    stop("在 obs 中未找到 'celltype' 列，请检查列名")
  }
  
  # 筛选目标细胞（注意大小写和空格完全匹配）
  cells_keep <- rownames(meta)[meta$celltype == "matri-prod CAF"]
  n_cells <- length(cells_keep)
  cat("  找到 ", n_cells, " 个 matri-prod CAF 细胞\n")
  
  if (n_cells == 0) {
    warning("未找到 matri-prod CAF 细胞，将跳过保存")
    next
  }
  
  # ----- 根据对象类型进行子集提取 -----
  if (inherits(obj, "Seurat")) {
    obj_sub <- subset(obj, cells = cells_keep)
  } else if (inherits(obj, "SingleCellExperiment")) {
    obj_sub <- obj[, cells_keep]
  } else if (is.data.frame(obj)) {
    obj_sub <- obj[cells_keep, , drop = FALSE]
  } else {
    stop("无法对该对象类型进行子集提取")
  }
  
  # 保存结果
  saveRDS(obj_sub, file = outfile)
  cat("  已保存: ", outfile, "\n")
}

cat("\n处理完成！\n")
