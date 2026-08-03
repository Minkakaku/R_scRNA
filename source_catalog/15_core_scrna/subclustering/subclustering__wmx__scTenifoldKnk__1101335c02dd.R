#!/usr/bin/env Rscript
# =============================================================================
# scTenifoldKnk 虚拟基因敲除分析
# 输入: Seurat RDS 对象 + 目标基因名
# 输出: 差异调控基因表 + KO网络图 + 富集分析
# =============================================================================
suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(scTenifoldKnk)
  library(igraph)
  library(ggplot2)
  set.seed(1234)
})

# =============================================================================
# 命令行参数
# =============================================================================
option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = NULL, help = "输入Seurat RDS文件 [必需]"),
  make_option(c("-o", "--outdir"), type = "character", default = "./knk_results", help = "输出目录"),
  make_option(c("-g", "--ko_gene"), type = "character", default = NULL, help = "目标敲除基因名 [必需]"),
  make_option(c("--celltype_col"), type = "character", default = NULL, help = "细胞类型列名(可选,用于筛选子集)"),
  make_option(c("--celltype"), type = "character", default = NULL, help = "目标细胞类型,逗号分隔多个,如'T cell,B cell'"),
  make_option(c("--max_cells"), type = "integer", default = 10000, help = "最大细胞数,超出随机降采样(默认2000)"),
  make_option(c("--n_hvg"), type = "integer", default = 2000, help = "高变基因数量(默认2000,设为0使用全部基因)"),
  make_option(c("--min_cells"), type = "integer", default = 10, help = "基因最少表达细胞数(默认10)"),
  make_option(c("--nComp"), type = "integer", default = 3, help = "TICA主成分数(默认3)"),
  make_option(c("--fdr_threshold"), type = "double", default = 0.05, help = "FDR显著性阈值"),
  make_option(c("--width"), type = "double", default = 10, help = "图片宽度(英寸)"),
  make_option(c("--height"), type = "double", default = 8, help = "图片高度(英寸)"),
  make_option(c("--dpi"), type = "integer", default = 300, help = "图片分辨率")
)

opt_parser <- OptionParser(option_list = option_list, usage = "%prog -i <input.rds> -g <gene> -o <output_dir> [options]")
opt <- parse_args(opt_parser)

if (is.null(opt$input)) { print_help(opt_parser); stop("错误: 必须指定 --input", call. = FALSE) }
if (is.null(opt$ko_gene)) { print_help(opt_parser); stop("错误: 必须指定 --ko_gene", call. = FALSE) }
if (!file.exists(opt$input)) stop(paste("错误: 文件不存在:", opt$input), call. = FALSE)
if (!dir.exists(opt$outdir)) dir.create(opt$outdir, recursive = TRUE)

cat("=== scTenifoldKnk 虚拟基因敲除分析 ===\n")
cat("输入:", opt$input, "\n")
cat("敲除基因:", opt$ko_gene, "\n")
cat("输出目录:", opt$outdir, "\n")

# =============================================================================
# 0. 兼容各版本的seurat的读取函数
# =============================================================================

get_counts_matrix <- function(seu, assay = "RNA") {
  
  DefaultAssay(seu) <- assay
  
  if ("Assay5" %in% class(seu[[assay]])) {
    cat("检测到 Seurat v5 Assay5，执行 JoinLayers...\n")
    seu <- JoinLayers(seu, assay = assay)
  }
  
  mat <- tryCatch({
    GetAssayData(seu, assay = assay, layer = "counts")
  }, error = function(e1) {
    tryCatch({
      GetAssayData(seu, assay = assay, slot = "counts")
    }, error = function(e2) {
      stop("无法提取 counts 矩阵，请检查 RNA assay 是否包含 counts/layer counts")
    })
  })
  
  return(as.matrix(mat))
}

# =============================================================================
# 1. 读取数据
# =============================================================================
cat("\n[1/5] 读取 Seurat 对象...\n")
file_ext <- tolower(tools::file_ext(opt$input))
if (file_ext == "h5ad") {
  cat("检测到 H5AD 格式，使用 schard::h5ad2seurat 转换...\n")
  seu <- schard::h5ad2seurat(opt$input)
} else {
  seu <- readRDS(opt$input)
}
cat("原始数据:", ncol(seu), "cells x", nrow(seu), "genes\n")

# 可选:按细胞类型筛选
if (!is.null(opt$celltype_col) && !is.null(opt$celltype)) {
  if (!opt$celltype_col %in% colnames(seu@meta.data)) {
    stop(paste("错误: metadata中不存在列:", opt$celltype_col), call. = FALSE)
  }
  types <- trimws(unlist(strsplit(opt$celltype, ",")))
  seu <- subset(seu, cells = colnames(seu)[seu@meta.data[[opt$celltype_col]] %in% types])
  cat("筛选后:", ncol(seu), "cells (", paste(types, collapse = ", "), ")\n")
}

# 降采样
if (ncol(seu) > opt$max_cells) {
  set.seed(1234)
  cells_keep <- sample(colnames(seu), opt$max_cells)
  seu <- subset(seu, cells = cells_keep)
  cat("降采样至:", ncol(seu), "cells\n")
}

# 检查目标基因是否存在
if (!opt$ko_gene %in% rownames(seu)) {
  available_genes <- rownames(seu)
  similar <- agrep(opt$ko_gene, available_genes, max.distance = 0.2, value = TRUE)
  msg <- paste("错误: 基因", opt$ko_gene, "不在表达矩阵中。")
  if (length(similar) > 0) msg <- paste(msg, "相近的基因名:", paste(similar[1:min(5, length(similar))], collapse = ", "))
  stop(msg, call. = FALSE)
}

# --- HVG 筛选 (性能优化) ---
if (opt$n_hvg > 0 && nrow(seu) > opt$n_hvg) {
  cat("筛选高变基因 (n_hvg =", opt$n_hvg, ")...\n")

  seu <- NormalizeData(seu, verbose = FALSE)
  seu <- FindVariableFeatures(seu, nfeatures = opt$n_hvg, verbose = FALSE)
  hvg <- VariableFeatures(seu)

  # 确保 KO 目标基因一定包含在 HVG 中
  if (!opt$ko_gene %in% hvg) {
    hvg <- c(opt$ko_gene, hvg)
    cat("注意: KO基因", opt$ko_gene, "不在 HVG 中,已强制加入\n")
  }

  seu <- subset(seu, features = hvg)
  cat("HVG 筛选后:", nrow(seu), "genes x", ncol(seu), "cells\n")
} else {
  cat("跳过 HVG 筛选,使用全部", nrow(seu), "基因\n")
}

# 提取 counts 矩阵 (genes x cells)
countMatrix <- get_counts_matrix(seu, assay = "RNA")
cat("表达矩阵:", nrow(countMatrix), "genes x", ncol(countMatrix), "cells\n")

# --- 过滤低表达基因 (防止 lm.fit 0 non-NA cases 错误) ---
gene_cell_count <- rowSums(countMatrix > 0)
keep_genes <- gene_cell_count >= opt$min_cells
# 确保 KO 基因不被过滤
if (opt$ko_gene %in% rownames(countMatrix)) {
  keep_genes[opt$ko_gene] <- TRUE
}
if (sum(keep_genes) < nrow(countMatrix)) {
  removed <- nrow(countMatrix) - sum(keep_genes)
  countMatrix <- countMatrix[keep_genes, ]
  cat("过滤低表达基因:", removed, "个基因在 <", opt$min_cells, "个细胞中表达,剩余", nrow(countMatrix), "个\n")
}

# =============================================================================
# 2. 运行 scTenifoldKnk
# =============================================================================
cat("\n[2/5] 运行 scTenifoldKnk 虚拟敲除...\n")
cat("这可能需要几分钟,取决于数据量...\n")

knk_result <- scTenifoldKnk(
  countMatrix = countMatrix,
  gKO = opt$ko_gene,
  qc = FALSE,
  nc_nComp = opt$nComp
)

cat("虚拟敲除完成!\n")

# =============================================================================
# 3. 差异调控分析结果
# =============================================================================
cat("\n[3/5] 提取差异调控基因...\n")

diff_reg <- knk_result$diffRegulation
diff_reg <- diff_reg[order(diff_reg$p.adj), ]

sig_genes <- diff_reg[diff_reg$p.adj < opt$fdr_threshold & diff_reg$distance > 1e-10, ]
cat("显著差异调控基因数量 (FDR <", opt$fdr_threshold, "):", nrow(sig_genes), "\n")

# 保存完整表格
write.csv(diff_reg, file.path(opt$outdir, "diff_regulation_all.csv"), row.names = FALSE)
cat("完整结果:", file.path(opt$outdir, "diff_regulation_all.csv"), "\n")

# 保存显著基因表
if (nrow(sig_genes) > 0) {
  write.csv(sig_genes, file.path(opt$outdir, "diff_regulation_significant.csv"), row.names = FALSE)
  cat("显著基因:", file.path(opt$outdir, "diff_regulation_significant.csv"), "\n")
}

# =============================================================================
# 4. KO 网络可视化
# =============================================================================
cat("\n[4/5] 生成 KO 网络图...\n")

plot_ko_network <- function(knk_result, gKO, outdir, width, height, dpi, fdr_thresh = 0.05, q = 0.99) {
  dr <- knk_result$diffRegulation
  gList <- unique(c(gKO, dr$gene[dr$distance > 1e-10 & dr$p.adj < fdr_thresh]))

  if (length(gList) < 2) {
    cat("警告: 显著基因过少,跳过网络图\n")
    return(invisible(NULL))
  }

  # 从WT网络提取子网
  wt_net <- as.matrix(knk_result$tensorNetworks$WT)
  common_genes <- intersect(gList, rownames(wt_net))
  if (length(common_genes) < 2) {
    cat("警告: 网络中匹配基因过少,跳过网络图\n")
    return(invisible(NULL))
  }

  sCluster <- wt_net[common_genes, common_genes]
  koInfo <- sCluster[gKO, ]

  # 过滤弱边
  threshold <- quantile(abs(sCluster), q)
  sCluster[abs(sCluster) <= threshold] <- 0
  sCluster[gKO, ] <- koInfo
  diag(sCluster) <- 0

  # 转为边列表
  edges <- reshape2::melt(sCluster)
  colnames(edges) <- c("from", "to", "W")
  edges <- edges[edges$W != 0, ]

  if (nrow(edges) == 0) {
    cat("警告: 无有效边,跳过网络图\n")
    return(invisible(NULL))
  }

  # 构建igraph
  g <- graph_from_data_frame(edges, directed = TRUE)
  deg <- degree(g, mode = "all")
  node_size <- 8 + (deg / max(deg)) * 15
  edge_color <- ifelse(E(g)$W > 0, "#E74C3C", "#3498DB")
  v_color <- ifelse(names(V(g)) == gKO, "#E74C3C80", "#2ECC7180")

  png(file.path(outdir, "ko_network.png"), width = width, height = height, units = "in", res = dpi)
  par(mar = c(2, 2, 3, 2))
  set.seed(42)
  plot(g,
    layout = layout_with_fr(g),
    edge.arrow.size = 0.3,
    edge.color = edge_color,
    edge.curved = 0.1,
    vertex.size = node_size,
    vertex.color = v_color,
    vertex.frame.color = NA,
    vertex.label.color = "black",
    vertex.label.cex = 0.7,
    vertex.label.family = "sans",
    main = paste0("Virtual Knockout: ", gKO)
  )
  dev.off()
  cat("网络图:", file.path(outdir, "ko_network.png"), "\n")
}

tryCatch(
  plot_ko_network(knk_result, opt$ko_gene, opt$outdir, opt$width, opt$height, opt$dpi, opt$fdr_threshold),
  error = function(e) cat("网络图生成失败:", conditionMessage(e), "\n")
)


# =============================================================================
# 5. 差异调控基因统计可视化
# =============================================================================
cat("\n[5/5] 生成统计可视化...\n")

# 5a. FC 排名图
tryCatch({
  # 排除敲除基因自身
  sig_genes_plot <- sig_genes[sig_genes$gene != opt$ko_gene, ]
  top_n <- min(20, nrow(sig_genes_plot))
  if (top_n > 0) {
    plot_data <- head(sig_genes_plot, top_n)
    plot_data$gene <- factor(plot_data$gene, levels = rev(plot_data$gene))
    
    p <- ggplot(plot_data, aes(x = gene, y = FC)) +
      geom_bar(stat = "identity", width = 0.7, fill = "#3498DB") +
      coord_flip() +
      scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
      labs(
        title = paste0("Top ", top_n, " Perturbed Genes (", opt$ko_gene, " KO)"),
        subtitle = paste0("FDR < ", opt$fdr_threshold),
        x = NULL, y = "Fold Change"
      ) +
      theme_minimal(base_size = 13) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
        plot.subtitle = element_text(hjust = 0.5, color = "black", size = 11),
        axis.text.y = element_text(face = "bold.italic", size = 12, margin = margin(r = 2)),
        axis.text.x = element_text(face = "bold", size = 12),
        axis.line = element_line(color = "black", linewidth = 0.5),
        panel.grid = element_blank(),
        legend.position = "none",
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA),
        plot.margin = margin(10, 15, 10, 10)
      )
    ggsave(file.path(opt$outdir, "top_perturbed_genes.png"), p,
           width = opt$width, height = opt$height, dpi = opt$dpi
    )
    cat("排名图:", file.path(opt$outdir, "top_perturbed_genes.png"), "\n")
  }
}, error = function(e) cat("排名图生成失败:", conditionMessage(e), "\n"))

# 5b. Z vs -log10(p-value) 散点图
tryCatch({
  df_vol <- knk_result$diffRegulation
  df_vol <- df_vol[df_vol$gene != opt$ko_gene, ]
  df_vol$log_pval <- -log10(df_vol$p.adj)
  label_genes <- subset(df_vol, abs(Z) > 2 & p.adj < 0.01)
  
  p2 <- ggplot(df_vol, aes(x = Z, y = log_pval)) +
    geom_point(alpha = 0.5) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
    geom_vline(xintercept = c(2), linetype = "dashed", color = "blue") +
    ggrepel::geom_text_repel(data = label_genes, aes(label = gene),
                             size = 3, max.overlaps = 50) +
    labs(
      title = "Z vs -log10(p-value)",
      x = "Z-score", y = "-log10(p-value)"
    ) +
    theme_classic(base_size = 13) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
      plot.background = element_rect(fill = "white", color = NA)
    )
  ggsave(file.path(opt$outdir, "volcano_plot.png"), p2,
         width = opt$width, height = opt$height, dpi = opt$dpi
  )
  cat("火山图:", file.path(opt$outdir, "volcano_plot.png"), "\n")
}, error = function(e) cat("火山图生成失败:", conditionMessage(e), "\n"))

# =============================================================================
# 完成
# =============================================================================
cat("\n=== 分析完成 ===\n")
cat("输出文件:\n")
for (f in list.files(opt$outdir, full.names = FALSE)) {
  cat("  ", f, "\n")
}
