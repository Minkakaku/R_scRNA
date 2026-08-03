# ==============================================================================
# Seurat + Harmony 自动化分析流程
# ==============================================================================
#
# 依赖安装:
#   install.packages(c("Seurat", "harmony", "ggplot2", "patchwork"))
#
# ---- 使用示例 ----
#
# source("seurat_harmony_pipeline.R")
#
# # 示例 1: 运行流程 (自动选PC + Harmony + 多resolution展示)
# res <- run_seurat_harmony(epi)
# # >>> 查看 clustering_results/dimplot_all_resolutions.png，选择合适的 resolution
#
# # 示例 2: 确定最终聚类结果
# obj_final <- finalize_clusters(res$obj, resolution = 0.5)
# DimPlot(obj_final, label = TRUE)
#
# # 示例 3: 自定义参数
# res <- run_seurat_harmony(
#   epi,
#   group_by = "sample",              # Harmony 批次变量
#   resolution_range = seq(0.2, 1.5, by = 0.1),  # DimPlot 输出范围
#   max_pcs = 50,                     # 最大 PC 数
#   min_pcs = 10,                     # 最少 PC 数
#   output_dir = "my_clustering_results",  # PNG 输出目录
#   plot_width = 10,                  # PNG 宽度 (英寸)
#   plot_height = 8                   # PNG 高度 (英寸)
# )
#
# # 示例 4: 手动指定 resolution (跳过多resolution展示)
# res <- run_seurat_harmony(epi, resolution = 0.7)
# DimPlot(res$obj)
#
# # 示例 5: 不跑 Harmony (单样本无需批次校正)
# res <- run_seurat_harmony(epi, run_harmony = FALSE)
#
# # ---- 输出文件说明 ----
# # clustering_results/
# #   ├── elbow_plot.png                    # Elbow Plot (标注选定PC)
# #   └── dimplot_all_resolutions.png       # 所有 resolution 拼合 UMAP
#
# ==============================================================================

#' Seurat + Harmony 自动化分析流程
#'
#' 自动完成 NormalizeData -> FindVariableFeatures -> ScaleData -> RunPCA ->
#' RunHarmony -> 自动选PC数 -> RunUMAP -> FindNeighbors ->
#' 多 resolution 聚类展示
#'
#' @param obj Seurat 对象
#' @param group_by Harmony 批次校正的分组变量，默认 "sample"
#' @param harmony_max_iter Harmony 最大迭代次数，默认 20
#' @param resolution 若指定数值则直接用该 resolution 聚类；
#'   默认 "auto" 表示展示多个 resolution 的结果供用户选择
#' @param resolution_range 多 resolution 展示的范围，默认 seq(0.1, 1.5, by = 0.1)
#' @param max_pcs 考虑的最大 PC 数量，默认 50
#' @param min_pcs 最少使用的 PC 数量，默认 5
#' @param reduction_name 使用的降维名称，默认 "harmony"；如不需要 Harmony，设为 "pca"
#' @param run_harmony 是否运行 Harmony 批次校正，默认 TRUE
#' @param plot_elbow 是否绘制 Elbow Plot 并标注选定的 PC 数，默认 TRUE
#' @param output_dir 输出 PNG 图片的目录，默认 "clustering_results"
#' @param plot_width PNG 图片宽度（英寸），默认 8
#' @param plot_height PNG 图片高度（英寸），默认 6
#' @return 一个 list，包含:
#'   - obj: 完成降维和建邻域图的 Seurat 对象 (尚未设定最终聚类)
#'   - optimal_pc: 自动选定的 PC 数
#'   - resolution_range: 展示过的 resolution 范围 (auto 模式)
#'
run_seurat_harmony <- function(obj,
                                group_by = "sample",
                                harmony_max_iter = 20,
                                resolution = "auto",
                                resolution_range = seq(0.1, 1.5, by = 0.1),
                                max_pcs = 50,
                                min_pcs = 5,
                                reduction_name = "harmony",
                                run_harmony = TRUE,
                                plot_elbow = TRUE,
                                output_dir = "clustering_results",
                                plot_width = 8,
                                plot_height = 6) {
  
  # 创建输出目录
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    message(">>> 创建输出目录: ", output_dir)
  }
  
  # ---- 1. 标准预处理 ----
  message(">>> Step 1: NormalizeData -> FindVariableFeatures -> ScaleData -> RunPCA")
  obj <- NormalizeData(obj) %>%
    FindVariableFeatures() %>%
    ScaleData() %>%
    RunPCA(npcs = max_pcs, verbose = FALSE)
  
  # ---- 2. Harmony 批次校正 (可选) ----
  if (run_harmony) {
    message(">>> Step 2: RunHarmony (group.by = '", group_by, "')")
    obj <- RunHarmony(obj, group.by.vars = group_by, max_iter = harmony_max_iter)
  } else {
    reduction_name <- "pca"
    message(">>> Step 2: 跳过 Harmony，使用 PCA 降维")
  }
  
  # ---- 3. 自动选择 PC 数 (Elbow 拐点法) ----
  message(">>> Step 3: 自动确定最佳 PC 数量...")
  sdev <- Stdev(obj, reduction = reduction_name)
  n_pcs <- length(sdev)
  
  var_explained <- sdev^2
  cumulative_var <- cumsum(var_explained) / sum(var_explained)
  delta_sdev <- -diff(sdev)
  
  # 最大曲率法
  x <- seq_len(n_pcs)
  y <- sdev
  x1 <- x[1]; y1 <- y[1]
  x2 <- x[n_pcs]; y2 <- y[n_pcs]
  line_len <- sqrt((x2 - x1)^2 + (y2 - y1)^2)
  distances <- abs((y2 - y1) * x - (x2 - x1) * y + x2 * y1 - y2 * x1) / line_len
  optimal_pc <- which.max(distances)
  optimal_pc <- max(min_pcs, optimal_pc)
  optimal_pc <- min(optimal_pc, n_pcs)
  
  message(">>> 自动选定 PC 数: ", optimal_pc,
          " (累计方差贡献: ", round(cumulative_var[optimal_pc] * 100, 1), "%)")
  
  # ---- 4. Elbow Plot 可视化 ----
  if (plot_elbow) {
    p <- ElbowPlot(obj, ndims = max_pcs, reduction = reduction_name) +
      geom_vline(xintercept = optimal_pc, linetype = "dashed", color = "red", linewidth = 0.8) +
      annotate("text", x = optimal_pc + 1.5, y = max(sdev) * 0.9,
               label = paste0("PC = ", optimal_pc),
               color = "red", fontface = "bold", hjust = 0) +
      ggtitle(paste0("Elbow Plot (auto-selected PC = ", optimal_pc, ")"))
    
    elbow_path <- file.path(output_dir, "elbow_plot.png")
    ggsave(elbow_path, p, width = plot_width, height = plot_height, dpi = 300)
    print(p)
    message(">>> Elbow Plot 已保存: ", elbow_path)
  }
  
  # ---- 5. UMAP 降维 ----
  message(">>> Step 4: RunUMAP (dims = 1:", optimal_pc, ")")
  obj <- RunUMAP(obj, reduction = reduction_name, dims = 1:optimal_pc)
  
  # ---- 6. FindNeighbors ----
  message(">>> Step 5: FindNeighbors (dims = 1:", optimal_pc, ")")
  obj <- FindNeighbors(obj, reduction = reduction_name, dims = 1:optimal_pc)
  
  # ---- 7. 聚类 ----
  if (is.character(resolution) && resolution == "auto") {
    
    # === 多 resolution 展示模式 ===
    message("\n>>> ==========================================")
    message(">>> 生成各 resolution 的聚类 DimPlot 拼合图...")
    message(">>> ==========================================")
    
    plot_list <- list()
    for (res in resolution_range) {
      obj_tmp <- FindClusters(obj, resolution = res, verbose = FALSE)
      p_dim <- DimPlot(obj_tmp, reduction = "umap", label = TRUE, repel = TRUE) +
        ggtitle(paste0("Res = ", res, 
                       " (k=", length(levels(Idents(obj_tmp))), ")")) +
        theme(plot.title = element_text(face = "bold", size = 10),
              legend.position = "none")
      plot_list[[length(plot_list) + 1]] <- p_dim
    }
    
    # 自动计算布局: 尽量接近正方形
    n_plots <- length(plot_list)
    n_col <- ceiling(sqrt(n_plots))
    n_row <- ceiling(n_plots / n_col)
    
    combined_dimplot <- patchwork::wrap_plots(plot_list, ncol = n_col) +
      patchwork::plot_annotation(
        title = "DimPlot across Resolutions",
        theme = theme(plot.title = element_text(face = "bold", size = 14))
      )
    
    combined_path <- file.path(output_dir, "dimplot_all_resolutions.png")
    ggsave(combined_path, combined_dimplot,
           width = plot_width * n_col * 0.6,
           height = plot_height * n_row * 0.6,
           dpi = 300, limitsize = FALSE)
    print(combined_dimplot)
    message(">>> 所有 resolution DimPlot 拼合图已保存: ", combined_path)
    
    message("\n>>> ==========================================")
    message(">>> 请查看拼合图，选择合适的 resolution 后调用:")
    message(">>>   obj_final <- finalize_clusters(res$obj, resolution = <你选择的值>)")
    message(">>> ==========================================")
    
  } else {
    # === 手动指定 resolution ===
    message(">>> Step 6: FindClusters (resolution = ", resolution, ")")
    obj <- FindClusters(obj, resolution = resolution)
    
    # 保存 DimPlot
    p_dim <- DimPlot(obj, reduction = "umap", label = TRUE, repel = TRUE) +
      ggtitle(paste0("Resolution = ", resolution, 
                     " (clusters = ", length(levels(Idents(obj))), ")")) +
      theme(plot.title = element_text(face = "bold", size = 14))
    
    png_path <- file.path(output_dir, paste0("dimplot_res_", gsub("\\.", "_", as.character(resolution)), ".png"))
    ggsave(png_path, p_dim, width = plot_width, height = plot_height, dpi = 300)
    
    message(">>> 完成! 共检测到 ", length(levels(Idents(obj))), " 个 clusters")
  }
  
  # ---- 输出总结 ----
  message("\n>>> ==========================================")
  message(">>> 分析完成! 结果总结:")
  message(">>> ==========================================")
  message(">>> 自动选定 PC 数: ", optimal_pc)
  message(">>> 所有 PNG 图片保存在: ", normalizePath(output_dir, mustWork = FALSE))
  message(">>> ==========================================")
  
  return(list(
    obj = obj,
    optimal_pc = optimal_pc,
    resolution_range = if (is.character(resolution) && resolution == "auto") resolution_range else resolution
  ))
}


# ======================================================================
# finalize_clusters: 用户确认 resolution 后设定最终聚类
# ======================================================================

#' 确定最终聚类结果
#'
#' 在查看 run_seurat_harmony() 输出的多 resolution DimPlot 后，
#' 调用此函数以选定的 resolution 进行最终聚类。
#'
#' @param obj run_seurat_harmony() 返回的 Seurat 对象 (res$obj)
#' @param resolution 选定的 resolution 值
#' @param output_dir 输出目录，默认 "clustering_results"
#' @param plot_width PNG 宽度，默认 8
#' @param plot_height PNG 高度，默认 6
#' @return 聚类完成的 Seurat 对象
#'
finalize_clusters <- function(obj,
                               resolution,
                               output_dir = "clustering_results",
                               plot_width = 8,
                               plot_height = 6) {
  
  message(">>> 使用 resolution = ", resolution, " 进行最终聚类...")
  obj <- FindClusters(obj, resolution = resolution)
  
  n_clusters <- length(levels(Idents(obj)))
  message(">>> 最终聚类完成! 共 ", n_clusters, " 个 clusters")
  
  # 保存最终 DimPlot
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  p_final <- DimPlot(obj, reduction = "umap", label = TRUE, repel = TRUE) +
    ggtitle(paste0("Final Clustering (resolution = ", resolution,
                   ", clusters = ", n_clusters, ")")) +
    theme(plot.title = element_text(face = "bold", size = 14))
  
  final_path <- file.path(output_dir, "dimplot_final.png")
  ggsave(final_path, p_final, width = plot_width, height = plot_height, dpi = 300)
  print(p_final)
  message(">>> 最终聚类 DimPlot 已保存: ", final_path)
  
  return(obj)
}
