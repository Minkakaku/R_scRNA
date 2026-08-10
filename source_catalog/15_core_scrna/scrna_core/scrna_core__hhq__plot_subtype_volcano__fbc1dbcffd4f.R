library(ggplot2)
library(ggrepel)
library(dplyr)

#' 绘制细胞亚类火山图（Volcano Plot）
#'
#' 基于 FindAllMarkers 结果，按细胞亚类分面绘制 pct_diff vs log2FC 的火山图。
#'
#' @param subtype_marker FindAllMarkers 输出的 data.frame，需包含
#'   gene, avg_log2FC, p_val_adj, pct.1, pct.2, cluster 列。
#' @param log2fc_thresh log2FC 显著性阈值，默认 1。
#' @param pval_thresh  调整后 p 值阈值，默认 0.05。
#' @param max_overlaps  ggrepel 标签最大重叠数，默认 20。
#' @param label_size    基因标签字号，默认 2.5。
#' @param facet_nrow    分面行数，默认 1（所有亚类排成一行）。
#' @param save_path     保存路径。若为 NULL 则不保存。
#' @param width         保存图片宽度（英寸），默认 16。
#' @param height        保存图片高度（英寸），默认 8。
#' @param dpi           保存分辨率，默认 300。
#'
#' @return 返回 ggplot 对象（不可见）。
plot_subtype_volcano <- function(subtype_marker,
                                  log2fc_thresh = 1,
                                  pval_thresh   = 0.05,
                                  max_overlaps  = 20,
                                  label_size    = 2.5,
                                  facet_nrow    = 1,
                                  save_path     = NULL,
                                  width          = 16,
                                  height         = 8,
                                  dpi            = 300) {

  # ---------- 数据预处理 ----------
  markers <- subtype_marker %>%
    mutate(
      pct_diff = pct.1 - pct.2,
      significance = case_when(
        abs(avg_log2FC) > log2fc_thresh & p_val_adj < pval_thresh ~ "Significant",
        TRUE ~ "Not Significant"
      ),
      direction = case_when(
        avg_log2FC >  log2fc_thresh & p_val_adj < pval_thresh ~ "Up-regulated",
        avg_log2FC < -log2fc_thresh & p_val_adj < pval_thresh ~ "Down-regulated",
        TRUE ~ "Not Significant"
      ),
      is_significant = abs(avg_log2FC) > log2fc_thresh & p_val_adj < pval_thresh
    )

  markers_to_label <- markers %>%
    filter(is_significant)

  # ---------- 亚类颜色 ----------
  clusters   <- unique(markers$cluster)
  n_clusters <- length(clusters)

  if (n_clusters <= 8) {
    cluster_colors <- RColorBrewer::brewer.pal(max(3, n_clusters), "Set2")[1:n_clusters]
  } else if (n_clusters <= 12) {
    cluster_colors <- RColorBrewer::brewer.pal(n_clusters, "Set3")
  } else {
    cluster_colors <- scales::hue_pal()(n_clusters)
  }
  names(cluster_colors) <- clusters

  # ---------- 绑图 ----------
  p <- ggplot(markers, aes(x = pct_diff, y = avg_log2FC)) +
    # 非显著点（灰色）
    geom_point(
      data  = markers %>% filter(!is_significant),
      color = "#BDC3C7", alpha = 0.3, size = 1
    ) +
    # 显著点（按亚类着色）
    geom_point(
      data = markers %>% filter(is_significant),
      aes(color = cluster), alpha = 0.8, size = 2
    ) +
    scale_color_manual(values = cluster_colors, name = "Cell Subtype") +
    # 参考线
    geom_hline(yintercept = c(-log2fc_thresh, log2fc_thresh),
               linetype = "dashed", color = "#7F8C8D", linewidth = 0.5) +
    geom_hline(yintercept = 0, linetype = "solid", color = "#2C3E50", linewidth = 0.3) +
    geom_vline(xintercept = 0, linetype = "solid", color = "#2C3E50", linewidth = 0.3) +
    # 基因标签
    geom_text_repel(
      data               = markers_to_label,
      aes(label = gene),
      size               = label_size,
      max.overlaps       = max_overlaps,
      box.padding        = 0.3,
      point.padding      = 0.2,
      segment.color      = "#7F8C8D",
      segment.size       = 0.3,
      min.segment.length = 0.2,
      force              = 2,
      force_pull         = 1
    ) +
    # 分面
    facet_wrap(~ cluster, scales = "free", nrow = facet_nrow) +
    labs(
      x = "Percentage Difference (pct.1 - pct.2)",
      y = expression(Log[2] ~ Fold ~ Change)
    ) +
    theme_bw() +
    theme(
      plot.title        = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle     = element_text(size = 10, hjust = 0.5, color = "#7F8C8D"),
      strip.text        = element_text(size = 10, face = "bold"),
      strip.background  = element_blank(),
      axis.title        = element_text(size = 10, face = "bold"),
      axis.text         = element_text(size = 8),
      legend.position   = "none",
      legend.title      = element_text(size = 9, face = "bold"),
      legend.text       = element_text(size = 8),
      panel.grid.minor  = element_blank(),
      panel.grid.major  = element_blank()
    )

  # ---------- 显示 & 保存 ----------
  print(p)

  if (!is.null(save_path)) {
    ggsave(
      filename = save_path,
      plot     = p,
      width    = width,
      height   = height,
      dpi      = dpi
    )
    message("Plot saved to: ", save_path)
  }

  invisible(p)
}


# ============================================================
# 使用示例
# ============================================================
# plot_subtype_volcano(subtype_marker)
#
# plot_subtype_volcano(
#   subtype_marker,
#   log2fc_thresh = 0.5,
#   save_path     = "CD8T_volcano_plot.png",
#   facet_nrow    = 2,
#   width         = 20,
#   height        = 12
# )
