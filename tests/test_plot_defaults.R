# 当前 ggplot2 的离散配色选项要求返回 scale；同时核对 >10 类仍来自 NPG。
suppressPackageStartupMessages(library(ggplot2))
source("src/R/00_runtime.R")
d <- data.frame(label = factor(seq_len(15)), value = 1)
p <- ggplot(d, aes(label, value, fill = label, color = label)) + geom_col()
b <- ggplot_build(p)
expected <- unname(npg_palette(levels(d$label)))
stopifnot(identical(b$data[[1]]$fill, expected), identical(b$data[[1]]$colour, expected),
  identical(expected[1:10], HONGF_NPG), !anyDuplicated(expected))
cat("PLOT_DEFAULTS_OK: default fill/color; 15 NPG categories\n")
source("src/R/03_quality_control.R")
for (classes in list(c("singlet", "doublet", "singlet"), rep("singlet", 3))) {
  md <- data.frame(orig.ident = "sample", doublet_class = classes)
  built <- ggplot_build(plot_qc_doublet_ratio(md))
  stopifnot(identical(built$plot$scales$get_scales("fill")$map(c("singlet", "doublet")),
    c("#4DBBD5", "#E64B35")))
  if (all(classes == "singlet")) stopifnot(all(built$data[[1]]$fill == "#4DBBD5"))
}
cat("QC_CLASS_COLORS_OK: singlet stays blue before and after filtering\n")

# 微型对象仅用于布局/数值回归；不是生物学分析数据。
suppressPackageStartupMessages(library(Seurat))
source("src/R/13_marker_dotplot.R")
counts <- matrix(c(0, 2, 1, 3, 0, 1, 5, 2, 0, 4, 3, 2, 1, 0, 3, 0, 4, 2),
  nrow = 3, dimnames = list(c("GeneA", "GeneB", "GeneC"), paste0("Cell", 1:6)))
obj <- NormalizeData(CreateSeuratObject(counts), verbose = FALSE)
obj$annotation <- factor(c("B", "A", "B", "A", "B", "B"), levels = c("B", "A"))
ref <- DotPlot(obj, features = c("GeneC", "GeneA"), group.by = "annotation")$data
panel <- hongf_marker_dotplot(obj, list(one = c("GeneC", "GeneA"), two = "GeneC"), group.by = "annotation")
stopifnot(inherits(panel, "patchwork"),
  isTRUE(all.equal(panel[[1]]$data$pct.exp, ref$pct.exp)),
  isTRUE(all.equal(panel[[1]]$data$avg.exp.scaled, ref$avg.exp.scaled)),
  identical(panel[[2]]$data$n_cells, c(4L, 2L)),
  identical(panel[[1]]$scales$get_scales("y")$limits, panel[[2]]$scales$get_scales("y")$limits),
  inherits(panel[[1]]$theme$panel.grid, "element_blank"),
  identical(panel[[1]]$theme$legend.position, "bottom"),
  identical(panel[[1]]$scales$get_scales("x")$limits, c("GeneC", "GeneA")),
  inherits(panel[[1]]$facet, "FacetNull"))
# 强制排版，检验底部图例和两侧面板能同时渲染。
stopifnot(inherits(patchwork::patchworkGrob(panel), "gtable"))
cat("REFERENCE_DOTPLOT_OK: Seurat values preserved; counts aligned; one panel; bottom legends\n")

# 核对连续图元的默认色阶端点，显式 scale 仍由调用者控制。
continuous <- data.frame(x = 1:2, y = 1:2, value = c(0, 1))
continuous_plot <- ggplot(continuous, aes(x, y, colour = value, fill = value)) + geom_point(shape = 21)
continuous_built <- ggplot_build(continuous_plot)$data[[1]]
stopifnot(identical(continuous_built$colour, HONGF_NPG_SEQUENTIAL[c(1, 4)]),
          identical(continuous_built$fill, HONGF_NPG_SEQUENTIAL[c(1, 4)]))
