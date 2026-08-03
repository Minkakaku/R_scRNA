#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
case <- if (length(args) >= 1) args[1] else stop("need case")
control <- if (length(args) >= 2) args[2] else stop("need control")
output_dir <- if (length(args) >= 3) args[3] else getwd()

report_dir <- file.path(output_dir, "Report", paste0(case, "_vs_", control))
dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)

deg_dir <- file.path(output_dir, "DE_Analysis", paste0(case, "_vs_", control, "_DEG"))
enrich_dir <- file.path(output_dir, "Enrichment", paste0(case, "_vs_", control, "_Enrichment"))
gsea_dir <- file.path(output_dir, "Enrichment", paste0(case, "_vs_", control, "_GSEA"))

deg_file <- file.path(deg_dir, paste0(case, "vs", control, "_DEG_inEachCelltype.xlsx"))
deg_num_file <- file.path(deg_dir, "DEG_Number.csv")
enrich_rds <- file.path(enrich_dir, "RDS", "All_Results.rds")

report_rmd <- file.path(report_dir, "report.Rmd")
report_html <- file.path(report_dir, "report.html")

lines <- c(
  "---",
  paste0("title: \"", case, " vs ", control, " DEG/Enrichment Report\""),
  "output:",
  "  html_document:",
  "    toc: true",
  "    toc_depth: 3",
  "    theme: cosmo",
  "    number_sections: true",
  "---",
  "",
  "```{r setup, include=FALSE}",
  "knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE)",
  "suppressPackageStartupMessages({",
  "  library(openxlsx)",
  "  library(dplyr)",
  "  library(knitr)",
  "})",
  paste0("case <- '", case, "'"),
  paste0("control <- '", control, "'"),
  paste0("deg_dir <- '", gsub("\\\\", "/", deg_dir), "'"),
  paste0("enrich_dir <- '", gsub("\\\\", "/", enrich_dir), "'"),
  paste0("gsea_dir <- '", gsub("\\\\", "/", gsea_dir), "'"),
  paste0("deg_file <- '", gsub("\\\\", "/", deg_file), "'"),
  paste0("deg_num_file <- '", gsub("\\\\", "/", deg_num_file), "'"),
  paste0("enrich_rds <- '", gsub("\\\\", "/", enrich_rds), "'"),
  "deg_df <- if (file.exists(deg_file)) openxlsx::read.xlsx(deg_file) else data.frame()",
  "deg_num <- if (file.exists(deg_num_file)) read.csv(deg_num_file) else data.frame()",
  "all_enrich <- if (file.exists(enrich_rds)) readRDS(enrich_rds) else NULL",
  "```",
  "",
  "# 1. 项目概述",
  "",
  paste0("本报告汇总 `", case, " vs ", control, "` 的差异表达、富集分析与 GSEA 结果。主分析逻辑沿用原始脚本，报告由独立脚本自动整理生成。"),
  "",
  "# 2. DEG 概览",
  "",
  "```{r}",
  "if (nrow(deg_num) > 0) knitr::kable(deg_num, caption = '各细胞类型上下调基因数量统计')",
  "```",
  "",
  "```{r, results='asis'}",
  "img_files <- c(",
  "  file.path(deg_dir, 'all_cluster_volcano.png'),",
  "  file.path(deg_dir, 'all_cluster_volcano_pct.png'),",
  "  file.path(deg_dir, 'Number_of_DEGs.png'),",
  "  file.path(deg_dir, 'DEG_top5_Heatmap.png')",
  ")",
  "for (img in img_files) if (file.exists(img)) cat(sprintf('![](%s)\\n\\n', img))",
  "```",
  "",
  "# 3. 富集分析",
  "",
  "```{r}",
  "if (!is.null(all_enrich)) {",
  "  ct_names <- names(all_enrich)",
  "  enrich_summary <- data.frame(Celltype = character(), Database = character(), Terms = integer())",
  "  for (ct in ct_names) {",
  "    for (db in names(all_enrich[[ct]])) {",
  "      obj <- all_enrich[[ct]][[db]]",
  "      if (!is.null(obj)) enrich_summary <- rbind(enrich_summary, data.frame(Celltype = ct, Database = db, Terms = nrow(as.data.frame(obj))))",
  "    }",
  "  }",
  "  if (nrow(enrich_summary) > 0) knitr::kable(enrich_summary, caption = '富集结果条目数汇总')",
  "}",
  "```",
  "",
  "```{r, results='asis'}",
  "plot_dir <- file.path(enrich_dir, 'Enrichment_Plots')",
  "plot_files <- c(",
  "  file.path(plot_dir, 'Network_Plot.png'),",
  "  file.path(plot_dir, 'Comparison.png'),",
  "  file.path(plot_dir, 'Lollipop.png'),",
  "  file.path(plot_dir, 'Complex_combined_bar.png')",
  ")",
  "for (img in plot_files) if (file.exists(img)) cat(sprintf('![](%s)\\n\\n', img))",
  "```",
  "",
  "# 4. GSEA",
  "",
  "```{r, results='asis'}",
  "gsea_plot_root <- file.path(gsea_dir, 'GSEA_Plots', 'KEGG')",
  "if (dir.exists(gsea_plot_root)) {",
  "  ct_dirs <- list.dirs(gsea_plot_root, recursive = FALSE, full.names = TRUE)",
  "  for (ct_dir in ct_dirs) {",
  "    cat(sprintf('## %s\\n\\n', basename(ct_dir)))",
  "    imgs <- file.path(ct_dir, c('1_Barplot.png', '2_RidgePlot.png', '3_GseaVis.png', '4_GSEAPlot_Combined.png'))",
  "    for (img in imgs) if (file.exists(img)) cat(sprintf('![](%s)\\n\\n', img))",
  "  }",
  "}",
  "```",
  "",
  "# 5. 结论建议",
  "",
  "1. 建议优先关注差异基因数量较多且富集方向一致的细胞类型。",
  "2. 建议结合 KEGG、GO 与 GSEA 共同支持的通路筛选关键机制。",
  "3. 建议从高变化幅度基因及核心通路交集基因中挑选后续验证候选。"
)

writeLines(lines, report_rmd, useBytes = TRUE)
rmarkdown::render(report_rmd, output_file = basename(report_html), output_dir = dirname(report_html), quiet = TRUE)
