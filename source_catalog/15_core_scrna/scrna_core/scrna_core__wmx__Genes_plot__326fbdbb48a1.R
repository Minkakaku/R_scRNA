library('igraph', lib.loc = '/mnt/NexData_file/home/wmx/envs/scRNA/lib/R/library/')
.libPaths("/mnt/NexData_file/share/scRNA-tools/envs/pip_monocle2/lib/R/library/")

suppressPackageStartupMessages({
  library(monocle)
  library(ggplot2)
  library(cowplot)
  library(dplyr)
  library(reshape2)
  library(getopt)
})

options(stringsAsFactors = FALSE, warn = -1)

# 绘制特定基因的变化图
# T
marker_list <- list(
  Stemness = c(
    "CCR7", "IL7R", "TCF7", "SELL"
  ),
  
  Exhaustion = c(
    "TOX", "CXCL13", "PDCD1", "LAG3", "CD101", "HAVCR2"
  ),
  
  cytotoxic = c(
    "PRF1", "IFNG", "GNLY", "GZMB", "GZMA", "GZMK", "GZMH"
  ),
  
  Inhibition = c(
    "PDCD1", "CTLA4", "LAG3", "HAVCR2", "TIGIT", "VSIR"
  )
)

# B/Plasma
marker_list_B <- list(
activation_genes = c(
  "CD27", "CD40", "EBI3", "MIF", "BATF"
),
csr_mechanism_genes = c(
  "APEX1", "APEX2", "XRCC5", "XRCC6", "POLD2", "POLE3"
),
csr_interactor_genes = c(
  "NPM1", "SERBP1"
))


# =========================
# functions
# =========================

safe_name <- function(x) {
  x <- gsub("[/\\\\:;,*?\"<>| ]+", "_", x)
  x <- gsub("_+", "_", x)
  x
}

save_plot <- function(p, prefix, width = 6.5, height = 4.8) {
  ggsave(
    filename = file.path(outdir, paste0(prefix, ".pdf")),
    plot = p,
    width = width,
    height = height,
    limitsize = FALSE
  )
  
  ggsave(
    filename = file.path(outdir, paste0(prefix, ".png")),
    plot = p,
    width = width,
    height = height,
    dpi = 300,
    limitsize = FALSE
  )
}


# =========================
# 作图函数
# =========================

plot_one_set <- function(cds, genes, set_name,outdir) {
  
  genes_use <- intersect(genes, rownames(cds))
  genes_miss <- setdiff(genes, rownames(cds))
  
  message("========================================")
  message("[", set_name, "]")
  message("Genes used: ", paste(genes_use, collapse = ", "))
  
  if (length(genes_miss) > 0) {
    message("Genes missing: ", paste(genes_miss, collapse = ", "))
  }
  
  if (length(genes_use) == 0) {
    message("[Skip] No genes found for ", set_name)
    return(NULL)
  }
  
  cds_sub <- cds[genes_use, ]
  
  n_gene <- length(genes_use)
  ncol_use <- ifelse(n_gene == 1, 1, 2)
  nrow_use <- ceiling(n_gene / ncol_use)
  
  # 不分组版本
  p1 <- plot_genes_in_pseudotime(
    cds_sub,
    color_by = "Pseudotime",
    min_expr = 0.1,
    nrow = nrow_use,
    ncol = ncol_use
  ) +
    ggtitle(set_name)
  
  ggsave(
    file.path(outdir, paste0(set_name, ".genes_in_pseudotime.pdf")),
    p1,
    width = 9,
    height = max(5, 2.5 * ceiling(length(genes_use) / 2)),
    limitsize = FALSE,
    create.dir = TRUE
  )
  
  ggsave(
    file.path(outdir, paste0(set_name, ".genes_in_pseudotime.png")),
    p1,
    width = 9,
    height = max(5, 2.5 * ceiling(length(genes_use) / 2)),
    dpi = 300,
    limitsize = FALSE,
    create.dir = TRUE
  )
  
}

# =========================
# 循环绘制每个 gene set
# =========================

# T cells
cds_files = c(
  "/mnt/NexData_file/home/wmx/Project/AE2511046228-02/shouhou12/Monocle_result/N_CD8T/Pseudotime.monocle.reRoot.rds",
  "/mnt/NexData_file/home/wmx/Project/AE2511046228-02/shouhou12/Monocle_result/N_CD4T/Pseudotime.monocle.reRoot.rds",
  "/mnt/NexData_file/home/wmx/Project/AE2511046228-02/shouhou12/Monocle_result/LARC_CD8T/Pseudotime.monocle.reRoot.rds",
  "/mnt/NexData_file/home/wmx/Project/AE2511046228-02/shouhou12/Monocle_result/LARC_CD4T/Pseudotime.monocle.reRoot.rds"
)

for (cds_file in cds_files){
  message("[Read CDS] ", cds_file)
  cds <- readRDS(cds_file)
  outdir <- gsub(x = cds_file,pattern = "Pseudotime.monocle.reRoot.rds",replacement = "genePlot/")
  meta <- pData(cds)
  dir.create(outdir, recursive = TRUE, showWarnings = TRUE)
  if (!"Pseudotime" %in% colnames(meta)) {
    stop("No Pseudotime column found in pData(cds). Please run orderCells first.")
  }

  for (set_name in names(marker_list)) {
    plot_one_set(
      cds = cds,
      genes = marker_list[[set_name]],
      set_name = set_name,
      outdir = outdir
   )
  }
  message("========================================")
  message("[Done] Results saved in: ", outdir)
  message("========================================")
}

# B and plasma cells
cds_files = c(
  "/mnt/NexData_file/home/wmx/Project/AE2511046228-02/shouhou12/Monocle_result/LARC_BandPlasma/Pseudotime.monocle.reRoot.rds",
  "/mnt/NexData_file/home/wmx/Project/AE2511046228-02/shouhou12/Monocle_result/N_BandPlasma/Pseudotime.monocle.reRoot.rds"
)

for (cds_file in cds_files){
  message("[Read CDS] ", cds_file)
  cds <- readRDS(cds_file)
  outdir <- gsub(x = cds_file,pattern = "Pseudotime.monocle.reRoot.rds",replacement = "/genePlot/")
  dir.create(outdir,showWarnings = F)
  meta <- pData(cds)
  
  if (!"Pseudotime" %in% colnames(meta)) {
    stop("No Pseudotime column found in pData(cds). Please run orderCells first.")
  }
  
  for (set_name in names(marker_list_B)) {
    plot_one_set(
      cds = cds,
      genes = marker_list_B[[set_name]],
      set_name = set_name,
      outdir = outdir
    )
  }
  message("========================================")
  message("[Done] Results saved in: ", outdir)
  message("========================================")
}
