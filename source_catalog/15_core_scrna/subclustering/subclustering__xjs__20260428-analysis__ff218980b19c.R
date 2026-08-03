# 20260428 ROI -------------

roi <- read.csv("ROI.csv")

roi$cell_id_new <- roi$cell_id
roi$cell_id <- paste0("F2410395-A3_F2410396-A3-1_F2410396-A3-2_", roi$cell_id_new)
roi$cell_id_new <- NULL
colnames(roi) <- c("cell_id", "ROI")
roi$ROI[roi$ROI %in% c("Non_classic_1", "Non_classic_2")] <- "Non_classic"
roi <- roi %>% column_to_rownames("cell_id")

all_data <- AddMetaData(all_data, metadata = roi["ROI"])

all_data@meta.data$ROI[is.na(all_data@meta.data$ROI)] <- "Other"
table(all_data$ROI, useNA = "always")


saveRDS(all_data, "all_data_ROI.rds")
# all_data <- readRDS("all_data_ROI.rds")


# ROI 区域高表达基因 FindAllMarkers ---------
library(tidyverse)
library(Seurat)
library(patchwork)
library(stringr)
library("clusterProfiler")
library("org.Hs.eg.db")
library(RColorBrewer)
library(EnhancedVolcano)
library(dplyr)


marker_gene_outdir <- "20260428-FindAllMarkers"
dir.create(marker_gene_outdir)

all_data_filt <- all_data[, all_data@meta.data$ROI != "Other"]

Idents(all_data_filt) <- "ROI"
DefaultAssay(all_data_filt) <- "RNA"



### enrichment --------
library(tidyverse)
library(Seurat)
library(patchwork)
library(stringr)
library("clusterProfiler")
library("org.Hs.eg.db")
library(RColorBrewer)
library(EnhancedVolcano)
library(dplyr)

ROIs <- unique(all_data_filt@meta.data$ROI)


diff_genes <- FindAllMarkers(
  JoinLayers(all_data_filt),
  test.use = "wilcox", # 使用Wilcoxon秩和检验
  logfc.threshold = 0.25, # 设置log2折叠变化阈值
  min.pct = 0.1, # 基因至少在10%的细胞中表达
  only.pos = TRUE, # 包含上调和下调基因
  return.thresh = 0.05 # 只返回调整p值<0.05的基因
)

significant_genes <- diff_genes %>%
  group_by(gene) %>%
  filter(all(p_val_adj < 0.05)) %>%
  ungroup()

write.csv(significant_genes, file.path(marker_gene_outdir, "ROI diffmarker.csv"))



### top marker dotplot --------
all_top5_gene <-
  diffmarkers %>%
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  dplyr::filter(pct.1 > 0.25) %>%
  dplyr::filter(abs(avg_log2FC) > 1) %>% ##
  arrange(desc(avg_log2FC), .by_group = TRUE) %>%
  slice_head(n = 5) %>%
  with(structure(
    gene,
    names = as.character(cluster)
  ))
all_top5_gene_unique <- all_top5_gene[!duplicated(all_top5_gene)]

p_dot_top5 <-
  DotPlot(
    all_data_filt,
    features = all_top5_gene_unique,
    dot.scale = 10
  ) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

ggsave(
  p_dot_top5,
  filename = file.path(marker_gene_outdir, "top 5 gene dotplot.pdf"),
  height = 6,
  width = 12
)




clusters <- levels(significant_genes$cluster) # 使用levels获取因子水平

for (cluster_name in clusters) {
  cat("正在处理:", cluster_name, "\n")

  # 1. 提取当前cluster的基因
  current_genes <- significant_genes[significant_genes$cluster == cluster_name, ]$gene

  # 如果基因列表为空，跳过
  if (length(current_genes) == 0) {
    cat(cluster_name, "没有基因，跳过\n")
    next
  }

  # 创建当前cluster的输出子目录
  cluster_dir <- file.path(marker_gene_outdir, cluster_name)
  if (!dir.exists(cluster_dir)) {
    dir.create(cluster_dir, recursive = TRUE)
  }

  # 2. GO富集分析
  cat("  GO富集分析...\n")
  enrichGO_up_ls <- enrichGO(
    gene = current_genes,
    OrgDb = org.Hs.eg.db,
    ont = "ALL",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    keyType = "SYMBOL",
    minGSSize = 10,
    maxGSSize = 500
  )

  # 保存GO结果
  if (!is.null(enrichGO_up_ls) && nrow(enrichGO_up_ls@result) > 0) {
    write.csv(enrichGO_up_ls@result,
      file = file.path(cluster_dir, paste0(cluster_name, "_GO_enrichment_table.csv")),
      row.names = FALSE
    )

    # GO点图
    enrichGO_up_dot <- dotplot(
      enrichGO_up_ls,
      showCategory = 10,
      label_format = 50
    ) +
      ggtitle(paste(cluster_name, "GO Enrichment"))

    ggsave(
      plot = enrichGO_up_dot,
      filename = file.path(cluster_dir, paste0(cluster_name, "_GO_dotplot.pdf")),
      height = 6,
      width = 8
    )

    # GO条形图
    enrichGO_up_bar <- barplot(
      enrichGO_up_ls,
      showCategory = 10,
      label_format = 50
    ) +
      ggtitle(paste(cluster_name, "GO Enrichment"))

    ggsave(
      plot = enrichGO_up_bar,
      filename = file.path(cluster_dir, paste0(cluster_name, "_GO_barplot.pdf")),
      height = 6,
      width = 8
    )

    cat("  GO分析完成，找到", nrow(enrichGO_up_ls@result), "个富集项\n")
  } else {
    cat("  GO分析无显著结果\n")
  }

  # 3. KEGG富集分析
  cat("  KEGG富集分析...\n")

  # 转换基因ID
  de_up_genes_id_for_enrichment <- bitr(
    current_genes,
    OrgDb = "org.Hs.eg.db",
    fromType = "SYMBOL",
    toType = "ENTREZID"
  )

  if (nrow(de_up_genes_id_for_enrichment) > 0) {
    entrez_ids <- de_up_genes_id_for_enrichment$ENTREZID

    enrichKEGG_up_ls <- enrichKEGG(
      gene = entrez_ids,
      organism = "hsa",
      keyType = "ncbi-geneid",
      pvalueCutoff = 0.05,
      pAdjustMethod = "BH",
      minGSSize = 10,
      maxGSSize = 500,
      qvalueCutoff = 0.2,
      use_internal_data = FALSE
    )

    # 保存KEGG结果
    if (!is.null(enrichKEGG_up_ls) && nrow(enrichKEGG_up_ls@result) > 0) {
      # 转换为可读的基因名
      enrichKEGG_up_ls <- setReadable(enrichKEGG_up_ls, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

      write.csv(enrichKEGG_up_ls@result,
        file = file.path(cluster_dir, paste0(cluster_name, "_KEGG_enrichment_table.csv")),
        row.names = FALSE
      )

      # KEGG点图
      enrichKEGG_up_dot <- dotplot(
        enrichKEGG_up_ls,
        showCategory = 10,
        label_format = 50
      ) +
        ggtitle(paste(cluster_name, "KEGG Enrichment"))

      ggsave(
        plot = enrichKEGG_up_dot,
        filename = file.path(cluster_dir, paste0(cluster_name, "_KEGG_dotplot.pdf")),
        height = 6,
        width = 8
      )

      # KEGG条形图
      enrichKEGG_up_bar <- barplot(
        enrichKEGG_up_ls,
        showCategory = 10,
        label_format = 50
      ) +
        ggtitle(paste(cluster_name, "KEGG Enrichment"))

      ggsave(
        plot = enrichKEGG_up_bar,
        filename = file.path(cluster_dir, paste0(cluster_name, "_KEGG_barplot.pdf")),
        height = 6,
        width = 8
      )

      cat("  KEGG分析完成，找到", nrow(enrichKEGG_up_ls@result), "个富集项\n")
    } else {
      cat("  KEGG分析无显著结果\n")
    }
  } else {
    cat("  无法转换基因ID，跳过KEGG分析\n")
  }

  cat("完成处理:", cluster_name, "\n")
  cat("=", rep("=", 50), "\n", sep = "")
}

cat("所有cluster分析完成！\n")






### 三组富集分析 ----------
library(tidyverse)
library(Seurat)
library(patchwork)
library(stringr)
library("clusterProfiler")
library("org.Hs.eg.db")
library(RColorBrewer)
library(EnhancedVolcano)
library(dplyr)
library(clusterProfiler)
library(ggplot2)



diff_outdir <- file.path(marker_gene_outdir, "diff analysis outdir")
dir.create(diff_outdir)


all_data <- readRDS("all_data_ROI.rds")
all_data_filt <- all_data[, all_data@meta.data$ROI != "Other"]

Idents(all_data_filt) <- "ROI"
DefaultAssay(all_data_filt) <- "RNA"




diff_genes <- FindAllMarkers(
  JoinLayers(all_data_filt),
  test.use = "wilcox", # 使用Wilcoxon秩和检验
  logfc.threshold = 0.25, # 设置log2折叠变化阈值
  min.pct = 0.1, # 基因至少在10%的细胞中表达
  only.pos = TRUE, # 包含上调和下调基因
  return.thresh = 0.05 # 只返回调整p值<0.05的基因
)

significant_genes <- diff_genes %>%
  group_by(gene) %>%
  filter(all(p_val_adj < 0.05)) %>%
  ungroup()

write.csv(significant_genes, file.path(diff_outdir, "group diffmarker.csv"))


table(significant_genes$cluster)
ids <- bitr(significant_genes$gene, "SYMBOL", "ENTREZID", "org.Hs.eg.db")
significant_genes <- merge(significant_genes, ids, by.x = "gene", by.y = "SYMBOL")

# significant_genes$group[significant_genes$avg_log2FC > 0] <- "upregulated"
# significant_genes$group[significant_genes$avg_log2FC < 0] <- "downregulated"
# significant_genes$cluster_group <- paste0(significant_genes$group, "_", significant_genes$cluster)
# significant_genes$cluster_group <- factor(significant_genes$cluster_group,
#                                           levels = c("upregulated_Tumor",
#                                                      "downregulated_Tumor",
#                                                      "upregulated_adenoma",
#                                                      "downregulated_adenoma",
#                                                      "upregulated_Non_classic",
#                                                      "downregulated_Non_classic",
#                                                      "upregulated_Normal_cell",
#                                                      "downregulated_Normal_cell"))

gcSample <- split(significant_genes$ENTREZID, significant_genes$cluster)
gcSample

Sys.setenv(KEGG_REST_URL = "https://rest.kegg.jp")
options(timeout = 600)

enrich_KEGG <- compareCluster(
  gcSample,
  fun = "enrichKEGG",
  organism = "hsa",
  pvalueCutoff = 0.05
)


enrich_KEGG <- compareCluster(
  ENTREZID~cluster,
  data = significant_genes,
  fun = "enrichKEGG",
  organism = "hsa",
  pvalueCutoff = 0.05
)


p1 <- dotplot(enrich_KEGG, showCategory = 5, includeAll = T) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 30, vjust = 0.5, hjust = 0.5, size = 10, colour = "black"),
    axis.text.y = element_text(size = 10, colour = "black"),
    axis.title = element_text(size = 12, face = "bold", colour = "black"),
    panel.grid.major = element_line(colour = "gray90", linewidth = 0.2),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = NA)
  ) +
  # scale_fill_continuous(low = "#E26EB8", high = "#7F95D1", name = "p.adjust",
  #                       guide = guide_colorbar(reverse = T)) +
  labs(
    title = "KEGG Pathway Enrichment Analysis",
    x = "Group",
    y = "Pathway"
  )

ggsave(file.path(diff_outdir, "compareCluster-ROI.markers_KEGG.pdf"), p1, width = 6, height = 8)





enrich_GO <- compareCluster(
  gcSample,
  fun = "enrichGO",
  OrgDb = "org.Hs.eg.db",
  ont = "BP",
  pvalueCutoff = 0.01,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05
)

enrich_GO@compareClusterResult$Cluster <- factor(enrich_GO@compareClusterResult$Cluster, 
                                                 levels = c("Normal_cell", "adenoma", "Non_classic", "Tumor"))
write.csv(enrich_GO@compareClusterResult, file.path(diff_outdir, "compareCluster-ROI.markers_GOBP.csv"))


p2 <- dotplot(enrich_GO, showCategory = 5, includeAll = T) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 30, vjust = 0.5, hjust = 0.5, size = 10, colour = "black"),
    axis.text.y = element_text(size = 10, colour = "black"),
    axis.title = element_text(size = 12, face = "bold", colour = "black"),
    panel.grid.major = element_line(colour = "gray90", linewidth = 0.2),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = NA)
  ) +
  # scale_fill_continuous(low = "#E26EB8", high = "#7F95D1", name = "p.adjust",
  #                       guide = guide_colorbar(reverse = T)) +
  labs(
    title = "GO Pathway Enrichment Analysis",
    x = "Group",
    y = "Pathway"
  )

ggsave(file.path(diff_outdir, "compareCluster-ROI.markers_GOBP.pdf"), p2, width = 6, height = 8)



enrich_GO <- compareCluster(
  gcSample,
  fun = "enrichGO",
  OrgDb = "org.Hs.eg.db",
  ont = "CC",
  pvalueCutoff = 0.01,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05
)

enrich_GO@compareClusterResult$Cluster <- factor(enrich_GO@compareClusterResult$Cluster, 
                                                 levels = c("Normal_cell", "adenoma", "Non_classic", "Tumor"))
write.csv(enrich_GO@compareClusterResult, file.path(diff_outdir, "compareCluster-ROI.markers_GOCC.csv"))

p3 <- dotplot(enrich_GO, showCategory = 5, includeAll = T) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 30, vjust = 0.5, hjust = 0.5, size = 10, colour = "black"),
    axis.text.y = element_text(size = 10, colour = "black"),
    axis.title = element_text(size = 12, face = "bold", colour = "black"),
    panel.grid.major = element_line(colour = "gray90", linewidth = 0.2),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = NA)
  ) +
  # scale_fill_continuous(low = "#E26EB8", high = "#7F95D1", name = "p.adjust",
  #                       guide = guide_colorbar(reverse = T)) +
  labs(
    title = "GO Pathway Enrichment Analysis",
    x = "Group",
    y = "Pathway"
  )

ggsave(file.path(diff_outdir, "compareCluster-ROI.markers_GOCC.pdf"), p3, width = 6, height = 6)



enrich_GO <- compareCluster(
  gcSample,
  fun = "enrichGO",
  OrgDb = "org.Hs.eg.db",
  ont = "MF",
  pvalueCutoff = 0.01,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05
)

p4 <- dotplot(enrich_GO, showCategory = 5, includeAll = T) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 30, vjust = 0.5, hjust = 0.5, size = 10, colour = "black"),
    axis.text.y = element_text(size = 9, colour = "black"),
    axis.title = element_text(size = 12, face = "bold", colour = "black"),
    panel.grid.major = element_line(colour = "gray90", linewidth = 0.2),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(10, 10, 10, 20)
  ) +
  # scale_fill_continuous(low = "#E26EB8", high = "#7F95D1", name = "p.adjust",
  #                       guide = guide_colorbar(reverse = T)) +
  labs(
    title = "GO Pathway Enrichment Analysis",
    x = "Group",
    y = "Pathway"
  )

ggsave(file.path(diff_outdir, "compareCluster-ROI.markers_GOMF.pdf"), p4, width = 6, height = 7)










# 20260603 三组富集分析 modify-----------
library(tidyverse)
library(Seurat)
library(patchwork)
library(stringr)
library("clusterProfiler")
library("org.Hs.eg.db")
library(RColorBrewer)
library(EnhancedVolcano)
library(dplyr)
library(clusterProfiler)
library(ggplot2)



diff_outdir <- file.path("20260603-diff analysis modify")
dir.create(diff_outdir)


all_data <- readRDS("all_data_ROI.rds")
all_data_filt <- all_data[, all_data@meta.data$ROI != "Other"]

Idents(all_data_filt) <- "ROI"
DefaultAssay(all_data_filt) <- "RNA"




diff_genes <- FindAllMarkers(
  JoinLayers(all_data_filt),
  test.use = "wilcox", # 使用Wilcoxon秩和检验
  logfc.threshold = 0.25, # 设置log2折叠变化阈值
  min.pct = 0.1, # 基因至少在10%的细胞中表达
  only.pos = TRUE, # 包含上调和下调基因
  return.thresh = 0.05 # 只返回调整p值<0.05的基因
)

significant_genes <- diff_genes %>%
  group_by(gene) %>%
  filter(all(p_val_adj < 0.05)) %>%
  ungroup()

write.csv(significant_genes, file.path(diff_outdir, "group diffmarker.csv"))


table(significant_genes$cluster)
ids <- bitr(significant_genes$gene, "SYMBOL", "ENTREZID", "org.Hs.eg.db")
significant_genes <- merge(significant_genes, ids, by.x = "gene", by.y = "SYMBOL")

# significant_genes$group[significant_genes$avg_log2FC > 0] <- "upregulated"
# significant_genes$group[significant_genes$avg_log2FC < 0] <- "downregulated"
# significant_genes$cluster_group <- paste0(significant_genes$group, "_", significant_genes$cluster)
# significant_genes$cluster_group <- factor(significant_genes$cluster_group,
#                                           levels = c("upregulated_Tumor",
#                                                      "downregulated_Tumor",
#                                                      "upregulated_adenoma",
#                                                      "downregulated_adenoma",
#                                                      "upregulated_Non_classic",
#                                                      "downregulated_Non_classic",
#                                                      "upregulated_Normal_cell",
#                                                      "downregulated_Normal_cell"))

gcSample <- split(significant_genes$ENTREZID, significant_genes$cluster)
gcSample

Sys.setenv(KEGG_REST_URL = "https://rest.kegg.jp")
options(timeout = 600)

enrich_KEGG <- compareCluster(
  gcSample,
  fun = "enrichKEGG",
  organism = "hsa",
  pvalueCutoff = 0.05
)


# enrich_KEGG <- compareCluster(
#   ENTREZID~cluster,
#   data = significant_genes,
#   fun = "enrichKEGG",
#   organism = "hsa",
#   pvalueCutoff = 0.05
# )

enrich_KEGG@compareClusterResult$Cluster <- 
  factor(enrich_KEGG@compareClusterResult$Cluster, 
         levels = c("Normal_cell", "adenoma", "Non_classic", "Tumor"))

write.csv(enrich_KEGG@compareClusterResult, 
          file.path(diff_outdir, "compareCluster-ROI.markers_KEGG.csv"))


p1 <- dotplot(enrich_KEGG, showCategory = 5, includeAll = T) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 30, vjust = 0.5, hjust = 0.5, size = 10, colour = "black"),
    axis.text.y = element_text(size = 10, colour = "black"),
    axis.title = element_text(size = 12, face = "bold", colour = "black"),
    panel.grid.major = element_line(colour = "gray90", linewidth = 0.2),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = NA)
  ) +
  # scale_fill_continuous(low = "#E26EB8", high = "#7F95D1", name = "p.adjust",
  #                       guide = guide_colorbar(reverse = T)) +
  labs(
    title = "KEGG Pathway Enrichment Analysis",
    x = "Group",
    y = "Pathway"
  )

ggsave(file.path(diff_outdir, "compareCluster-ROI.markers_KEGG.pdf"), p1, width = 6, height = 8)





enrich_GO <- compareCluster(
  gcSample,
  fun = "enrichGO",
  OrgDb = "org.Hs.eg.db",
  ont = "BP",
  pvalueCutoff = 0.01,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05
)

enrich_GO@compareClusterResult$Cluster <- factor(enrich_GO@compareClusterResult$Cluster, 
                                                 levels = c("Normal_cell", "adenoma", "Non_classic", "Tumor"))
write.csv(enrich_GO@compareClusterResult, file.path(diff_outdir, "compareCluster-ROI.markers_GOBP.csv"))


p2 <- dotplot(enrich_GO, showCategory = 5, includeAll = T) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 30, vjust = 0.5, hjust = 0.5, size = 10, colour = "black"),
    axis.text.y = element_text(size = 10, colour = "black"),
    axis.title = element_text(size = 12, face = "bold", colour = "black"),
    panel.grid.major = element_line(colour = "gray90", linewidth = 0.2),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = NA)
  ) +
  # scale_fill_continuous(low = "#E26EB8", high = "#7F95D1", name = "p.adjust",
  #                       guide = guide_colorbar(reverse = T)) +
  labs(
    title = "GO Pathway Enrichment Analysis",
    x = "Group",
    y = "Pathway"
  )

ggsave(file.path(diff_outdir, "compareCluster-ROI.markers_GOBP.pdf"), p2, width = 6, height = 8)



enrich_GO <- compareCluster(
  gcSample,
  fun = "enrichGO",
  OrgDb = "org.Hs.eg.db",
  ont = "CC",
  pvalueCutoff = 0.01,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05
)

enrich_GO@compareClusterResult$Cluster <- 
  factor(enrich_GO@compareClusterResult$Cluster,
         levels = c("Normal_cell", "adenoma", "Non_classic", "Tumor"))

write.csv(enrich_GO@compareClusterResult, 
          file.path(diff_outdir, "compareCluster-ROI.markers_GOCC.csv"))

p3 <- dotplot(enrich_GO, showCategory = 5, includeAll = T) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 30, vjust = 0.5, hjust = 0.5, size = 10, colour = "black"),
    axis.text.y = element_text(size = 10, colour = "black"),
    axis.title = element_text(size = 12, face = "bold", colour = "black"),
    panel.grid.major = element_line(colour = "gray90", linewidth = 0.2),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = NA)
  ) +
  # scale_fill_continuous(low = "#E26EB8", high = "#7F95D1", name = "p.adjust",
  #                       guide = guide_colorbar(reverse = T)) +
  labs(
    title = "GO Pathway Enrichment Analysis",
    x = "Group",
    y = "Pathway"
  )

ggsave(file.path(diff_outdir, "compareCluster-ROI.markers_GOCC.pdf"), p3, width = 6, height = 6)



enrich_GO <- compareCluster(
  gcSample,
  fun = "enrichGO",
  OrgDb = "org.Hs.eg.db",
  ont = "MF",
  pvalueCutoff = 0.01,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05
)

enrich_GO@compareClusterResult$Cluster <- 
  factor(enrich_GO@compareClusterResult$Cluster,
         levels = c("Normal_cell", "adenoma", "Non_classic", "Tumor"))

write.csv(enrich_GO@compareClusterResult, 
          file.path(diff_outdir, "compareCluster-ROI.markers_GOMF.csv"))


p4 <- dotplot(enrich_GO, showCategory = 5, includeAll = T) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 30, vjust = 0.5, hjust = 0.5, size = 10, colour = "black"),
    axis.text.y = element_text(size = 9, colour = "black"),
    axis.title = element_text(size = 12, face = "bold", colour = "black"),
    panel.grid.major = element_line(colour = "gray90", linewidth = 0.2),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(10, 10, 10, 20)
  ) +
  # scale_fill_continuous(low = "#E26EB8", high = "#7F95D1", name = "p.adjust",
  #                       guide = guide_colorbar(reverse = T)) +
  labs(
    title = "GO Pathway Enrichment Analysis",
    x = "Group",
    y = "Pathway"
  )

ggsave(file.path(diff_outdir, "compareCluster-ROI.markers_GOMF.pdf"), p4, width = 6, height = 7)












