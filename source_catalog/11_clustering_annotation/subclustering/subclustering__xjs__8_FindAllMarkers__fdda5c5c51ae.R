
library(tidyverse)
library(Seurat)
library(patchwork)
library(stringr)
library('clusterProfiler')
library("org.Hs.eg.db")
library(RColorBrewer)
library(EnhancedVolcano)
library(dplyr) 



# 读取数据
# all_data <- readRDS("all_data_RCTD_subtype_1023.rds")
# 
# # enrichment--- FindAllMarkers ------------------
# unique(all_data$sample)
# unique(all_data$group)

Tumor_filt <- readRDS("RCTD_new/Tumor/Tumor_filt.rds")

tumor_PC <- subset(Tumor_filt, subset = group == "PC")

Idents(tumor_PC) <- "subtype"
DefaultAssay(tumor_PC) <- "RNA"

tumor_enrich_outdir <- "Tumor_enrich"
dir.create(tumor_enrich_outdir)



diffmarkers <- FindAllMarkers(
  JoinLayers(tumor_PC), 
  test.use = "wilcox", 
  logfc.threshold = 0.25, 
  min.pct = 0.1, 
  only.pos = TRUE, 
  
)

write.csv(diffmarkers, file.path(tumor_enrich_outdir, 'cluster diffmarker.csv'))


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
    tumor_PC, 
    features = all_top5_gene_unique
  ) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

ggsave(
  p_dot_top5,
  filename = file.path(tumor_enrich_outdir, "top 5 gene dotplot.pdf"),
  height = 6,
  width  = 18
)




### Tumor enrich -------
celltype_diffmarkers <- diffmarkers %>%
  filter(cluster == "C1", p_val_adj < 0.05)
write.csv(celltype_diffmarkers, file.path(tumor_enrich_outdir, 'cluster C1 diffmarker.csv'))


enrichGO_ls <- enrichGO(
  rownames(celltype_diffmarkers), 
  OrgDb = org.Hs.eg.db, 
  ont = "ALL", 
  pAdjustMethod = "BH", 
  pvalueCutoff = 0.05, 
  qvalueCutoff = 0.05, 
  keyType = "SYMBOL"
)
write.csv(enrichGO_ls@result, file = file.path(tumor_enrich_outdir, "C1 GO enrichment table.csv"), row.names = FALSE)

enrichGO_dot <- enrichplot::dotplot(
  enrichGO_ls, 
  showCategory = 10, 
  label_format = 50
)
ggsave(enrichGO_dot, filename = file.path(tumor_enrich_outdir, "C1 GO dotplot.pdf"), height = 6, width = 8)

enrichGO_bar <- barplot(
  enrichGO_ls, 
  showCategory = 10, 
  label_format = 50
)
ggsave(enrichGO_bar, filename = file.path(tumor_enrich_outdir, "C1 GO barplot.pdf"), height = 6, width = 8)


de_genes_id_for_enrichmenet <- 
  clusterProfiler::bitr(
    rownames(celltype_diffmarkers),
    OrgDb = "org.Hs.eg.db",
    fromType = "SYMBOL",
    toType = "ENTREZID"
  ) %>%
  .$ENTREZID

enrichKEGG_ls <- enrichKEGG(
  de_genes_id_for_enrichmenet,
  organism = "hsa", 
  keyType = "ncbi-geneid",
  # keyType = "kegg", 
  pvalueCutoff = 0.05, 
  pAdjustMethod = "BH", 
  minGSSize = 10,
  maxGSSize = 500, 
  qvalueCutoff = 0.2,
  use_internal_data = FALSE
  # use_internal_data = TRUE
)
enrichKEGG_ls <- setReadable(enrichKEGG_ls, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
write.csv(enrichKEGG_ls@result,file = file.path(tumor_enrich_outdir, 'C1 KEGG enriment table.csv'),row.names = FALSE)

enrichKEGG_dot <- enrichplot::dotplot(
  enrichKEGG_ls,
  showCategory = 10,
  label_format = 50)
enrichKEGG_bar <- barplot(
  enrichKEGG_ls,
  showCategory = 10,
  label_format = 50)
ggsave(
  enrichKEGG_dot,
  filename = file.path(tumor_enrich_outdir, 'C1 KEGG dotplot.pdf'),
  height = 6,
  width = 8)
ggsave(
  enrichKEGG_bar,
  filename = file.path(tumor_enrich_outdir, 'C1 KEGG barplot.pdf'),
  height = 6,
  width = 8)





celltype_diffmarkers <- diffmarkers %>%
  filter(cluster == "C4", p_val_adj < 0.05)
write.csv(celltype_diffmarkers, file.path(tumor_enrich_outdir, 'cluster C4 diffmarker.csv'))


enrichGO_ls <- enrichGO(
  rownames(celltype_diffmarkers), 
  OrgDb = org.Hs.eg.db, 
  ont = "ALL", 
  pAdjustMethod = "BH", 
  pvalueCutoff = 0.05, 
  qvalueCutoff = 0.05, 
  keyType = "SYMBOL"
)
write.csv(enrichGO_ls@result, file = file.path(tumor_enrich_outdir, "C4 GO enrichment table.csv"), row.names = FALSE)

enrichGO_dot <- enrichplot::dotplot(
  enrichGO_ls, 
  showCategory = 10, 
  label_format = 50
)
ggsave(enrichGO_dot, filename = file.path(tumor_enrich_outdir, "C4 GO dotplot.pdf"), height = 6, width = 8)

enrichGO_bar <- barplot(
  enrichGO_ls, 
  showCategory = 10, 
  label_format = 50
)
ggsave(enrichGO_bar, filename = file.path(tumor_enrich_outdir, "C4 GO barplot.pdf"), height = 6, width = 8)


de_genes_id_for_enrichmenet <- 
  clusterProfiler::bitr(
    rownames(celltype_diffmarkers),
    OrgDb = "org.Hs.eg.db",
    fromType = "SYMBOL",
    toType = "ENTREZID"
  ) %>%
  .$ENTREZID

enrichKEGG_ls <- enrichKEGG(
  de_genes_id_for_enrichmenet,
  organism = "hsa", 
  keyType = "ncbi-geneid",
  # keyType = "kegg", 
  pvalueCutoff = 0.05, 
  pAdjustMethod = "BH", 
  minGSSize = 10,
  maxGSSize = 500, 
  qvalueCutoff = 0.2,
  use_internal_data = FALSE
  # use_internal_data = TRUE
)
enrichKEGG_ls <- setReadable(enrichKEGG_ls, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
write.csv(enrichKEGG_ls@result,file = file.path(tumor_enrich_outdir, 'C4 KEGG enriment table.csv'),row.names = FALSE)

enrichKEGG_dot <- enrichplot::dotplot(
  enrichKEGG_ls,
  showCategory = 10,
  label_format = 50)
enrichKEGG_bar <- barplot(
  enrichKEGG_ls,
  showCategory = 10,
  label_format = 50)
ggsave(
  enrichKEGG_dot,
  filename = file.path(tumor_enrich_outdir, 'C4 KEGG dotplot.pdf'),
  height = 6,
  width = 8)
ggsave(
  enrichKEGG_bar,
  filename = file.path(tumor_enrich_outdir, 'C4 KEGG barplot.pdf'),
  height = 6,
  width = 8)









