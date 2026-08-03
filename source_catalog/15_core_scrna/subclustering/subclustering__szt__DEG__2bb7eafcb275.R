
library(Seurat)
library(ggplot2)
library(scRNAtoolVis)
library(openxlsx)

allData_harmony <- readRDS("../22-AUCell分析/alldata_harmony.rds")

table(allData_harmony$subType, allData_harmony$AUC_TLS_Group)

Idents(allData_harmony) = paste0('c_', allData_harmony$AUC_TLS_Group)
table(Idents(allData_harmony))
#-----------------------------------------------------------------------------
#进行细胞类型间的差异分析
degs <- data.frame()
for (i in unique(allData_harmony$subType)) {
  tmp <- FindMarkers(allData_harmony[,allData_harmony$subType==i], 
              avg_log2FC.threshold = 0.1,
              ident.1 = 'c_activate', 
              ident.2 = 'c_non_activate')
  tmp$cluster <- i
  tmp$gene <- rownames(tmp)
  
  degs <- rbind(degs, tmp)
}
write.xlsx(degs, file = "DEG.xlsx", asTable = T, overwrite = T)

col_Cell_subType <- c("CD8+ T cells" = "#49B1DDCC", "Tgd" = "#11A579CC", "Other CD4+ T cells" = "#3969ACCC",
                      "Th1" = "#3969ACCC", "Th2" = "#3969ACCC", "Th17" = "#3969ACCC", "Tfh" = "#3969ACCC",
                      "CD56DimNK" = "#FE76C8CC", "iDC" = "#FDFD16CC",
                      "B cells" = "#FEB705CC", "Mast cells" = "#FF0085CC",  "Type II pneumocytes" = "#FFBB78CC",
                      "Monocytes" = "#95C6C6CC",  "Alveolar macrophages" = "#7F3C8DCC", "Endothelial cells" = "#008695CC", 
                      "Fibroblasts" = "#00F5FFCC", "Ciliated cells" = "#000000CC", "Plasma cells" = "#8C564BCC",
                      "Type I pneumocytes" = "#FF7F0ECC")

p <- jjVolcano(diffData = degs, aesCol = c('#219EBCCC', '#E73F74CC'), tile.col = col_Cell_subType, 
               log2FC.cutoff = 0.1, topGeneN = 0,  pSize  = 1.5, size = 3,
               legend.position = c(0.8,0.2), flip = T)
p
ggsave(plot = p, filename = "DEG.tif", device = "tiff", dpi = 300, width = 7, height = 12)
ggsave(plot = p, filename = "DEG.pdf", device = "pdf", dpi = 300, width = 7, height = 12)
#-----------------------------------------------------------------------------
library(clusterProfiler)

degs <- read.xlsx("DEG.xlsx", sheet = 1)
degs.sel <- degs[degs$avg_log2FC > 0,]
 
#构建分组，转化genesymbol-gene ID
group <- data.frame(gene = degs.sel$gene,
                    group = degs.sel$cluster)

Gene_ID <- bitr(degs.sel$gene, 
                fromType="SYMBOL", 
                toType="ENTREZID", 
                OrgDb="org.Hs.eg.db")

#构建文件并分析
data  <- merge(Gene_ID, group, by.x='SYMBOL', by.y='gene')

go.res <- compareCluster(ENTREZID~group, 
  data=data, 
  fun="enrichGO", 
  OrgDb="org.Hs.eg.db",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05, 
  readable= TRUE)

data_GO_sim <- simplify(go.res, cutoff = 0.7, by = "p.adjust", select_fun = min)
data_GO_sim_fil <- data_GO_sim@compareClusterResult

p <- dotplot(data_GO_sim, showCategory = 5, includeAll = F, font.size = 9, label_format = 60, title = "Top 5 GO BP Enrichment Results")
p <- p + scale_fill_gradient(high ="#219EBC", low = "#E73F74")
p <- p + theme(axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5),
              axis.text = element_text(color = 'black', size = 10))
p
ggsave(p, filename = "enrich_gobp_dotplot.png", width = 10, height = 15)
ggsave(p, filename = "enrich_gobp_dotplot.pdf", width = 10, height = 15)
#-----------------------------------------------------------------------------
kegg.res <- compareCluster(geneClusters = ENTREZID~group, data = data, fun = "enrichKEGG", organism="hsa")
data_GO_sim_fil <- data_GO_sim@compareClusterResult

p <- dotplot(kegg.res, showCategory = 5, includeAll = F, font.size = 10, label_format = 50, title = "Top 10 KEGG")
p <- p + scale_fill_gradient(high ="#219EBC", low = "#E73F74")
p <- p + theme(axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5),
               axis.text = element_text(color = 'black', size = 10))
p
ggsave(p, filename = "enrich_kegg_dotplot.pdf", width = 10, height = 13)










