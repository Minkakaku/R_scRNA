
library(Seurat)
library(ggplot2)
library(scRNAtoolVis)
library(openxlsx)

allData_harmony <- readRDS("../22-AUCell分析/alldata_harmony.rds")

table(allData_harmony$subType, allData_harmony$AUC_TLS_Group)

Idents(allData_harmony) = paste0('c_', allData_harmony$AUC_TLS_Group)
table(Idents(allData_harmony))
#-------------------------------------------------------------------
#整体细胞的差异分析
degs_all <- FindMarkers(allData_harmony,
            avg_log2FC.threshold = 0,
            ident.1 = 'c_activate', 
            ident.2 = 'c_non_activate')
write.xlsx(degs_all, file = "DEG_all.xlsx", asTable = T, overwrite = T, rowNames = T)

diffSig <- degs_all[abs(degs_all$avg_log2FC) > 0.1,]
diffSig$Tag <- ifelse(diffSig$avg_log2FC > 0, yes = "UP", no = "DOWN")
write.xlsx(diffSig, file = "DEG_diffSig.xlsx", asTable = T, overwrite = T, rowNames = T)
#-------------------------------------------------------------------
library(tidyverse)
library(ggrepel)
library(ggfun)
library(grid)

y.sel <- "p_val_adj"
p.val <- 0.05
fc <- 2^0.1

allDiff_use <- degs_all
test <- which(allDiff_use[,y.sel] == 0)
allDiff_use[,y.sel][test] <- min(degs_all[,y.sel][degs_all[,y.sel]!=min(degs_all[,y.sel])])*(10^-1)

####----load data----####
df <- allDiff_use %>%
  rownames_to_column(var = "SYMBOL")
####将选择的基因进行可视化和添加名称####
####----------默认为差异基因前20------------####
df.sel <- df[df$SYMBOL %in% rownames(diffSig),] %>%
  tidyr::drop_na() %>%
  dplyr::arrange(desc(-log10(get(y.sel)))) %>%
  dplyr::slice(1:20)
####-------自己进行指定------------####
df.sel <- df %>%
  filter(SYMBOL %in% diffSig.sel$Gene)

####-----------图片绘制------------####
p <- ggplot(data = df)
p <- p + geom_point(aes(x = avg_log2FC, y = -log10(get(y.sel)), color = avg_log2FC, size = -log10(get(y.sel))))

p <- p + geom_point(data = df.sel, aes(x = avg_log2FC, y = -log10(get(y.sel)),
                                       size = -log10(get(y.sel))), shape = 21, show.legend = F, color = "#000000")

p <- p + geom_text_repel(data = df.sel %>%
                           dplyr::filter(avg_log2FC > 0),
                         aes(x = avg_log2FC, y = -log10(get(y.sel)), label = SYMBOL),
                         box.padding = 0.5,
                         max.overlaps = 20,
                         nudge_x = 0.5,
                         nudge_y = 0.2,
                         segment.curvature = -0.1,
                         segment.ncp = 3,
                         direction = "y",
                         hjust = "left")

p <- p + geom_text_repel(data =  df.sel %>%
                           dplyr::filter(avg_log2FC < 0),
                         aes(x = avg_log2FC, y = -log10(get(y.sel)), label = SYMBOL),
                         box.padding = 0.5,
                         max.overlaps = 20,
                         nudge_x = -0.2,
                         nudge_y = 0.2,
                         segment.curvature = -0.1,
                         segment.ncp = 3,
                         segment.angle = 20,
                         direction = "y",
                         hjust = "right")
####-------设置点的颜色----------####
p <- p + scale_color_gradientn(colours = c("#49B1DD", "#79C0D1","#E9CF91", "#E98DAF", "#E73F74"),
                               values = seq(0, 1, 0.2))
p <- p + scale_fill_gradientn(colours = c("#49B1DD", "#79C0D1","#E9CF91", "#E98DAF", "#E73F74"),
                              values = seq(0, 1, 0.2))
####-------设置坐标轴----------####
p <- p + xlim(c(min(df$avg_log2FC)-0.01, max(df$avg_log2FC)+0.01)) +
  ylim(c(0, max(-log10(df[,y.sel]))+1.5))
p <- p + geom_vline(xintercept = c(-log2(fc), log2(fc)), linetype = 2) +
  geom_vline(xintercept = c(0), linetype = "solid") +
  geom_hline(yintercept = -log10(p.val), linetype = 4) +
  scale_size(range = c(1,7), name = paste0("-log10(", y.sel, ")")) +
  ggtitle(label = paste0('TLS activated vs. non-activated Different Expression mRNAs'),
          subtitle = paste(paste(y.sel, "<", p.val ,"& |log2FC| >", log2(fc), ", UP DEG = "),
                           nrow(diffSig[diffSig$avg_log2FC > 0,]), ", DOWN DEG = ",
                           nrow(diffSig[diffSig$avg_log2FC < 0,]), sep = "")) +
  ylab(label = paste0("-log10(", y.sel, ")"))
####-------设置theme----------####
p <- p + theme_bw() +
  theme(panel.grid = element_blank(),
        legend.background = element_roundrect(color = "#808080", linetype = 1),
        axis.text = element_text(size = 13, color = "#000000"),
        axis.title = element_text(size = 15),
        plot.title = element_text(hjust = 0, size = 18),
        plot.subtitle = element_text(hjust = 0, size = 13)) +
  coord_cartesian(clip = "off")

####-------添加方向----------####
p <- p + annotation_custom(
  grob = grid::segmentsGrob(y0 = unit(-10, "pt"), y1 = unit(-10, "pt"),
                            arrow = arrow(angle = 45, length = unit(.2, "cm"), ends = "first"),
                            gp = grid::gpar(lwd = 3, col = "#219EBC")),
  xmin = min(df$avg_log2FC),
  xmax = -0.5,
  ymin = max(-log10(df[,y.sel]))+1,
  ymax = max(-log10(df[,y.sel]))+1)
p <- p + annotation_custom(
  grob = grid::textGrob(label = "Down",gp = grid::gpar(col = "#219EBC")),
  xmin = min(df$avg_log2FC),
  xmax = -0.5,
  ymin = max(-log10(df[,y.sel]))+1,
  ymax = max(-log10(df[,y.sel]))+1)
p <- p + annotation_custom(
  grob = grid::segmentsGrob(y0 = unit(-10, "pt"), y1 = unit(-10, "pt"),
                            arrow = arrow(angle = 45, length = unit(.2, "cm"), ends = "last"),
                            gp = grid::gpar(lwd = 3, col = "#E73F74")),
  xmin = max(df$avg_log2FC),
  xmax = 0.5,
  ymin = max(-log10(df[,y.sel]))+1,
  ymax = max(-log10(df[,y.sel]))+1)
p <- p + annotation_custom(
  grob = grid::textGrob(label = "Up", gp = grid::gpar(col = "#E73F74")),
  xmin = max(df$avg_log2FC),
  xmax = 0.5,
  ymin = max(-log10(df[,y.sel]))+1,
  ymax = max(-log10(df[,y.sel]))+1)
p
####----save result----####
ggsave(filename = "volcano_plot.pdf", height = 9, width = 9, dpi = 300)
ggsave(filename = "volcano_plot.tif", height = 9, width = 9, device = "tiff", dpi = 300)
#-------------------------------------------------------------------
## 差异基因的功能富集分析
library(openxlsx)
library(clusterProfiler)
library(R.utils)
setOption("clusterProfiler.download.method", "auto")
library(AnnotationDbi)
library(org.Hs.eg.db)
library(VennDiagram)

diffSig <- read.xlsx("DEG_diffSig.xlsx", sheet = 1)
diffSig.sel <- diffSig[diffSig$Tag %in% "DOWN",]

gene.sel <- "DOWN"
write.xlsx(diffSig.sel, file = paste0(gene.sel, "_DEG.xlsx"), overwrite = T, asTable = T)

#symbol ID转换为ENTREZ ID
DEG.gene_symbol <- as.character(diffSig.sel$row.names)
DEG.entrez_id <- mapIds(x = org.Hs.eg.db,
                        keys = DEG.gene_symbol,
                        keytype = "SYMBOL",
                        column = "ENTREZID")
DEG.entrez_id <- na.omit(DEG.entrez_id)

#GO分析
erich.go.BP <- enrichGO(gene = DEG.entrez_id,
                        OrgDb = org.Hs.eg.db,
                        keyType = "ENTREZID",
                        ont = "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff = 1,
                        qvalueCutoff = 1,
                        readable= TRUE)
erich.go.BP <- as.data.frame(erich.go.BP)
erich.go.BP <- erich.go.BP[erich.go.BP$pvalue < 0.05,]
write.xlsx(erich.go.BP, paste0("GO.bp.Sig.xlsx"), overwrite = T, asTable = T)

erich.go.CC <- enrichGO(gene = DEG.entrez_id,
                        OrgDb = org.Hs.eg.db,
                        keyType = "ENTREZID",
                        ont = "CC",
                        pAdjustMethod = "BH",
                        pvalueCutoff = 1,
                        qvalueCutoff = 1,
                        readable= TRUE)
erich.go.CC <- as.data.frame(erich.go.CC)
erich.go.CC <- erich.go.CC[erich.go.CC$pvalue < 0.05,]
write.xlsx(erich.go.CC, paste0("GO.cc.Sig.xlsx"), overwrite = T, asTable = T)

erich.go.MF <- enrichGO(gene = DEG.entrez_id,
                        OrgDb = org.Hs.eg.db,
                        keyType = "ENTREZID",
                        ont = "MF",
                        pAdjustMethod = "BH",
                        pvalueCutoff = 1,
                        qvalueCutoff = 1,
                        readable= TRUE)
erich.go.MF <- as.data.frame(erich.go.MF)
erich.go.MF <- erich.go.MF[erich.go.MF$pvalue < 0.05,]
write.xlsx(erich.go.MF, paste0("GO.mf.Sig.xlsx"), overwrite = T, asTable = T)

#KEGG分析
egoPathway <- enrichKEGG(gene = DEG.entrez_id,
                         organism = "hsa",
                         keyType = "kegg",
                         pAdjustMethod = "BH",
                         pvalueCutoff = 1,
                         qvalueCutoff = 1)
#将entrez重新转换为symbol：
egoPathway2 <- setReadable(egoPathway, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
#转换前：
egoPathway@result$geneID[1]
#转换后：
egoPathway2@result$geneID[1]
egoPathway2 <- as.data.frame(egoPathway2)
egoPathway2 <- egoPathway2[egoPathway2$pvalue < 0.05,]
egoPathway2$Description <- gsub(pattern = " \\- Mus musculus \\(house mouse\\)", replacement = "", x = egoPathway2$Description)
write.xlsx(egoPathway2, paste0("KEGG.Sig.xlsx"), overwrite = T, asTable = T)
#----------------------------------------------------------------
#将富集结果整合做长条形图
library(ggplot2)
library(stringr)
library(openxlsx)
library(dplyr)

GO.bp.sel <- read.xlsx("GO.bp.Sig.xlsx", sheet = 2)
GO.bp.sel$Category <- "GO BP"
GO.cc.sel <- read.xlsx("GO.cc.Sig.xlsx", sheet = 2)
GO.cc.sel$Category <- "GO CC"
GO.mf.sel <- read.xlsx("GO.mf.Sig.xlsx", sheet = 2)
GO.mf.sel$Category <- "GO MF"
kegg.sel <- read.xlsx("KEGG.Sig.xlsx", sheet = 2)
kegg.sel$Category <- "KEGG"

enrich <- rbind(GO.bp.sel, GO.cc.sel)
enrich <- rbind(enrich, GO.mf.sel)
enrich <- rbind(enrich, kegg.sel)
##计算Rich Factor（富集因子）：
enrich <- mutate(enrich, RichFactor = Count/as.numeric(sub("/\\d+", "", BgRatio)))
##计算Fold Enrichment（富集倍数）：
enrich$FoldEnrichment <- apply(enrich,1,function(x){
  GeneRatio <- eval(parse(text = x["GeneRatio"]))
  BgRatio <- eval(parse(text = x["BgRatio"]))
  foldEnrichment <- round(GeneRatio/BgRatio, 2)
  foldEnrichment})

#指定绘图顺序（转换为因子）：
enrich$Description <- factor(paste(enrich$Description), levels = rev(paste(enrich$Description)))
mytheme <- theme(axis.title = element_text(size = 13),
                 axis.text = element_text(size = 11),
                 plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
                 legend.title = element_text(size = 13),
                 legend.text = element_text(size = 11))

#先在自定义主题中隐去y轴文本标签显示
mytheme2 <- mytheme + theme(axis.text.y = element_blank())
enrich$text_x <- rep(0.005, nrow(enrich)) #新增一列重复数组，使绘图时文本标签能从固定位置开始；
p2 <- ggplot(data = enrich, aes(x = RichFactor, y = Description))
p2 <- p2 + geom_bar(aes(fill = Category), stat = "identity", width = 0.8, alpha = 0.7)
p2 <- p2 + scale_fill_manual(values = c("GO BP" = "#067467", "GO CC" = "#F5C35E",
                                        "GO MF" = "#E98DAF", "KEGG" = "#49B1DD"))
p2 <- p2 + labs(x = "RichFactor", y = "pathway", title = "Top GO & KEGG enrichment barplot")
p2 <- p2 + geom_text(aes(x = text_x, label = Description), hjust = 0) #hjust=0，左对齐
p2 <- p2 + theme_classic() + mytheme2
p2
ggsave(p2, filename = "enrich_barplot.pdf", width = 5, height = 14)
ggsave(p2, filename = "enrich_barplot.png", width = 5, height = 14)
#----------------------------------------------------------------
#分组富集分析
diffSig.up <- diffSig[diffSig$Tag %in% "UP",]
diffSig.down <- diffSig[diffSig$Tag %in% "DOWN",]

gc <- list(all = na.omit(mapIds(x = org.Hs.eg.db, keys = diffSig$row.names, keytype = "SYMBOL", column = "ENTREZID")), 
           UP = na.omit(mapIds(x = org.Hs.eg.db, keys = diffSig.up$row.names, keytype = "SYMBOL", column = "ENTREZID")), 
           DOWN = na.omit(mapIds(x = org.Hs.eg.db, keys = diffSig.down$row.names, keytype = "SYMBOL", column = "ENTREZID")))

go.res <- compareCluster(gc, fun = "enrichGO", OrgDb = 'org.Hs.eg.db', ont = "BP", readable= TRUE)
#使用语义分析去除冗余的GO条目
data_GO_sim <- simplify(go.res, cutoff = 0.7, by = "p.adjust", select_fun = min)
data_GO_sim_fil <- data_GO_sim@compareClusterResult
write.xlsx(data_GO_sim_fil, file = "GP_BP.xlsx", asTable = T, overwrite = T)

p <- dotplot(data_GO_sim, showCategory = 20, includeAll = F, font.size = 9, label_format = 50, title = "Top 10 GO BP")
p <- p + scale_fill_gradient(high ="#219EBC", low = "#E73F74")
p
ggsave(p, filename = "enrich_gobp_dotplot.pdf", width = 5.4, height = 14)
#--------------------------------------------------------------------
go.res <- compareCluster(gc, fun = "enrichGO", OrgDb = 'org.Hs.eg.db', ont = "CC", readable= TRUE)
#使用语义分析去除冗余的GO条目
data_GO_sim <- simplify(go.res, cutoff = 0.7, by = "p.adjust", select_fun = min)
data_GO_sim_fil <- data_GO_sim@compareClusterResult
write.xlsx(data_GO_sim_fil, file = "GP_CC.xlsx", asTable = T, overwrite = T)
#--------------------------------------------------------------------
go.res <- compareCluster(gc, fun = "enrichGO", OrgDb = 'org.Hs.eg.db', ont = "MF", readable= TRUE)
#使用语义分析去除冗余的GO条目
data_GO_sim <- simplify(go.res, cutoff = 0.7, by = "p.adjust", select_fun = min)
data_GO_sim_fil <- data_GO_sim@compareClusterResult
write.xlsx(data_GO_sim_fil, file = "GP_MF.xlsx", asTable = T, overwrite = T)
#----------------------------------------------------------------
kegg.res <- compareCluster(gc, fun = "enrichKEGG", organism="hsa")
#将entrez重新转换为symbol：
kegg.res <- setReadable(kegg.res, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
kegg.res_fil <- kegg.res@compareClusterResult
write.xlsx(kegg.res_fil, file = "KEGG.xlsx", asTable = T, overwrite = T)

p <- dotplot(kegg.res, showCategory = 20, includeAll = F, font.size = 9, label_format = 50, title = "Top 10 KEGG")
p <- p + scale_fill_gradient(high ="#219EBC", low = "#E73F74")
p
ggsave(p, filename = "enrich_kegg_dotplot.pdf", width = 5, height = 14)

