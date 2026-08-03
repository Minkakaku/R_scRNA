setwd('/mnt/NexData_file/home/zwl/Project/personal/demo/muscat/')
obj = readRDS('../anno.RDS')
DimPlot(obj,reduction = 'umap',group.by = c('orig.ident','group','cell_type'))

ALLDEG = NULL
for (i in unique(obj$cell_type)) {
  subobj = subset(obj,cell_type == i)
  DEG = FindAllMarkers(subobj,group.by = 'group',only.pos = FALSE, min.pct = 0.1, logfc.threshold = 0.25)
  DEG$cell_type = i
  ALLDEG = rbind(ALLDEG,DEG)
}
write.csv(file = './Muscat/4.FindAllMarkers_DEG.csv',ALLDEG)
library(UpSetR)
ALLDEG = ALLDEG %>% dplyr::filter(avg_log2FC >1)
de_gs_by_k <- split(ALLDEG$gene, ALLDEG$cell_type)  
p1 = upset(fromList(de_gs_by_k))
png(paste0('./Muscat/4.FindAllMarkers_DEG_upsetplot.png'),width = 10, height = 7, units = "in", res = 300)
print(p1) 
dev.off()

muscat = read.csv(file = './Muscat/2.GroupDEG_logFC1.csv',row.names = 1)
de_gs_by_k1 <- split(muscat$gene, muscat$cluster_id)  
p2 = upset(fromList(de_gs_by_k1))

library(ggvenn)
A = list(Muscat = muscat$gene, FindAllMarkers = ALLDEG$gene)
p = ggvenn(A,fill_color= c("#E41A1C", "#1E90FF"),fill_alpha = .7,#stroke_linetype = "longdash",
           set_name_size = 8,text_size=5)+theme(plot.background = element_rect('white')) 

library(clusterProfiler)
library(msigdbr)
library(org.Hs.eg.db)
gene_m = bitr(de_gs_by_k1$B, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
gene_f = bitr(de_gs_by_k$B, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
m_t2g <- msigdbr(species = "Homo sapiens", category = "H")%>% 
  dplyr::select(gs_name, entrez_gene)

em1 <- enricher(gene_m$ENTREZID, TERM2GENE=m_t2g,pvalueCutoff = 0.05)
em1 <- setReadable(em1, OrgDb = 'org.Hs.eg.db', keyType="ENTREZID")
write.csv(file = './Muscat/5.Bcell_Muscat_hallmark_pathway.csv',em1@result)
em2 <- enricher(gene_f$ENTREZID, TERM2GENE=m_t2g,pvalueCutoff = 0.05)
em2 <- setReadable(em2, OrgDb = 'org.Hs.eg.db', keyType="ENTREZID")
write.csv(file = './Muscat/5.Bcell_FindallMarker_hallmark_pathway.csv',em2@result)


library(enrichplot)
A = barplot(em1, showCategory=5) 
ggsave(filename = './Muscat/5.Bcell_Muscat_hallmark_pathway.png',A,w=6,h=4)
ggsave(filename = './Muscat/5.Bcell_Muscat_hallmark_pathway.pdf',A,w=6,h=4)
B = barplot(em2, showCategory=5) 
ggsave(filename = './Muscat/5.Bcell_FindallMarker_hallmark_pathway.png',B,w=6,h=4)
ggsave(filename = './Muscat/5.Bcell_FindallMarker_hallmark_pathway.pdf',B,w=6,h=4)
A+B
