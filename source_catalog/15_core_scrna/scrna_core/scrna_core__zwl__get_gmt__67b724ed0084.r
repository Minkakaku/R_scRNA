setwd('/mnt/NexData_file/home/zwl/Project/SC2507233889-01/hallmark/')
library(clusterProfiler)
library(org.Mm.eg.db)
library(msigdbr)
library(msigdbr)
rm(list = ls())
gc()

# 定义需要提取的肿瘤侵袭相关Hallmark通路（包含原有间充质转化相关通路+新增通路）
hallmark_paths <- c(
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  #"HALLMARK_MESENCHYMAL_EPITHELIAL_TRANSITION",
  "HALLMARK_TGF_BETA_SIGNALING",
  "HALLMARK_WNT_BETA_CATENIN_SIGNALING",
  "HALLMARK_INFLAMMATORY_RESPONSE"
 # "HALLMARK_IL6_JAK_STAT3_SIGNALING"  # 新增的IL6通路
)

# 获取人类Hallmark基因集数据库
mdb <- msigdbr(species = "Homo sapiens", category = "H")

# 循环处理每个通路，生成对应的gmt文件
for (path in hallmark_paths) {
  fgsea_sets <- mdb[mdb$gs_name == path, ] %>% 
    split(x = .$gene_symbol, f = .$gs_name)
  
  genes <- fgsea_sets[[1]]
  
  line <- paste(path, paste(genes, collapse = "\t"), sep = "\t")
  writeLines(line, paste0(path, ".gmt"))
}


##################### 
setwd('/mnt/NexData_file/home/zwl/Project/SC2507233889-01/monocle')
dir.create('common_up')
dir.create('common_down')
up = read.csv('./1.pseudotime_upgene_commen.csv')
down = read.csv('./1.pseudotime_downgene_commen.csv')

upg = bitr(up$gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
downg = bitr(down$x, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

m_t2g <- msigdbr(species = "Homo sapiens", category = "H")%>% 
  dplyr::select(gs_name, entrez_gene)

em1 <- enricher(upg$ENTREZID, TERM2GENE=m_t2g,pvalueCutoff = 0.05)
em1 <- setReadable(em1, OrgDb = 'org.Hs.eg.db', keyType="ENTREZID")
em2 <- enricher(downg$ENTREZID, TERM2GENE=m_t2g,pvalueCutoff = 0.05)
em2 <- setReadable(em2, OrgDb = 'org.Hs.eg.db', keyType="ENTREZID")
write.csv(file = './2.pseudotime_up_commongene_hallmark_pathway.csv',em1@result)
write.csv(file = './2.pseudotime_down_commongene_hallmark_pathway.csv',em2@result)


library(enrichplot)
A = barplot(em1, showCategory=10) 
ggsave(filename = './2.pseudotime_up_commongene_hallmark_pathway.png',A,w=6,h=4)
ggsave(filename = './2.pseudotime_up_commongene_hallmark_pathway.pdf',A,w=6,h=4)
B = barplot(em2, showCategory=10) 
ggsave(filename = './2.pseudotime_down_commongene_hallmark_pathway.png',B,w=6,h=4)
ggsave(filename = './2.pseudotime_down_commongene_hallmark_pathway.pdf',B,w=6,h=4)

