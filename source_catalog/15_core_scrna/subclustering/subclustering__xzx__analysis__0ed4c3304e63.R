setwd("/mnt/NexData_file/home/xzx/project/SC2507113630/shouhou2/")
library(Seurat)
library(ggplot2)

SMC <- readRDS("./recluster/Smooth_muscle_cells/h5_trans.rds")
GO_term <- read.csv("./gene_term.csv")
GO1 <- GO_term[1,3]
GO1 <- unlist(strsplit(GO1, split = ","))
GO1

genesets_addmodulescore <- function(gene_sets,gene_name,scRNA){
  genes_list <- list(gene_sets)
  
  scRNA_fin <- scRNA
  scRNA_fin <- AddModuleScore(
    object = scRNA_fin,
    features = genes_list,
    ctrl = 50,
    name = gene_name
  )
  # 提取需要的数据（细胞ID、分组、模块分数）
  plot_data <- data.frame(
    cell_type = scRNA_fin$sub_cell_type,  # 分组信息
    score = scRNA_fin@meta.data[[paste0(gene_name,"1")]]  # 模块分数
  )
  
  #plot_data$cell_type <- factor(plot_data$cell_type, levels = cell_type_order)
  p2 <- ggplot(plot_data, aes(x = cell_type, y = score, fill = cell_type)) +
    geom_violin(scale = "width") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = paste0(gene_name," Module Score"), x = "Cell Type", y = "Score")
  #+
  #  scale_fill_manual(values = cell_type_colors)
  p2
  ggsave(plot = p2,filename = paste0("/mnt/NexData_file/home/xzx/project/SC2507113630/shouhou2/addmodulescore/",gene_name,"_addmodulescore.png"),
        width = 16,height = 9,dpi = 300,bg = "white")
  ggsave(plot = p2,filename = paste0("/mnt/NexData_file/home/xzx/project/SC2507113630/shouhou2/addmodulescore/",gene_name,"_addmodulescore.pdf"),
        width = 16,height = 9,dpi = 300,bg = "white")
  
  # plot_data <- data.frame(
  #   cluster = scRNA_fin$cluster,  # 分组信息
  #   score = scRNA_fin@meta.data[[paste0(gene_name,"1")]]  # 模块分数
  # )
  # p3 <- ggplot(plot_data, aes(x = cluster, y = score, fill = cluster)) +
  #   geom_violin(scale = "width") +
  #   theme_classic() +
  #   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  #   labs(title = paste0(gene_name," Module Score"), x = "Cell Type", y = "Score")
  # ggsave(plot = p3,filename = paste0("/mnt/NexData_file/home/xzx/project/SC2507113630/shouhou2/",gene_name,"_addmodulescore_cluster.png"),
  #        width = 16,height = 9,dpi = 300,bg = "white")
  # ggsave(plot = p3,filename = paste0("/mnt/NexData_file/home/xzx/project/SC2507113630/shouhou2/",gene_name,"_addmodulescore_cluster.pdf"),
  #        width = 16,height = 9,dpi = 300,bg = "white")
  return(scRNA_fin)
}

SMC <- genesets_addmodulescore(gene_sets = GO1,gene_name = "smooth_muscle_contraction_v2",scRNA = SMC)


GO2 <- GO_term[2,3]
GO2 <- unlist(strsplit(GO2, split = ","))
GO2

SMC <- genesets_addmodulescore(gene_sets = GO2,gene_name = "regulation_of_smooth_muscle_contraction_v2",scRNA = SMC)

GO3 <- GO_term[3,3]
GO3 <- unlist(strsplit(GO3, split = ","))
GO3

SMC <- genesets_addmodulescore(gene_sets = GO3,gene_name = "negative_regulation_of_smooth_muscle_contraction_v2",scRNA = SMC)

###
mdb <- msigdbr(species = "Homo sapiens", category = "C5")## 定义基因集，选取C2
fgsea_sets = mdb [grep("^GOBP_RELAXATION_OF_MUSCLE",mdb $gs_name),] %>% split(x = .$gene_symbol, f = .$gs_name)
genes <- fgsea_sets[[1]]
SMC <- genesets_addmodulescore(gene_sets = genes,gene_name = "RELAXATION_OF_SMOOTH_MUSCLE",scRNA = SMC)



###KEGG\GO富集分析
library(grid)
library(ggplot2)
library(gridExtra)
library(clusterProfiler)
#library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(org.Hs.eg.db)
library(GOplot)
library(dplyr)
library(readxl)
markers <- read_delim("./analysis_file/cluster1_differentially_expressed_genes.xls")
gene = markers[markers$logfoldchanges >1,]$names

entrezIDs = bitr(gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb= "org.Hs.eg.db", drop = TRUE)
#> Warning in bitr(gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb =
#> "org.Hs.eg.db", : 7.1% of input gene IDs are fail to map...
gene<- entrezIDs$ENTREZID

## GO
#gcSample=split(marker_new$ID, marker_new$cluster) 
###参数可以更改，看看?compareCluster
#One of "groupGO", "enrichGO", "enrichKEGG", "enrichDO" or "enrichPathway" 
xx <- enrichGO(gene = gene,
               OrgDb = "org.Hs.eg.db",
               ont = "ALL",
               pAdjustMethod = "BH"
)

p <- dotplot(xx)
p
ggsave(plot = p,filename = paste0("/mnt/NexData_file/home/xzx/project/SC2507113630/shouhou2/addmodulescore/cluster1_enrichGO.png"),
       width = 6,height = 6,dpi = 300,bg = "white")
ggsave(plot = p,filename = paste0("/mnt/NexData_file/home/xzx/project/SC2507113630/shouhou2/addmodulescore/cluster1_enrichGO.pdf"),
       width = 6,height = 6,dpi = 300,bg = "white")
kegg <- enrichKEGG(gene = gene,
                   organism = 'hsa',
                   keyType = 'kegg')
p <- dotplot(kegg)
p
ggsave(plot = p,filename = paste0("/mnt/NexData_file/home/xzx/project/SC2507113630/shouhou2/addmodulescore/cluster1_enrichKEGG.png"),
       width = 6,height = 6,dpi = 300,bg = "white")
ggsave(plot = p,filename = paste0("/mnt/NexData_file/home/xzx/project/SC2507113630/shouhou2/addmodulescore/cluster1_enrichKEGG.pdf"),
       width = 6,height = 6,dpi = 300,bg = "white")



cajal <- c("KIT","ANO1","ETV1","PRKCQ", "EDN3")
sc <- genesets_addmodulescore(gene_sets = cajal,gene_name = "cajal",scRNA = sc)
neurons <- c("RBFOX1","ELAVL4","UCHL1", "DSCAM","SNAP25","SYT1","RBFOX3","MAP2","TUBB3")
sc <- genesets_addmodulescore(gene_sets = neurons,gene_name = "neurons",scRNA = sc)
Glia <- c("SOX10","PLP1","S100B")
sc <- genesets_addmodulescore(gene_sets = Glia,gene_name = "Glia",scRNA = sc)
Enteroendocrine <- c("CHGA","CHGB","NEUROD1","TPH1","TAC1","EPCAM","KRT18","KRT8")
sc <- genesets_addmodulescore(gene_sets = Enteroendocrine,gene_name = "Enteroendocrine",scRNA = sc)

table(sc$cell_type)
FeaturePlot(object = sc,features = c("RBFOX3","ELAVL4","TUBB3", "MAP2","SYN1","RBFOX1","NOS1","VIP","CHAT","TAC1","CALB2","HAND2","ASCL1","NEFH","NPY","SLC5A7","ACHE","TIA1","G3BP1","CALB1","PBX3","TBX3","PENK","GDNF","NRTN","PDGFRA"),pt.size = 0)
FeaturePlot(object = sc,features = c("CHGA","CHGB","NEUROD1","TPH1","TAC1","PDGFRA"),pt.size = 0)
FeaturePlot(object = sc,features = c("RBFOX1","ELAVL4","UCHL1", "DSCAM","SNAP25","SYT1"),pt.size = 0)
