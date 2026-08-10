setwd('/mnt/NexData_file/home/zwl/Project/AE2404160288/require2/myeloid')
rm(list = ls())
gc()

library(cowplot)
library(Seurat)
library(RColorBrewer)

# ###########################################exchange h5ad to rds############################################
library(sceasy)
library(Seurat)
library(reticulate)
#use_condaenv('myname')

# sceasy::convertFormat('../1177/celltype/anno.h5ad', from="anndata", to="seurat",outFile='1177_anno.rds')
obj <- readRDS(file = './raw_anno.rds')
obj <- FindVariableFeatures(obj,nfeatures = 2000)
Seurat::DimPlot(obj,reduction = 'umap',group.by = c('group','cell_type','leiden'))

############################################reorder##########################################################
library(monocle)
library(Seurat)
library(ggplot2)
library(argparse)
library(dplyr)
library(tidyverse)
library(ggsci)
library(patchwork)
library(stringr)


HSMM = readRDS(paste0('/mnt/NexData_file/home/zwl/Project/AE2404160288/require1/Macro_monocle/monocle.rds'))
deg <- read.csv(paste0('/mnt/NexData_file/home/zwl/Project/AE2404160288/require1/Macro_monocle/monocle.DEG.xls'),sep="\t")
rownames(deg) = deg$Genes
ordergene <- rownames(deg)[1:min(1000,nrow(deg))]
# 重新排序
plot_cell_trajectory(HSMM, color_by = "State")+
  facet_wrap(~cell_type,ncol = 4)
HSMM_myo <- orderCells(HSMM, root_state = 1)


# 重新画图
# setwd(sample[i])
getwd()
clustercolname = 'cell_type'

p = plot_cell_trajectory(HSMM_myo,color_by="Pseudotime",size=1,show_backbone=TRUE)
ggsave("monocle.trajectory.pseudotime.pdf",p,width = 7,height = 7.5)
ggsave("monocle.trajectory.pseudotime.png",p,width = 7,height = 7.5)

p = plot_complex_cell_trajectory(HSMM_myo,x=1,y=2,color_by = clustercolname)
ggsave("monocle.trajectory.tree.pdf",p,width = 7,height = 7.5)
ggsave("monocle.trajectory.tree.png",p,width = 7,height = 7.5)
ncols = ceiling(sqrt(length(unique(HSMM_myo@phenoData@data[,clustercolname]))))
nrows = ceiling(length(unique(HSMM_myo@phenoData@data[,clustercolname]))/ncols)
p = plot_complex_cell_trajectory(HSMM_myo,x=1,y=2,color_by = clustercolname)+facet_wrap(vars(!!sym(clustercolname)),nrow=nrows)
ggsave("monocle.trajectory.tree.faceted.pdf",p,width = 4*ncols,height = 4.2*nrows)
ggsave("monocle.trajectory.tree.faceted.png",p,width = 4*ncols,height = 4.2*nrows)

pseudotime_df = HSMM_myo@phenoData@data[,c("sample",clustercolname,"Size_Factor","num_genes_expressed","Pseudotime","State")]
df = data.frame(Cells=rownames(pseudotime_df),pseudotime_df)
write.table(df,"3.monocle.metadata_reorder.xls",sep = "\t",quote = F,row.names = F)

saveRDS(file = '3.obj_reorder.rds',HSMM_myo,compress = F)

# 逆时序山脊图
library(ggplot2)
library(ggridges)
df = HSMM_myo@phenoData@data %>% dplyr::select(cell_type,Pseudotime,group) 
median_df <- df %>%
  group_by(cell_type) %>%
  summarise(median_pseudo = median(Pseudotime, na.rm = TRUE))%>%
  mutate(y_pos = as.numeric(as.factor(cell_type)))

df$y_pos <- as.numeric(as.factor(df$cell_type))

# 3. 绘图：山脊图 + 仅在本行的中位数短虚线
p <- ggplot(df, aes(x = Pseudotime, y = cell_type, fill = cell_type)) +
  geom_density_ridges(rel_min_height = 0.005, alpha = 0.6) +
  
  # 🔥 只画在山脊内部的短虚线（关键！）
  geom_segment(
    data = median_df,
    aes(x = median_pseudo, xend = median_pseudo,
        y = y_pos, yend = y_pos + 1), # 只画在当前山脊高度内
    color = "black",
    linetype = "dashed",
    linewidth = 0.7
  ) +
  
  theme_ridges() +
  theme(plot.background = element_rect(fill = 'white'))

p
ggsave(filename = "3.monocle.pseudotime_celltype_group.png",p,h=7,w=10,dpi = 150)
ggsave(filename = "3.monocle.pseudotime_celltype_group.pdf",p,h=7,w=10,dpi = 150)


# p = plot_cell_trajectory(HSMM_myo,color_by=clustercolname,size=1,show_backbone=TRUE)+
#   facet_grid(cell_type~group)

# 拟时相关的差异gene,ordergene做回归
diff_test_res <- differentialGeneTest(HSMM_myo[ordergene,],
                                      fullModelFormulaStr = "~sm.ns(Pseudotime)")
sig_gene_names <- diff_test_res %>% dplyr::filter(num_cells_expressed >10,qval<0.05)%>% #防止小类干扰
  dplyr::arrange(qval) %>% row.names() %>% head(50)
write.csv(file = '3.obj_Pseudotime_ordergene.csv',sig_gene_names)
p = plot_pseudotime_heatmap(HSMM_myo[sig_gene_names,],
                            num_clusters = 5,
                            cores = 1,
                            show_rownames = T,return_heatmap = TRUE)
ggsave(filename = "monocle.top50gene_cortime.png",p,h=8,w=6,dpi = 150)
ggsave(filename = "monocle.top50gene_cortime.pdf",p,h=8,w=6,dpi = 150)


#############################BEAM分析###########################################
cds = HSMM_myo
outpath = '/mnt/NexData_file/home/zwl/Project/AE2404160288/require2/myeloid'
# for(branch_point in 1:((length(levels(cds@phenoData@data[,"State"]))-1)/2)){
branch_point = 2
outbeam = paste0(outpath,"/BEAM_analysis.Node",branch_point)
print(paste0("Analysis Beam Node ",branch_point))
if(!dir.exists(outbeam)){
  dir.create(outbeam)
}

BEAM_res = BEAM(cds,branch_point = branch_point,cores = 4, progenitor_method = "duplicate")
BEAM_res = BEAM_res[,c("gene_short_name","pval","qval")]
BEAM_res = BEAM_res[order(BEAM_res$qval),]
write.csv(BEAM_res,paste0(outbeam,"/BEAM_res.csv"),row.names = F)
selectgenes = row.names(subset(BEAM_res,qval<0.05))

p=plot_genes_branched_pseudotime(cds[selectgenes[1:min(20,length(selectgenes))],],branch_point = branch_point,color_by = "State",ncol = 4)
ggsave(paste0(outbeam,"/BEAM_pseudotime_top20.png"),p,width=15,height = 12)
ggsave(paste0(outbeam,"/BEAM_pseudotime_top20.pdf"),p,width=15,height = 12)

p = plot_genes_branched_heatmap(cds[selectgenes,],
                                branch_point = branch_point,
                                num_clusters = 4,
                                cores=1,
                                use_gene_short_name = T,
                                show_rownames = T,return_heatmap = T)
ggsave(paste0(outbeam,"/BEAM_heatmap.png"),p$ph_res,width=7,height = 8)
ggsave(paste0(outbeam,"/BEAM_heatmap.pdf"),p$ph_res,width=7,height = 8)

write.csv(p$annotation_row,file = '/mnt/NexData_file/home/zwl/Project/SC2507233885/require4/BEAM_analysis.Node1/BEAM_heatmap_row.csv')



# }

library(clusterProfiler)
data = p$annotation_row
for (i in unique(p$annotation_row$Cluster)) {
  gene_set = data %>% dplyr::filter(Cluster == i) %>% rownames()
  GO_ALL <- enrichGO(gene = gene_set,  #基因列表(转换的ID)
                     keyType = "SYMBOL",  #指定的基因ID类型，默认为ENTREZID
                     OrgDb=org.Mm.eg.db,  #物种对应的org包
                     ont = "ALL",   #CC细胞组件，MF分子功能，BP生物学过程
                     pvalueCutoff = 0.05,  #p值阈值
                     pAdjustMethod = "BH",  #多重假设检验校正方式
                     minGSSize = 3,   #注释的最小基因集，默认为10
                     maxGSSize = 500,  #注释的最大基因集，默认为500
                     qvalueCutoff = 1,  #q值阈值
                     readable = TRUE)  #基因ID转换为基因名
  write.table(GO_ALL, paste0(outbeam,"/BEAM_cluster",i, "_enrich.go.csv"), sep = ',', row.names = FALSE, quote = TRUE)
  
  if(!is.null(GO_ALL) && nrow(GO_ALL) > 0){
    dotplot(GO_ALL,  #GO富集分析结果
            x = "GeneRatio",  #横坐标,默认GeneRation,也可以为Count
            color = "p.adjust",  #右纵坐标,默认p.adjust,也可以为pvalue和qvalue
            showCategory = 10,  #展示前20个点，默认为10个
            size = NULL,  #点的大小
            title = "GO_dotplot",  #设置图片的标题
            split='ONTOLOGY',
            label_format = 80)+
      facet_grid(ONTOLOGY~.,scale="free")
    
    ggsave(paste0(outbeam,"/BEAM_cluster",i, ".GO_dotplot.pdf"),width = 12,height = 9)
    ggsave(paste0(outbeam,"/BEAM_cluster",i, ".GO_dotplot.png"),width = 12,height = 9)
    
    barplot(GO_ALL,  #GO富集分析结果
            x = "Count",  #横坐标,默认Count,也可以为GeneRation
            color = "p.adjust",  #右纵坐标,默认p.adjust,也可以为pvalue和qvalue
            showCategory = 10,  #展示前20个，默认为10个
            size = NULL,  
            title = "GO_barplot",  #设置图片的标题
            split='ONTOLOGY',
            label_format = 80)+
      facet_grid(ONTOLOGY~.,scale="free")
    
    ggsave(paste0(outbeam,"/BEAM_cluster",i, ".GO_barplot.pdf"),width = 12,height = 9)
    ggsave(paste0(outbeam,"/BEAM_cluster",i, ".GO_barplot.png"),width = 12,height = 9)
  }
  
  gene_symbol <- bitr(geneID = gene_set,  #感兴趣的基因集
                      fromType="SYMBOL",   #输入ID的类型
                      toType=c("ENTREZID"),   #输出ID的类型，可为多个
                      OrgDb="org.Mm.eg.db")  #物种注释数据库
  gene_KEGG <- gene_symbol[,2]
  
  # KEGG<- enrichKEGG(gene = gene_KEGG,   #基因列表(同GO) 
  #                   organism = "mmu",  #物种
  #                   keyType = "kegg",  #指定的基因ID类型，默认为kegg
  #                   minGSSize = 3, 
  #                   maxGSSize = 500,
  #                   pvalueCutoff = 0.05,  
  #                   pAdjustMethod = "BH",
  #                   qvalueCutoff = 1)
  # 
  # id_map <- setNames(gene_symbol$SYMBOL, gene_symbol$ENTREZID)
  # KEGG@result$geneID <- sapply(KEGG@result$geneID,
  #                              function(x) paste(id_map[strsplit(x, "/")[[1]]], collapse = "/"))
  # write.table(KEGG, paste0(outbeam,"/BEAM_cluster",i,"_enrich.KEGG.csv"), sep = ',', row.names = FALSE, quote = TRUE)
  # 
  # KEGG@result$Description <- sub(" - Mus musculus \\(house mouse\\)$", "", KEGG@result$Description)
  # if(!is.null(KEGG) && nrow(KEGG) > 0){
  #   p <- dotplot(KEGG,  #GO富集分析结果
  #                x = "GeneRatio",  #横坐标,默认GeneRation,也可以为Count
  #                color = "p.adjust",  #右纵坐标,默认p.adjust,也可以为pvalue和qvalue
  #                showCategory = 10,  #展示前20个点，默认为10个
  #                size = NULL,  #点的大小
  #                title = "KEGG_dotplot",  #设置图片的标题
  #                label_format = 100)
  #   p + guides(colour = "Count")
  #   ggsave(paste0(outbeam,"/BEAM_cluster",i, ".KEGG_dotplot.pdf"),width = 12,height = 5)
  #   ggsave(paste0(outbeam,"/BEAM_cluster",i, ".KEGG_dotplot.png"),width = 12,height = 5)
  #   
  #   barplot(KEGG,  #GO富集分析结果
  #           x = "Count",  #横坐标,默认Count,也可以为GeneRation
  #           color = "p.adjust",  #右纵坐标,默认p.adjust,也可以为pvalue和qvalue
  #           showCategory = 10,  #展示前20个，默认为10个
  #           size = NULL,  
  #           title = "KEGG_barplot",  #设置图片的标题
  #           label_format = 80)
  #   
  #   ggsave(paste0(outbeam,"/BEAM_cluster",i, ".KEGG_barplot.pdf"),width = 12,height = 5)
  #   ggsave(paste0(outbeam,"/BEAM_cluster",i, ".KEGG_barplot.png"),width = 12,height = 5)
  # }
}


