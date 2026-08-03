setwd('/mnt/NexData_file/home/zwl/Project/FWSC20242216/require2/')
# 绘制堆叠柱形图
library(ggplot2)
library(RColorBrewer)
library(scales)
library(dplyr)
# #基础色
# brewer_palette <- brewer.pal(n = 8, name = "Set3")
# #扩展色
# final_palette = colorRampPalette(brewer_palette)(length(unique(proportion_df$cell_type)))
rm(list = ls())
gc()
##############################组间差异柱形图#########################
ALLmetadata = read.csv('./1.metadata.csv',row.names = 1) 
Tmetadata = read.csv('1.T_metadata.csv',row.names = 1)
CD8TeXmetadata = read.csv('1.CD8TEX_metadata.csv',row.names = 1)
CD8TeXmetadata$cell_type = paste0('C',as.character(CD8TeXmetadata$leiden))
metadata = list(ALL = ALLmetadata,
                Tcell = Tmetadata,
                CD8TEX = CD8TeXmetadata)
for (i in names(metadata)) {
  count_df <- table(metadata[[i]]$sample, metadata[[i]]$cell_type)
  write.csv(file = paste0(i,'.cell_number.csv'),count_df)
  proportion_df <- prop.table(count_df, margin = 1)
  proportion_df <- as.data.frame(proportion_df)
  colnames(proportion_df) <- c("sample", "cell_type", "proportion")
  
  p0 = ggplot(proportion_df, aes(x = sample, y = proportion, fill = cell_type)) +
    geom_col(color = "black", size = 0.2) +  # 添加边框，颜色为黑色，线宽为 0.2
    geom_text(aes(label = scales::percent(proportion, accuracy = 0.1)), 
              position = position_stack(vjust = 0.5), size = 3) + 
    scale_y_continuous(labels = scales::percent_format()) +
    #scale_fill_manual(values = final_palette)+
    theme_classic()
  ggsave(paste0(i,'.cell_proportion_barchart.png'),p0,w=6,h=4,dpi=200)
  ggsave(paste0(i,'.cell_proportion_barchart.pdf'),p0,w=6,h=4,dpi=200)
  ##############################组间差异箱线图#########################
  
  data = proportion_df %>% merge(unique(metadata[[i]][,c("sample","group")]),by="sample") 
  
  p_value = data%>%
    group_by(cell_type) %>% 
    rstatix::t_test(formula = proportion~group)%>%
    add_significance() %>% add_xy_position()#x = "cell_type"
  ##新罗马没办法加粗
  p = ggboxplot(
    data, x = "group", y = "proportion", 
    color = "group", palette = "jco",
    facet.by = "cell_type",add = "jitter"
  ) + 
    stat_pvalue_manual(p_value, label = "p",hide.ns = F)
  ggsave(paste0(i,'.cell_proportion_boxplot.png'),p,w=4,h=6,dpi=200)
  ggsave(paste0(i,'.cell_proportion_boxplot.pdf'),p,w=4,h=6,dpi=200)
}

##############################组间差异火山图#########################








##############################CD8Tmonocle3###############################
#############h5ad转变成seurat对象####################
setwd('./CD8 Tex_monocle')
sceasy::convertFormat('/mnt/NexData_file/home/zwl/Project/FWSC20242216/require2/CD8 Tex/result.h5ad', from="anndata", to="seurat",
                      outFile='CD8Tex_result.rds')

################导入注释信息####################

obj <- readRDS(file = '/mnt/NexData_file/home/zwl/Project/FWSC20242216/require2/CD8 Tex_monocle/CD8Tex_result.rds')
obj <- FindVariableFeatures(obj,nfeatures = 2000)
Seurat::DimPlot(obj,reduction = 'umap',group.by = c('cell_type','leiden'))

library(monocle3)
f_data = data.frame(gene_short_name = rownames(obj),
                    row.names = rownames(obj))
expression_matrix = as.matrix(LayerData(obj, assay="RNA", layer='counts'))
cds <- new_cell_data_set(expression_matrix,
                         cell_metadata = obj@meta.data,
                         gene_metadata =  f_data)
cds <- preprocess_cds(cds, num_dim = 100)
plot_pc_variance_explained(cds)
cds <- reduce_dimension(cds)

##从seurat导入整合过的umap坐标
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(obj, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed
plot_cells(cds, color_cells_by="cell_type")

## Monocle3聚类分区
cds <- cluster_cells(cds)
p1 <- plot_cells(cds,color_cells_by = "cell_type", show_trajectory_graph = FALSE) + 
  ggtitle("label by clusterID")
p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE) + 
  ggtitle("label by partitionID")
wrap_plots(p1, p2)
#用注释后的cluster
library(do)
obj$cluster1 = 1
cds@clusters@listData[["UMAP"]][["partitions"]] = obj$
plot_cells(cds,color_cells_by = "partition", show_trajectory_graph = FALSE) + 
  ggtitle("label by newclusterID")

## 识别轨迹
cds <- learn_graph(cds,close_loop=FALSE, #分支不画圈
                   use_partition = FALSE)# 所有细胞看称整体画分支
                   #learn_graph_control=list(minimal_branch_len = 6))#保留一些比较短的分支
p = plot_cells(cds, color_cells_by = "leiden",label_groups_by_cluster = FALSE,
               label_leaves = FALSE,group_label_size  = 5, 
               label_branch_points = FALSE)

cds <- order_cells(cds)
p2 = plot_cells(cds,
                color_cells_by = "pseudotime",
                label_cell_groups=FALSE,
                label_leaves=FALSE,
                label_branch_points=FALSE,
                graph_label_size=3)
ggsave(filename = 'monocle3_trajectory_pseudotime_umap.png',p+p2,w=10,h=5)
save(file = 'CD8T_order.Rdata',cds)

###################################################T cell组间GSVA############################################### 
setwd('/mnt/NexData_file/home/zwl/Project/FWSC20242216/require2/')
rm(list = ls())
gc()
library(clusterProfiler)
library(org.Hs.eg.db)

cell = c('CD4 Tem_Tcm','CD4 Tn','CD8 GZMK','CD8 Tem','NK','Treg')
for (i in cell) {
  DEG = read.csv(file = paste0('./3.T_enrichment/',i,'/S_vs_F/all_cluster_marker.xls'),sep = '\t')
  DEG <- DEG %>% mutate(logfoldchanges = ifelse(cluster == "S", -logfoldchanges, logfoldchanges))
  
  eg = bitr(DEG$Symbol, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
  #降序排列
  DEG_ALL = merge(DEG,eg,by.x='Symbol',by.y='SYMBOL')%>%arrange(.,desc(logfoldchanges))
  geneList <- DEG_ALL$logfoldchanges
  names(geneList) <- DEG_ALL$ENTREZID
  
  ego_BP <- gseGO(geneList     = geneList,
                  OrgDb        = org.Hs.eg.db,
                  ont          = "BP",
                  pvalueCutoff = 0.05,
                  verbose      = FALSE)
  ego_BP <- setReadable(ego_BP, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
  res = ego_BP@result
  write.csv(file = paste0('./3.T_enrichment/',i,'/1.GSEA_F_vs_S_GO_BP_pathway_res.csv'),res)
  
  ego_CC <- gseGO(geneList     = geneList,
                  OrgDb        = org.Hs.eg.db,
                  ont          = "CC",
                  pvalueCutoff = 0.05,
                  verbose      = FALSE)
  ego_CC <- setReadable(ego_CC, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
  res = ego_CC@result
  write.csv(file = paste0('./3.T_enrichment/',i,'/1.GSEA_F_vs_S_GO_CC_pathway_res.csv'),res)
  
  ego_MF <- gseGO(geneList     = geneList,
                  OrgDb        = org.Hs.eg.db,
                  ont          = "MF",
                  pvalueCutoff = 0.05,
                  verbose      = FALSE)
  ego_MF <- setReadable(ego_MF, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
  res = ego_MF@result
  write.csv(file = paste0('./3.T_enrichment/',i,'/1.GSEA_F_vs_S_GO_MF_pathway_res.csv'),res)
  
  kk2 <- gseKEGG(geneList     = geneList,
                 organism     = 'human',
                 pvalueCutoff = 0.5,
                 verbose      = FALSE)
  kk2 <- setReadable(kk2, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
  write.csv(file = paste0('./3.T_enrichment/',i,'/1.GSEA_F_vs_S_KEGG_pathway_res.csv'),kk2@result)
  
  
  library(DOSE)
  library(enrichplot)
  library(scales)
  library(RColorBrewer)
  if (dim(ego_BP@result)[1]>0) {
    if (dim(ego_BP@result)[1]<5) {
      n=dim(ego_BP@result)[1]
    }else{n=5}
    p = gseaplot2(ego_BP, geneSetID = 1:n)
    ggsave(filename = paste0('./3.T_enrichment/',i,'/2.GO_BPpathway_group_F_vs_S.png'),w=6,h=5,dpi = 200,p)
    ggsave(filename = paste0('./3.T_enrichment/',i,'/2.GO_BPpathway_group_F_vs_S.png'),w=6,h=5,dpi = 200,p)
  }
  if (dim(ego_CC@result)[1]>0) {
    if (dim(ego_CC@result)[1]<5) {
      n=dim(ego_CC@result)[1]
    }else{n=5}
    p2 = gseaplot2(ego_CC, geneSetID = 1:n)
    ggsave(filename = paste0('./3.T_enrichment/',i,'/2.GO_CCpathway_group_F_vs_S.png'),w=6,h=4,dpi = 200,p2)
    ggsave(filename = paste0('./3.T_enrichment/',i,'/2.GO_CCpathway_group_F_vs_S.png'),w=6,h=4,dpi = 200,p2)
  }
  if (dim(ego_MF@result)[1]>0) {
    if (dim(ego_MF@result)[1]<5) {
      n=dim(ego_MF@result)[1]
    }else{n=5}
    p3 = gseaplot2(ego_MF, geneSetID = 1:n)
    ggsave(filename = paste0('./3.T_enrichment/',i,'/2.GO_MFpathway_group_F_vs_S.png'),w=6,h=4,dpi = 200,p3)
    ggsave(filename = paste0('./3.T_enrichment/',i,'/2.GO_MFpathway_group_F_vs_S.png'),w=6,h=4,dpi = 200,p3)
  }
  if (dim(kk2@result)[1]>0) {
    if (dim(kk2@result)[1]<5) {
      n=dim(kk2@result)[1]
    }else{n=5}
    p4 = gseaplot2(kk2, geneSetID = 1:n)
    ggsave(filename = paste0('./3.T_enrichment/',i,'/2.KEGG_pathway_group_F_vs_S.png'),w=6,h=4,dpi = 200,p4)
    ggsave(filename = paste0('./3.T_enrichment/',i,'/2.KEGG_pathway_group_F_vs_S.png'),w=6,h=4,dpi = 200,p4)
  }
}
