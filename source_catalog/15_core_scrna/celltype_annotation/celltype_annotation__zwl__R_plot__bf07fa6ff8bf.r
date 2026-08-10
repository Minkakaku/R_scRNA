setwd('/mnt/NexData_file/home/zwl/Project/FW2024-225-G/08.plot/')
library(reshape2)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(dplyr)
rm(list=ls())
gc()
#设置字体#######################################
library(showtext)
font_add(family = "TimesNewRoman", regular = "/home/zwl/.fonts/times.ttf",) 
#加载字体
showtext_auto()
#exchange h5ad to rds############################################
library(sceasy)
library(Seurat)
library(reticulate)
sceasy::convertFormat('/home/zwl/Project/FW2024-225-G/04.scanpy_sub/2_HSPC_celltype.h5ad', from="anndata", to="seurat",
                      outFile='HSPC_celltype.rds')
HSPC <- readRDS(file = './HSPC_celltype.rds')
HSPC <- FindVariableFeatures(HSPC,nfeatures = 2000)

Seurat::DimPlot(HSPC,reduction = 'umap',group.by = c('celltype'))
#fig2.B IL6 score小提琴图#######################################
forplot1 = read.csv('./Fig2.b3_IL6_score_violin_plots.csv',row.names = 1)


pdf(file = 'Fig2.b3_IL6_score_violin_plots.pdf')
for (i in unique(forplot1$main_celltype)){
  forplot = forplot1[forplot1$main_celltype==i,]
  
  p_value = forplot%>%
    group_by(score) %>% 
    rstatix::t_test(formula = value~group)%>%
    add_significance() %>% add_xy_position(x = 'score')
  
  gg = ggplot(forplot, aes(score, value, colour = group)) +
    geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8), size=0.5,stroke=0.5, shape=16) +
    geom_boxplot(position = position_dodge(width = 0.8), outlier.size = 0, width = 0.3, show.legend = FALSE) +
    geom_violin(position = position_dodge(width = 0.8), scale = "width", alpha = 0.5, width = 0.6, size = 0.8) + 
    scale_color_manual(values = c('iMCD'="#ff7f0e", 'NC'="#1f77b4"))+
    stat_pvalue_manual(p_value, label = "p.signif",hide.ns = T)+
    theme_classic()+labs(y=i,x=NULL)+
    theme(text = element_text(family = "TimesNewRoman"))
  print(gg)
}
dev.off()
# 亚群小提琴图
pdf(file = 'HSPC_subcelltype_IL6_score_violin_plots.pdf')
forplot11 = forplot1[forplot1$main_celltype=="HSC",]
for (i in unique(forplot11$sub_celltype2)){
  forplot = forplot11[forplot11$sub_celltype2== i,]
  p_value = forplot%>%
    group_by(score) %>% 
    rstatix::t_test(formula = value~group)%>%
    add_significance() %>% add_xy_position(x = 'score')
  
  gg = ggplot(forplot, aes(score, value, colour = group)) +
    geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8), size=0.5,stroke=0.5, shape=16) +
    geom_boxplot(position = position_dodge(width = 0.8), outlier.size = 0, width = 0.3, show.legend = FALSE) +
    geom_violin(position = position_dodge(width = 0.8), scale = "width", alpha = 0.5, width = 0.6, size = 0.8) + 
    scale_color_manual(values = c('iMCD'="#ff7f0e", 'NC'="#1f77b4"))+
    stat_pvalue_manual(p_value, label = "p.signif",hide.ns = T)+
    theme_classic()+labs(y=i,x=NULL)+
    theme(text = element_text(family = "TimesNewRoman"))
  print(gg)
}
dev.off()

pdf(file = 'Myeloid_subcelltype_IL6_score_violin_plots2.pdf')
forplot11 = forplot1[forplot1$main_celltype=="Erythrocyte",]
for (i in unique(forplot11$sub_celltype2)){
  forplot = forplot11[forplot11$sub_celltype2== i,]
  p_value = forplot%>%
    group_by(score) %>% 
    rstatix::t_test(formula = value~group)%>%
    add_significance() %>% add_xy_position(x = 'score')
  
  gg = ggplot(forplot, aes(score, value, colour = group)) +
    geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8), size=0.5,stroke=0.5, shape=16) +
    geom_boxplot(position = position_dodge(width = 0.8), outlier.size = 0, width = 0.3, show.legend = FALSE) +
    geom_violin(position = position_dodge(width = 0.8), scale = "width", alpha = 0.5, width = 0.6, size = 0.8) + 
    scale_color_manual(values = c('iMCD'="#ff7f0e", 'NC'="#1f77b4"))+
    stat_pvalue_manual(p_value, label = "p.signif",hide.ns = T)+
    theme_classic()+labs(y=i,x=NULL)+
    theme(text = element_text(family = "TimesNewRoman"))
  print(gg)
}
dev.off()

##fig2.B IL6 score亚群热图##############
IL6_metadata = read.csv(file = './IL6_Signaling_Score.csv')
average_expression <- IL6_metadata %>%
  dplyr::select(sub_celltype2,disease,IL6_Signaling_Score,BIOCARTA_IL6_PATHWAY,KEGG_PATHWAY) %>%
  group_by(sub_celltype2, disease) %>%
  summarise(across(all_of(setdiff(c('IL6_Signaling_Score','BIOCARTA_IL6_PATHWAY','KEGG_PATHWAY'),
                                  c("sub_celltype2", "disease"))),mean))

average_expression = data.frame(average_expression)
#rownames(average_expression)= paste0(average_expression$disease,'_',average_expression$sub_celltype2)
average_expression1 = tidyr::pivot_wider(average_expression,names_from = 'disease',
                   values_from = c('IL6_Signaling_Score','BIOCARTA_IL6_PATHWAY','KEGG_PATHWAY'))%>%
  tibble::column_to_rownames('sub_celltype2')

Sexcolor <- c("#ff7f0e","#1f77b4") 
names(Sexcolor) <- c("iMCD","NC") #类型颜色
ann_colors <- list( Group=Sexcolor) #颜色设置

p1 = pheatmap::pheatmap(average_expression1[c(2,4,6,1,3,5)],
                        annotation_col = data.frame(Group = rep(c('NC','iMCD'),3),
                                                    row.names = colnames(average_expression1)),
                        #color = colorRampPalette(c("navy","white","firebrick3"))(100),
                        annotation_colors = ann_colors,
                        scale = 'none',gaps_col = 3,
                        cluster_cols = F,cluster_rows = F,show_rownames =T,
                        fontfamily = "serif")
ggsave(filename = './Fig2.IL6_score_heatmap.png',p1,w=6,h=7,dpi = 200)
ggsave(filename = './Fig2.IL6_score_heatmap.pdf',p1,w=8,h=8,dpi = 200)

#fig3.1:细胞比例图#######################################
metadata = read.csv("../04.scanpy_sub/3.ALL_subtype_metadata.csv")
list_g = list(Meyloid_HSPC = c("Myeloid","Megakaryocyte","Erythrocyte","HSC"))
              # HSPC = c("HSC"),
              # T_NK = c("T","NK"),
              # B_Plasma = c("B","Plasma"))
for (ci in 1:length(list_g)) {
  stat_cell = metadata %>% dplyr::filter(main_celltype %in% list_g[[ci]]) %>%
    dplyr::select("sample","sub_celltype2") %>%table()
  percent = t(stat_cell)/colSums(stat_cell)*100
  
  data = data.frame(percent) %>% merge(unique(metadata[,c("sample","disease")]),by="sample") 
  
  p_value = data%>%
    group_by(sub_celltype2) %>% 
    rstatix::t_test(formula = Freq~disease)%>%
    add_significance() %>% add_xy_position()#x = "cell_type"
  ##新罗马没办法加粗
  p = ggboxplot(
    data, x = "disease", y = "Freq", 
    color = "disease", palette = "jco",
    facet.by = "sub_celltype2",add = "jitter"
  ) + 
    stat_pvalue_manual(p_value, label = "p.signif",hide.ns = T)+
    theme(text = element_text(family = "TimesNewRoman"))
  ggsave(filename = paste0("Fig3.a1_",names(list_g)[ci],"_celltype_percent.pdf"),w=7,h=6,p)
}


##fig3.4：MK DEG火山图+富集分析#######################################
rm(list = ls())
gc()
MK_DEG = read.csv(file = './Fig3.c_MK_groupDEG.csv',row.names = 1)
gene_exp = read.csv(file = './Fig3.c_MK_ALLgene_meanexp.csv')

plot_data = merge(MK_DEG,gene_exp,by.y = 'X',by.x='names',all = TRUE) %>%
  dplyr::mutate(sig = ifelse(is.na(logfoldchanges),'Stable',ifelse(logfoldchanges>0,'Up in iMCD','Down in iMCD')))
colnames(plot_data)

top_genes = plot_data %>% na.omit() %>% arrange(logfoldchanges) %>%
  slice(1:5, (n() - 4):n())

gg = ggplot(plot_data, aes(NC, iMCD, colour = sig)) +
  geom_jitter(size=1,stroke=0.5) +
  scale_color_manual(values = c('Up in iMCD'="firebrick", 'Down in iMCD'="navy",'Stable'='grey'))+
  theme_classic()+labs(y='Average expression in iMCD',x='Average expression in NC')+
  ggrepel::geom_text_repel(data = top_genes, aes(label = names), size = 4, max.overlaps = Inf) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black")
  theme(text = element_text(family = "TimesNewRoman"))
ggsave(filename = paste0('./Fig3.c1_MK_iMCD_DEG_Vcano.png'),gg,w=4,h=3,dpi = 200)
ggsave(filename = paste0('./Fig3.c1_MK_iMCD_DEG_Vcano.pdf'),gg,w=5,h=4,dpi = 200)
##################################DEG富集
library(clusterProfiler)
library(org.Hs.eg.db)
#geneID转换
gene.df <- bitr(MK_DEG$names,fromType="SYMBOL",toType="ENTREZID", OrgDb = org.Hs.eg.db)[,'ENTREZID']
ego_KEGG <- enrichKEGG(gene          = gene.df,
                       organism     = 'hsa',
                       # keyType = 'SYMBOL',
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.01,
                       qvalueCutoff  = 0.05)
p = barplot(ego_KEGG)  #富集柱形图
ggsave(filename = paste0('./Fig3.c2_Cluster_MK_iMCD.KEGG_enrichment.png'),p,w=5,h=3,dpi = 200)
ggsave(filename = paste0('./Fig3.c2_Cluster_MK_iMCD.KEGG_enrichment.pdf'),p,w=5,h=3,dpi = 200)
write.csv(file = paste0('./Fig3.c2_Cluster_MK_iMCD.KEGG_enrichment.csv'),ego_KEGG@result)
#GO
ego_CC <- enrichGO(gene          = MK_DEG$names,
                   OrgDb         = org.Hs.eg.db,
                   keyType = 'SYMBOL',
                   ont           = "CC",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05,
                   readable      = TRUE)
write.csv(file = paste0('./Fig3.c2_Cluster_MK_iMCD.GO_CC_enrichment.csv'),ego_CC@result)

ego_MF <- enrichGO(gene          = MK_DEG$names,
                   OrgDb         = org.Hs.eg.db,
                   keyType = 'SYMBOL',
                   ont           = "MF",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05,
                   readable      = TRUE)
write.csv(file = paste0('./Fig3.c2_Cluster_MK_iMCD.GO_MF_enrichment.csv'),ego_MF@result)

ego_BP <- enrichGO(gene          = MK_DEG$names,
                   OrgDb         = org.Hs.eg.db,
                   keyType = 'SYMBOL',
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05,
                   readable      = TRUE)
write.csv(file = paste0('./Fig3.c2_Cluster_MK_iMCD.GO_BP_enrichment.csv'),ego_BP@result)

ALL_data = rbind(slice_min(ego_BP@result,order_by = p.adjust,n = 5)%>% dplyr::mutate(category = 'BP'),
                 slice_min(ego_CC@result,order_by = p.adjust,n = 5)%>% dplyr::mutate(category = 'CC'),
                 slice_min(ego_MF@result,order_by = p.adjust,n = 5)%>% dplyr::mutate(category = 'MF'))%>% 
  arrange(category,Count) %>% ungroup()
ALL_data$Description1=factor(ALL_data$Description,levels =ALL_data$Description, ordered = TRUE)
COLS <- c("#66C3A5", "#8DA1CB", "#FD8D62")
p1 = ggplot(data=ALL_data, aes(x=Count,y=Description1, fill=category)) + 
  geom_bar(stat="identity", width=0.8) + 
  scale_fill_manual(values = COLS) + 
  theme_bw() + 
  labs(title = "MK iMCD Enriched GO Terms",y = "GO term",x= "Num of Genes")+ 
  theme(axis.text.x=element_text(face = "bold" ),
        text = element_text(family = "TimesNewRoman"),
        panel.grid=element_blank())

ggsave(filename = paste0('./Fig3.c2_MK_iMCD_GOenrich.png'),p1,w=3,h=3,dpi = 200)
ggsave(filename = paste0('./Fig3.c2_MK_iMCD_GOenrich.pdf'),p1,w=4,h=4,dpi = 200)

#Fig4群体差异基因#################################
rm(list = ls())
gc()
ALL_DEG = data.frame()
for (ci in c('Erythroid progenitor','HSC_MPP','Megakaryocyte','Erythrocyte','GMP_LMPP','MkP','MEP')) {
  data = read.csv(paste0('./Fig4.groupDEG',ci,'.csv'),row.names = 1) %>% dplyr::mutate(celltype = ci)
  ALL_DEG = rbind(ALL_DEG,data)
}

# 依据数值对应颜色进行区分文字
ALL_DEG$celltype = gsub('Megakaryocyte','MK',gsub('Erythroid progenitor','EP',ALL_DEG$celltype))
ALL_DEG$threshold = factor(ifelse(ALL_DEG$logfoldchanges>0,'Up','Down'),levels = c('Up','Down'))
# 调整横轴顺序
ALL_DEG$celltype = factor(ALL_DEG$celltype,levels = c("HSC_MPP","GMP_LMPP","MK","EP","Erythrocyte"))
# 对齐
jitter_position <- position_jitter(width = 0.2, height = 0, seed = 123)
set.seed(123)  
plot_bar = ALL_DEG %>% data.frame() %>%group_by(celltype)  %>%
  dplyr::mutate(up = max(logfoldchanges)+0.3,down = min(logfoldchanges)-0.3)%>%
  dplyr::select(celltype,up,down) %>% unique()

p=ggplot(ALL_DEG,aes(x=celltype,y=logfoldchanges))+  
  geom_point(aes(x=celltype,color=threshold), 
  # position="jitter",                                                 #添加微小扰动点
  position = jitter_position,pch=16, cex=1)+                            #点大小
  scale_color_manual(values=c("#DC143C",                                        #设置theshold第一种level up红色
                              "#00008B"))+                                        #设置theshold第二种level down蓝色
  theme_bw()+
  geom_hline(yintercept = 0,lty=2,                                                        
             col="black",lwd=0.5)+
  geom_col(data = plot_bar,aes(x=celltype,y=up),alpha = 0.2)+
  geom_col(data = plot_bar,aes(x=celltype,y=down),alpha = 0.2)+
  theme_classic()+                                                          #添加线的粗细
  theme(axis.text.x=element_text(hjust=1,vjust=0.5,angle=90),
        text = element_text(family = "TimesNewRoman"))
p

ggsave(paste0("Fig4.a1_ALL_celltype.png"),p,width=4,height=3,dpi = 200)
ggsave(paste0("Fig4.a1_ALL_celltype.pdf"),p,width=6,height=5,dpi = 200)

##Fig4 cellchat#################################
library(CellChat)
load('/home/zwl/Project/FW2024-225-G/cellchat/red_M/cellchat_object.RData')
gg1 <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = NULL, 
               targets.use = NULL, stacked = T, do.stat = TRUE)+
  theme(text = element_text(family = "TimesNewRoman"))
ggsave(filename = './Fig4.c1.relative_information_flow_weight_stacked.png',gg1,w=3,h=3,dpi = 200)
ggsave(filename = './Fig4.c1.relative_information_flow_weight_stacked.pdf',gg1,w=5,h=6,dpi = 200)

gg2 <- rankNet(cellchat, mode = "comparison", measure = "count", sources.use = NULL, 
               targets.use = NULL, stacked = T, do.stat = TRUE)+
  theme(text = element_text(family = "TimesNewRoman"))
ggsave(filename = './Fig4.c1.relative_information_flow_count_stacked.png',gg2,w=3,h=3,dpi = 200)
ggsave(filename = './Fig4.c1.relative_information_flow_count_stacked.pdf',gg2,w=5,h=6,dpi = 200)

#############################Fig4 CXCL12-CXCR4
load('/home/zwl/Project/FW2024-225-G/cellchat/red_M/cellchat_object_list.RData')
pathways.show = c("CXCL")
pairLR.CXCL <- extractEnrichedLR(object.list[["iMCD"]], signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.CXCL[1,] 
# Hierarchy plot
vertex.receiver = seq(1,4) # a numeric vector
pdf(file = 'Fig4.e4.CXCL_signaling_circle.pdf')
netVisual_individual(object.list[["iMCD"]], signaling = pathways.show,  
                     pairLR.use = LR.show, vertex.receiver = vertex.receiver)
dev.off()

#############################Fig4 imcd vs nc 分7个细胞群体 
# 画不出来他想要的
# gg9 = netVisual_bubble(cellchat,  
#                        signaling = c("CCL","CXCL","TGFb"),
#                        comparison = c(1, 2),
#                        #title.name = paste0('iMCD singnal communication'),
#                        angle.x = 90)+
#   theme(legend.position = 'top',
#         legend.key.width = unit(1,'cm'),
#         legend.margin = margin(0.5,0.5,0,0,'cm'),
#         plot.margin=unit(c(2.5, 2.5, 2.5, 2.5),'cm'),
#         text = element_text(family = "TimesNewRoman"))+
#   coord_cartesian(clip = 'off')
# ggsave(filename = '3.CCL_CXCL_TGFb_Singnal_communication_iMCD.png',w=6,h=6,gg9)

##CCL和CXCLgene表达##############################################
obj <- readRDS(file = '/home/zwl/Project/FW2024-225-G/cellchat/red_M//HSPC_redsc.rds')
CellChatDB = CellChatDB.human
data = CellChatDB[["interaction"]]
target = data %>% dplyr::filter(pathway_name %in% c("CCL",'CXCL')) %>% 
  pull(ligand,receptor) %>% unique()
f_target = target[target %in% rownames(obj)]
# 筛选出特定基因集的表达数据
subset_expr <- FetchData(obj,vars = c("subcelltype", "disease",f_target),layer = 'count')
average_expression <- subset_expr %>%
  group_by(subcelltype, disease) %>%
  summarise(across(all_of(setdiff(f_target, c("subcelltype", "disease"))), mean))
#setdiff 函数用于筛选出需要处理的基因列，而 all_of 函数用于确保这些列名在数据框中存在
expression_percentage <- subset_expr %>%
  group_by(subcelltype, disease) %>%
  summarise(across(all_of(setdiff(f_target, c("subcelltype", "disease"))), ~ sum(.x > 0) / length(.x) * 100))
result <- merge(reshape2::melt(average_expression,variable.name = 'gene',value.name = 'Average_express'),
                reshape2::melt(expression_percentage,variable.name = 'gene',value.name = 'Percent_Express'),
                by = c('subcelltype','disease','gene')) %>%
  dplyr::mutate(categary = ifelse(grepl("^CCL", gene), "CCL", 'CXCL'))
result$disease = factor(result$disease,levels = c('iMCD','NC'))
# 输出结果
p = ggplot(result,aes(x=gene,y=subcelltype,size = Percent_Express,fill=Average_express))+
  geom_point(stroke=0.3,shape = 21)+
  scale_fill_gradientn(colors = c("#abddff",'#ffc556','red'))+
  coord_flip()+  theme_bw()+
  theme(panel.grid = element_blank(),axis.text = element_text(color='black'))+
  scale_size_continuous(range = c(0,5))+
  facet_grid(categary~disease, scales = "free", space = "free")+
  theme(axis.text.x = element_text(angle=45,vjust=1,hjust=1), 
        text = element_text(family = "TimesNewRoman",size = 18))+
  labs(x='',y=NULL)+scale_color_manual(values=c('black','black'))
ggsave(filename = './Fig4.d1.CCL_CXCL_dotplot.png',p,w=5,h=6,dpi = 200)
ggsave(filename = './Fig4.d1.CCL_CXCL_dotplot.pdf',p,w=7,h=8,dpi = 200)
################画热图:无占比
average_expression1 = data.frame(average_expression)
rownames(average_expression1)= paste0(average_expression1$disease,'_',average_expression1$subcelltype)

library(RColorBrewer)
groupcolor <- RColorBrewer::brewer.pal(7,'Set2')
names(groupcolor) <- unique(average_expression1$subcelltype) #类型颜色
Sexcolor <- c("#ff7f0e","#1f77b4") 
names(Sexcolor) <- c("iMCD","NC") #类型颜色
ann_colors <- list(subcelltype=groupcolor, disease=Sexcolor) #颜色设置

p1 = pheatmap::pheatmap(average_expression1[c(3:6,8:dim(average_expression1)[2])],
                   annotation_row = average_expression1[1:2],
                   #color = colorRampPalette(c("navy","white","firebrick3"))(100),
                   annotation_colors = ann_colors,
                   scale = 'column',gaps_col = 15,
                   cluster_cols = F,cluster_rows = F,show_rownames =F,
                   fontfamily = "serif")
ggsave(filename = './Fig4.d1.CCL_CXCL_scale_heatmap.png',p1,w=6,h=4,dpi = 200)
ggsave(filename = './Fig4.d1.CCL_CXCL_scale_heatmap.pdf',p1,w=8,h=6,dpi = 200)
#gg.gap()分割y轴

#Fig5##############################################
#"Erythroid progenitor","HSC/MPP","Megakaryocyte","Erythrocyte","GMP/LMPP" ,"MkP","MEP" 
obj<-readRDS('/home/zwl/Project/FW2024-225-G/cellchat/red_M/HSPC_redsc.rds')
Seurat::DimPlot(obj,group.by = c('subcelltype','disease'),reduction = 'umap')
Seurat::DotPlot(obj,features = c('TGFB1'),group.by = 'disease')

dotpltg_cellgp <- function(geneset,filename1,with,high){
  gene_set = c(geneset, c("subcelltype", "disease"))
  subset_expr <- FetchData(obj,vars = gene_set,layer = 'count')
  average_expression <- subset_expr %>%
    group_by(subcelltype, disease) %>%
    summarise(across(all_of(setdiff(gene_set, c("subcelltype", "disease"))), mean))
  #setdiff 函数用于筛选出需要处理的基因列，而 all_of 函数用于确保这些列名在数据框中存在
  expression_percentage <- subset_expr %>%
    group_by(subcelltype, disease) %>%
    summarise(across(all_of(setdiff(gene_set, c("subcelltype", "disease"))), ~ sum(.x > 0) / length(.x) * 100))
  result <- merge(reshape2::melt(average_expression,variable.name = 'gene',value.name = 'Average_express'),
                  reshape2::melt(expression_percentage,variable.name = 'gene',value.name = 'Percent_Express'),
                  by = c('subcelltype','disease','gene'))
  result$disease = factor(result$disease,levels = c('iMCD','NC'))
  
  p = ggplot(result,aes(x=gene,y=subcelltype,size = Percent_Express,fill=Average_express))+
    geom_point(stroke=0.3,shape = 21)+
    scale_fill_gradientn(colors = c("#abddff",'#ffc556','red'))+
    coord_flip()+  theme_bw()+
    theme(panel.grid = element_blank(),axis.text = element_text(color='black'))+
    scale_size_continuous(range = c(0,5))+
    facet_grid(~disease, scales = "free", space = "free")+
    theme(axis.text.x = element_text(angle=45,vjust=1,hjust=1), 
          text = element_text(family = "TimesNewRoman",size = 18))+
    labs(x='',y=NULL)+scale_color_manual(values=c('black','black'))
  ggsave(filename = paste0(filename1,'.png'),p,w=with-2,h=high-1,dpi = 200)
  ggsave(filename = paste0(filename1,'.pdf'),p,w=with,h=high,dpi = 200)
  
}
##1 TGFgene###############
gene_set = c('TGFB1')
obj = HSPC
obj$subcelltype = obj$celltype
dotpltg_cellgp(gene_set,'./Fig5.1.TGFB1_HSPC_dotplot',7,3)

##2 TGF受配体##############
CellChatDB = CellChatDB.human
data = CellChatDB[["interaction"]]
target = data %>% dplyr::filter(pathway_name %in% c("TGFb")) %>% 
  pull(ligand,receptor) %>% unique()
f_target = target[target %in% rownames(obj)]
dotpltg_cellgp(f_target,'./Fig5.2.TGFb_receptor_ligen_dotplot',7,3)

#画热图:无占比
average_expression1 = data.frame(average_expression)
rownames(average_expression1)= paste0(average_expression1$disease,'_',average_expression1$subcelltype)

library(RColorBrewer)
groupcolor <- RColorBrewer::brewer.pal(7,'Set2')
names(groupcolor) <- unique(average_expression1$subcelltype) #类型颜色
Sexcolor <- c("#ff7f0e","#1f77b4") 
names(Sexcolor) <- c("iMCD","NC") #类型颜色
ann_colors <- list(subcelltype=groupcolor, disease=Sexcolor) #颜色设置

p1 = pheatmap::pheatmap(average_expression1[c(3:dim(average_expression1)[2])],
                        annotation_row = average_expression1[1:2],
                        #color = colorRampPalette(c("navy","white","firebrick3"))(100),
                        annotation_colors = ann_colors,
                        scale = 'column',
                        cluster_cols = F,cluster_rows = F,show_rownames =F,
                        fontfamily = "serif")
ggsave(filename = './Fig5.2.TGFb_receptor_ligen_heatmap.png',p1,w=4,h=4,dpi = 200)
ggsave(filename = './Fig5.2.TGFb_receptor_ligen_heatmap.pdf',p1,w=5,h=5,dpi = 200)
##3筛选出SMAD集的表达数据##############
gene_set = c('FOS','DDX5','JUN','TRIM33','DROSHA','ZC3H3','PPM1A','SMURF1')
dotpltg_cellgp(gene_set,'./Fig5.3.SMAD_dotplot',7,5)

##5细胞周期基因###########
gene_set1 = c('CDK6','CDKN1A','CDKN2A')
dotpltg_cellgp(gene_set1,'./Fig5.5.cycling_gene_dotplot',with = 7,high = 3)

##4 7个群体R-SMAD通路富集GSVA###########
library(clusterProfiler)
library(org.Hs.eg.db)
HSC_MPP <-read.csv('/mnt/NexData_file/home/zwl/Project/FW2024-225-G/04.scanpy_sub/HSPC_group/clusterHSC/MPP_iMCD/differentially_expressed_genes.xls',
                   sep = '\t')
EP <- read.csv('/mnt/NexData_file/home/zwl/Project/FW2024-225-G/03.scanpy/group/clusterErythroid progenitor_iMCD/differentially_expressed_genes.xls',
               sep = '\t')

geneHSC = bitr(HSC_MPP$names, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
HSC_data = merge(geneHSC,HSC_MPP,by.y = 'names' ,by.x='SYMBOL')%>%
  arrange(desc(logfoldchanges))
geneList = HSC_data$logfoldchanges
names(geneList) = HSC_data$ENTREZID
ego3 <- gseGO(geneList     = geneList,
              OrgDb        = org.Hs.eg.db,
              ont          = "MF",
              minGSSize    = 100,
              maxGSSize    = 500,
              pvalueCutoff = 0.1,
              verbose      = FALSE)

geneEP = bitr(EP$names, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
EP_data = merge(geneEP,EP,by.y = 'names' ,by.x='SYMBOL')%>%
  arrange(desc(logfoldchanges))
geneList = EP_data$logfoldchanges
names(geneList) = EP_data$SYMBOL
ego3 <- gseGO(geneList     = geneList,
              OrgDb        = org.Hs.eg.db,
              ont          = "MF",
              keyType      = 'SYMBOL',
              pvalueCutoff = 0.1,
              verbose      = FALSE)

gseaplot2(ego3, geneSetID = 1, title = ego3$Description[1])
R_SMAD =c("PARP1","LDLRAD4","DDX5","FOS","JUN","SMAD2","SMAD3","SMAD4","SMAD6","MEN1","PAX6","PPARG","PPM1A","AXIN1","RANBP3", 
          "FOXH1","STUB1","ZC3H3","ANKRD1","RGCC","DROSHA","TRIM33","PMEPA1","SMURF1","MYOCD") 
obj <-AddModuleScore(obj,features = list(R_SMAD),name = 'R-SMAD Binding' )
FeaturePlot(obj,features = 'R-SMAD Binding1' )
Seurat::DotPlot(obj,features = 'R-SMAD Binding1',group.by = c('subcelltype'))

gene_set = 'R-SMAD Binding1'
dotpltg_cellgp(gene_set,'./Fig5.4.R-SMAD Binding',8,3)



























