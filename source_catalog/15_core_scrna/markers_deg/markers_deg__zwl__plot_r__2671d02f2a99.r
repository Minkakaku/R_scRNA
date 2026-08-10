setwd('/mnt/NexData_file/home/zwl/Project/FW2024-629/require/')
rm(list = ls())
gc()
# ###########################################exchange h5ad to rds############################################
library(sceasy)
library(Seurat)
library(reticulate)

sceasy::convertFormat('/mnt/NexData_file/home/zwl/Project/FW2024-629/Cardio_raw.h5ad', from="anndata", to="seurat",
                      outFile='Cardio_result.rds')
obj <- readRDS(file = '../Cardio_result.rds')
obj$cell_type = paste0('C',obj$leiden)
Seurat::DimPlot(obj,reduction = 'umap',group.by = c('cell_type','sample'))
#需求1############################################
genelist = list('Parthanatos' = c('Afg3l2','Mif','Parp1'),
                'Death of iron' = c('Acsl4','Ptgs2','Hmox1'),
                'MPT_necrosis' = c('Aifm1','Diablo','Ppif' ),
                'Necroptosis' = c('Ripk1','Ripk3','MlKl'),
                'pyroptosis' = c( 'Gsdmd','Pycard','Nlrp3'),
                'Apoptosis' = c('Apaf1','Casp9','Casp8','Casp3'))
p = Seurat::DotPlot(obj,features = genelist,group.by = 'cell_type')+
  scale_color_gradientn(colours = c("navy",'white',"firebrick")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        plot.background = element_rect(fill = "white"))
ggsave(filename = '1.geneset_dotplot.png',w=12,h=4,p,dpi = 200)
ggsave(filename = '1.geneset_dotplot.pdf',w=12,h=4,p,dpi = 200)
#需求2############################################
library(stringr)
for (i in unique(obj$leiden)) {
  GO_BP = read.csv(file = paste0('/mnt/NexData_file/home/zwl/Project/FW2024-629/Cardiomyocytes_0.1/cluster',i,
                                 '/GO_Biological_Process_2023.Mouse.enrichr.reports.txt'),sep = '\t')
  top5_GO = GO_BP %>% slice_min(Adjusted.P.value,n = 5) %>% head(5) %>%
    dplyr::mutate(nCount = as.numeric(str_split_fixed(Overlap, "/", 2)[, 1]))
  p = ggplot(top5_GO,aes(x=nCount,y=Term))+
    geom_bar(stat = "identity",fill = "skyblue")+theme_classic()
  ggsave(filename = paste0('2.Cardiomyocytes_cluster',i,'_barplot.png'),p,w=8,h=4,dpi = 200)
  ggsave(filename = paste0('2.Cardiomyocytes_cluster',i,'_barplot.pdf'),p,w=8,h=4,dpi = 200)
}

#需求3############################################
library(monocle3)
mono3_obj = readRDS(file = '../monocle3/Cardiomyocytes_0.1/allClusterCDS.rds')
cytotrace = read.csv(file = '../cytotrace/Cardiomyocytes_0.1/CytoTRACE.value.txt',sep = '\t',row.names = 1)
mono3_obj@colData@listData$cytotrace = cytotrace[mono3_obj@colData@rownames,]
p1 <- plot_cells(mono3_obj,color_cells_by = "cytotrace",
                 show_trajectory_graph = T,label_leaves=F,
                 label_branch_points=FALSE) + 
  ggtitle("cytotrace score")
ggsave(filename = '3.monocle_cytotrace_trajectory.png',w=5,h=5,dpi = 200,p1)
ggsave(filename = '3.monocle_cytotrace_trajectory.pdf',w=5,h=5,dpi = 200,p1)

##需求3.2############################################
library(ggplot2)
mono3_obj <- order_cells(mono3_obj)
#获取pseudotime
mono3_pseudo = data.frame(mono3_obj@principal_graph_aux[["UMAP"]]$pseudotime)
geneset = read.csv(file = '../3.2_geneset_score.csv',row.names = 1)
geneset_n = colnames(geneset)
geneset$mono3_pseudo = mono3_pseudo[rownames(geneset),]

for (i in geneset_n) {
  p2 = ggplot(geneset, aes(x = mono3_pseudo, y = !!sym(i))) +
    geom_smooth(method = "loess",fill = "grey",alpha = 0.5) + theme_classic()
  ggsave(filename = paste0('3.2_',i,'_geneset_pseudotime_smooth.png'),w=4,h=4,p2,dpi = 200)
  ggsave(filename = paste0('3.2_',i,'_geneset_pseudotime_smooth.pdf'),w=4,h=4,p2,dpi = 200)
}

#需求4############################################
shrink_g = c('Tnnc1','Mybpc3','Tnni3','Myh6','Myl2','Myl3','Atp2a2','Ryr2','Hrc',
             'Casq2','Cacna2d1','Cacna2d3','Kcnb1','Kcnd3','Scn5a')
metabolism = c('Acsl1','Ogdh','Sdhb','Sdhc','Mdh2','Sucla2','Dlst','Ogdhl','Aco2','Dlat')
obj <- ScaleData(obj, features = rownames(obj))
p41 = Seurat::DoHeatmap(obj,features = shrink_g,group.by = 'cell_type')
ggsave(filename = paste0('4.shrink_gene_heatmap.png'),w=6,h=4,p41,dpi = 200)
ggsave(filename = paste0('4.shrink_gene_heatmap.pdf'),w=6,h=4,p41,dpi = 200)
p42 = Seurat::DoHeatmap(obj,features = metabolism,group.by = 'cell_type')
ggsave(filename = paste0('4.metabolism_gene_heatmap.png'),w=6,h=4,p42,dpi = 200)
ggsave(filename = paste0('4.metabolism_gene_heatmap.pdf'),w=6,h=4,p42,dpi = 200)

#需求5############################################
library(dplyr)
DEG = read.csv(file = '../5.C1_DEG_plot.csv',row.names = 1)
DEG$threshold = factor(ifelse(DEG$pvals_adj>0.05|abs(DEG$logfoldchanges) <1,'Nosign',
                                 ifelse(DEG$logfoldchanges>1,'Up','Down')),
                          levels = c('Up','Nosign','Down'))
top_genes = DEG %>% dplyr::filter(threshold != 'Nosign') %>% 
  arrange(logfoldchanges) %>%
  slice(1:10, (n() - 9):n())
DEG$y_value <- -log10(DEG$pvals_adj)
max_y <- max(DEG$y_value[is.finite(DEG$y_value)]) * 1.1  # 设置 y 轴上限为最大值的 1.1 倍

gg = ggplot(DEG, aes(x = logfoldchanges, -log10(pvals_adj), colour = threshold)) +
  geom_jitter(size=0.8,stroke=0.5, alpha = 0.6) +
  scale_color_manual(values = c('Up'="firebrick", 'Down'="navy",'Nosign'='grey'))+
  theme_classic()+labs(x='logFC C1/C0',y='-log10 pvalue_adj')+
  coord_cartesian(ylim = c(0, max_y)) +
  ggrepel::geom_text_repel(data = top_genes, aes(label = names), size = 4, max.overlaps = Inf,box.padding = 0.5 ) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black")+
  geom_vline(xintercept = c(-1,1), linetype = "dashed", color = "black")
 
ggsave(filename = paste0('./5.C1_C0_group_DEG_Vcano.png'),gg,w=6,h=4,dpi = 200)
ggsave(filename = paste0('./5.C1_C0_group_DEG_Vcano.pdf'),gg,w=6,h=4,dpi = 200)
#需求6############################################
top5_gene = DEG %>% dplyr::filter(threshold != 'Nosign') %>%  slice_max(order_by = logfoldchanges,n=5)
cds_subset <- mono3_obj[rowData(mono3_obj)$gene_short_name %in% top5_gene$names,]
gg6 = plot_genes_violin(cds_subset, group_cells_by="leiden") +
  theme(axis.text.x=element_text(angle=45, hjust=1))
ggsave(filename = paste0('./6.C1_C0_top5DEG_cluster.png'),gg6,w=6,h=7,dpi = 200)
ggsave(filename = paste0('./6.C1_C0_top5DEG_cluster.pdf'),gg6,w=6,h=7,dpi = 200)

#差异基因中和拟时间相关的基因
DEG = read.csv(file = '5.C0_C1/all_cluster_marker.xls',sep = '\t')
#colData(mono3_obj)$pseudotime = data.frame(pseudotime(mono3_obj))[,1]
cds_subset <- mono3_obj[rowData(mono3_obj)$gene_short_name %in% DEG$Symbol,]
gene_fits <- monocle3::fit_models(cds_subset, model_formula_str = "~pseudotime")
fit_coefs <- monocle3::coefficient_table(gene_fits)
emb_time_terms <- fit_coefs %>% filter(term == "pseudotime")%>% filter (q_value < 0.05) %>%
  select(gene_short_name, term, q_value, estimate) %>% arrange(q_value)

gg6.2 = plot_genes_in_pseudotime(mono3_obj[head(emb_time_terms$gene_short_name,n=5),],
                         color_cells_by="leiden",
                         min_expr=0.5)
ggsave(filename = paste0('./6.C1_C0_top5_pseudotime_DEG_cluster.png'),gg6.2,w=6,h=7,dpi = 200)
ggsave(filename = paste0('./6.C1_C0_top5_pseudotime_DEG_cluster.pdf'),gg6.2,w=6,h=7,dpi = 200)

#需求7############################################
topgene = c(head(emb_time_terms$gene_short_name,n=5),top5_gene$names) %>% unique()
gene_module_df <- find_gene_modules(mono3_obj[topgene,], resolution=0.001)
g7 = plot_cells(mono3_obj,genes = topgene,label_cell_groups=FALSE,show_trajectory_graph=T,
           label_leaves=F,label_branch_points=FALSE)
ggsave(filename = paste0('./7.C1_C0_pseudotime_DEG_UMAP.png'),g7,w=8,h=8,dpi = 200)
ggsave(filename = paste0('./7.C1_C0_pseudotime_DEG_UMAP.pdf'),g7,w=9,h=9,dpi = 200)
