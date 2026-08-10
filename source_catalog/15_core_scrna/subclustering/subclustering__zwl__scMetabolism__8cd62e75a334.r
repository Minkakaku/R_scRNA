setwd('/mnt/NexData_file/home/zwl/Project/FWSC20242216/require4/scMetabolism')
library(scMetabolism)
library(ggplot2)
library(rsvd)
library(Seurat)
library(dplyr)
library(rstatix)
rm(list = ls())
gc()


#h5ad转变成seurat对象####################
library(sceasy)
library(reticulate)
sceasy::convertFormat('/mnt/NexData_file/home/zwl/Project/FWSC20242216/require4/Tcell_result1.h5ad', from="anndata", to="seurat",
                      outFile='T_cell.rds')
sceasy::convertFormat('/mnt/NexData_file/home/zwl/Project/FWSC20242216/require4/NK_raw.h5ad', from="anndata", to="seurat",
                      outFile='NK_cell.rds')
obj = readRDS('/mnt/NexData_file/home/zwl/Project/FWSC20242216/require4/T_cell.rds')
#obj = readRDS('/mnt/NexData_file/home/zwl/Project/FWSC20242216/require4/NK_cell.rds')
obj = FindVariableFeatures(obj, nfeatures = 2000)
DimPlot(obj,reduction = 'umap',group.by = c('cell_type','group'))

#scMetabolism###############################
countexp.Seurat<-sc.metabolism.Seurat(obj = obj, 
                                      method = "AUCell",
                                      imputation = F, 
                                      ncores = 2,
                                      metabolism.type = "KEGG")
saveRDS(file = './scMetabolism_KEGG_Tcell.rds',compress = F,countexp.Seurat)


##提取结果,计算组间差异
res = cbind(t(countexp.Seurat@assays[["METABOLISM"]][["score"]]),
            countexp.Seurat@meta.data[,c('cell_type','sample','group')])
write.csv(file = '0.METABOLISM_AUC.csv',res)

avg_ALL <- res %>% 
  group_by(sample,group, cell_type) %>% 
  summarise(across(where(is.numeric), mean), .groups = 'drop')

plot_avg = tidyr::pivot_longer(avg_ALL,cols = !c('sample','cell_type','group'),
                               names_to = 'Pathway',
                               values_to = 'Auc_score')
p_value = plot_avg%>%
  group_by(Pathway,cell_type) %>% 
  rstatix::t_test(formula = Auc_score~group)%>%
  add_significance() %>% add_xy_position()
write.csv(file = '1.Tcell_METABOLISM_celltype_meanAUC.csv',avg_ALL)
write.csv(file = '1.Tcell_METABOLISM_celltype_group_P.csv',p_value[,1:10])


##plot diff pathway#######################
pathway = unique(p_value[p_value$p<0.05,]$Pathway)

library(RColorBrewer)
groupcolor <- RColorBrewer::brewer.pal(length(unique(avg_ALL$cell_type)),'Set3')
names(groupcolor) <- unique(avg_ALL$cell_type) #类型颜色
agecolor <- c("#1f77b4",'firebrick3') 
names(agecolor) <- c('S','F') #类型颜色
ann_colors <- list(cell_type=groupcolor, group=agecolor) #颜色设置

avg_res1 = data.frame(avg_ALL[,c(pathway,'group','sample','cell_type')])
rownames(avg_res1) = paste0(avg_res1$sample,'_',avg_res1$cell_type)
avg_res1 <- avg_res1[order(avg_res1$cell_type), ]
p1 = pheatmap::pheatmap(t(avg_res1[1:(dim(avg_res1)[2]-3)]),
                        annotation_col = avg_res1[,c('group','cell_type')],
                        #color = colorRampPalette(c("navy","white","firebrick3"))(100),
                        annotation_colors = ann_colors,
                        scale = 'row',#gaps_row = c(8,16),
                        fontsize = 12,
                        cluster_cols = F,cluster_rows = T,show_colnames =T)
ggsave(filename = './1.Tcell_scMetabolism_scale_heatmap.png',p1,w=18,h=6,dpi = 200)
ggsave(filename = './1.Tcell_scMetabolism_scale_heatmap.pdf',p1,w=20,h=6,dpi = 200)

library(clusterProfiler)
library(org.Hs.eg.db)
p_value = read.csv('1.Tcell_METABOLISM_celltype_group_P.csv',row.names = 1)
obj <-readRDS('../T_cell.rds')

pathway = p_value %>% dplyr::filter(p<0.05) %>%pull(Pathway) %>% unique()
hsa_kegg <- clusterProfiler::download_KEGG("hsa")

data_kegg =  merge(hsa_kegg[["KEGGPATHID2EXTID"]], hsa_kegg[["KEGGPATHID2NAME"]], by="from")%>%
  rename(ID = from,ENTREZID = to.x,Term = to.y)%>% 
  dplyr::filter(Term %in% pathway)
gene_symbols <- bitr(data_kegg$ENTREZID, fromType="ENTREZID", toType="SYMBOL", OrgDb="org.Hs.eg.db")
data_kegg1 = merge(data_kegg,gene_symbols,by = 'ENTREZID')
write.csv(file = '2.p0.05_pathway_genelist.csv',data_kegg1)

pdf(file = './2.T cells_diff_pathway_genelist.pdf',width = 12,height = 6)
for (i in pathway) {
  genelist = data_kegg1 %>% dplyr::filter(Term == i)%>%pull(SYMBOL)
  p = DotPlot(obj,features = genelist,group.by = 'group')+
    theme(axis.text.x = element_text(angle = 30,hjust = 1))+labs(y='T cells',x=i)
  print(p)
}
dev.off()

objTex = subset(obj,cell_type == 'CD8 Tex')
pdf(file = './2.CD8Tex cells_diff_pathway_genelist.pdf',width = 12,height = 6)
for (i in c("Phosphonate and phosphinate metabolism","Propanoate metabolism")) {
  genelist = data_kegg1 %>% dplyr::filter(Term == i)%>%pull(SYMBOL)
  p = DotPlot(objTex,features = genelist,group.by = 'group')+
    theme(axis.text.x = element_text(angle = 30,hjust = 1))+labs(y='CD8Tex cells',x=i)
  print(p)
}
dev.off()



#TEX###############################
avg_res <- res %>% dplyr::filter(cell_type == 'CD8 Tex')%>%
  group_by(sample,group, cell_type) %>%
  summarise(across(where(is.numeric), mean), .groups = 'drop')

plot_avg = tidyr::pivot_longer(avg_res,cols = !c('cell_type','sample','group'),
                               names_to = 'Pathway',
                               values_to = 'Auc_score')
p_value = plot_avg%>%
  group_by(Pathway,cell_type) %>%
  rstatix::t_test(formula = Auc_score~group)%>%
  add_significance() %>% add_xy_position()
write.csv(file = '1.CD8TEX_METABOLISM_celltype_meanAUC.csv',avg_res)
write.csv(file = '1.CD8TEX_METABOLISM_celltype_group_P.csv',p_value[,1:10])


#plot diff pathway
pathway = p_value[p_value$p<0.05,]$Pathway

library(RColorBrewer)
agecolor <- c("#1f77b4",'firebrick3') 
names(agecolor) <- c('S','F') #类型颜色
ann_colors <- list(group=agecolor) #颜色设置


avg_res1 = data.frame(avg_res[,c(pathway,'group','sample')])
rownames(avg_res1) = avg_res$sample

anno_g = data.frame(group=avg_res1[,'group'])
rownames(anno_g) = avg_res1[,'sample']
p1 = pheatmap::pheatmap(t(avg_res1[1:(dim(avg_res1)[2]-2)]),
                        annotation_col = anno_g,
                        #color = colorRampPalette(c("navy","white","firebrick3"))(100),
                        annotation_colors = ann_colors,
                        scale = 'row',#gaps_row = c(8,16),
                        fontsize = 12,
                        cluster_cols = F,cluster_rows = T,show_colnames =T)
ggsave(filename = './1.CD8TEX_scMetabolism_scale_heatmap.png',p1,w=10,h=4,dpi = 200)
ggsave(filename = './1.CD8TEX_scMetabolism_scale_heatmap.pdf',p1,w=12,h=4,dpi = 200)


# 要求pathway gene dotplot####################################
obj <-readRDS('scMetabolism_KEGG_Tcell.rds')
subobj = subset(obj,cell_type == 'NK')
NK_cell = read.csv(file = '/mnt/NexData_file/home/zwl/Project/FWSC20242216/require5/cell_type.csv',row.names = 1)
subobj$cell_type = NK_cell[colnames(subobj),'cell_type'] 
DimPlot(subobj,reduction = 'umap',group.by = 'cell_type')
subobj$group_cell = paste0(subobj$cell_type,'_',subobj$group)
  
pathway = c('Glycolysis / Gluconeogenesis',
            'Citrate cycle (TCA cycle)',
            'Pentose phosphate pathway',
            'Pentose and glucuronate interconversions',
            'Fructose and mannose metabolism',
            'Galactose metabolism',
            'Pyruvate metabolism',
            'Inositol phosphate metabolism',
            'Oxidative phosphorylation',
            'Glycosylphosphatidylinositol (GPI)-anchor biosynthesis',
            'N-Glycan biosynthesis',
            'Various types of N-glycan biosynthesis',
            'Mucin type O-glycan biosynthesis',
            'Mannose type O-glycan biosynthesis',
            'Other types of O-glycan biosynthesis',
            
            'Fatty acid biosynthesis',
            'Fatty acid elongation',
            'Fatty acid degradation',
            'Glycerolipid metabolism',
            'Glycerophospholipid metabolism',
            'Sphingolipid metabolism',
            'Arachidonic acid metabolism',
            'Linoleic acid metabolism',
            'alpha-Linolenic acid metabolism',
            'Biosynthesis of unsaturated fatty acids')


hsa_kegg <- clusterProfiler::download_KEGG("hsa")
data_kegg =  merge(hsa_kegg[["KEGGPATHID2EXTID"]], hsa_kegg[["KEGGPATHID2NAME"]], by="from")%>%
  rename('from'='ID','to.x'='ENTREZID','to.y'='Term')%>% 
  dplyr::filter(Term %in% pathway)
gene_symbols <- bitr(data_kegg$ENTREZID, fromType="ENTREZID", toType="SYMBOL", OrgDb="org.Hs.eg.db")
data_kegg1 = merge(data_kegg,gene_symbols,by = 'ENTREZID')

pdf(file = '../../require5/2.T cells_special_pathway_genelist.pdf',width = 35,height = 6)
for (i in pathway) {
  genelist = data_kegg1 %>% dplyr::filter(Term == i)%>%pull(SYMBOL)
  p = DotPlot(subobj,features = genelist,group.by = 'group_cell')+
    theme(axis.text.x = element_text(angle = 30,hjust = 1))+labs(x=i,y=NULL)+
    scale_color_gradientn(colours = c('navy','white','firebrick3'))
  print(p)
}
dev.off()

## pathway pvalue ##########################################################
NK_meta = read.csv(file = '../../require6/metadata.csv',row.names = 1)

rownames(NK_meta) %in% rownames(obj@meta.data)
obj$cell_type <- as.character(obj$cell_type)
obj@meta.data[rownames(NK_meta),'cell_type'] = NK_meta$cell_type
unique(obj$cell_type)

res = cbind(t(obj@assays[["METABOLISM"]][["score"]][pathway,]),
            obj@meta.data[,c('cell_type','sample','group')])
res <- res[, colSums(is.na(res)) != nrow(res)]

plot_avg = tidyr::pivot_longer(res,cols = !c('sample','cell_type','group'),
                               names_to = 'Pathway',
                               values_to = 'Auc_score')
p_value = plot_avg%>%
  group_by(Pathway,cell_type) %>% 
  rstatix::t_test(formula = Auc_score~group)%>%
  add_significance() %>% add_xy_position()
write.csv(file = '0.T-NK_METABOLISM_celltype_cellAUC.csv',res)
write.csv(file = '0.T-NK_METABOLISM_celltype_cellAUC_group_P.csv',p_value[,1:11])


avg_ALL <- res %>% 
  group_by(sample,group, cell_type) %>% 
  summarise(across(where(is.numeric), mean), .groups = 'drop')

plot_avg = tidyr::pivot_longer(avg_ALL,cols = !c('sample','cell_type','group'),
                               names_to = 'Pathway',
                               values_to = 'Auc_score')

p_value = plot_avg%>%
  group_by(Pathway,cell_type) %>% 
  rstatix::t_test(formula = Auc_score~group)%>%
  add_significance() %>% add_xy_position()
write.csv(file = '0.T-NK_METABOLISM_celltype_meanAUC.csv',avg_ALL)
write.csv(file = '0.T-NK_METABOLISM_celltype_meanAUC_group_P.csv',p_value1[,1:11])

############获取 genelist 的表达矩阵#####################################
for (i in pathway) {
  genelist = data_kegg1 %>% dplyr::filter(Term == i)%>%pull(SYMBOL)
  
  # 1. 创建组合分组列
  obj$combo_group <- paste(obj$sample, obj$cell_type, obj$group, sep = "|")
  
  # 2. 计算平均值
  avg_expr <- AverageExpression(
    obj,
    assays = "RNA",
    features = genelist,  # 所有基因或指定基因列表
    #group.by = 'combo_group',
    slot = "data"  # 使用标准化数据
  )$RNA %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column("gene")
  write.csv(file = paste0('1.',gsub('/','-',i),'_genlist_meanby_sample_celltype.csv'),avg_expr)
  
  avg_expl = tidyr::pivot_longer(avg_expr,cols = !c('gene'),
                                 names_to = 'Pathway',
                                 values_to = 'gene_mean')%>% 
    separate(col = 'Pathway', into = c("sample", "cell_type", "group"), sep = "\\|")
  p_valueg = avg_expl%>%
    group_by(gene,cell_type) %>% 
    rstatix::t_test(formula = gene_mean~group)%>%
    add_significance() %>% add_xy_position()
}

for (i in pathway) {
  genelist = data_kegg1 %>% dplyr::filter(Term == i,SYMBOL %in% rownames(obj))%>%pull(SYMBOL)
  gene_exp = GetAssayData(obj,
                          assays = "RNA",
                          layer =  "counts"  # 使用标准化数据
                          )[genelist,] %>% 
    as.matrix() %>% t()%>%  
    as.data.frame()
  
  # 3. 合并元数据
  combined_data <- cbind(
    gene_exp,
    obj@meta.data[, c("sample", "cell_type", "group")]
  )
  write.csv(file = paste0('2.',gsub('/','-',i),'_genlist_cellexp.csv'),combined_data)
  
  avg_expl = tidyr::pivot_longer(combined_data,cols = !c("sample", "cell_type", "group"),
                                 names_to = 'Pathway_gene',
                                 values_to = 'gene_exp')
  p_valueg1 = avg_expl%>%
    group_by(Pathway_gene,cell_type) %>% 
    rstatix::t_test(formula = gene_exp~group)%>%
    add_significance() %>% add_xy_position() 
  write.csv(file = paste0('2.',gsub('/','-',i),'_genlist_cellexp_p.csv'),p_valueg1[,1:11])
}

# #计算每个细胞类型组间差异(cell)
# pathways <- setdiff(colnames(data.frame(res)), c("cell_type", "group"))
# 
# # 对每列循环计算
# p_values_list <- lapply(pathways, function(p) {
#   data.frame(res) %>%
#     group_by(cell_type) %>%
#     t_test(reformulate("group", response = p)) %>%
#     add_significance() %>%
#     mutate(pathway = p)  # 添加列名标识
# })
# 
# # 合并结果
# p_values <- bind_rows(p_values_list)
# write.csv(file = '2.METABOLISM_celltype_group_cellDEG_P.csv',p_values)

# #计算每个细胞类型组间差异(sample)
# avg_res <- res %>% 
#   group_by(sample,group, cell_type) %>% 
#   summarise(across(where(is.numeric), mean), .groups = 'drop')
# 
# plot_avg = tidyr::pivot_longer(avg_res,cols = !c('cell_type','sample','group'),
#                                names_to = 'Pathway',
#                                values_to = 'Auc_score')
# p_value = plot_avg%>%
#   group_by(Pathway,cell_type) %>% 
#   rstatix::t_test(formula = Auc_score~group)%>%
#   add_significance() %>% add_xy_position()

