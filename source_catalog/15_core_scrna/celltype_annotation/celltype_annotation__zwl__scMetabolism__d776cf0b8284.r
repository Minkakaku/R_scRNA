setwd('/mnt/NexData_file/home/zwl/Project/FWAE20250109/require3/scMetabolism/')
library(scMetabolism)
library(ggplot2)
library(rsvd)
library(Seurat)
rm(list = ls())
gc()
HSC = readRDS('/mnt/NexData_file/home/zwl/Project/FWAE20250109/require2/HSPC_celltype.rds')
DimPlot(HSC,reduction = 'umap',group.by = 'cell_type')

countexp.Seurat<-sc.metabolism.Seurat(obj = HSC, 
                                      method = "AUCell",
                                      imputation = F, 
                                      ncores = 2,
                                      metabolism.type = "KEGG")
saveRDS(file = './scMetabolism_KEGG_HSC.rds',compress = F,countexp.Seurat)
DimPlot.metabolism(obj = countexp.Seurat, 
                   pathway = pathway,
                   dimention.reduction.type = "umap", 
                   dimention.reduction.run = F, 
                   size = 1)

pathway = c('Glycolysis / Gluconeogenesis',
            'Pentose phosphate pathway',
            'Pyruvate metabolism',
            
            'Oxidative phosphorylation',
            "Parkinson's disease",
            "Huntington's disease",
            
            #"Fatty acid biosynthesis",
            'Fatty acid degradation',
            #'Steroid biosynthesis',
            'Glycerolipid metabolism',
            'Glycerophospholipid metabolism',
            
            'Arginine biosynthesis',
            'Alanine, aspartate, and glutamate metabolism',
            'Glycine, serine, and threonine metabolism',
            'Cysteine and methionine metabolism',
            'Valine, leucine, and isoleucine degradation',
            'Arginine and proline metabolism',
            'Histidine metabolism',
            'Tyrosine metabolism',
            'Phenylalanine metabolism',
            'Phenylalanine, tyrosine, and tryptophan biosynthesis')
pathway = pathway[pathway %in% rownames(countexp.Seurat@assays[["METABOLISM"]][["score"]])]
DotPlot.metabolism(obj = countexp.Seurat, pathway = pathway, phenotype = "cell_type", norm = "y")
#提取结果##################
res = cbind(t(countexp.Seurat@assays[["METABOLISM"]][["score"]][pathway,]),
            countexp.Seurat@meta.data[,c('cell_type','age')])
avg_res = res %>% group_by(age,cell_type) %>% 
  summarise(across(all_of(pathway), mean))

# plot_avg = tidyr::pivot_longer(avg_res,cols = !c('cell_type','age'),
#                                names_to = 'Pathway',
#                                values_to = 'Auc_score')
# plot_avg$age = factor(plot_avg$age,levels = c('Young','Middle','Elderly'))

library(RColorBrewer)
groupcolor <- RColorBrewer::brewer.pal(8,'Set2')
names(groupcolor) <- unique(avg_res$cell_type) #类型颜色
agecolor <- c("#1f77b4","#ff7f0e",'firebrick3') 
names(agecolor) <- c('Young','Middle','Elderly') #类型颜色
ann_colors <- list(cell_type=groupcolor, age=agecolor) #颜色设置

avg_res = data.frame(avg_res)
rownames(avg_res) = paste0(avg_res$age,'_',avg_res$cell_type)
p1 = pheatmap::pheatmap(avg_res[3:dim(avg_res)[2]],
                        annotation_row = avg_res[1:2],
                        #color = colorRampPalette(c("navy","white","firebrick3"))(100),
                        annotation_colors = ann_colors,
                        scale = 'column',gaps_row = c(8,16),
                        fontsize = 12,
                        cluster_cols = F,cluster_rows = F,show_rownames =F)
ggsave(filename = './scMetabolism_scale_heatmap.png',p1,w=10,h=10,dpi = 200)
ggsave(filename = './scMetabolism_scale_heatmap.pdf',p1,w=12,h=12,dpi = 200)



res = cbind(t(countexp.Seurat@assays[["METABOLISM"]][["score"]][pathway,]),
            countexp.Seurat@meta.data[,c('cell_type','age')])
avg_res = res %>% group_by(age) %>% 
  summarise(across(all_of(pathway), mean))
agecolor <- c("#1f77b4","#ff7f0e",'firebrick3') 
names(agecolor) <- c('Young','Middle','Elderly') #类型颜色
ann_colors <- list(age=agecolor) #颜色设置

avg_res = data.frame(avg_res)
rownames(avg_res) = avg_res$age
p1 = pheatmap::pheatmap(avg_res[2:dim(avg_res)[2]],
                        annotation_row = avg_res[1],
                        #color = colorRampPalette(c("navy","white","firebrick3"))(100),
                        annotation_colors = ann_colors,
                        scale = 'column',
                        fontsize = 12,
                        cluster_cols = F,cluster_rows = F,show_rownames =F)
ggsave(filename = './scMetabolism_scale_heatmap_age.png',p1,w=10,h=5,dpi = 200)
ggsave(filename = './scMetabolism_scale_heatmap_age.pdf',p1,w=12,h=6,dpi = 200)


res = cbind(t(countexp.Seurat@assays[["METABOLISM"]][["score"]][pathway,]),
            countexp.Seurat@meta.data[,c('cell_type','age')])
avg_res = res %>% group_by(age) %>% 
  summarise(across(all_of(pathway), mean))


res = cbind(t(countexp.Seurat@assays[["METABOLISM"]][["score"]][pathway,]),
            countexp.Seurat@meta.data[,c('type','age')])
avg_res = res %>% group_by(age,type) %>% 
  summarise(across(all_of(pathway), mean))
library(RColorBrewer)
groupcolor <- RColorBrewer::brewer.pal(8,'Set2')[1:2]
names(groupcolor) <- unique(avg_res$type) #类型颜色
agecolor <- c("#1f77b4","#ff7f0e",'firebrick3') 
names(agecolor) <- c('Young','Middle','Elderly') #类型颜色
ann_colors <- list(type=groupcolor, age=agecolor) #颜色设置

avg_res = data.frame(avg_res)
rownames(avg_res) = paste0(avg_res$age,'_',avg_res$type)
p1 = pheatmap::pheatmap(avg_res[3:dim(avg_res)[2]],
                        annotation_row = avg_res[1:2],
                        #color = colorRampPalette(c("navy","white","firebrick3"))(100),
                        annotation_colors = ann_colors,
                        scale = 'column',gaps_row = c(2,4),
                        fontsize = 12,
                        cluster_cols = F,cluster_rows = F,show_rownames =F)
ggsave(filename = './scMetabolism_scale_heatmap_age_type.png',p1,w=10,h=7,dpi = 200)
ggsave(filename = './scMetabolism_scale_heatmap_age_type.pdf',p1,w=12,h=8,dpi = 200)
