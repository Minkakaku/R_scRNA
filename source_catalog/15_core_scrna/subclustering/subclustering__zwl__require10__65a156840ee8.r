setwd('/mnt/NexData_file/home/zwl/Project/SC2511076348/require10')
rm(list=ls())
gc()
library(monocle)
HSMM = readRDS(paste0('/mnt/NexData_file/home/zwl/Project/SC2511076348/require8/3.Myeloid_monocle-raw/monocle.rds'))
deg <- read.csv(paste0('/mnt/NexData_file/home/zwl/Project/SC2511076348/require8/3.Myeloid_monocle-raw/monocle.metadata.xls'),sep="\t")
rownames(deg) = deg$Genes
ordergene <- rownames(deg)[1:min(1000,nrow(deg))]
# 重新排序
plot_cell_trajectory(HSMM, color_by = "State")+
  facet_wrap(~cell_type,ncol = 4)
HSMM_myo <- orderCells(HSMM, root_state = 3)

# 重新画图
setwd('/mnt/NexData_file/home/zwl/Project/SC2511076348/require8/3.Myeloid_monocle-raw/')
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
#########################gene
gs = c('FABP4')
cds_subset = HSMM_myo[gs,]
p = plot_genes_in_pseudotime(cds_subset, color_by = "cell_type")
ggsave("3.FABP4_monocle.trajectory.tree.faceted.pdf",p,width = 6,height = 4)
ggsave("3.FABP4_monocle.trajectory.tree.faceted.png",p,width = 6,height = 4)

library(ggsci)
pData(HSMM_myo)$FABP4 =exprs(HSMM_myo)['FABP4',]#+1)
p1=plot_cell_trajectory(HSMM_myo, color_by = "FABP4") + scale_color_gsea()
ggsave("3.FABP4_monocle.trajectory.pdf",p1,width = 5,height = 4)
ggsave("3.FABP4_monocle.trajectory.png",p1,width = 5,height = 4)


pseudotime_df = HSMM_myo@phenoData@data[,c("sample",'group','cell_type',"FABP4","Pseudotime","State")]
df = data.frame(Cells=rownames(pseudotime_df),pseudotime_df)
write.table(df,"4.FABP4_monocle.metadata_reorder.xls",sep = "\t",quote = F,row.names = F)

Blank_data = pseudotime_df %>% dplyr::filter(State %in% c(2,3,4,1))
cor_res <- cor.test(Blank_data$FABP4, Blank_data[['Pseudotime']])

gg <- ggplot(Blank_data, aes(y = FABP4, x = .data[['Pseudotime']])) +
  geom_jitter(size = 0.5, stroke = 0.5, shape = 16, alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE, color = "red3", formula = y ~ x) +
  annotate("text", x = Inf, y = Inf,
           label = sprintf("r = %.2f\np = %.2g", 
                           cor_res$estimate, 
                           cor_res$p.value),
           hjust = 1.1, vjust = 1.1, size = 5) +
  labs(y = "FABP4", x = 'Pseudotime') +
  theme_classic(base_size = 12)
gg
ggsave("4.FABP4_Pseudotime_cor_state1.pdf",gg,width = 5,height = 4)
ggsave("4.FABP4_Pseudotime_cor_state1.png",gg,width = 5,height = 4)

Blank_data = pseudotime_df %>% dplyr::filter(State %in% c(2,3,4,5))
cor_res <- cor.test(Blank_data$FABP4, Blank_data[['Pseudotime']])

gg <- ggplot(Blank_data, aes(y = FABP4, x = .data[['Pseudotime']])) +
  geom_jitter(size = 0.5, stroke = 0.5, shape = 16, alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE, color = "red3", formula = y ~ x) +
  annotate("text", x = Inf, y = Inf,
           label = sprintf("r = %.2f\np = %.2g", 
                           cor_res$estimate, 
                           cor_res$p.value),
           hjust = 1.1, vjust = 1.1, size = 5) +
  labs(y = "FABP4", x = 'Pseudotime') +
  theme_classic(base_size = 12)
gg
ggsave("4.FABP4_Pseudotime_cor_state5.pdf",gg,width = 5,height = 4)
ggsave("4.FABP4_Pseudotime_cor_state5.png",gg,width = 5,height = 4)
#####################################Fibro########################################


HSMM = readRDS(paste0('/mnt/NexData_file/home/zwl/Project/SC2511076348/require8/4.Fibro_monocle-raw/monocle.rds'))
deg <- read.csv(paste0('/mnt/NexData_file/home/zwl/Project/SC2511076348/require8/4.Fibro_monocle-raw/monocle.metadata.xls'),sep="\t")
rownames(deg) = deg$Genes
ordergene <- rownames(deg)[1:min(1000,nrow(deg))]
# 重新排序
plot_cell_trajectory(HSMM, color_by = "State")+
  facet_wrap(~cell_type,ncol = 4)
HSMM_myo <- orderCells(HSMM, root_state = 1)

# 重新画图
setwd('/mnt/NexData_file/home/zwl/Project/SC2511076348/require8/4.Fibro_monocle-raw/')
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

gs = c('FABP4')
cds_subset = HSMM_myo[gs,]
p = plot_genes_in_pseudotime(cds_subset, color_by = "cell_type")
p
ggsave("3.FABP4_monocle.trajectory.tree.faceted.pdf",p,width = 6,height = 4)
ggsave("3.FABP4_monocle.trajectory.tree.faceted.png",p,width = 6,height = 4)

library(ggsci)
pData(HSMM_myo)$FABP4 =exprs(HSMM_myo)['FABP4',]#+1)
p1=plot_cell_trajectory(HSMM_myo, color_by = "FABP4") + scale_color_gsea()
ggsave("3.FABP4_monocle.trajectory.pdf",p1,width = 5,height = 4)
ggsave("3.FABP4_monocle.trajectory.png",p1,width = 5,height = 4)

pseudotime_df = HSMM_myo@phenoData@data[,c("sample",'group','cell_type',"FABP4","Pseudotime","State")]
df = data.frame(Cells=rownames(pseudotime_df),pseudotime_df)
write.table(df,"3.FABP4_monocle.metadata_reorder.xls",sep = "\t",quote = F,row.names = F)

Blank_data = pseudotime_df %>% dplyr::filter(State %in% c(1,2,3,5))
cor_res <- cor.test(Blank_data$FABP4, Blank_data[['Pseudotime']])

gg <- ggplot(Blank_data, aes(y = FABP4, x = .data[['Pseudotime']])) +
  geom_jitter(size = 0.5, stroke = 0.5, shape = 16, alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE, color = "red3", formula = y ~ x) +
  annotate("text", x = Inf, y = Inf,
           label = sprintf("r = %.2f\np = %.2g", 
                           cor_res$estimate, 
                           cor_res$p.value),
           hjust = 1.1, vjust = 1.1, size = 5) +
  labs(y = "FABP4", x = 'Pseudotime') +
  theme_classic(base_size = 12)
gg
ggsave("4.FABP4_Pseudotime_cor_state1235.pdf",gg,width = 5,height = 4)
ggsave("4.FABP4_Pseudotime_cor_state1235.png",gg,width = 5,height = 4)

Blank_data = pseudotime_df %>% dplyr::filter(State %in% c(1,7))
cor_res <- cor.test(Blank_data$FABP4, Blank_data[['Pseudotime']])

gg <- ggplot(Blank_data, aes(y = FABP4, x = .data[['Pseudotime']])) +
  geom_jitter(size = 0.5, stroke = 0.5, shape = 16, alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE, color = "red3", formula = y ~ x) +
  annotate("text", x = Inf, y = Inf,
           label = sprintf("r = %.2f\np = %.2g", 
                           cor_res$estimate, 
                           cor_res$p.value),
           hjust = 1.1, vjust = 1.1, size = 5) +
  labs(y = "FABP4", x = 'Pseudotime') +
  theme_classic(base_size = 12)
gg
ggsave("4.FABP4_Pseudotime_cor_state17.pdf",gg,width = 5,height = 4)
ggsave("4.FABP4_Pseudotime_cor_state17.png",gg,width = 5,height = 4)
########################################FABP4+—#################################
data1 = read.csv(file = '/mnt/NexData_file/home/zwl/Project/SC2511076348/require10/Myeloids/FABP4_exp.csv',row.names = 1)
target_score = c('HALLMARK_INFLAMMATORY_RESPONSE',"GOBP_OSTEOCLAST_DIFFERENTIATION","GOBP_REGULATION_OF_OSTEOCLAST_DIFFERENTIATION")

Blank_data = data1 %>% dplyr::filter(cell_type_new %in% c('FABP4+ Macro_Mono'))
setwd('/mnt/NexData_file/home/zwl/Project/SC2511076348/require10/Myeloids')

for(gene in target_score){
  
  # 1. 计算相关性
  cor_res <- cor.test(Blank_data$FABP4, Blank_data[[gene]])
  
  # 2. 绘图
  gg <- ggplot(Blank_data, aes(x = FABP4, y = .data[[gene]])) +
    geom_jitter(size = 0.5, stroke = 0.5, shape = 16, alpha = 0.6) +
    geom_smooth(method = "lm", se = TRUE, color = "red3", formula = y ~ x) +
    annotate("text", x = Inf, y = Inf,
             label = sprintf("r = %.2f\np = %.2g", 
                             cor_res$estimate, 
                             cor_res$p.value),
             hjust = 1.1, vjust = 1.1, size = 5) +
    labs(x = "FABP4", y = gene) +
    theme_classic(base_size = 12)
  
  # 3. 保存图片（自动命名）
  ggsave(filename = paste0('2.FABP4+Macro_Mono_', gene, '_cor_scatterplot.png'),
         dpi = 200, w=5, h=5, gg)
  ggsave(filename = paste0('2.FABP4+Macro_Mono_', gene, '_cor_scatterplot.pdf'),
         dpi = 200, w=5, h=5, gg)
}

################################################################################
##############相关性图2############################
################################################################################
library(readxl)
obj = readRDS('/mnt/NexData_file/home/zwl/Project/SC2511076348/require7/Myeloids/annoraw.rds')
score = read.csv('/mnt/NexData_file/home/zwl/Project/SC2511076348/require8/Fib_Mye_pathway_score.csv',row.names = 1)
target_score = c('HALLMARK_INFLAMMATORY_RESPONSE',"GOBP_OSTEOCLAST_DIFFERENTIATION","GOBP_REGULATION_OF_OSTEOCLAST_DIFFERENTIATION")

subdata = subset(obj,cell_type %in% 'Macro1')
target_genes = c('FABP4',"IL1B","TNF","NLRP3","NFKBIA")
expr_matrix <- Seurat::GetAssayData(object = subdata, assay = "RNA", slot = "data")[target_genes, ]
expr_df <- as.data.frame(t(expr_matrix))
final_df <- cbind(expr_df, subdata@meta.data[, c("cell_type", "group")],score[rownames(expr_df),target_score])

Blank_data = final_df %>% dplyr::filter(group =='DM',FABP4 >0)
setwd('/mnt/NexData_file/home/zwl/Project/SC2511076348/require10')

for(gene in c(target_genes[-1],target_score)){
  
  Blank_data1 = Blank_data[Blank_data[[gene]]>0,]
  # 1. 计算相关性
  cor_res <- cor.test(Blank_data1$FABP4, Blank_data1[[gene]])
  
  # 2. 绘图
  gg <- ggplot(Blank_data, aes(x = FABP4, y = .data[[gene]])) +
    geom_jitter(size = 0.5, stroke = 0.5, shape = 16, alpha = 0.6) +
    geom_smooth(method = "lm", se = TRUE, color = "red3", formula = y ~ x) +
    annotate("text", x = Inf, y = Inf,
             label = sprintf("r = %.2f\np = %.2g", 
                             cor_res$estimate, 
                             cor_res$p.value),
             hjust = 1.1, vjust = 1.1, size = 5) +
    labs(x = "FABP4", y = gene) +
    theme_classic(base_size = 12)
  
  # 3. 保存图片（自动命名）
  ggsave(filename = paste0('16.Macro1_DM_FABP4_vs_', gene, '_cor_scatterplot.png'),
         dpi = 200, w=5, h=5, gg)
  ggsave(filename = paste0('16.Macro1_DM_FABP4_vs_', gene, '_cor_scatterplot.pdf'),
         dpi = 200, w=5, h=5, gg)
}

subdata = subset(obj,cell_type %in% 'Macro2')
target_genes = c('FABP4',"IL1B","TNF","NLRP3","NFKBIA")
expr_matrix <- Seurat::GetAssayData(object = subdata, assay = "RNA", slot = "data")[target_genes, ]
expr_df <- as.data.frame(t(expr_matrix))
final_df <- cbind(expr_df, subdata@meta.data[, c("cell_type", "group")],score[rownames(expr_df),target_score])

Blank_data = final_df %>% dplyr::filter(group =='DM')
setwd('/mnt/NexData_file/home/zwl/Project/SC2511076348/require9')

for(gene in c(target_genes[-1],target_score)){
  
  # 1. 计算相关性
  cor_res <- cor.test(Blank_data$FABP4, Blank_data[[gene]])
  
  # 2. 绘图
  gg <- ggplot(Blank_data, aes(x = FABP4, y = .data[[gene]])) +
    geom_jitter(size = 0.5, stroke = 0.5, shape = 16, alpha = 0.6) +
    geom_smooth(method = "lm", se = TRUE, color = "red3", formula = y ~ x) +
    annotate("text", x = Inf, y = Inf,
             label = sprintf("r = %.2f\np = %.2g", 
                             cor_res$estimate, 
                             cor_res$p.value),
             hjust = 1.1, vjust = 1.1, size = 5) +
    labs(x = "FABP4", y = gene) +
    theme_classic(base_size = 12)
  
  # 3. 保存图片（自动命名）
  ggsave(filename = paste0('16.Macro2_DM_FABP4_vs_', gene, '_cor_scatterplot.png'),
         dpi = 200, w=5, h=5, gg)
  ggsave(filename = paste0('16.Macro2_DM_FABP4_vs_', gene, '_cor_scatterplot.pdf'),
         dpi = 200, w=5, h=5, gg)
}