setwd('/mnt/NexData_file/home/zwl/Project/SC2507233885/require7')
down_gene1 = read.table(file = '/mnt/NexData_file/home/zwl/Project/SC2507233885/require6/scenic/L2-3IT-MT/res/0.down_regulon_csv.txt')
down_gene2 = read.table(file = '/mnt/NexData_file/home/zwl/Project/SC2507233885/require6/scenic/L5-6-NP2/res/0.down_regulon_csv.txt')
intersect(down_gene1$V1,down_gene2$V1)
up_gene1 = read.table(file = '/mnt/NexData_file/home/zwl/Project/SC2507233885/require6/scenic/L2-3IT-MT/res/0.up_regulon_csv.txt')
up_gene2 = read.table(file = '/mnt/NexData_file/home/zwl/Project/SC2507233885/require6/scenic/L5-6-NP2/res/0.up_regulon_csv.txt')
intersect(up_gene1$V1,up_gene2$V1)


sceasy::convertFormat('/mnt/NexData_file/home/zwl/Project/SC2507233885/require7/raw_Neu.h5ad', from="anndata", to="seurat",
                      outFile='/mnt/NexData_file/home/zwl/Project/SC2507233885/require7/raw_Neu.rds')
rm(list=ls()) 
gc()

scRNA=readRDS("./raw_Neu.rds")
unique(scRNA$cell_type)
scRNA$g_celltype = paste0(scRNA$cell_type,'_',scRNA$group)

up_gene = read.table('/mnt/NexData_file/home/zwl/Project/SC2507233885/require7/0.commen_up_regulon_csv.txt',sep = '\t')
down_gene = read.table('/mnt/NexData_file/home/zwl/Project/SC2507233885/require7/0.commen_down_regulon_csv.txt',sep = '\t')


for (gi in 1:length(down_gene$V1)) {
  genelist = unlist(strsplit(down_gene[gi,'V2'], split = ","))
  scRNA <- AddModuleScore(  object = scRNA,  features = list(genelist),  
                                  name = paste0(down_gene[gi,'V1'],' target gene score'))
  plot = scRNA@meta.data# %>% dplyr::filter(cell_type %in% c('L5/6 NP2','L2-3 IT-MT'))

  p = ggplot(plot, aes(x = group,
                       y = .data[[paste0(down_gene[gi,'V1'],' target gene score1')]])) +
    # 小提琴图
    geom_violin(aes(fill = group), alpha = 1, width = 0.8) +
    # 抖动点
    geom_jitter(
      color = 'black',
      size = 0.1,
      alpha = 0.5
    ) +
    # 内层箱线图（窄一点更好看）
    geom_boxplot(width = 0.2,  color = "black", outlier.shape = NA) +
    # P 值放在正中间
    stat_compare_means(
      method = "t.test",
      label.x = 1.5,
      label = "p",
      size = 4.5,
      vjust = 0.5
    ) +
    labs(
      x = "Group",
      y = paste0(down_gene[gi,'V1'],' target gene score'),
      fill = "Group"
    ) +
    theme_classic()
  p
  # ggsave(paste0(down_gene[gi,'V1'],' target gene score_volcano.pdf'), plot = p, width = 4, height = 6)
  # ggsave(paste0(down_gene[gi,'V1'],' target gene score_volcano.png'), plot = p, width = 4, height = 6)
  
  p2 = ggplot(plot, aes(x = sample,
                       y = .data[[paste0(down_gene[gi,'V1'],' target gene score1')]],
                       group = sample)) +
    geom_violin(aes(fill = group)) +
    geom_jitter(
      aes(group = sample),
      color = 'black',
      size = 0.1,
      alpha = 0.8
    )+
    geom_boxplot(width = 0.2,  color = "black", outlier.shape = NA) +
    labs(
      x = "sample",
      y = paste0(down_gene[gi,'V1'],' target gene score'),
      color = "Group"
    ) +
    theme_classic()
  p2
  # ggsave(paste0(down_gene[gi,'V1'],' target gene score_sample_volcano.pdf'), plot = p2, width = 6, height = 6)
  # ggsave(paste0(down_gene[gi,'V1'],' target gene score_sample_volcano.png'), plot = p2, width = 6, height = 6)
  
  
}

for (gi in 1:length(up_gene$V1)) {
  gi=1
  genelist = unlist(strsplit(up_gene[gi,'V2'], split = ","))
  print(up_gene[gi,'V1'])
  scRNA <- AddModuleScore(  object = scRNA,  features = list(genelist),  
                            name = paste0(up_gene[gi,'V1'],' target gene score'))
  plot = scRNA@meta.data %>% dplyr::filter(cell_type %in% c('L5/6 NP2','L2-3 IT-MT'))
  mean_df <- plot %>%
    group_by(group) %>%  # 按分组求平均
    summarise(
      mean_score = mean(.data[[paste0(up_gene[gi,'V1'],' target gene score1')]], na.rm = TRUE),  # 平均值
      sd_score = sd(.data[[paste0(up_gene[gi,'V1'],' target gene score1')]], na.rm = TRUE),      # 可选：标准差
      .groups = "drop"
    )
  paste0("均值更高组：", mean_df$group[which.max(mean_df$mean_score)])
  p = ggplot(plot, aes(x = group,
                       y = .data[[paste0(up_gene[gi,'V1'],' target gene score1')]])) +
    # 小提琴图
    geom_violin(aes(fill = group), alpha = 1, width = 0.8) +
    # 抖动点
    geom_jitter(
      color = 'black',
      size = 0.1,
      alpha = 0.5
    ) +
    # 内层箱线图（窄一点更好看）
    geom_boxplot(width = 0.2,  color = "black", outlier.shape = NA) +
    # P 值放在正中间
    stat_compare_means(
      method = "wilcox.test",
      label.x = 1.5,
      label = "p",
      size = 4.5,
      vjust = 0.5
    ) +
    labs(
      x = "Group",
      y = paste0(up_gene[gi,'V1'],' target gene score'),
      fill = "Group"
    ) +
    theme_classic()
  p
  # ggsave(paste0(up_gene[gi,'V1'],' target gene score_volcano.pdf'), plot = p, width = 4, height = 6)
  # ggsave(paste0(up_gene[gi,'V1'],' target gene score_volcano.png'), plot = p, width = 4, height = 6)
  
  p2 = ggplot(plot, aes(x = sample,
                        y = .data[[paste0(up_gene[gi,'V1'],' target gene score1')]],
                        group = sample)) +
    geom_violin(aes(fill = group)) +
    geom_jitter(
      aes(group = sample),
      color = 'black',
      size = 0.1,
      alpha = 0.8
    )+
    geom_boxplot(width = 0.2,  color = "black", outlier.shape = NA) +
    labs(
      x = "sample",
      y = paste0(up_gene[gi,'V1'],' target gene score'),
      color = "Group"
    ) +
    theme_classic()
  p2
  # ggsave(paste0(up_gene[gi,'V1'],' target gene score_sample_volcano.pdf'), plot = p2, width = 6, height = 6)
  # ggsave(paste0(up_gene[gi,'V1'],' target gene score_sample_volcano.png'), plot = p2, width = 6, height = 6)
  
  
}


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
obj <- readRDS(file = './raw_anno.rds')
obj <- FindVariableFeatures(obj,nfeatures = 2000)
Seurat::DimPlot(obj,reduction = 'umap',group.by = c('group','cell_type','leiden'))

HSMM = readRDS(paste0('/mnt/NexData_file/home/zwl/Project/SC2507233885/require3/monocle/monocle.rds'))
deg <- read.csv(paste0('/mnt/NexData_file/home/zwl/Project/SC2507233885/require3/monocle/monocle.DEG.xls'),sep="\t")
rownames(deg) = deg$Genes
ordergene <- rownames(deg)[1:min(1000,nrow(deg))]
# 重新排序
plot_cell_trajectory(HSMM, color_by = "State")+
  facet_wrap(~cell_type,ncol = 4)
HSMM_myo <- orderCells(HSMM, root_state = 3)

setwd('/mnt/NexData_file/home/zwl/Project/SC2507233885/require3/monocle/')
# 重新画图
setwd(sample[i])
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

p = ggplot(df, aes(x = Pseudotime, y = cell_type,fill =  cell_type)) +
  geom_density_ridges(rel_min_height = 0.005,alpha = 0.6)+
  theme_ridges()+ theme(plot.background = element_rect('white'))
ggsave(filename = "3.monocle.pseudotime_celltype_group.png",p,h=7,w=10,dpi = 150)
ggsave(filename = "3.monocle.pseudotime_celltype_group.pdf",p,h=7,w=10,dpi = 150)


plot_df <- data.frame(
  cell = rownames(pData(HSMM_myo)),
  Comp1 = HSMM_myo@reducedDimS[1, ],
  Comp2 = HSMM_myo@reducedDimS[2, ],
  group = HSMM_myo$group,
  sample = HSMM_myo$sample,
  cell_type= HSMM_myo$cell_type
)

library(ggplot2)
library(ggpointdensity)

p = ggplot(plot_df, aes(Comp1, Comp2)) +
  geom_pointdensity(size = 0.6,shape=19) +
  scale_color_viridis_c() +
  facet_wrap(~group) +
  theme_classic() +
  labs(color = "Density")
p
ggsave(filename = "3.monocle.pseudotime_density_group.png",p,h=5,w=10,dpi = 150)
ggsave(filename = "3.monocle.pseudotime_density_group.pdf",p,h=5,w=10,dpi = 150)

p = ggplot(plot_df, aes(Comp1, Comp2)) +
  geom_pointdensity(size = 0.6,shape=19) +
  scale_color_viridis_c() +
  facet_wrap(~sample) +
  theme_classic() +
  labs(color = "Density")
p
ggsave(filename = "3.monocle.pseudotime_density_sample.png",p,h=8,w=12,dpi = 150)
ggsave(filename = "3.monocle.pseudotime_density_sample.pdf",p,h=8,w=12,dpi = 150)
