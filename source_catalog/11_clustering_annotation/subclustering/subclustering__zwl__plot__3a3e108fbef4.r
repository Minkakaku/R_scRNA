setwd('/mnt/NexData_file/home/zwl/Project/FW2023-179/subcluster/')
rm(list = ls())
gc()
#set cell_type color########################
color_cell = c('#ea7070', '#fdc4b6', '#e59572', '#2694ab', '#96ceb4', '#ffeead')
names(color_cell) = c('CCR8+ Treg', 'CCR8- Treg', 'CD4T_Tem_GZMA-p', 'CD4T_Tem_GZMA-w',
                     'CD4T_Th1-like_CXCL13', 'CD4T_Tm_ANXA1')

metadata = read.csv("./CD4T_metadata.csv")
stat_cell = metadata %>% dplyr::select("sample","cell_type") %>%table()
percent = t(stat_cell)/colSums(stat_cell)*100
rowSums(percent) #检查

data = data.frame(percent) %>% merge(unique(metadata[,c("sample","group")]),by="sample") 

p_value = data%>%
  group_by(cell_type) %>% 
  rstatix::t_test(formula = Freq~group)%>%
  add_significance() %>% add_xy_position()#x = "cell_type"

p = ggboxplot(
  data, x = "group", y = "Freq", 
  color = "group", palette = "jco",
  facet.by = "cell_type",add = "jitter"
) + 
  stat_pvalue_manual(p_value, label = "p",hide.ns = F)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
ggsave(filename = "./3.CD4T_celltype_percent_addp.png",w=6,h=10,p)
ggsave(filename = "./3.CD4T_celltype_percent_addp.pdf",w=6,h=10,p)

##############################CD8Tmonocle3###############################
#############h5ad转变成seurat对象####################
setwd('./Monocle3')
sceasy::convertFormat('/mnt/NexData_file/home/zwl/Project/FW2023-179/CD4T_raw.h5ad', from="anndata", to="seurat",
                      outFile='CD4T_result.rds')

################导入注释信息####################

obj <- readRDS(file = '/mnt/NexData_file/home/zwl/Project/FW2023-179/subcluster/CD4T_result.rds')
metadata = read.csv(file = './CD4T_metadata.csv',row.names = 1)
obj@meta.data = metadata
obj <- FindVariableFeatures(obj,nfeatures = 2000)
Seurat::DimPlot(obj,reduction = 'umap',group.by = c('cell_type','leiden','group'))

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
p1 <- plot_cells(cds,color_cells_by = "cell_type", show_trajectory_graph = FALSE,
                 group_cells_by = 'cell_type',group_label_size  = 5 )
p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE,group_label_size  = 5) + 
  ggtitle("label by partitionID")
ggsave(filename = './CD4T_monocle3/1.monocle3_cell_type_umap.png',p1,w=6,h=5)
ggsave(filename = './CD4T_monocle3/1.monocle3_cell_type_umap.pdf',p1,w=6,h=5)
library(patchwork)
wrap_plots(p1, p2)
#用注释后的cluster
library(do)
obj$cluster1 = 1
monocle3::partitions(cds) = obj$cell_type

## 识别轨迹
cds <- learn_graph(cds,close_loop=FALSE, #分支不画圈
                   use_partition = FALSE)# 所有细胞看称整体画分支
#learn_graph_control=list(minimal_branch_len = 6))#保留一些比较短的分支
p = plot_cells(cds, color_cells_by = "group",label_groups_by_cluster = FALSE,
               label_leaves = FALSE,group_label_size  = 5, 
               label_branch_points = FALSE)
ggsave(filename = './CD4T_monocle3/2.monocle3_group_umap.png',p,w=6,h=5)
ggsave(filename = './CD4T_monocle3/2.monocle3_group_umap.pdf',p,w=6,h=5)
cds <- order_cells(cds)
p2 = plot_cells(cds,
                color_cells_by = "pseudotime",
                label_cell_groups=FALSE,
                label_leaves=FALSE,
                label_branch_points=FALSE,
                graph_label_size=3)
ggsave(filename = './CD4T_monocle3/3.monocle3_trajectory_pseudotime_umap.png',p2,w=6,h=5)
ggsave(filename = './CD4T_monocle3/3.monocle3_trajectory_pseudotime_umap.pdf',p2,w=6,h=5)
save(file = './CD4T_monocle3/CD4T_order.Rdata',cds)

pse = data.frame(monocle3::pseudotime(cds))
obj$pseudotime = pse[colnames(obj),]

obj$cell_type = factor(obj$cell_type,levels = c('CD4T_Th1-like_CXCL13', 'CCR8- Treg', 'CD4T_Tem_GZMA-w',
                                                'CD4T_Tm_ANXA1', 'CD4T_Tem_GZMA-p','CCR8+ Treg'))
p2 = ggplot(obj@meta.data,aes(x = cell_type,y=pseudotime,fill = cell_type))+
  geom_boxplot()+theme_classic()+
  scale_fill_manual(values = color_cell)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(filename = './4.monocle3_celltype_pseudotime_umap.png',p2,w=8,h=5)
ggsave(filename = './4.monocle3_celltype_pseudotime_umap.pdf',p2,w=8,h=5)


library(magick)
library(grid)
image1 = magick::image_read('/mnt/NexData_file/home/zwl/Project/FW2023-179/subcluster/CD4T_cytotrace/cytotrace2_results/cytotrace2_potency_score_by_phenotype.png')
magick::image_write(image = image1,'/mnt/NexData_file/home/zwl/Project/FW2023-179/subcluster/CD4T_cytotrace/cytotrace2_results/cytotrace2_potency_score_by_phenotype.pdf')

setwd('/mnt/NexData_file/home/zwl/Project/FW2023-179/subcluster/CD4T_cytotrace/cytotrace2_results/')
cyto = read.csv(file = './cytotrace2_results.txt',sep = '\t',row.names = 1)
meta = read.csv(file = '/mnt/NexData_file/home/zwl/Project/FW2023-179/subcluster/CD4T_cytotrace/CD4T_celltype.txt',sep = '\t',row.names = 1)
cyto$cell_type = gsub('T_','_',meta[rownames(cyto),'cell_type'])

color_cell = c('#ea7070', '#fdc4b6', '#e59572', '#2694ab', '#96ceb4', '#ffeead')
names(color_cell) = gsub('T_','_',c('CCR8+ Treg', 'CCR8- Treg', 'CD4T_Tem_GZMA-p', 'CD4T_Tem_GZMA-w',
                                    'CD4T_Th1-like_CXCL13', 'CD4T_Tm_ANXA1'))
group = 'cell_type'
potencyBoxplot_byPheno <- ggplot(cyto, aes(x = reorder(!!sym(group), -CytoTRACE2_Score), y = CytoTRACE2_Score)) +
  geom_boxplot(aes(fill = !!sym(group)), width = 0.8, alpha = 0.8, outlier.shape = NA) +
  scale_fill_manual(values = color_cell) +
  geom_jitter(width = 0.1, size = 0.05) +  # 添加散点
  theme_classic() +
  scale_y_continuous(
    breaks = seq(0, 1, by = 0.2),
    limits = c(0, 0.4),
    sec.axis = sec_axis(
      transform = ~.,
      breaks = seq(0, 1, by = 1/12),
      labels = c("", "Differentiated", "", "Unipotent", "", "Oligopotent", 
                 "", "Multipotent", "", "Pluripotent", "", "Totipotent", "")
    )
    # sec.axis = sec_axis(
    #   transform = ~.,
    #   breaks = seq(0, 1, by = 1/12),
    #   labels = c("", "Differentiated", "", "Unipotent", "", "Oligopotent", 
    #              "", "Multipotent", "", "Pluripotent", "", "Totipotent", "")
    # )
  ) +
  labs(x = group, y = "Potency score") +
  ggtitle(sprintf("Developmental potential by %s", group)) +
  theme(
    panel.grid.major.y = element_line(color = "grey90"),  # 添加灰色网格
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),  # x轴标签旋转90度
    legend.position = "none",
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5, margin = margin(b = 20)),
    axis.ticks.y.right = element_line(color = c("black", NA, "black", NA, "black", NA, 
                                                "black", NA, "black", NA, "black", NA, "black")),
    aspect.ratio = 0.8,
    #axis.ticks.length.y.right = unit(0.3, "cm")
  )
ggsave(filename = './cytotrace2_potency_score_by_phenotype2.png',potencyBoxplot_byPheno,w = 4,h=4,dpi = 200)
ggsave(filename = './cytotrace2_potency_score_by_phenotype2.pdf',potencyBoxplot_byPheno,w = 4,h=4,dpi = 200)


setwd('/mnt/NexData_file/home/zwl/Project/FW2023-179/')
res = read.csv('../all_cells_features.xls',sep = '\t',row.names = 1)
res$group = do::Replace(res$group,pattern = c('LC_18A:Aged mouse','LC_6W:Young mouse'))
library(dplyr)
library(ggplot2)

# 统计每个样本中细胞类型的数量
stat_cell = res %>% dplyr::select("sample","cell_type") %>%table()
percent = t(stat_cell)/colSums(stat_cell)*100
rowSums(percent) #检查

data = data.frame(percent) %>% merge(unique(res[,c("sample","group")]),by="sample")
data_summary = data %>%group_by(group,cell_type) %>% 
  mutate(upper =  quantile(Freq, 0.75),
         lower = quantile(Freq, 0.25),
         mean = mean(Freq),
         median = median(Freq))



dodge_scale_val <- 0.4
seg_width <- 0.15
data_summary2 <-
  data_summary %>%
  ungroup() %>%
  mutate(
    xx1 = as.numeric(factor(cell_type, levels = sort(unique(cell_type)))),
    xx2 = as.numeric(factor(group, levels = sort(unique(group)))) %>% {. - median(.)}
  ) %>%
  mutate(
    x_statrt = xx1 + dodge_scale_val * xx2 - seg_width/2,
    x_end = xx1 + dodge_scale_val * xx2 + seg_width/2
  )


set.seed(1234)
p_add_segment <- ggplot(data, aes(x = cell_type, y = Freq)) +
  # 箱线图元素：展示中位值(线)、四分位范围(箱体)
  geom_errorbar(
    data = data_summary,
    aes(ymin = lower, ymax = upper, group = group),
    width = 0.15,
    color = "grey60", 
    position = position_dodge(0.8)
  ) +
  
  geom_segment(
    data = data_summary2,
    aes(x = x_statrt, xend = x_end, y = median),
    color = "grey60",
    inherit.aes = FALSE
  ) +
  
  # 分组散点：空心圆/方块
  geom_point(
    aes(shape = group,colour = group),
    size = 3,
    stroke = 1.2, 
    fill = NA, 
    position = position_jitterdodge(jitter.width = 0.3,dodge.width = 0.75)
  ) +
  scale_shape_manual(values = c(21, 22)) +  # 21=空心圆，22=空心方块
  scale_color_manual(values = c('Aged mouse'="#ff7f0e", 'Young mouse'="#1f77b4"))+
  # stat_pvalue_manual(p_value, label = "p",hide.ns = F)+
  #geom_pwc(aes(group = group), method = "wilcox.test") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x=NULL,y=NULL)

p_add_segment

ggsave(filename = "./1.celltype_percent_addp.png",w=12,h=5,p_add_segment)
ggsave(filename = "./1.celltype_percent_addp.pdf",w=12,h=5,p_add_segment)
ggsave(filename = "./1.celltype_percent.png",w=12,h=5,p_add_segment)
ggsave(filename = "./1.celltype_percent.pdf",w=12,h=5,p_add_segment)


#########################B细胞#################################
setwd('/mnt/NexData_file/home/zwl/Project/FW2023-179/require4/')
sceasy::convertFormat('./adata_BP.h5ad', from="anndata", to="seurat",outFile='adata_BP.rds')
obj <- readRDS(file = './adata_BP.rds')
obj <- FindVariableFeatures(obj,nfeatures = 2000)
Seurat::DimPlot(obj,reduction = 'umap',group.by = c('group','cell_type','leiden'))


color_cell = c('#ea7070', '#fdc4b6', '#e59572', '#2694ab', '#96ceb4', '#ffeead')
names(color_cell) = unique(obj$cell_type)

metadata = obj@meta.data
stat_cell = metadata %>% dplyr::select("sample","cell_type") %>%table()
percent = t(stat_cell)/colSums(stat_cell)*100
rowSums(percent) #检查

data = data.frame(percent) %>% merge(unique(metadata[,c("sample","group")]),by="sample") 

p_value = data%>%
  group_by(cell_type) %>% 
  rstatix::wilcox_test(formula = Freq~group)%>%
  add_significance() %>% add_xy_position()#x = "cell_type"

p = ggboxplot(
  data, x = "group", y = "Freq", 
  color = "group", palette = "jco",
  facet.by = "cell_type",add = "jitter"
) + 
  stat_pvalue_manual(p_value, label = "p",hide.ns = T)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
ggsave(filename = "./3.B_celltype_percent.png",w=6,h=10,p)
ggsave(filename = "./3.B_celltype_percent.pdf",w=6,h=10,p)
