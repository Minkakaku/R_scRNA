setwd('/mnt/NexData_file/home/zwl/Project/FWSC20251792/require6')
##########################################细胞比例数量堆叠图############################################
metadata = read.csv(file = '/mnt/NexData_file/home/zwl/Project/FWSC20251792/require6/1.ALL_data_counts.csv')
color_df = read.csv(file = '/mnt/NexData_file/home/zwl/Project/FWSC20251792/require6/1.ALL_data_counts_color.csv',row.names = 1)
metadata = read.csv(file = '/mnt/NexData_file/home/zwl/Project/FWSC20251792/require6/1.Myeloids_data_counts.csv')
color_df = read.csv(file = '/mnt/NexData_file/home/zwl/Project/FWSC20251792/require6/1.Myeloids_data_counts_color.csv',row.names = 1)
metadata = read.csv(file = '/mnt/NexData_file/home/zwl/Project/FWSC20251792/require6/1.Tcells_data_counts.csv')
color_df = read.csv(file = '/mnt/NexData_file/home/zwl/Project/FWSC20251792/require6/1.Tcells_data_counts_color.csv',row.names = 1)

group_color <-setNames(c("#1f77b4","#ff7f0e"),c("Control1","EDA1"))
my_colors <- setNames(color_df$color, rownames(color_df))
metadata = metadata[,-1]
colnames(metadata) = c("sample" ,"celltype","count")
library(dplyr)
# 比例
metadata_prop <- metadata %>%
  group_by(celltype) %>%
  mutate(prop = count / sum(count))

# 总数
metadata_sum <- metadata %>%
  group_by(celltype) %>%
  summarise(total = sum(count)) %>%
  arrange(desc(total))

# 保持顺序一致
metadata_prop$celltype <- factor(metadata_prop$celltype, levels = metadata_sum$celltype)
metadata_sum$celltype <- factor(metadata_sum$celltype, levels = metadata_sum$celltype)

library(ggplot2)

p1 <- ggplot(metadata_prop, aes(x = prop, y = celltype, fill = sample)) +
  geom_bar(stat = "identity") +
  scale_x_continuous(labels = scales::percent) +
  scale_fill_manual(values = group_color) + 
  labs(x = "Proportion", y = NULL) +
  theme_classic() +
  theme(
    legend.position = "top"
  )
p1
p2 <- ggplot(metadata_sum, aes(x = total, y = celltype, fill = celltype)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = total), hjust = -0.1) +
  expand_limits(x = max(metadata_sum$total) * 1.2) +  # ⭐关键
  labs(x = "Count", y = NULL) +
  scale_fill_manual(values = my_colors) + 
  theme_classic() +
  theme(legend.position = "none",
        axis.text.y = element_blank(),   # ❌ 去掉文字
        axis.ticks.y = element_blank()   # ❌ 去掉刻度（可选）
        )
p2
library(patchwork)

p = p1 + p2 + plot_layout(widths = c(1, 1))
p
ggsave(filename = '/mnt/NexData_file/home/zwl/Project/FWSC20251792/require6/1.ALL_data_counts.png',p,w=6,h=4,dpi=150)
ggsave(filename = '/mnt/NexData_file/home/zwl/Project/FWSC20251792/require6/1.ALL_data_counts.pdf',p,w=6,h=4,dpi=150)

ggsave(filename = '/mnt/NexData_file/home/zwl/Project/FWSC20251792/require6/1.Myeloids_data_counts.png',p,w=6,h=4,dpi=150)
ggsave(filename = '/mnt/NexData_file/home/zwl/Project/FWSC20251792/require6/1.Myeloids_data_counts.pdf',p,w=6,h=4,dpi=150)

ggsave(filename = '/mnt/NexData_file/home/zwl/Project/FWSC20251792/require6/1.Tcells_data_counts.png',p,w=6,h=4,dpi=150)
ggsave(filename = '/mnt/NexData_file/home/zwl/Project/FWSC20251792/require6/1.Tcells_data_counts.pdf',p,w=6,h=4,dpi=150)

##########################################分组气泡图############################################
library(dplyr)
library(tidyr)
library(stringr)

scRNA = readRDS('/mnt/NexData_file/home/zwl/Project/FWSC20251792/require6/T_anno.rds')
scRNA$g_celltype = paste0(scRNA$cell_type,' ',scRNA$group)
data_use = Seurat::DotPlot(scRNA, 
                           features = c("Gzma", "Gzmb", "Gzmk", "Prf1", "Ifng", "Fasl"),
                           group.by  = c('g_celltype'))$data
data_use[c('cell_type', 'group')] <- str_split_fixed(data_use$id , ' ', 2)

p <- ggplot(data_use, aes(x = group, y = features.plot)) +
  geom_point(aes(size = pct.exp, color = avg.exp.scaled)) +
  scale_color_gradient2(low = "#2166ac", mid = "white", high = "#b2182b") +
  facet_grid(. ~ cell_type) +
  scale_size(range = c(1, 5)) +
  
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_blank(),
  )
p
ggsave(filename = '/mnt/NexData_file/home/zwl/Project/FWSC20251792/require6/2.Tcells_effect_exp.png',p,w=10,h=4,dpi=150)
ggsave(filename = '/mnt/NexData_file/home/zwl/Project/FWSC20251792/require6/2.Tcells_effect_exp.pdf',p,w=10,h=4,dpi=150)

data_use = Seurat::DotPlot(scRNA, 
                           features = c("Pdcd1", "Ctla4", "Lag3", "Tigit", "Havcr2", "Entpd1"),
                           group.by  = c('g_celltype'))$data
data_use[c('cell_type', 'group')] <- str_split_fixed(data_use$id , ' ', 2)

p <- ggplot(data_use, aes(x = group, y = features.plot)) +
  geom_point(aes(size = pct.exp, color = avg.exp.scaled)) +
  scale_color_gradient2(low = "#2166ac", mid = "white", high = "#b2182b") +
  facet_grid(. ~ cell_type) +
  scale_size(range = c(1, 5)) +
  
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_blank(),
  )
p
ggsave(filename = '/mnt/NexData_file/home/zwl/Project/FWSC20251792/require6/2.Tcells_exhuast_exp.png',p,w=10,h=4,dpi=150)
ggsave(filename = '/mnt/NexData_file/home/zwl/Project/FWSC20251792/require6/2.Tcells_exhuast_exp.pdf',p,w=10,h=4,dpi=150)

####################################################TCR箱线图#################################################
TCR_res= read.table('/mnt/NexData_file/home/zwl/Project/FWSC20251792/require2/TCR/Result.xls',sep = '\t',header = T,row.names = 1)
TCR_res$size_group <- cut(
  TCR_res$airr.cc_aa_alignment_size,
  breaks = c(1, 5, 15, 25, 35, Inf),
  right = FALSE,
  labels = c("[1-5)", "[5-15)", "[15-25)", "[25-35)", ">=35")
)
df_plot <- TCR_res %>%
  group_by(gex.group, gex.cell_type, size_group) %>%
  summarise(count = n(), .groups = "drop")
colnames(df_plot) = c('group','cell_type', 'size_group','count')

color_df = read.csv(file = '/mnt/NexData_file/home/zwl/Project/FWSC20251792/require6/1.Tcells_data_counts_color.csv',row.names = 1)
group_color <-setNames(c("#1f77b4","#ff7f0e"),c("Control1","EDA1"))
my_colors <- setNames(color_df$color, rownames(color_df))

p = ggplot(df_plot, aes(x = size_group, y = count, color = group)) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA,alpha = 0.7) +
  geom_jitter(
    aes(fill = cell_type,group = group),
    color = 'black',
    shape = 21,
    position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8),
    size = 2,
    alpha = 0.8
  ) +
  scale_color_manual(values = group_color) +
  # scale_color_manual(values = my_colors) + 
  labs(
    x = "Size Group",
    y = "Count",
    fill = "Cell Type",
    color = "Group"
  ) +
  theme_classic()
ggsave(filename = '/mnt/NexData_file/home/zwl/Project/FWSC20251792/require6/3.TCR_size_group.png',p,w=6,h=4,dpi=150)
ggsave(filename = '/mnt/NexData_file/home/zwl/Project/FWSC20251792/require6/3.TCR_size_group.pdf',p,w=6,h=4,dpi=150)



##################################################cellchat气泡图修改#######################################################
library(tidyverse)
library(CellChat)
library(RColorBrewer)
load("/mnt/NexData_file/home/zwl/Project/FWSC20251792/require3/cellchat/EDA1_Control1/tmp_object.RData")
pathway = c( "Tnf  - Tnfrsf1b","Thbs1  - Cd47",
             "Sirpb1a - Cd47","Sirpb1c - Cd47",
             "Lgals9  - Cd45","Lgals9  - Havcr2", 
             "Il18  - (Il18r1+Il18rap)",
             # "H2-t23  - (Klrd1+Klrc1)","H2-k1  - Cd8a",
             "Cxcl16  - Cxcr6",
             # "Cd86  - Cd28","Cd86  - Ctla4", 
             # "Cd80  - Cd28","Cd80  - Ctla4",
             "Cd274  - Pdcd1")

p1 = object.list %>% 
  map(subsetCommunication) %>%
  map(function(x){
    x %>% 
      dplyr::filter(
        source %in% 'Mono_Vcan',
        target  %in% c("CD4_Cd5","CD4_Lef1","CD8_Ccl4","CD8_Ccl5","CD8_Ccl6","CD8_Cxcr4","CD8_Gzmf","CD8_Mki67","CD8_Xcl1" ),
        interaction_name_2 %in% pathway,
        # target  %in% c('KC Trem2+', 'HSC'),
        # pathway_name == "TGFb",
        pval < 0.05
      )  %>%
      mutate(cell_pair = str_c(source, target, sep = " -> \n")) 
  }) %>%
  enframe('group') %>%
  unnest(value) %>%
  
  
  ggplot(aes(x = group , y = interaction_name_2)) +
  geom_point(aes(color = prob, size = pval)) +
  facet_wrap(~ cell_pair, nrow = 1) +
  theme_classic() +
  theme(strip.background = element_rect(colour = NA),
        axis.text.x = element_text(size = 10, face = "bold"),
        axis.text.y = element_text(size = 10, face = "bold"),
        axis.title.x = element_text(size = 12, face = "bold"),
        axis.title.y = element_text(size = 12, face = "bold"),
        strip.text = element_text(size = 10, face = "bold")) +
  scale_color_gradientn(colours  = rev(RColorBrewer::brewer.pal(Inf, "RdYlBu"))) +
  scale_size(transform = "reverse", range = c(1, 3)) +
  Seurat::RotatedAxis()
p1
ggsave(filename = '/mnt/NexData_file/home/zwl/Project/FWSC20251792/require6/4.cellchat_mono_Vcan_dotplot.png',p1,w=14,h=8,dpi=150)
ggsave(filename = '/mnt/NexData_file/home/zwl/Project/FWSC20251792/require6/4.cellchat_mono_Vcan_dotplot.pdf',p1,w=14,h=8,dpi=150)

ggsave(filename = '/mnt/NexData_file/home/zwl/Project/FWSC20251792/require9/1.cellchat_mono_Vcan_dotplot2.png',p1,w=12,h=6,dpi=150)
ggsave(filename = '/mnt/NexData_file/home/zwl/Project/FWSC20251792/require9/1.cellchat_mono_Vcan_dotplot2.pdf',p1,w=12,h=6,dpi=150)
################################火山图#######################################
library(ggplot2)
library(ggrepel)
file_path = '/mnt/NexData_file/home/zwl/Project/FWSC20251792/require7/DEG_plot'
path = list.files(
  path = file_path,
  pattern = "^0.*All_expressed_genes\\.xls$",  # 严格匹配 0xxxAll_expressed_genes.xls
  recursive = TRUE,  # 进入所有子文件夹
  full.names = TRUE   # 返回完整路径（方便直接读取）
)
path
for(n in path){
  
  data = read.table(file = n,sep = "\t",head=TRUE)
  
  colnames(data) <- c("Genes", "scores", "avg_log2FC","pvals", "p_val_adj","pct_nz_group" , "pct_nz_reference", "cluster")
  data$log10pvals = -log10(data$pvals)
  
  #data$log10p_val_adj[is.infinite(data$log10p_val_adj)] = 300
  cut_off_PValue =0.05 #设置FDR的阈值
  cut_off_log2FC =1 #设置log2FC的阈值
  data$Sig = ifelse(data$pvals < cut_off_PValue &    #根据阈值筛选差异显著的上下调基因，与差异不显著的基因
                      abs(data$avg_log2FC) >= cut_off_log2FC,  #abs绝对值
                    ifelse(data$avg_log2FC > cut_off_log2FC ,'Up','Down'),'no')
  data = data.frame(data)
  table(data$Sig)
  # data = data %>%
  #   filter(pct_nz_group >= 0.2) %>%
  #   filter(pct_nz_reference <= 0.6) %>%
  #   arrange(desc(abs(log2FC)),FDR)
  xlim_min<-min(data$avg_log2FC)-0.3
  xlim_max<-max(data$avg_log2FC)+0.3
  ylim_max <- min(ceiling(max(data$log10pvals)),300)# * 50
  ###绘图——基础火山图###
  colmap = data.frame(type = c("Down","no","Up"),color=c("#546de5","#d2dae2","#F08080"))
  coloruse = colmap[colmap$type %in% data$Sig,]
  
  p1 <- ggplot(data, aes(x =avg_log2FC, y=log10pvals, colour=Sig)) + #x、y轴取值限制，颜色根据"Sig"
    geom_point(alpha=0.8, size=1) +  #点的透明度、大小
    scale_color_manual(values=coloruse$color) + 
    # xlim(-7,7) +  #调整点的颜色和x轴的取值范围
    geom_vline(xintercept=c(-cut_off_log2FC,cut_off_log2FC),lty=4,col="black",lwd=0.8) + 
    geom_hline(yintercept = -log10(cut_off_PValue), lty=4,col="black",lwd=0.8) + 
    labs(x="log2FC", y="-log10pvals") +  #x、y轴标签
    ggtitle("Volcano Plot") + #标题
    scale_x_continuous(limits = c(xlim_min, xlim_max)) +
    scale_y_continuous(limits = c(0, ylim_max)) +#, breaks = seq(0, ylim_max, 50)
    theme_bw() + # 主题，help(theme)查找其他个性化设置
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(size = 20),
          axis.text.y = element_text(size = 20),
          axis.title.x = element_text(size = 20),
          axis.title.y = element_text(size = 20),
          legend.text = element_text(size = 15),
          legend.position="right", 
          legend.title = element_blank()
    ) 
  p1
  
  Up_top_10 =(     #筛选差异显著上调的前10个Gene
    data %>%
      filter(Sig == 'Up') %>%
      arrange(desc(abs(avg_log2FC)),pvals) %>%
      head(5)
  )
  # Up_top_10 %>% gt() #数据制成表
  Down_top_10 = (       #筛选差异显著下调的前10个Gene
    data %>%
      filter(Sig == 'Down') %>%
      arrange( desc(abs(avg_log2FC)),pvals) %>%
      head(5)
  )
  # Down_top_10 %>% gt() #数据制成表
  
  
  p2 <- p1 +
    geom_label_repel(data = Up_top_10,
                     aes(avg_log2FC, log10pvals, label = Genes),
                     size = 3, fill="#FFCCCC",
                     alpha = 0.65, color = "black")+
    geom_label_repel(data = Down_top_10,
                     aes(avg_log2FC, log10pvals, label = Genes),
                     size = 3, fill="#CCFFFF",
                     alpha = 0.65, color = "black",max.overlaps=10)
  
  p2 #出图
  output_dir <- dirname(normalizePath(n, winslash = "/", mustWork = TRUE))
  output_prefix <- file.path(output_dir, tools::file_path_sans_ext(basename(n)))
  ggsave(paste0(output_prefix, "_volcano.pdf"), plot = p2, width = 8, height = 6)
  ggsave(paste0(output_prefix, "_volcano.png"), plot = p2, width = 8, height = 6)
  
  
  
}

######################################基因小提琴图#################################
setwd('/mnt/NexData_file/home/zwl/Project/FWSC20251792/require9')
metadata = read.csv('/mnt/NexData_file/home/zwl/Project/FWSC20251792/require9/gene_exp.csv',row.names = 1)
metadata_f = metadata# %>%
# dplyr::filter(group %in% c("TST","CON"))#"O_TST",
metadata_f$group = factor(metadata_f$group,levels = c("EDA1","Control1" ))
# metadata_f$celltype = 'ALL cells'
library(rstatix)
library(dplyr)
library(ggpubr)
for (i in colnames(metadata_f)[4:7]) {
  print(i)
  metadata_sub = metadata_f[,c("sample","group","cell_type",i)] 
  colnames(metadata_sub) = c("sample","group","cell_type",'exp')
  cell_counts = table(metadata_sub$cell_type)
  metadata_sub = metadata_sub[metadata_sub$cell_type %in% names(cell_counts)[cell_counts >= 2], ]
  
  
  p1 = ggplot(metadata_sub, aes(x = cell_type, y = exp, fill = cell_type)) +
    geom_violin(position = position_dodge(0.9), width=0.7,alpha = 0.7,scale = 'width') +  # 必加
    geom_jitter(position = position_jitterdodge(jitter.width=0.1, dodge.width=0.9), 
                color='black', size=0.1, alpha=0.5) +
    geom_boxplot(position = position_dodge(0.9), width=0.2, 
                 color="black", outlier.shape=NA) + # 必加
    # stat_pvalue_manual(p_value1, label = "p.adj.signif", hide.ns = TRUE) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size=10))+
    labs(y=i)
  p1
  ggsave(filename = paste0('2.celltype_',i,".png"),w=7,h=4,p1)
  ggsave(filename = paste0('2.celltype_',i,".pdf"),w=7,h=4,p1)
  
  p_value = metadata_sub%>%
    group_by(cell_type) %>% 
    rstatix::wilcox_test(formula = exp~group)%>%
    add_significance() %>% add_xy_position(x = 'cell_type',step.increase = 0.3)
  
  p = ggplot(metadata_sub, aes(x = cell_type, y = exp, fill = group)) +
    geom_violin(position = position_dodge(0.9), width=0.7,alpha = 0.7,scale = 'width')  +  # 必加
    geom_jitter(position = position_jitterdodge(jitter.width=0.1, dodge.width=0.9), 
                color='black', size=0.1, alpha=0.5) +
    geom_boxplot(position = position_dodge(0.9), width=0.2, 
                 color="black", outlier.shape=NA) + # 必加
    stat_pvalue_manual(p_value, label = "p.signif", hide.ns = TRUE) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1, size=10))+
    labs(y=i)
  p
  ggsave(filename = paste0('2.celltype_group_',i,".png"),w=12,h=6,p)
  ggsave(filename = paste0('2.celltype_group_',i,".pdf"),w=12,h=6,p)
  
}
######################################克隆数箱线图#################################
setwd('/mnt/NexData_file/home/zwl/Project/FWSC20251792/require8')
rm(list=ls())
gc()
TCR_res= read.table('/mnt/NexData_file/home/zwl/Project/FWSC20251792/require2/TCR/Result.xls',sep = '\t',header = T,row.names = 1)
library(dplyr)
library(ggplot2)

#--------------------------------------------------
# 1. 克隆数分组
#--------------------------------------------------
TCR_res2 <- TCR_res %>%
  mutate(
    size_group = cut(
      airr.cc_aa_alignment_size,
      breaks = c(1, 5, 15, 25, 35, Inf),
      right = FALSE,
      labels = c("1-5", "5-15", "15-25", "25-35", ">35")
    )
  )

#--------------------------------------------------
# 2. 每个样本每个区间计算比例
#--------------------------------------------------
plot_df <- TCR_res2 %>%
  group_by(gex.sample, gex.group, size_group) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(gex.sample, gex.group) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

#--------------------------------------------------
# 3. 绘制组间差异箱线图
#--------------------------------------------------
group_color <-setNames(c("#1f77b4","#ff7f0e"),c("Control1","EDA1"))
p = ggplot(plot_df, aes(x = size_group, y = prop, fill = gex.group)) +
  geom_bar(
    stat = "identity",
    position = position_dodge(width = 0.75),
    width = 0.65
  ) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(
    x = "expanded clones size",
    y = "Proportion of expanded clones in different tissue groups",
    fill = "tissue groups"
  ) +
  scale_fill_manual(values = group_color) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black")
  )
ggsave(filename = '/mnt/NexData_file/home/zwl/Project/FWSC20251792/require8/4.TCR_size_group.png',p,w=6,h=4,dpi=150)
ggsave(filename = '/mnt/NexData_file/home/zwl/Project/FWSC20251792/require8/4.TCR_size_group.pdf',p,w=6,h=4,dpi=150)

