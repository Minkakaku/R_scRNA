setwd("/mnt/NexData_file/home/xzx/project/SC2507294017/shouhou14/PPT_2/")

rm(list = ls())
library(igraph)
library(monocle)
library(Seurat)
library(ggplot2)
#library(argparse)
library(dplyr)
library(tidyverse)
library(ggsci)
library(patchwork)
library(stringr)
library(getopt)

load("/mnt/NexData_file/home/xzx/project/SC2507294017/shouhou4/Epithelial_all/monocle/mono.Rdata")

#### 添加新的celltype分类
cnv = read.csv(file = '/mnt/NexData_file/home/xzx/project/SC2507294017/shouhou15/meta.csv',row.names = 1)
cds@phenoData@data$cell_type = cnv[rownames(cds@phenoData@data),'cell_type']#
cds@phenoData@data$inferCNV = cnv[rownames(cds@phenoData@data),'inferCNV']
table(cds@phenoData@data$cell_type)

plot_cell_trajectory(cds,x=1,y=2,color_by = 'cell_type',cell_size = .5)

cell_type_colors <- c(
  "BGN+ malignant Epi" = "#1f77b4",
  "CCNL2+ malignant Epi" = "#ff7f0e",
  "KRT13+ malignant Epi" = "#2ca02c",
  "SPRR1B+ malignant Epi" = "#d62728",
  "Apoptosis Epithelium" = "#9467bd",
  "MUC5B+ Epithelium" = "#8c564b"
)


p <- plot_cell_trajectory(cds,x=1,y=2,color_by = 'cell_type',cell_size = .5)+
  facet_wrap(~cell_type)+
  scale_color_manual(values = cell_type_colors)
p
ggsave(filename = "./PPT_2/pseudotime_celltype_facet.png",p,h=7,w=10,dpi = 300)
ggsave(filename = "./PPT_2/pseudotime_celltype_facet.pdf",p,h=7,w=10)



p <- plot_cell_trajectory(cds,x=1,y=2,color_by = 'cell_type',cell_size = .5)+
  scale_color_manual(values = cell_type_colors)
p
ggsave(filename = "./PPT_2/pseudotime_celltype.png",p,h=7,w=10,dpi = 300)
ggsave(filename = "./PPT_2/pseudotime_celltype.pdf",p,h=7,w=10)


################################################################################
sc <- readRDS("./adata_epi_anno.rds")
DimPlot(sc,group.by = 'cell_type')


# 1. 提取 Monocle 对象的 State 信息
# pData(cds) 提取的是 metadata 数据框
monocle_states <- pData(cds)$State

# 2. 为提取出的 State 命名（确保能和 Seurat 匹配）
names(monocle_states) <- colnames(cds)

# 3. 将 State 添加到 Seurat 对象的 metadata 中
# Seurat 会自动根据 names(monocle_states) 匹配细胞 ID
sc[["Monocle_State"]] <- monocle_states[colnames(sc)]

# 4. 验证是否成功
DimPlot(sc, group.by = "Monocle_State", label = TRUE)

library(msigdbr)
# 提取人类 Hallmark 集合中的 EMT 基因
emt_genes <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::filter(gs_name == "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION")

emt_genes <- emt_genes$human_gene_symbol

genesets_addmodulescore <- function(gene_sets,gene_name,scRNA,group_name,out_dir){
  genes_list <- list(gene_sets)
  
  scRNA_fin <- scRNA
  scRNA_fin <- AddModuleScore(
    object = scRNA_fin,
    features = genes_list,
    ctrl = 100,
    nbin = 5,
    name = gene_name
  )
  # 提取需要的数据（细胞ID、分组、模块分数）
  plot_data <- data.frame(
    Monocle_State = sc@meta.data[[group_name]],  # 分组信息
    score = scRNA_fin@meta.data[[paste0(gene_name,"1")]]  # 模块分数
  )
  
  #plot_data$cell_type <- factor(plot_data$cell_type, levels = cell_type_order)
  p2 <- ggplot(plot_data, aes(x = Monocle_State, y = score, fill = Monocle_State)) +
    geom_violin(scale = "area") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = paste0(gene_name," Module Score"), x = NULL, y = "Score")
  #+
  #  scale_fill_manual(values = cell_type_colors)
  p2
  ggsave(plot = p2,filename = paste0(out_dir,gene_name,"_",group_name,"_addmodulescore.png"),
         width = 6,height = 6,dpi = 300,bg = "white")
  ggsave(plot = p2,filename = paste0(out_dir,gene_name,"_",group_name,"_addmodulescore.pdf"),
         width = 6,height = 6,dpi = 300,bg = "white")
  
  return(scRNA_fin)
}


##
#sc <- readRDS('./B/celltype/adata_anno.rds')

EMT_Gene <- emt_genes
sc <- genesets_addmodulescore(gene_sets = EMT_Gene,gene_name = "EMT_Gene",
                              scRNA = sc,group_name = "Monocle_State",out_dir = "./")

################################################################################

library(getopt)
library(dplyr)
library(ggplot2)
library(ggsci)
library(rstatix)
library(reshape2)
library(RColorBrewer)
library(ggpubr)
library(pheatmap)
options(stringsAsFactors = FALSE,warn = -1)
command = matrix(c(
  # 参数名称   缩写  需要的值个数 数据类型   # 参数注释
  'MetaData',     'm', 1, "character", # 细胞元数据文件路径
  'Celltype',     'c', 1, "character", # 细胞类型列名称
  'Group',        'g', 1, "character", # 分组信息列名称
  'Output',       'o', 1, "character", # 输出文件路径
  'CelltypeColor','a', 2, "character", # 细胞类型对应的颜色映射（键值对）
  'GroupColor',   'b', 2, "character"  # 组别对应的颜色映射（键值对）
), byrow=TRUE, ncol=4)
args=getopt(command)

args$MetaData = '/mnt/NexData_file/home/xzx/project/SC2507294017/shouhou15/meta.xls'
args$Celltype = 'cell_type'
args$Group = 'group'
args$Output = "./"
args$GroupColor = '/mnt/NexData_file/home/xzx/project/SC2507294017/shouhou15/group_color.xls'
args$CelltypeColor = '/mnt/NexData_file/home/xzx/project/SC2507294017/shouhou15/cell_type_color.xls'
# rm(list = ls())

print(args)
source('/mnt/NexData_file/home/xzx/scripts/Proption/Rscript/function.R')
outdir = paste(args$Output,"Proportion Plot",sep = "/")
dir.create(outdir)
setwd(outdir)

#setwd("/mnt/NexData_file/home/xzx/project/SC2507294017/shouhou14/PPT_5/")

data <- args$MetaData %>% read.csv(sep = "\t",header = T,row.names = 1,check.names = F)
prop = data %>% select(args$Celltype,args$Group) %>% table() %>% prop.table() %>%
  as.data.frame() %>% rename_with(~c("Celltype","Group", "Proportion"))
prop = prop %>% group_by(Group) %>% mutate(Proportion = Proportion / sum(Proportion))
grouporder <- getOrderColor(args$GroupColor)[[1]]
groupcolor <- getOrderColor(args$GroupColor)[[2]]
celltypeorder <- getOrderColor(args$CelltypeColor)[[1]]
celltypecolor <- getOrderColor(args$CelltypeColor)[[2]]
prop1 = prop %>% mutate(Proportion = Proportion*100)
write.table(prop1,paste0("./celltype_proportion.xls"),sep = "\t",quote = F,row.names = F)

ifelse(!is.null(grouporder),(prop$Group = factor(prop$Group,levels= grouporder[,1])),(prop$Group = factor(prop$Group)))
ifelse(!is.null(celltypeorder),(prop$Celltype = factor(prop$Celltype,levels= celltypeorder[,1])),(prop$Celltype = factor(prop$Celltype)))

celltypeusecolor <- getcolor(celltypecolor,length(unique(prop$Celltype))) 
groupusecolor <- getcolor(groupcolor,length(unique(prop$Group))) 
# print(groupusecolor)
# print(celltypeusecolor)
# 绘制细胞类型比例图
p1 = ggplot(prop,aes(x=Group, y=Proportion, fill=Celltype))+
  geom_bar(stat="identity",position="fill")+
  theme_classic()+
  scale_fill_manual(values=celltypeusecolor[,1])+
  ggtitle("Proportional Barplot of celltype")+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),
        axis.text.x=element_text(color='black', angle=45, hjust=1))
p1
p2 <- ggplot(prop,aes(x=Celltype, y=Proportion, fill=Group))+
  geom_bar(stat="identity",position="fill")+ ##（position="fill"强制堆叠为百分比）
  theme_classic()+
  #scale_fill_manual(values=groupusecolor[,1])+
  ggtitle("Proportional Barplot of group")+
  theme(panel.border = element_rect(fill=NA, color="black", size=0.5, linetype="solid"),
        axis.text.x=element_text(color='black', angle=90, hjust=1,vjust = 0.5)) +
  guides(fill=guide_legend(title=NULL))+geom_hline(yintercept  = c(0.25, 0.5, 0.75), linetype = "dashed", color = "blue", size = 1)
p2
# 绘制圆圈比例图
p3 <- ggplot(prop, aes(x = "", y = Proportion, fill = Celltype)) +
  # geom_bar(stat = "identity", width = 1,position = "fill") + 
  geom_bar(stat = "identity", width = 1) + 
  coord_polar(theta = "y", direction = -1, clip = "off") +  
  theme_classic() +
  scale_fill_manual(values = celltypeusecolor[, 1]) +  # Apply custom colors
  ggtitle("Proportional Circle Plot of celltype") +  # Add a title
  facet_wrap(~Group, ncol = 4) +  # Facet by 'Group' with 4 columns
  theme(
    axis.text.x = element_blank(),  
    axis.line = element_blank(),
    axis.title.x = element_blank(), 
    axis.title.y = element_blank(), 
    axis.ticks = element_blank()
  )
p3
p4 <- ggplot(prop, aes(x = "", y = Proportion, fill = Group)) +
  geom_bar(stat = "identity", width = 1,position = "fill") + 
  coord_polar(theta = "y", direction = -1, clip = "off") +  
  theme_classic() +
  #scale_fill_manual(values = groupusecolor[, 1]) + 
  ggtitle("Proportional Circle Plot by Group") +  
  facet_wrap(~Celltype, ncol = 3) +  
  theme(
    axis.text.x = element_blank(),  
    axis.line = element_blank(),    
    axis.title.x = element_blank(), 
    axis.title.y = element_blank(), 
    axis.ticks = element_blank()    
  )
p4
savePlot(p1,outdir,"celltype_barplot_proportion",unit_width = 2, unit_height = 2)
savePlot(p2,outdir,"group_barplot_proportion",unit_width = 1, unit_height = 2)
savePlot(p3,outdir,"celltype_circleplot_proportion",unit_width = 4, unit_height = 4)
savePlot(p4,outdir,"group_circleplot_proportion",unit_width = 4, unit_height = 4)

outdir = paste(args$Output,"ROE Analysis",sep = "/")
dir.create(outdir)
setwd(outdir)
if ('orig.ident' %in% colnames(data) ) {
  data$sample = data$orig.ident
}

df = data %>% select(sample,args$Group,args$Celltype) %>% rename_with(~c("sample","group",'celltype'))
df = df %>% group_by(sample,celltype) %>% summarise(counts = n(),group=group) %>% unique()
df1 = calROE(df)
df1 = df1 %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance() %>%
  mutate(class=ifelse((`Ro/e`>1.5) & (p.value.adj <0.05),"Enrichment",ifelse((p.value.adj <0.05)&(`Ro/e`<0.5),"Depletion","No Change")))

write.table(df1,file = 'RoeResult.xls',sep = "\t",quote = F,row.names = F)
# 修改顺序
ifelse(!is.null(grouporder),(df1$group = factor(df1$group,levels= grouporder[,1])),(df1$group = factor(df1$group)))
ifelse(!is.null(celltypeorder),(df1$celltype = factor(df1$celltype,levels= celltypeorder[,1])),(df1$celltype = factor(df1$celltype)))
# 预设参数
ang = 60 ;bubblelowcolor = '#4a6fe3';bubblehighcolor = '#8e063b';bubblemidcolor =  'grey';bubblesizelength = 5
# unique(df1$class)
df1$class = factor(df1$class,levels = c("Enrichment","No Change","Depletion"))

p5 = ggplot(df1,aes(x=group,y=celltype))+
  # color显示颜色数值大小
  geom_point(aes(color=class,size=`Ro/e`))+
  # 颜色设置针对上面的color参数low为低值对应颜色，high为高值对应颜色
  scale_color_manual(values = rev(c(bubblelowcolor,bubblemidcolor,bubblehighcolor)))+
  # scale_color_manual(values = rev(c(bubblelowcolor,bubblemidcolor)))+
  # 调整圈大小
  scale_size(range = c(2,8),
             limits = c(floor(min(df1$`Ro/e`)),ceiling(max(df1$`Ro/e`))),
             breaks = seq(floor(min(df1$`Ro/e`)),floor(max(df1$`Ro/e`)),length.out = bubblesizelength))+
  # 横纵轴转换
  coord_flip()+
  labs(x = "Group",y = "Celltype",color="",size="Ro/e")+bubble_theme
p5
savePlot(p5,outdir,"ROE_Bubble",unit_width = 7, unit_height = 1)


# install.packages('gtable')
# BiocManager::install('gtable')
# library(gtable)
# outdir = paste(args$Output,"ball stick",sep = "\t")
p6 = ggplot(df1,aes(x=group,y = `Ro/e`))+
  geom_segment(aes(xend = group, yend = 0),
               size = 2, color = 'grey') +
  geom_point(size = 6,shape=20, mapping = aes(color = group)) +
  geom_text(aes(label = round(`Ro/e`,2)),
            size = 5,
            vjust = -.75)+                     #添加Roe值到圆圈中
  # geom_hline(yintercept =1,linetype = 2)+
  geom_hline(yintercept =1.5,linetype = 2,color=bubblehighcolor)+
  geom_hline(yintercept =0.5,linetype = 2,color=bubblelowcolor)+
  scale_color_manual(values = groupusecolor[,1])+
  facet_wrap(~celltype, ncol = 3)+
  labs(y = "Ro/e",x="")+
  ylim(c(0,max(df1$`Ro/e`)+0.5))+
  bubble_theme
p6
savePlot(p6,outdir,"ROE_ball_stick",unit_width = 5, unit_height = 2)
df2 = df1 %>% select(group,celltype,`Ro/e`) %>% dcast(celltype~group) 
rownames(df2) = df2[,1]; if(!is.null(celltypeorder)){df2 = df2[celltypeorder[,1],]} ; ifelse(!is.null(grouporder),(df2 = df2[,grouporder[,1]]),(df2=df2[,-1]))
#p7 = pheatmap::
p7=pheatmap(df2, 
            color = colorRampPalette(brewer.pal(n = 7, name = "Reds"))(10),
            cluster_rows = F, 
            cluster_cols = F, 
            display_numbers = TRUE,
            number_color = "black",
            cellwidth = 30, cellheight = 16, fontsize = 10,                    # fontsize字号大小
            border_color = '#ffffff',
            angle_col='45');dev.off()
p7
# p7$gtable$heights[2] = unit(15,"bigpts");p7$gtable[[2]][p7$gtable[[2]]$name =="legend",1:4] = c(4,5,4,6)
# p7$gtable <- gtable_add_grob(p7$gtable, textGrob("Ro/e",x = unit(0.32,"npc"),y = unit(0.7,"npc")), t = 2, l = 5, r= 5, b=2)
pdf("Roe_heatmap.pdf",width = ncol(df2)*1.5,height = nrow(df2)*0.7);p7;dev.off()
png("Roe_heatmap.png",width = ncol(df2)*130,height = nrow(df2)*40,pointsize = 1);p7;dev.off()
# 绘制组间差异
library(ggpubr)
library(rstatix)
df3 = data %>% select(sample,args$Celltype) %>% rename_with(~c("sample",'celltype'))
df3 = df3 %>% group_by(sample,celltype) %>% table() %>% prop.table() %>% as.data.frame() %>% rename_with(~c("sample","celltype","Proportion"))
df4 = data %>% select(sample,args$Group) %>% unique()
df3 = merge(df3,df4)
df3 = df3 %>% group_by(group) %>% mutate(Proportion = Proportion / sum(Proportion))
ifelse(!is.null(grouporder),(df3$group = factor(df3$group,levels= grouporder[,1])),(df3$group = factor(df3$group)))
ifelse(!is.null(celltypeorder),(df3$celltype = factor(df3$celltype,levels= celltypeorder[,1])),(df3$celltype = factor(df3$celltype)))
outdir = paste(args$Output,"Proportion Plot",sep = "/")
colnames(df3)
group = unique(as.character(df3$group))
y_positions <- df3 %>%
  group_by(celltype) %>%
  summarize(y.position = max(Proportion) * 1.1) 
p8 <- ggplot(df3, aes(x = group, y = Proportion)) +
  geom_boxplot(aes(fill=group)) +
  #scale_fill_manual(values = groupusecolor[,1]) +
  theme_classic() +
  stat_compare_means(method = "wilcox.test", comparisons = combn(group, 2, simplify = FALSE), label = "p.format") +
  facet_wrap(~ celltype, scales = 'free', nrow = 2) +
  theme(axis.text.x=element_text(color='black', angle=45, hjust=1))
p8
savePlot(p8,outdir,'wilcox', unit_width = 4, unit_height = 2)
p9 <- ggplot(df3, aes(x = group, y = Proportion)) +
  geom_boxplot(aes(fill=group)) +
  #scale_fill_manual(values = groupusecolor[,1]) +
  theme_classic() +
  stat_compare_means(method = "t.test", comparisons = combn(group, 2, simplify = FALSE), label = "p.format",paired = TRUE) +
  facet_wrap(~ celltype, scales = 'free', nrow = 2) +
  theme(axis.text.x=element_text(color='black', angle=45, hjust=1))
p9
savePlot(p9,outdir,'t-test', unit_width = 4, unit_height = 2)


################################################################################
# ===================== 1. 加载必需包 =====================
setwd("/mnt/NexData_file/home/xzx/project/SC2507294017/shouhou14/PPT_2/")
library(tidyverse)  # 包含数据处理(dplyr)、绘图(ggplot2)、文件读取(readr)

# ===================== 2. 设置文件路径（仅需修改这里！） =====================
# 核心路径：定位到 recluster 文件夹（cluster0~5的上级目录）
main_path <- "/mnt/NexData_file/home/xzx/project/SC2507294017/shouhou2/recluster_fin/Epithelium/"
# 文件名（所有cluster下的文件名一致）
file_name <- "GO_Biological_Process_2023.Human.enrichr.reports.txt"

# ===================== 3. 批量读取6个cluster的GO结果 =====================
# 创建空列表存储数据
all_data <- list()

# 循环读取 cluster0 ~ cluster5
for (i in 0:5) {
  # 拼接每个文件的完整路径
  cluster_path <- file.path(main_path, paste0("cluster", i), file_name)
  
  # 读取txt文件（enrichr输出的制表符分隔文件）
  df <- read.delim(cluster_path, sep = "\t", stringsAsFactors = FALSE)
  
  # 添加分组列（标记cluster编号）
  df$cluster <- paste0("cluster", i)
  
  # 存入列表
  all_data[[i+1]] <- df
}

# 合并所有数据
combined_data <- bind_rows(all_data)

combined_data <- combined_data %>%
  mutate(cluster = recode(cluster,
                          "cluster0" = "BGN+ Epithelium",
                          "cluster1" = "CCNL2+ Epithelium",
                          "cluster2" = "KRT13+ Epithelium",
                          "cluster3" = "SPRR1B+ Epithelium",
                          "cluster4" = "Apoptosis Epithelium",
                          "cluster5" = "MUC5B+ Epithelium"))
# ===================== 4. 数据筛选：Adjusted P-value <0.05 + 每组TOP5 =====================
# 关键列名说明（enrichr标准输出）：
# Adjusted.P.value = 校正后P值 | Odds.Ratio = 优势比 | Term = GO功能条目
final_data <- combined_data %>%
  # 1. 筛选显著条目：Adjusted P-value < 0.05，每组取TOP5（按P值升序，最显著在前）
  dplyr::filter(Adjusted.P.value < 0.05) %>%
  group_by(cluster) %>%
  arrange(Adjusted.P.value, .by_group = TRUE) %>%
  slice_head(n = 5) %>%
  ungroup() %>%
  
  # 2. 【核心1】彻底删除Term中的(GO:xxxx)编号
  mutate(
    Term = str_remove(Term, " \\(GO:\\d+\\)"), # 正则匹配并删除" (GO:数字)"
    Term = str_trim(Term) # 清除多余空格
  )

# 3. 【核心2】按横轴细胞类型顺序排列纵轴Term
# 先按cluster顺序（normal→tumor1→tumor2→tumor3→tumor4→tumor5）+ P值升序排序
term_order <- final_data %>%
  mutate(cluster = factor(cluster, levels = c("MUC5B+ Epithelium", "BGN+ Epithelium", "CCNL2+ Epithelium", "KRT13+ Epithelium", "SPRR1B+ Epithelium", "Apoptosis Epithelium")
                          )) %>%
  arrange(cluster, Adjusted.P.value) %>%
  pull(Term) %>%
  unique() # 去重，避免重复term

# 将Term设置为因子，固定顺序（反转后normal在最上方，符合阅读习惯）
final_data <- final_data %>%
  mutate(
    Term = factor(Term, levels = rev(term_order)),
    cluster = factor(cluster, levels = c("MUC5B+ Epithelium", "BGN+ Epithelium", "CCNL2+ Epithelium", "KRT13+ Epithelium", "SPRR1B+ Epithelium", "Apoptosis Epithelium")
                     ) # 固定横轴顺序
  )

# ===================== 5. 科研风美化气泡图绘制 =====================
p <- ggplot(final_data, 
            aes(x = cluster, 
                y = Term, 
                size = Odds.Ratio, 
                fill = Adjusted.P.value)) + # 填充色映射P值
  # 气泡样式：带黑色边框+半透明，避免重叠，更清晰
  geom_point(
    shape = 21, # 带边框的圆形
    alpha = 0.9, 
    stroke = 0.8, # 边框粗细
    color = "black" # 边框统一黑色
  ) +
  
  # 颜色渐变：P值越小（越显著）→ 颜色越深（红→蓝，符合科研习惯）
  scale_fill_gradient(
    low = "#D73027",   # 低P值（显著）：深红色
    high = "#3274A1",  # 高P值（接近0.05）：蓝色
    limits = c(0, 0.05), # 固定颜色标尺，统一标准
    name = "Adjusted P-value",
    guide = guide_colorbar(reverse = TRUE) # 反转图例，0.00在最上方，更直观
  ) +
  
  # 气泡大小：映射Odds Ratio，调整范围让差异更明显
  scale_size_continuous(
    range = c(4, 12), # 最小4，最大12（可根据数据调整）
    name = "Odds Ratio",
    breaks = c(20, 40, 60, 80) # 对应图例刻度，和原图一致
  ) +
  
  # 纵轴标签自动换行，避免过长拥挤
  scale_y_discrete(labels = function(x) str_wrap(x, width = 45)) +
  
  # 图表标签与标题
  labs(
    x = "Cell Cluster",
    y = "GO Biological Process",
    title = "GO Biological Process Enrichment Analysis"
    #subtitle = "Top 5 Significant Terms (Adjusted P-value < 0.05)"
  ) +
  
  # 主题美化：极简科研风，符合期刊要求
  theme_minimal() +
  theme(
    # 字体统一为Arial（绝大多数期刊要求）
    text = element_text(size = 10,face = 'bold'),
    # 标题样式：居中、加粗、放大
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5, margin = margin(b = 5)),
    plot.subtitle = element_text(size = 11, hjust = 0.5, color = "gray40", margin = margin(b = 10)),
    # 坐标轴标签样式
    axis.title = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(size = 10, face = "bold"),
    axis.text.y = element_text(size = 9, hjust = 1),
    # 网格线：仅保留纵向网格，横向删除，更简洁
    panel.grid.major.x = element_line(color = "gray90", linewidth = 0.5),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    # 图例样式：右侧垂直排列，调整大小
    legend.position = "right",
    legend.box = "vertical",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9),
    legend.key.size = unit(0.8, "cm"),
    # 边距调整，避免标签被截断
    plot.margin = margin(1, 1.5, 1, 1.5, "cm")
  ) +
  
  # 【可选】替换为Nature子刊配色（更符合高分期刊风格）
  # scale_fill_npg(guide = guide_colorbar(reverse = TRUE), name = "Adjusted P-value") +
  
  NULL

# 显示图表
print(p)

# ===================== 6. 保存高清图片 =====================
ggsave("GOBP_dotplot.pdf", p, width = 16, height = 10, dpi = 300)
ggsave("GOBP_dotplot.png", p, width = 16, height = 10, dpi = 300,bg = 'white')

################################################################################

gene_list <- c(  "PHGDH", 
                 "SLC1A4", "MTHFD1", "DNMT1", "TYMS", 
                 "AHCY", "ATIC", "DNMT3A")

p <- VlnPlot(sc, 
             pt.size=0, # 设置点的大小，但是如果堆叠，则强制不显示点
             stack = T, # 是否堆叠
             flip = T, # 翻转
             group.by = 'cell_type',
             split.by = 'group',
             adjust=3, # 调整粗细
             fill.by = "ident", 
             cols = c('#3274A1','#D73027'), 
             features = gene_list) + 
  labs(x = "") + 
  theme(#legend.position = 'none', 
        axis.text.x = element_text(size = 12, color = "black",angle = 0,hjust = 0.5),
        axis.title.y = element_blank(),
        strip.text.y = element_text(hjust = 0))
p

ggsave(p, file = "../Fig4_ocm_vlnplot.pdf", width = 16, height = 9)
ggsave(p, file = "../Fig4_ocm_vlnplot.png", width = 16, height = 9, dpi = 300,bg = 'white')


################################################################################
sc <- readRDS("../adata_anno_epi.rds")
DimPlot(sc,group.by = 'cell_type')

gene_list <- c("BAX", "CAPN1", "CAPN2", "IRAK1", "NFKBIA","PRKACA", "PRKAR1A", "TNFRSF1A")

p <- VlnPlot(sc, 
             pt.size=0, # 设置点的大小，但是如果堆叠，则强制不显示点
             stack = T, # 是否堆叠
             flip = T, # 翻转
             group.by = 'cell_type',
             #split.by = 'group',
             adjust=3, # 调整粗细
             fill.by = "ident", 
             #cols = c('#3274A1','#D73027'), 
             features = gene_list) + 
  labs(x = "") + 
  theme(#legend.position = 'none', 
    axis.text.x = element_text(size = 12, color = "black",angle = 0,hjust = 0.5),
    axis.title.y = element_blank(),
    strip.text.y = element_text(hjust = 0))
p
