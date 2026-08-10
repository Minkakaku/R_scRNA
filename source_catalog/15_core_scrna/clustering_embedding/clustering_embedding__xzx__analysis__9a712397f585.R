rm(list = ls())

gc()


setwd("/mnt/NexData_file/home/xzx/project/AE2509195142/shouhou3/")

library(Seurat)
library(ggplot2)
color = read.delim("./stromal/cell_type_color.xls")
celltype_color <- color$cell_type.1
names(celltype_color) <- color$cell_type
celltype_color
genesets_addmodulescore <- function(gene_sets,gene_name,scRNA,group_name,out_dir){
  genes_list <- list(gene_sets)
  
  scRNA_fin <- scRNA
  scRNA_fin <- AddModuleScore(
    object = scRNA_fin,
    features = genes_list,
    ctrl = 200,
    #nbin = 5,
    name = gene_name
  )
  # 提取需要的数据（细胞ID、分组、模块分数）
  plot_data <- data.frame(
    cell_type = sc@meta.data[[group_name]],  # 分组信息
    score = scRNA_fin@meta.data[[paste0(gene_name,"1")]]  # 模块分数
  )
  
  #plot_data$cell_type <- factor(plot_data$cell_type, levels = cell_type_order)
  p2 <- ggplot(plot_data, aes(x = cell_type, y = score, fill = cell_type)) +
    geom_violin(scale = "width") +
    scale_fill_manual(
      values = celltype_color
    )+
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.background = element_rect(fill = 'white'),
          panel.grid = element_blank()) +
    labs(title = paste0(gene_name," Module Score"), x = '', y = "Score",fill=NULL)
  #+
  #  scale_fill_manual(values = cell_type_colors)
  p2
  ggsave(plot = p2,filename = paste0(out_dir,gene_name,"_",group_name,"_addmodulescore.png"),
         width = 9,height = 6,dpi = 300,bg = "white")
  ggsave(plot = p2,filename = paste0(out_dir,gene_name,"_",group_name,"_addmodulescore.pdf"),
         width = 9,height = 6,dpi = 300,bg = "white")
  
  return(scRNA_fin)
}

sc <- readRDS("./stromal/anno_stomal.rds")

Pericyte <- c('PDGFRB', 'RGS5', 'CSPG4','MCAM','ABCC9','KCNJ8','NOTCH3','DES','ANGPT1')
sc <- genesets_addmodulescore(gene_sets = Pericyte,gene_name = "Pericyte",scRNA = sc,group_name = "leiden",
                              out_dir = "./stromal/")
sc <- genesets_addmodulescore(gene_sets = Pericyte,gene_name = "Pericyte",scRNA = sc,group_name = "cell_type",
                              out_dir = "./stromal/")

SMC <- c('ACTA2','TAGLN','MYH11','CNN1','TPM2','LMOD1','MYL9','CALD1','ACTG2','SMTN')
sc <- genesets_addmodulescore(gene_sets = SMC,gene_name = "SMC",scRNA = sc,group_name = "leiden",
                              out_dir = "./stromal/")
sc <- genesets_addmodulescore(gene_sets = SMC,gene_name = "SMC",scRNA = sc,group_name = "cell_type",
                              out_dir = "./stromal/")

################################################################################
### 客户需求

#### 关键marker在各群的平均表达与阳性比例（明阳性判定阈值）：
markers <- c('PDGFRB','RGS5','CSPG4','MCAM','NOTCH3','ABCC9','DES','KCNJ8','ANGPT1', 
             'PECAM1','VWF','KDR',
             'COL1A1','COL1A2','DCN','LUM','PI16','C7','POSTN','PDGFRA','COL3A1','CXCL14','MFAP4','FBLN1',
             'ACTA2','TAGLN','MYH11','CNN1','TPM2','LMOD1','MYL9','CALD1','ACTG2','SMTN',
             'VIM','MMP2','MMP14','ITGA5','ITGB1','FN1','VCAN','SPP1','CXCL8','CCL2')

# 平均表达
avg <- AverageExpression(
  sc,
  features = markers,
  group.by = "leiden",
  assays = "RNA"
)$RNA
avg <- as.data.frame(avg)
colnames(avg) <- paste0("leiden_",colnames(avg))
write.csv(avg,"./stromal/3.marker_leiden_avg.csv")

# 阳性比例
expr <- GetAssayData(sc, slot = "data")
pct <- sapply(markers, function(g){
  tapply(expr[g, ] > 0, sc$leiden, mean)
})
pct <- as.data.frame(pct)
rownames(pct) <- paste0("leiden_",rownames(pct))
write.csv(pct,"./stromal/3.exp_pos_pct.csv")
#pct <- t(pct)



#### Pericyte score与SMC score在cluster层面的统计:
df <- sc@meta.data

dt <- df %>%
  group_by(leiden) %>%
  summarise(
    Pericyte_mean = mean(Pericyte1),
    Pericyte_sd   = sd(Pericyte1),
    Pericyte_med  = median(Pericyte1),
    Pericyte_IQR  = IQR(Pericyte1),
    SMC_mean      = mean(SMC1),
    SMC_sd        = sd(SMC1),
    SMC_med       = median(SMC1),
    SMC_IQR       = IQR(SMC1)
  )

write.csv(dt,"./stromal/2.summary_ScoreinLeiden.csv",row.names = FALSE)



#### leiden亚群与原始cell_type大类的交叉组成比例表：

tab <- table(sc$leiden, sc$cell_type)

# 行内比例
prop <- prop.table(tab, margin = 1)



################################################################################
rm(list = ls())
#setwd("/mnt/NexData_file/home/xzx/project/SC2510165703/shouhou1/")
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

args$MetaData = '/mnt/NexData_file/home/xzx/project/AE2509195142/shouhou3/stromal/meta.xls'
args$Celltype = 'cell_type'
args$Group = 'group'
args$Output = "/mnt/NexData_file/home/xzx/project/AE2509195142/shouhou3/stromal/celltype/"
args$GroupColor = '/mnt/NexData_file/home/xzx/project/AE2509195142/shouhou3/stromal/group_color.xls'
args$CelltypeColor = '/mnt/NexData_file/home/xzx/project/AE2509195142/shouhou3/stromal/cell_type_color.xls'
# rm(list = ls())

print(args)
source('/mnt/NexData_file/home/xzx/scripts/Proption/Rscript/function.R')
outdir = paste(args$Output,"Proportion Plot",sep = "/")
dir.create(outdir)
setwd(outdir)

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
  geom_bar(stat="identity",position="fill")+
  theme_classic()+
  scale_fill_manual(values=groupusecolor[,1])+
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
  scale_fill_manual(values = groupusecolor[, 1]) + 
  ggtitle("Proportional Circle Plot by Group") +  
  facet_wrap(~Celltype, ncol = 5) +  
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
  scale_fill_manual(values = groupusecolor[,1]) +
  theme_classic() +
  stat_compare_means(method = "wilcox.test", comparisons = combn(group, 2, simplify = FALSE), label = "p.format") +
  facet_wrap(~ celltype, scales = 'free', nrow = 2) +
  theme(axis.text.x=element_text(color='black', angle=45, hjust=1))
p8
savePlot(p8,outdir,'wilcox', unit_width = 4, unit_height = 2)
p9 <- ggplot(df3, aes(x = group, y = Proportion)) +
  geom_boxplot(aes(fill=group)) +
  scale_fill_manual(values = groupusecolor[,1]) +
  theme_classic() +
  stat_compare_means(method = "t.test", comparisons = combn(group, 2, simplify = FALSE), label = "p.format",paired = TRUE) +
  facet_wrap(~ celltype, scales = 'free', nrow = 2) +
  theme(axis.text.x=element_text(color='black', angle=45, hjust=1))
p9
savePlot(p9,outdir,'t-test', unit_width = 4, unit_height = 2)

