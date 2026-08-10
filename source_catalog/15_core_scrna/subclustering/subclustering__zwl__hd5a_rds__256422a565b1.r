setwd('/mnt/NexData_file/home/zwl/Project/FWSC20251318/monocle')
rm(list = ls())
gc()

library(cowplot)
library(Seurat)
library(RColorBrewer)

# ###########################################exchange h5ad to rds############################################
library(sceasy)
library(Seurat)
library(reticulate)
#use_condaenv('myname')

sceasy::convertFormat('./Myeloid_anno.h5ad', from="anndata", to="seurat",outFile='Myeloid_anno.rds')
sceasy::convertFormat('./Tcell_anno.h5ad', from="anndata", to="seurat",outFile='Tcell_anno.rds')

obj <- readRDS(file = './Myeloid_anno.rds')
obj <- FindVariableFeatures(obj,nfeatures = 2000)
Seurat::DimPlot(obj,reduction = 'umap',group.by = c('group','cell_type','leiden'))
cellchat.list = SplitObject(obj,split.by = "group")
saveRDS(file = './Myeloid_Pre.rds',cellchat.list[["pre"]])
saveRDS(file = './Myeloid_Post.rds',cellchat.list[["post"]])



obj2 <- readRDS(file = './Tcell_anno.rds')
obj2 <- FindVariableFeatures(obj2,nfeatures = 2000)
Seurat::DimPlot(obj2,reduction = 'umap',group.by = c('group','cell_type','leiden_batch08'))
FeaturePlot(obj2,features = c('CD4','CD8A'))

CD4T = subset(obj2,cell_type %in% c('CD4 Tn','CD4 Tm','CD8 Tem','Cycling'))
CD8T = subset(obj2,cell_type %in% c('CD8 Tn','CD8 Tc','CD8 Teff-early','CD8 Teff','CD8 Tex'))
  
cellchat.list = SplitObject(CD4T,split.by = "group")
saveRDS(file = './CD4Tcell_Pre.rds',cellchat.list[["pre"]])
saveRDS(file = './CD4Tcell_Post.rds',cellchat.list[["post"]])
cellchat.list = SplitObject(CD8T,split.by = "group")
saveRDS(file = './CD8Tcell_Pre.rds',cellchat.list[["pre"]])
saveRDS(file = './CD8Tcell_Post.rds',cellchat.list[["post"]])

#################################################ROE####################################
setwd('/mnt/NexData_file/home/zwl/Project/FWSC20251318/subcelltype/')
outdir = '/mnt/NexData_file/home/zwl/Project/FWSC20251318/subcelltype/'
metadata = read.csv(file = '/mnt/NexData_file/home/zwl/Project/FWSC20251318/subcelltype/T cells/bbknn_celltype/metadata.csv',row.names = 1) 
cell_name = 'T cells'

source('/mnt/NexData_file/home/zwl/Project/FWSC20240744-02/function.R')
df = metadata %>% select(sample,group,cell_type) %>% rename_with(~c("sample","group",'celltype'))
df = df %>% group_by(sample,celltype) %>% summarise(counts = n(),group=group) %>% unique()
df1 = calROE(df)
df1 = df1 %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance() %>%
  mutate(class=ifelse((`Ro/e`>1.5) & (p.value.adj <0.05),"Enrichment",ifelse((p.value.adj <0.05)&(`Ro/e`<0.5),"Depletion","No Change")))

write.table(df1,file = paste0(cell_name,'_RoeResult.xls'),sep = "\t",quote = F,row.names = F)
# 修改顺序
# ifelse(!is.null(grouporder),(df1$group = factor(df1$group,levels= grouporder[,1])),(df1$group = factor(df1$group)))
# ifelse(!is.null(celltypeorder),(df1$celltype = factor(df1$celltype,levels= celltypeorder[,1])),(df1$celltype = factor(df1$celltype)))
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
savePlot(p5,outdir,paste0(cell_name,"_ROE_Bubble"),unit_width = 3, unit_height = 1)










