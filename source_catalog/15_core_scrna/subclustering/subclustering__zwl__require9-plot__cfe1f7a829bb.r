umap = read.csv('/mnt/NexData_file/home/zwl/Project/SC2601160389/售后3/pyscenic/umap_coordinates.csv',row.names = 1)
meta = read.csv('/mnt/NexData_file/home/zwl/Project/SC2601160389/售后3/pyscenic/metadata.csv',row.names = 1)
meta[, 'cell_type'] = gsub('/','_',meta [, 'cell_type'])

loom = open_loom('/mnt/NexData_file/home/zwl/Project/SC2601160389/售后3/pyscenic/aue_result.loom')
# Read information from loom file:
exprMat = get_dgem(loom)
exprMat_log = log2(exprMat + 1) # Better if it is logged/normalized
regulons_incidMat = get_regulons(loom, column.attr.name="Regulons")    
regulons = regulonsToGeneLists(regulons_incidMat)                     #调控子和对应的基因向量
regulonAUC = get_regulons_AUC(loom,column.attr.name = "RegulonsAUC")
regulonAucThresholds = get_regulon_thresholds(loom)                   #用于二值化 regulonAUC，判断细胞是否激活某调控子
embeddings = get_embeddings(loom)
close_loom(loom)


#每个cell type AUC热图 ###########################################################
# 划分cell_type
regulonAUC
regulonAUC  =  regulonAUC [ onlyNonDuplicatedExtended ( rownames ( regulonAUC )),] 

rss  <-  calcRSS ( AUC = getAUC ( regulonAUC ),
                   cellAnnotation = meta [ colnames ( regulonAUC ),  'Cluster' ])
target_sc = rownames(rss)[order(-rss[,'sC6'])][1:13]
target_express = getAUC( regulonAUC )[target_sc,]
write.csv(file = '/mnt/NexData_file/home/zwl/Project/SC2601160389/售后3/pyscenic/sC6_top13_express.csv',t(target_express))

obj = readRDS('/mnt/NexData_file/home/zwl/Project/SC2601160389/售后2/Fibro_raw.rds')
Seurat::DimPlot(obj,reduction = 'umap',group.by = c('cell_type','Cluster','group'))

obj@meta.data = cbind(obj@meta.data,t(target_express)[rownames(obj@meta.data),])
subobj = subset(obj,Cluster %in% 'sC6')

color = c('#1F577B','#E069A6','#FCBC10','#368650')
names(color) = c('Con','Con2','Stent','Stent2')

p = VlnPlot(subobj,features = target_sc,group.by = 'group',pt.size = 0,cols =color,ncol =4) 
ggsave(filename = paste0('/mnt/NexData_file/home/zwl/Project/SC2601160389/售后3/pyscenic/12.sC6_top13_Vlnplot.png'),p,w=10,h=8)
ggsave(filename = paste0('/mnt/NexData_file/home/zwl/Project/SC2601160389/售后3/pyscenic/12.sC6_top13_Vlnplot.pdf'),p,w=10,h=8)

# 直接运行这一句
obj$new_Cluster <- ifelse(obj$Cluster == "sC6", "sC6", "Other")
p = VlnPlot(obj,features = target_sc,group.by = 'new_Cluster',pt.size = 0,ncol =4) 
ggsave(filename = paste0('/mnt/NexData_file/home/zwl/Project/SC2601160389/售后3/pyscenic/12.sC6-other_top13_Vlnplot.png'),p,w=6,h=8)
ggsave(filename = paste0('/mnt/NexData_file/home/zwl/Project/SC2601160389/售后3/pyscenic/12.sC6-other_top13_Vlnplot.pdf'),p,w=6,h=8)


################################################################################
##################cellchat
################################################################################
setwd('/mnt/NexData_file/home/zwl/Project/SC2601160389/售后3/cellchat/')
library(CellChat)
library(ggplot2)
library(ggalluvial)
library(svglite)
library(patchwork)
library(ComplexHeatmap)
library(argparser)
# library(pdftools)
library(stringr)
library(jpeg)
library(cowplot)
library(Seurat)
library(RColorBrewer)
rm(list = ls())
gc()

load("/mnt/NexData_file/home/zwl/Project/SC2601160389/售后3/cellchat/Stent_Con/tmp_object.RData")
cellchat <- mergeCellChat(object.list,add.names = names(object.list))
outdir<-"/mnt/NexData_file/home/zwl/Project/SC2601160389/售后3/cellchat/Stent_Con/"

p1 = object.list %>% 
  map(subsetCommunication) %>%
  map(function(x){
    x %>% 
      dplyr::filter(
        source %in%  c("Macro1","Osteoclast","Macro2","BMSC","Osteoblast"),
        target  %in% c("Macro1","Osteoclast","Macro2","BMSC","Osteoblast"),
        # interaction_name_2 %in% pathway,
        # target  %in% c('KC Trem2+', 'HSC'),
        # pathway_name %in% c('APP','LAIR1'),
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
ggsave(filename = '/mnt/NexData_file/home/zwl/Project/SC2511076348/require9/12.cellchat/pathway/cellchat_mono_dotplot.png',p1,w=30,h=30,dpi=150)
ggsave(filename = '/mnt/NexData_file/home/zwl/Project/SC2511076348/require9/12.cellchat/pathway/cellchat_mono_dotplot.pdf',p1,w=30,h=30,dpi=150)


################################################################################
#######cellphoneDB######
################################################################################











