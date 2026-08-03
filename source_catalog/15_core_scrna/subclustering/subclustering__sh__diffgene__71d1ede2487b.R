# Sys.setenv(JAVA_HOME='C:/Program Files/Java/jre-9.0.4')
.libPaths()
library(Seurat)
library(dior)
# library(getopt)
library(hdf5r)
library(rhdf5)

rdsdir <- 'AE2508274589-01/01.Mission'
prefix <- 'Main_mapping'

setwd(paste0('D:/项目/',rdsdir))
#setwd(paste0('D:/服务器环境/',rdsdir))
getwd()

a <- dior::read_h5(paste0('D:/项目/',rdsdir,'/',prefix,'.h5'), target.object = 'seurat')
Idents(a) = a$sub_cell_type
a = FindVariableFeatures(a,assay = 'RNA')
saveRDS(a, paste0('D:/项目/',rdsdir,'/',prefix,'.rds'))
# DimPlot(a,reduction='umap',group.by='cluster',raster=FALSE,label = TRUE)
DimPlot(a,reduction='umap',group.by='sub_cell_type',raster=FALSE)

library(dplyr)
library(tidyverse)
library(reshape2)
# adata<-readRDS("result.rds")
adata<-a
adata$group<-factor(adata$group)
adata$Group<-adata$group
levels(adata$Group)

outdir <- "./diffgene-sub_cell_type/"
dir.create(outdir, recursive = TRUE)
Idents(adata) <- "Group"
diffgroup <- c("4cyc-pCP/non-pCR")
diffgroup <- strsplit(diffgroup, "/")
for(n in 1:length(levels(adata$sub_cell_type))){
  # n=11
  cluster_diff_sub<-levels(adata$sub_cell_type)[[n]]
  barcodes<-rownames(adata@meta.data[which(adata$sub_cell_type %in% cluster_diff_sub),])
  SeuratObject_splitlist<-subset(adata,cells = barcodes)
  SeuratObject_splitlist$sub_cell_type <- droplevels(SeuratObject_splitlist$sub_cell_type)
  x <- SeuratObject_splitlist
  tryCatch(
    {
      for(m in 1:length(diffgroup)){
        # 计算差异表达基因
        diff1<-diffgroup[[m]][1]
        diff2<-diffgroup[[m]][2]
        markers <- FindMarkers(x, ident.1 =  diff1, ident.2 = diff2,logfc.threshold = 0.05,min.pct = 0.01)
        markers$gene <- rownames(markers)
        markers <- select(markers,6,everything())
        # 写入文件，确保文件名唯一
        file_path <- paste0(outdir,cluster_diff_sub,"_",diff1,"-",diff2, "_diffgenes.csv")
        write.table(markers, file = file_path, sep = ",", append = TRUE, row.names = FALSE)
        markers <- markers %>% filter(p_val_adj < 0.05, avg_log2FC > 0.25, pct.1 > pct.2, pct.1 > 0.1)
        file_path <- paste0(outdir,cluster_diff_sub,"_",diff1,"-",diff2, "_diffgenes_significant.csv")
        write.table(markers, file = file_path, sep = ",", append = TRUE, row.names = FALSE)
        message(paste("成功处理并保存了对象", cluster_diff_sub))
      }
      return(x)
    },
    error = function(e) {
      message(paste("处理对象", cluster_diff_sub, "时出错:", e$message))
      return(NULL)
    }
  )
}
