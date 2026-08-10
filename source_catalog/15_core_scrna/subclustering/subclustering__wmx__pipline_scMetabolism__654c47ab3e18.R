# install.packages(c("devtools", "AUCell", "GSEABase", "GSVA", "ggplot2","rsvd"))
# devtools::install_github("YosefLab/VISION")
# devtools::install_github("wu-yc/scMetabolism")
# devtools::install_github("saeyslab/nichenetr")
# library(gtable,lib.loc = "/mnt/NexData_file/home/wmx/envs/scMetabolism/lib/R/library") # 
# library(ggplot2,lib.loc = "/mnt/NexData_file/home/wmx/envs/scRNA//lib/R/library")
# library(scales,lib.loc = "/mnt/NexData_file/home/wmx/envs/scRNA//lib/R/library")
library(matrixStats,  lib.loc = "/mnt/NexData_file/share/scRNA-tools/Proption/envs/scRNA/lib/R/library")
.libPaths("/mnt/NexData_file/home/wmx/envs/scMetabolism/lib/R/library")
library(scMetabolism)
library(Seurat)
library(rsvd)
library(pheatmap)
library(nichenetr)  
library(dplyr)
library(homologene)
library(magrittr)
library(tidyverse)
rm(list = ls())

Mouse.sc.metabolism <- function(obj, metabolism.type=c("KEGG","REACTOME")){
  #将鼠的基因名转化为人的,会有小问题，nichenetr要的R版本很低
  gene_trans = rownames(obj) %>% convert_mouse_to_human_symbols() %>% as.data.frame()
  gene_mouse <- as.data.frame(rownames(obj))
  gene_use <- cbind(gene_trans, gene_mouse)
  gene_use <- na.omit(gene_use)
  colnames(gene_use) <- c('human','mouse')
  mouse_data_trans <- subset(obj,features=gene_use$mouse)
  #函数参考
  #https://www.jianshu.com/p/6495706bac53
  RenameGenesSeurat <- function(obj,newnames,gene.use=NULL,de.assay) {
    # Replace gene names in different slots of a Seurat object. Run this before integration. Run this before integration.
    # It only changes obj@assays$RNA@counts, @data and @scale.data.
    print("Run this before integration. It only changes obj@assays$*@counts, @data and @scale.data, @var.features,@reductions$pca@feature.loadings")
    # lassays <- Assays(obj)
    lassays <-names(obj@assays)
    assay.use <- obj@reductions$pca@assay.used
    DefaultAssay(obj) <- de.assay
    if (is.null(gene.use)) {
      all_genenames <- rownames(obj)
    }else{
      all_genenames <- gene.use
      obj <- subset(obj,features=gene.use)
    }
    
    order_name <- function(v1,v2,ref){
      v2 <- make.names(v2,unique=T)
      df1 <- data.frame(v1,v2)
      rownames(df1) <- df1$v1
      df1 <- df1[ref,]
      return(df1)
    }
    
    df1 <- order_name(v1=all_genenames,v2=newnames,ref=rownames(obj))
    all_genenames <- df1$v1
    newnames <- df1$v2
    
    if ('SCT' %in% lassays) {
      if ('SCTModel.list' %in%  slotNames(obj@assays$SCT)) {
        obj@assays$SCT@SCTModel.list$model1@feature.attributes <- obj@assays$SCT@SCTModel.list$model1@feature.attributes[all_genenames,]
        rownames(obj@assays$SCT@SCTModel.list$model1@feature.attributes) <- newnames
      }
    }
    change_assay <- function(a1=de.assay,obj,newnames=NULL,all_genenames=NULL){
      RNA <- obj@assays[a1][[1]]
      if (nrow(RNA) == length(newnames)) {
        if (length(RNA@counts)) RNA@counts@Dimnames[[1]]            <- newnames
        if (length(RNA@data)) RNA@data@Dimnames[[1]]                <- newnames
        if (length(RNA@var.features)) {
          df1 <- order_name(v1=all_genenames,v2=newnames,ref=RNA@var.features)
          all_genenames1 <- df1$v1
          newnames1 <- df1$v2
          RNA@var.features            <- newnames1
        }
        if (length(RNA@scale.data)){
          df1 <- order_name(v1=all_genenames,v2=newnames,ref=rownames(RNA@scale.data))
          all_genenames1 <- df1$v1
          newnames1 <- df1$v2
          rownames(RNA@scale.data)    <- newnames1
        }
        
      } else {"Unequal gene sets: nrow(RNA) != nrow(newnames)"}
      obj@assays[a1][[1]] <- RNA
      return(obj)
    }
    
    for (a in lassays) {
      DefaultAssay(obj) <- a
      df1 <- order_name(v1=all_genenames,v2=newnames,ref=rownames(obj))
      all_genenames1 <- df1$v1
      newnames1 <- df1$v2
      obj <- change_assay(obj=obj,a1=a,newnames=newnames1,all_genenames=all_genenames1)
    }
    hvg <- VariableFeatures(obj,assay=assay.use)
    if (length(obj@reductions$pca)){
      df1 <- order_name(v1=all_genenames,v2=newnames,ref=hvg)
      df1 <- df1[rownames(obj@reductions$pca@feature.loadings),]
      all_genenames1 <- df1$v1
      newnames1 <- df1$v2
      rownames(obj@reductions$pca@feature.loadings) <- newnames1
    }
    try(obj[[de.assay]]@meta.features <- data.frame(row.names = rownames(obj[[de.assay]])))
    return(obj)
  }
  #转化
  mouse_data_trans <- RenameGenesSeurat(mouse_data_trans, 
                                        newnames = gene_use$human,
                                        gene.use = gene_use$mouse,
                                        de.assay = 'RNA')
  
  mouse_metabolism <- sc.metabolism.Seurat(obj = mouse_data_trans,
                                           method = "AUCell", 
                                           imputation =F, 
                                           ncores = 2, 
                                           metabolism.type = metabolism.type)
  return(mouse_metabolism)
}

## analysis
setwd('//mnt/NexData_file/home/wmx/Project/SC2510235895/shouhou7/scMetabolism/')
rdsfiles = list.files("/mnt/NexData_file/home/wmx/Project/SC2510235895/shouhou7/h5ad_to_rds/","final.rds",full.name =T)
for(rds in rdsfiles){
  scdata = readRDS(rds)
  mouse_results <- Mouse.sc.metabolism(scdata, metabolism.type = 'KEGG')
  save(mouse_results,file = paste0(outpath,'/','Kegg_metabolism.Rdata'))
  # reactome
  mouse_results <- Mouse.sc.metabolism(scdata, metabolism.type = "REACTOME")
  save(mouse_results,file = paste0(outpath,'/','Reactome_metabolism.Rdata'))
}





