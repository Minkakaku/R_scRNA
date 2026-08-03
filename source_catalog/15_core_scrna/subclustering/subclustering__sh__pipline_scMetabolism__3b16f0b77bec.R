# install.packages(c("devtools", "AUCell", "GSEABase", "GSVA", "ggplot2","rsvd"))
# devtools::install_github("YosefLab/VISION")
# devtools::install_github("wu-yc/scMetabolism")
# devtools::install_github("saeyslab/nichenetr")

library(scMetabolism)
library(Seurat)
library(ggplot2)
library(rsvd)
library(pheatmap)
library(nichenetr)
library(dplyr)
library(homologene)
library(magrittr)
library(tidyverse)

rm(list = ls())
setwd('/mnt/NexData_file/home/sh/Projects/FWSC2024-192/18.Mission/result_scMetabolism/')
scdata = readRDS("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/12.Mission/sc_mono_anno.rds")
# scdata <- subset(scdata, cell_type=='Microglia')

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

mouse_results <- Mouse.sc.metabolism(scdata, metabolism.type = 'KEGG')

# 取小胶质细胞
# Neu_metabolism <- subset(mouse_results, cell_type=='Microglia')
df = data.frame(t(mouse_results@assays[["METABOLISM"]][["score"]]),check.names = F)
# mouse_results$sample = mouse_results@active.ident ！！！应该是对细胞类型进行绘图

# names(mouse_results$sample)
# rownames(df) <- gsub(".", "-", rownames(df), fixed = TRUE)
# df = df[names(mouse_results$sample),]

# df$sample <- mouse_results$sample
df$group <- mouse_results$group

avg_df <- df %>%
  group_by(group) %>%
  summarise(across(everything(), mean), .groups = "drop") %>%
  column_to_rownames("group")

df$cell_type <- mouse_results$cell_type
write.table(df, file = 'metabolism_matrix.xls', sep='\t')

avg_df <- as.data.frame(t(avg_df))
# avg_df <- avg_df %>% filter(LAHS > NLAHS & NLAHS > NC) %>% arrange(-LAHS)
write.table(avg_df, file = 'metabolism_heatmap.xls', sep='\t')
avg_df_sec <- avg_df

pheatmap(avg_df_sec, 
         show_colnames = T,
         scale='row', 
         cluster_rows = T,
         color=colorRampPalette(c('#1A5592','white',"#B83D3D"))(100),
         cluster_cols = F,
         filename = 'metabolism_heatmap.png',
         width = 12,
         height = 14)
pheatmap(avg_df_sec, 
         show_colnames = T,
         scale='row', 
         cluster_rows = T,
         color=colorRampPalette(c('#1A5592','white',"#B83D3D"))(100),
         cluster_cols = F,
         filename = 'metabolism_heatmap.pdf',
         width = 12,
         height = 14)

result <- avg_df %>% arrange(-LAHS)
DotPlot.metabolism(obj = mouse_results, pathway = rownames(result), 
                   phenotype = "group", norm = "y")
ggsave('metabolism_dotplot.png', height = 12)
ggsave('metabolism_dotplot.pdf', height = 12)


avg_df <- avg_df %>% filter(LAHS > NLAHS & NLAHS > NC) %>% arrange(-LAHS)
avg_df_sec <- avg_df
pheatmap(avg_df_sec, 
         show_colnames = T,
         scale='row', 
         cluster_rows = F,
         color=colorRampPalette(c('#1A5592','white',"#B83D3D"))(100),
         cluster_cols = F,
         filename = 'metabolism_heatmap_filter.png',
         width = 6.5,
         height = 7)
pheatmap(avg_df_sec, 
         show_colnames = T,
         scale='row', 
         cluster_rows = F,
         color=colorRampPalette(c('#1A5592','white',"#B83D3D"))(100),
         cluster_cols = F,
         filename = 'metabolism_heatmap_filter.pdf',
         width = 6.5,
         height = 7)
DotPlot.metabolism(obj = mouse_results, pathway = rownames(avg_df_sec), 
                   phenotype = "group", norm = "y")
ggsave('metabolism_dotplot_filter.png', height = 7)
ggsave('metabolism_dotplot_filter.pdf', height = 7)



DotPlot.metabolism = function (obj, pathway, phenotype, norm = "y"){
    input.norm = norm
    input.pathway <- pathway
    input.parameter <- phenotype
    metadata <- obj@meta.data
    metabolism.matrix <- obj@assays$METABOLISM$score
    cat("\nPlease Cite: \nYingcheng Wu, Qiang Gao, et al. Cancer Discovery. 2021. \nhttps://pubmed.ncbi.nlm.nih.gov/34417225/   \n\n")
    metadata[, input.parameter] <- as.character(metadata[, input.parameter])
    metabolism.matrix_sub <- t(metabolism.matrix[input.pathway,
        ])
    gg_table <- c()
    for (i in 1:length(input.pathway)) {
        gg_table <- rbind(gg_table, cbind(metadata[, input.parameter],
            input.pathway[i], metabolism.matrix_sub[, i]))
    }
    gg_table <- data.frame(gg_table)
    gg_table_median <- c()
    input.group.x <- unique(as.character(gg_table[, 1]))
    input.group.y <- unique(as.character(gg_table[, 2]))
    for (x in 1:length(input.group.x)) {
        for (y in 1:length(input.group.y)) {
            gg_table_sub <- subset(gg_table, gg_table[, 1] ==
                input.group.x[x] & gg_table[, 2] == input.group.y[y])
            gg_table_median <- rbind(gg_table_median, cbind(input.group.x[x],
                input.group.y[y], median(as.numeric(as.character(gg_table_sub[,
                  3])))))
        }
    }
    gg_table_median <- data.frame(gg_table_median)
    gg_table_median[, 3] <- as.numeric(as.character(gg_table_median[,
        3]))
    gg_table_median_norm <- c()
    input.group.x <- unique(as.character(gg_table[, 1]))
    input.group.y <- unique(as.character(gg_table[, 2]))
    range01 <- function(x) {
        (x - min(x))/(max(x) - min(x))
    }
    if (input.norm == "y")
        for (y in 1:length(input.group.y)) {
            gg_table_median_sub <- subset(gg_table_median, gg_table_median[,
                2] == input.group.y[y])
            norm_value <- range01(as.numeric(as.character(gg_table_median_sub[,
                3])))
            gg_table_median_sub[, 3] <- norm_value
            gg_table_median_norm <- rbind(gg_table_median_norm,
                gg_table_median_sub)
        }
    if (input.norm == "x")
        for (x in 1:length(input.group.x)) {
            gg_table_median_sub <- subset(gg_table_median, gg_table_median[,
                1] == input.group.x[x])
            norm_value <- range01(as.numeric(as.character(gg_table_median_sub[,
                3])))
            gg_table_median_sub[, 3] <- norm_value
            gg_table_median_norm <- rbind(gg_table_median_norm,
                gg_table_median_sub)
        }
    if (input.norm == "na")
        gg_table_median_norm <- gg_table_median
    gg_table_median_norm <- data.frame(gg_table_median_norm)
    gg_table_median_norm[, 3] <- as.numeric(as.character(gg_table_median_norm[,
        3]))
    library(wesanderson)
    pal <- wes_palette("Zissou1", 100, type = "continuous")
		# gg_table_median_norm$X2<-factor(gg_table_median_norm$X2,levels=pathway)
		gg_table_median_norm$X1<-factor(gg_table_median_norm$X1,levels=c("LAHS","NLAHS","NC"))
    ggplot(data = gg_table_median_norm, aes(x = gg_table_median_norm[,
        1], y = gg_table_median_norm[, 2], color = gg_table_median_norm[,
        3])) + geom_point(data = gg_table_median_norm, aes(size = gg_table_median_norm[,
        3])) + ylab("Metabolic Pathway") + xlab(input.parameter) +
        theme_bw() + theme(axis.text.x = element_text(angle = 45,
        hjust = 1), panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
        scale_color_gradientn(colours = pal) + labs(color = "Value",
        size = "Value") + NULL
}


