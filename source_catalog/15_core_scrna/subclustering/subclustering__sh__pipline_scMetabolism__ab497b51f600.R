library(Seurat)
# library(harmony)
library(ggplot2)
library(magrittr)
library(dplyr)
# library(clusterProfiler)
library(msigdbr)
library(UCell)
library(AUCell)

setwd("/mnt/NexData_file/home/sh/Projects/SC2510235895/13.Mission/Main/")
S_subf<-readRDS("/mnt/NexData_file/home/sh/Projects/SC2510235895/10.Mission/cellchat_sub_cell_type.rds")
# DimPlot(S_subf,reduction='umap',label = TRUE)
# geneset<-read.table("/mnt/NexData_file/home/sh/Projects/SC2510235895/13.Mission/script/KEGG_metabolism_nc.gmt",sep="\t",header = T)
lines <- readLines("/mnt/NexData_file/home/sh/Projects/SC2510235895/13.Mission/script/KEGG_metabolism_nc.gmt")
geneset <- list()
for (line in lines) {
  parts <- strsplit(line, "\t")[[1]]
  if (length(parts) >= 3) {
    set_name <- parts[1]
    genes <- parts[3:length(parts)]
    geneset[[set_name]] <- genes
  }
}
lines <- readLines("/mnt/NexData_file/home/sh/Projects/FWAE20242338/02.Mission/script/mmu.pathway.gmt")
geneset_all <- list()
for (line in lines) {
  parts <- strsplit(line, "\t")[[1]]
  if (length(parts) >= 3) {
    set_name <- parts[1]
    genes <- parts[3:length(parts)]
    geneset_all[[set_name]] <- genes
  }
}

geneset_matched <- geneset_all[names(geneset)]
lines <- readLines("/mnt/NexData_file/home/sh/Projects/SC2510235895/13.Mission/script/add.gmt")
geneset_add <- list()
for (line in lines) {
  parts <- strsplit(line, "\t")[[1]]
  if (length(parts) >= 3) {
    set_name <- parts[1]
    genes <- parts[3:length(parts)]
    geneset_add[[set_name]] <- genes
  }
}
geneset<- c(geneset_matched, geneset_add)
geneset <- geneset[!is.na(names(geneset))]
names(geneset) <- gsub(" / ", "_", names(geneset))
write.table(geneset, file = "/mnt/NexData_file/home/sh/Projects/SC2510235895/13.Mission/script/all.gmt")
# geneset<-geneset$gene
# gene<-list(geneset)
adata<-AddModuleScore_UCell(S_subf,features=geneset)
colnames(adata@meta.data)

# 假设基因集名称存储在向量中
colnames(adata@meta.data) <- gsub(" / ", "_", colnames(adata@meta.data))
gene_set_names <- colnames(adata@meta.data)[16:100]  # 85个基因集的名称

# 创建输出目录
dir.create("./violin_plots", recursive = TRUE, showWarnings = FALSE)
dir.create("./umap_plots", recursive = TRUE, showWarnings = FALSE)
dir.create("./dot_plots", recursive = TRUE, showWarnings = FALSE)
# 循环处理每个基因集
adata$group <- factor(adata$group,levels = c("asham","aday1","aday3","aday7","aday14"))
for (i in seq_along(gene_set_names)) {
  score_col <- gene_set_names[i]
  if (!score_col %in% colnames(adata@meta.data)) {
    warning(paste("Column not found:", score_col, "- skipping"))
    next
  }
  cat("Processing", i, "/", length(gene_set_names), ":", score_col, "\n")
  # 绘制小提琴图
  p <- VlnPlot(adata, features = score_col, pt.size = 0, group.by = 'group') + NoLegend()+geom_boxplot(width=0.1,col="black",position = position_dodge(width = 0.7),outlier.shape = NA)
  ggsave(paste0("./violin_plots/", score_col, ".VlnPlot.pdf"), 
         width = 10, height = 10, bg = 'white')
  ggsave(paste0("./violin_plots/", score_col, ".VlnPlot.png"), 
         width = 10, height = 10, bg = 'white')
  # 绘制 UMAP 图
  p <- FeaturePlot(adata, features = score_col, reduction = 'umap', 
                   label = TRUE, max.cutoff = 0.99)
  ggsave(paste0("./umap_plots/", score_col, "_umap.png"), 
         width = 7, height = 6, bg = 'white')
  ggsave(paste0("./umap_plots/", score_col, "_umap.pdf"), 
         width = 7, height = 6, bg = 'white')
  p <- DotPlot(adata, features = score_col, group.by = 'group', cols = c("lightgrey","red")) + coord_flip() +
    theme(axis.text.y = element_blank(), plot.title = element_text(hjust = 0.5,size = 10)) +
    ggtitle(score_col) 
  ggsave(paste0("./dot_plots/", score_col, "_dotplot_group.png"), 
         width = 7, height = 4, bg = 'white')
  ggsave(paste0("./dot_plots/", score_col, "_dotplot_group.pdf"), 
         width = 7, height = 4, bg = 'white')
  p <- DotPlot(adata, features = score_col, group.by = 'cell_type', cols = c("lightgrey","red")) + coord_flip() +
    theme(axis.text.y = element_blank(), plot.title = element_text(hjust = 0.5,size = 10)) +
    ggtitle(score_col) 
  ggsave(paste0("./dot_plots/", score_col, "_dotplot_cell_type.png"), 
         width = 30, height = 4, bg = 'white')
  ggsave(paste0("./dot_plots/", score_col, "_dotplot_cell_type.pdf"), 
         width = 30, height = 4, bg = 'white')
}



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
rm(list = ls())
setwd("/mnt/NexData_file/home/sh/Projects/SC2510235895/13.Mission/Main/")
scdata = readRDS("/mnt/NexData_file/home/sh/Projects/SC2510235895/10.Mission/cellchat_sub_cell_type.rds")
# scdata <- subset(scdata, cell_type=='Astrocytes')

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
  # mouse_data_trans[["RNA"]]$data <- NULL
  mouse_metabolism <- sc.metabolism.Seurat(obj = mouse_data_trans,
                                           method = "VISION", 
                                           imputation =F, 
                                           ncores = 8, 
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
df$cell_type <- mouse_results$cell_type
write.table(df, file = 'metabolism_matrix.xls', sep='\t')

avg_df =aggregate(df[,1:ncol(df)-1],list(df$group),mean)
rownames(avg_df) = avg_df$Group.1
avg_df=avg_df[,-1]

avg_df <- as.data.frame(t(avg_df))
write.table(avg_df, file = 'metabolism_heatmap.xls', sep='\t')
avg_df_sec <- avg_df[1:30,]

pheatmap(avg_df_sec, 
         show_colnames = T,
         scale='none', 
         cluster_rows = T,
         color=colorRampPalette(c('#1A5592','white',"#B83D3D"))(100),
         cluster_cols = F,
         filename = 'metabolism_heatmap.png',
         width = 6.5,
         height = 7)
pheatmap(avg_df_sec, 
         show_colnames = T,
         scale='none', 
         cluster_rows = T,
         color=colorRampPalette(c('#1A5592','white',"#B83D3D"))(100),
         cluster_cols = F,
         filename = 'metabolism_heatmap.pdf',
         width = 6.5,
         height = 7)
DotPlot.metabolism(obj = mouse_results, pathway = rownames(mouse_results@assays[["METABOLISM"]][["score"]]), 
                   phenotype = "group", norm = "y")
ggsave('metabolism_dotplot.png', height = 12)
ggsave('metabolism_dotplot.pdf', height = 12)



