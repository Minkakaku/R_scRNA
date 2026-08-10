#group DEG muscat##################################################################
# demo data: GSE181919 
# pathway: /mnt/NexData_file/home/zwl/Project/FWSC20241362/normal_data/
# function: 计算组间差异基因;如果样本细胞亚群数量太少(<50)不适用;分组名不包含-

rm(list = ls())
gc()

library(dplyr)
library(scater)
library(ggplot2)
library(limma)
library('getopt')
library(muscat)
library(limma)
library(purrr)
library(Seurat)
library(SingleCellExperiment)
library(UpSetR)
options(stringsAsFactors = FALSE,warn = -1)
command = matrix(c('outpath','o',1,"character",       #输出地址
                   'seusetinput','r',1,"character", #rds输入地址
                   'celltype','c',1,"character",      #celltype列名
                   'group','g',1,"character",         #group列名
                   'grouplevel','l','1','character', #group比较顺序 test2:test1:control
                   'sample','s',1,"character")        #sample列名
                 ,byrow=TRUE, ncol=4)
args=getopt(command)
if (is.null(args$celltype)) args$celltype = "cell_type"
if (is.null(args$group)) args$group = "group"
if (is.null(args$grouplevel)) print("Need input group compare levels, eg: test2:test1:control")
if (is.null(args$sample)) args$sample = "sample" 

# args$outpath = tmpdir
# args$seusetinput = paste0(tmpdir,"/output/aue_result.loom")
# args$group = "tissue.type"
# args$celltype = "cell.type"
# args$sample = "sample.id"

setwd(args$outpath)
if (!dir.exists('Muscat')) {
  dir.create('Muscat')
}
setwd('./Muscat')

obj <- readRDS(args$seusetinput)
gl = unlist(stringr::str_split(args$grouplevel,':'))
obj@meta.data[[args$group]] = factor(obj@meta.data[[args$group]],levels = gl)
# Filter #####################################################
# obj should filter cell
# remove undetected or lowly expressed genes
obj <- obj[rowSums(GetAssayData(obj,layer = 'counts') > 0) > 10, ]
sce <- as.SingleCellExperiment(obj)

# compute sum-factors & normalize
# sce <- scater::computeLibraryFactors(sce)
# sce <- logNormCounts(sce)

# Data preparation #############################################################
(sce <- prepSCE(sce, 
                kid = args$celltype,    # subpopulation assignments
                gid = args$group,       # group IDs (ctrl/stim)
                sid = args$sample,      # sample IDs (ctrl/stim.1234)
                drop = TRUE))           # drop all other colData columns

nk <- length(kids <- levels(sce$cluster_id))
ns <- length(sids <- levels(sce$sample_id))
names(kids) <- kids; names(sids) <- sids

## Data overview ###############################################################
t(table(sce$cluster_id, sce$sample_id)) #样本细胞数小于50会被忽略
sce <- runUMAP(sce, pca = 30)

# wrapper to prettify reduced dimension plots
.plot_dr <- function(sce, dr, col)
  plotReducedDim(sce, dimred = dr, colour_by = col) +
  guides(fill = guide_legend(override.aes = list(alpha = 1, size = 3))) +
  theme_minimal() + theme(aspect.ratio = 1)

# downsample to max. 100 cells per cluster随机取样
cs_by_k <- split(colnames(sce), sce$cluster_id)
cs100 <- unlist(sapply(cs_by_k, function(u) 
  sample(u, min(length(u), 100))))

# plot t-SNE & UMAP colored by cluster & group ID
pdf(file = '0.Dataoverview_downsample100.pdf')
for (dr in c( "UMAP"))
  for (col in c("cluster_id", "group_id"))
    p = .plot_dr(sce[, cs100], dr, col)
    print(p)
dev.off()
# Begin to DS ##################################################################
# aggregate count to bin -----> matrix: genenum * cellnum(sample*celltype)
pb <- muscat::aggregateData(sce,
                    assay = "counts", fun = "sum",
                    by = c("cluster_id", "sample_id"))
if (dim(colData(pb))[2]==0){
  stop("不存在group_id,样本分组存在问题")
}
# one sheet per subpopulation
assayNames(pb)
t(head(assay(pb)))   #统计每个样本基因的bin
(pb_mds <- pbMDS(pb))
ggsave(filename = '1.data_aggr_celltype_sample.png',print(pb_mds),w=7,h=7)
ggsave(filename = '1.data_aggr_celltype_sample.pdf',print(pb_mds),w=7,h=7)

## run DS analysis##############################################################
ei <- metadata(sce)$experiment_info
g_n = levels(ei$group_id)
if (length(g_n) >2) {
  # construct design & contrast matrix
  mm <- model.matrix(~ 0 + ei$group_id)
  dimnames(mm) <- list(ei$sample_id, levels(ei$group_id))
  for ( gi1 in g_n[-length(g_n)] ) {                               #循环组间差异
    for (gi0 in g_n[(which(g_n == gi1)+1):length(g_n)] ) {
      print(paste0(gi1,'-',gi0,'_DEG:'))
      contrast <- makeContrasts(paste0(gi1,'-',gi0), levels = mm)
      res = pbDS(pb, design = mm, contrast = contrast)             #key function
      group_ALL <- bind_rows(res$table[[1]])
      write.csv(file = paste0('2.',gi1,'-',gi0,'_allgene_nofilter.csv'),group_ALL)
      
      tbl_fil = group_ALL %>% 
        dplyr::filter(p_adj.loc < 0.05, abs(logFC) > 1) %>%
        dplyr::arrange(cluster_id,desc(logFC)) %>% ungroup()
      write.csv(file = paste0('2.',gi1,'-',gi0,'_GroupDEG_logFC1.csv'),tbl_fil)
      
      ##############################visualizing
      
      de_gs_by_k <- split(tbl_fil$gene, tbl_fil$cluster_id)
      p1 = upset(fromList(de_gs_by_k))

      png(paste0('3.', gi1, '-', gi0, '_celltype_sample_DEG_upsetplot.png'),width = 10, height = 7, units = "in", res = 300)
      print(p1) 
      dev.off()
      
      pdf(paste0('3.', gi1, '-', gi0, '_celltype_sample_DEG_upsetplot.pdf'),width = 10, height = 7)
      print(p1) 
      dev.off()
    }
  }
}else if (length(g_n) ==2){
  mm <- model.matrix(~ 0 + ei$group_id)
  dimnames(mm) <- list(ei$sample_id, levels(ei$group_id))
  contrast <- makeContrasts(paste0(g_n[1],'-',g_n[2]), levels = mm)
  # run DS analysis
  res <- pbDS(pb, verbose = FALSE,contrast = contrast,filter = 'none')
  group_ALL <- bind_rows(res$table[[1]])
  write.csv(file = paste0('2.allgene_nofilter.csv'),group_ALL)
  
  tbl_fil = group_ALL %>% 
    dplyr::filter(p_adj.loc < 0.05, abs(logFC) > 1) %>%
    dplyr::arrange(cluster_id,desc(logFC)) %>% ungroup()
  write.csv(file = paste0('2.GroupDEG_logFC1.csv'),tbl_fil)
  
  ##############################visualizing
  
  de_gs_by_k <- split(tbl_fil$gene, tbl_fil$cluster_id)
  p1 = upset(fromList(de_gs_by_k))
  png(paste0('3.', g_n[1], '-', g_n[2], '_celltype_sample_DEG_upsetplot.png'),width = 10, height = 7, units = "in", res = 300)
  print(p1) 
  dev.off()
  
  pdf(paste0('3.', g_n[1], '-', g_n[2], '_celltype_sample_DEG_upsetplot.pdf'),width = 10, height = 7)
  print(p1) 
  dev.off()
}else{
  print('Group levels <2')
}

save(file = 'Aggr_data_info.Rdata',pb,ei)

# # filter list
# tbl_fil <- lapply(tbl, function(u) {
#   u <- dplyr::filter(u, p_adj.loc < 0.05, abs(logFC) > 1)
#   dplyr::arrange(u, p_adj.loc)
# })

