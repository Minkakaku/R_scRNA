library(Clonalscope)
library(tidyverse)
setwd("./Clonalscope-main/")

size=read.table("data-raw/sizes.cellranger-GRCh38-1.0.0.txt", stringsAsFactors = F)
size

# List of cyclegenes retrieved from the "CopyKAT"package (https://github.com/navinlabcode/copykat)
cyclegenes=readRDS("data-raw/cyclegenes.rds")
cyclegenes %>% head()

# bed file indicating gene positions (hg19 and GRCh38 are provided.)
bed=read.table("data-raw/hg38.genes.bed", sep='\t', header = T)
bed %>% head()




# subclone detection --------------
## Step0. Load the input files ------
dir_path <- "./samples/BC_ductal2/ST/output/"
# files for scRNA-seq data
mtx=readMM("data-raw/BC_ductal2/ST/filtered_feature_bc_matrix/matrix.mtx.gz")
barcodes=read.table("data-raw/BC_ductal2/ST/filtered_feature_bc_matrix/barcodes.tsv.gz", stringsAsFactors = F, sep='\t', header=F)
features=read.table("data-raw/BC_ductal2/ST/filtered_feature_bc_matrix/features.tsv.gz", stringsAsFactors = F, sep='\t', header=F)
celltype0=readRDS("data-raw/BC_ductal2/ST/celltype_all.rds")
celltype0 %>% head()


## Step1. Segment the genome using matched WGS/WES data (Optional) ------------
### without WGS/WES data --------------
# Generate segmentation table for each chromosome arm.
chrarm=read.table("data-raw/cytoarm_table_hg38.txt", stringsAsFactors = F, sep='\t', header=T)
chrarm=chrarm[order(as.numeric(chrarm[,1]),as.numeric(chrarm[,3])),]
chrarm %>% head()
bin_bed=chrarm[,-2]
bin_bed %>% head()

seg_table_filtered=data.frame(
  "chr"   =bin_bed[,1], 
  'start' =as.numeric(bin_bed[,2]),
  'end'   =as.numeric(bin_bed[,3]), 
  'states'=1,
  'length'=as.numeric(bin_bed[,3])-as.numeric(bin_bed[,2]),
  'mean'  =0,
  'var'   =0, 
  'Var1'  =1:nrow(bin_bed),
  'Freq'  =50000,
  'chrr'  =paste0(bin_bed[,1],":", bin_bed[,2]),
  stringsAsFactors = F
)
seg_table_filtered %>% head()


### with WGS/WES data --------------
# files for matched WGS/WES data
WGSt <- readRDS(paste0("data-raw/P5931/WGS/raw_counts.rds"))
WGSn <- readRDS(paste0("data-raw/P5931/WGS/ref_counts.rds"))

# your_list <- 
#   rownames(WGSt) %>% 
#   str_split('-')
# 
# df <- do.call(rbind, lapply(your_list, function(x) {
#   data.frame(
#     chr = x[1],
#     start = as.numeric(x[2]),
#     end = as.numeric(x[3]),
#     stringsAsFactors = FALSE
#   )
# }))
# df %>% head()
# 
# df %>% write_csv('bed.csv')



# First, create an object for segmentation
Obj_filtered <- Createobj_bulk(raw_counts=WGSt, # from matched DNA sequencing (bulk/single)
                               ref_counts=WGSn, # from matched DNA sequencing (bulk/single)
                               samplename= "P5931",
                               genome_assembly="GRCh38", dir_path=dir_path, size=size, assay='WGS')

# Segment the genome on matched WGS data
Obj_filtered <- Segmentation_bulk(Obj_filtered=Obj_filtered,
                                  plot_seg = TRUE, hmm_states = c(0.8, 1.1, 1.2))

Obj_filtered$seg_table %>% head()


## Step2. Subclone detection in spatial transcriptomic data ---------------
celltype <- celltype0
celltype$celltype %>% table()
# celltype[which(!grepl('Epithelial', celltype[,2])),2] <- 'normal' # 非上皮细胞定为normal
# clustering_barcodes <- celltype[which(grepl("Epithelial",celltype[,2])),1] # 只保留上皮细胞做后续分析
clustering_barcodes <- celltype[, 1]

### Filter input files -------
Input_filtered=FilterFeatures(mtx=mtx, barcodes=barcodes, features=features, cyclegenes=cyclegenes)

# Remove raw inputs
rm(mtx); rm(barcodes); rm(features)


### Subclone detection ------------
set.seed(2022)
Cov_obj=RunCovCluster(mtx=Input_filtered$mtx, barcodes=Input_filtered$barcodes, 
                      features=Input_filtered$features, bed=bed, 
                      celltype0=celltype0, var_pt=0.99, var_pt_ctrl=0.99, include='all',
                      alpha_source='all', ctrl_region=NULL, 
                      seg_table_filtered=seg_table_filtered, size=size,
                      dir_path=dir_path, breaks=50, prep_mode = 'intersect', seed=200) 

# # with WGS/WES
# Cov_obj <- RunCovCluster(mtx=Input_filtered$mtx, barcodes=Input_filtered$barcodes,
#                          features=Input_filtered$features, bed=bed,
#                          celltype0=celltype, var_pt=0.99, var_pt_ctrl=0.99, include='all',
#                          alpha_source='all', ctrl_region=NULL,
#                          seg_table_filtered=Obj_filtered$seg_table, size=size,
#                          dir_path=dir_path, breaks=50, prep_mode = 'intersect', est_cap = 2,
#                          clust_mode='cna_only',clustering_barcodes=clustering_barcodes)
# saveRDS(Cov_obj,paste0(dir_path,"/Cov_obj_WGS.rds"))

# save the object
saveRDS(Cov_obj,paste0(dir_path,"/Cov_obj_new.rds"))


## Step3. visiualization --------------
clustering= Cov_obj$result_final$clustering
clustering2= Cov_obj$result_final$clustering2
result=Cov_obj$result_final$result
Zest=result$Zest
table(result$Zest)

# 调整resolution，mincell越大cluster数越少
#result=AssignCluster(clustering2, mincell = 100)
#Zest=result$Zest
#table(result$Zest)

PlotClusters(df = clustering$data, celltype = celltype0, Assign_obj =result, mode = "genome",  fontsize = 7, lab_mode=FALSE)
# # with WGS
# PlotClusters(df = Cov_obj$result_final$df_obj$df, celltype = celltype0, Assign_obj =result, mode = "genome",  fontsize = 7, lab_mode='annot')

# Plot the UMAP and save the coordinates
set.seed(2022)
emb=PlotUMAP(df = clustering$data, celltype = celltype0, Assign_obj =result, mode = "Zest")
# # with WGS
# emb <- PlotUMAP(df = Cov_obj$result_final$df_obj$df, celltype = celltype0, Assign_obj =result, mode = "Zest")


## Step4. Integrate with Spatial information (Optional) ------------------
library(Seurat)
vs <- Load10X_Spatial("data-raw/BC_ductal2/ST/", 
                      filename = "V1_Breast_Cancer_Block_A_Section_2_filtered_feature_bc_matrix.h5")

vs$orig.ident %>% table()

clusters=Zest[match(names(vs$orig.ident), names(Zest))]; clusters[is.na(clusters)]=0
Idents(vs)=factor(clusters, levels=sort(as.numeric(unique(clusters))))
library(gplots)
# cols=c("#ea7070", "#fdc4b6", "#e59572", "#2694ab", "#96ceb4", "#ffeead", 
#        "#d9534f", "#ffad60", "#57D1C9", "#ED5485", "#7F95D1", "#FFE869", "#de4307", 
#        "#f29c2b", "#f6d04d", "#8bc24c", "#4695d6", "#fed95c", "#fa6e57", "#f69e53", 
#        "#248888", "#999999", "#E7475E", "#F0D879", "#29A2C6", "#FFCB18", "#73B66B", 
#        "#FF6D31", "#f1ac9d", "#f06966", "#dee2d1", "#6abe83", "#FCF4D9", "#FFB85F", 
#        "#FF7A5A", "#8ED2C9", "#a5dff9", "#ef5285", "#60c5ba", "#feee7d", "#e8a0b8", 
#        "#ffc300", "#bccf3d", "#02c9c9", "#fbf579", "#005995", "#fa625f", "#600473", 
#        "#f17d80", "#737495")[1:length(Idents(vs) %>% levels())]
cols <- c(colors()[c(1, 609, 536, 62, 652, 611, 463, 498, 71, 258, 84, 56, 26, 154, 59, 134, 78, 116, 85, 20, 259)])[1:length(Idents(vs) %>% levels())]
names(cols) <- Idents(vs) %>% levels()
SpatialDimPlot(vs, label = F, label.size = 3, stroke = 0, pt.size.factor = 2, cols=cols)




# subclone tracing --------------
dir_path <- "./samples/BC_ductal1/ST/output/"; dir.create(dir_path) # set up output directory


## Step0. Load the input files -----------
# files for scRNA-seq data
mtx <- readMM("data-raw/BC_ductal1/ST/filtered_feature_bc_matrix/matrix.mtx.gz")
barcodes <- read.table("data-raw/BC_ductal1/ST/filtered_feature_bc_matrix/barcodes.tsv.gz", stringsAsFactors = F, sep='\t', header=F)
features <- read.table("data-raw/BC_ductal1/ST/filtered_feature_bc_matrix/features.tsv.gz", stringsAsFactors = F, sep='\t', header=F)
celltype0 <- readRDS("data-raw/BC_ductal1/ST/celltype_all_c3.rds")
celltype0 %>% head()


## Step1. Segment the genome using matched WGS/WES data (Optional) -------------
### without WGS/WES data --------------
# Generate segmentation table for each chromosome arm.
chrarm <- read.table("data-raw/cytoarm_table_hg38.txt", stringsAsFactors = F, sep='\t', header=T)
chrarm <- chrarm[order(as.numeric(chrarm[,1]),as.numeric(chrarm[,3])),]
bin_bed <- chrarm[,-2]

seg_table_filtered <- data.frame("chr"=bin_bed[,1], 'start'=as.numeric(bin_bed[,2]),
                                 'end'=as.numeric(bin_bed[,3]), 'states'=1, 'length'=as.numeric(bin_bed[,3])-as.numeric(bin_bed[,2]),
                                 'mean'=0, 'var'=0, 'Var1'=1:nrow(bin_bed),'Freq'=50000,
                                 'chrr'=paste0(bin_bed[,1],":", bin_bed[,2]), stringsAsFactors = F)
seg_table_filtered %>% head()


# ### with WGS/WES data --------------
# # files for matched WGS/WES data
# WGSt <- readRDS(paste0("data-raw/P5931/WGS/raw_counts.rds"))
# WGSn <- readRDS(paste0("data-raw/P5931/WGS/ref_counts.rds"))
# 
# 
# # First, create an object for segmentation
# Obj_filtered <- Createobj_bulk(raw_counts=WGSt, # from matched DNA sequencing (bulk/single)
#                                ref_counts=WGSn, # from matched DNA sequencing (bulk/single)
#                                samplename= "P5931",
#                                genome_assembly="GRCh38", dir_path=dir_path, size=size, assay='WGS')
# 
# # Segment the genome on matched WGS data
# Obj_filtered <- Segmentation_bulk(Obj_filtered=Obj_filtered,
#                                   plot_seg = TRUE, hmm_states = c(0.8, 1.1, 1.2))
# 
# Obj_filtered$seg_table %>% head()


## Step2. Subclone tracing and novel subclone detection in spatial transcriptomic data -----------
celltype <- celltype0
celltype$celltype %>% table()
# celltype[which(!grepl('Epithelial', celltype[,2])),2] <- 'normal' # 非上皮细胞定为normal
# clustering_barcodes <- celltype[which(grepl("Epithelial",celltype[,2])),1] # 只保留上皮细胞做后续分析
clustering_barcodes <- celltype[, 1]

# Filter input files
Input_filtered <- FilterFeatures(mtx=mtx, barcodes=barcodes, features=features, cyclegenes=cyclegenes)

# Remove raw inputs
rm(mtx); rm(barcodes); rm(features)

# Subclone tracing in ST data
# Read in prior object from another associated dataset
Cov_obj0 <- readRDS("samples/BC_ductal2/ST/output/Cov_obj_WGS.rds") # <------------
clustering0 <- Cov_obj0$result_final$clustering
result0 <- Cov_obj0$result_final$result

# Trace subclone in the new dataset
undebug(RunCovCluster)

seg_table_filtered_validation <- Obj_filtered$seg_table
seg_table_filtered_validation$mean <- 0
seg_table_filtered_validation$var <- 0


set.seed(2022)
Cov_obj <- RunCovCluster(mtx=Input_filtered$mtx, barcodes=Input_filtered$barcodes, 
                         features=Input_filtered$features, bed=bed, 
                         celltype0=celltype0, var_pt=0.99, var_pt_ctrl=0.99, include='all',
                         alpha_source='all', ctrl_region=NULL, 
                         seg_table_filtered=seg_table_filtered_validation, size=size,
                         dir_path=dir_path, breaks=50, prep_mode = 'intersect', 
                         clustering0=clustering0, result0=result0, 
                         # est_cap = 2, clust_mode='cna_only',clustering_barcodes=clustering_barcodes, 
                         seed = 200
                         )
# save the object
saveRDS(Cov_obj,paste0(dir_path,"/Cov_obj_WGS.rds"))


## Step3. Visualization --------------
clustering <- Cov_obj$result_final$clustering
clustering2 <- Cov_obj$result_final$clustering2
result <- Cov_obj$result_final$result
Zest <- result$Zest
table(result$Zest)

# result <- AssignCluster(clustering2, mincell = 100)
# Zest <- result$Zest
# table(result$Zest)

PlotClusters(df = clustering$data, celltype = celltype0, Assign_obj =result, mode = "genome",  fontsize = 7, lab_mode=FALSE)
# Plot the UMAP and save the coordinates
set.seed(2022)
emb <- PlotUMAP(df = clustering$data, celltype = celltype0, Assign_obj =result, mode = "Zest")


## Step4. Integrate with Spatial information (Optional) ----------
library(Seurat)
vs <- Load10X_Spatial("data-raw/BC_ductal1/ST/", 
                      filename <- "V1_Breast_Cancer_Block_A_Section_1_filtered_feature_bc_matrix.h5")

library(gplots)
clusters <- Zest[match(names(vs$orig.ident), names(Zest))]; clusters[is.na(clusters)]=0
cols0 <- clusters
U0 <- result$U0
new_indc <- which(! names(table(cols0)) %in% as.character(1:max(which(!is.na(U0[,1])))))
col_use <- c(colors()[c(1, 609, 536, 62, 652, 611, 463, 498, 71, 258, 84, 56, 26, 154, 59, 134, 78, 116, 85, 20, 259)])[1:(length(new_indc)+length(which(!is.na(U0[,1]))))]
names(col_use) <- sort(as.numeric(c(which(!is.na(U0[,1])), names(table(cols0))[new_indc])))

Idents(vs) <- factor(clusters, levels=sort(as.numeric(unique(clusters))))
SpatialDimPlot(vs, label = F, label.size = 3, stroke = 0, pt.size.factor = 2, cols=col_use)
