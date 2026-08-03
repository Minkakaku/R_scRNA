library(Clonalscope)
library(tidyverse)

library(future)
plan(multisession, workers = 10)
options(future.globals.maxSize = 2000 * 1024 ^ 3) # 10Gb

dir_path <- "clonalscope_result/NEPC_subclone_tracing"; dir.create(dir_path) # set up output directory

size=read.table("Clonalscope-main/data-raw/sizes.cellranger-GRCh38-1.0.0.txt", stringsAsFactors = F)
size

# List of cyclegenes retrieved from the "CopyKAT"package (https://github.com/navinlabcode/copykat)
cyclegenes=readRDS("Clonalscope-main/data-raw/cyclegenes.rds")
cyclegenes %>% head()

# bed file indicating gene positions (hg19 and GRCh38 are provided.)
bed=read.table("Clonalscope-main/data-raw/hg38.genes.bed", sep='\t', header = T)
bed %>% head()

## Step0. Load the input files -----------
# files for scRNA-seq data
mtx <- readMM("/data/outputdata/outputdata-sp/FSWP20202434/NEPC/outs/binned_outputs/square_008um/filtered_feature_bc_matrix/matrix.mtx.gz")
barcodes <- read.table("/data/outputdata/outputdata-sp/FSWP20202434/NEPC/outs/binned_outputs/square_008um/filtered_feature_bc_matrix/barcodes.tsv.gz", stringsAsFactors = F, sep='\t', header=F)
features <- read.table("/data/outputdata/outputdata-sp/FSWP20202434/NEPC/outs/binned_outputs/square_008um/filtered_feature_bc_matrix/features.tsv.gz", stringsAsFactors = F, sep='\t', header=F)
celltype0 <- readRDS("NEPC_celltype.rds")
celltype0 %>% head()
barcodes$V1 <- str_c('NEPC_', barcodes$V1)
barcodes %>% head()




## Step1. Segment the genome using matched WGS/WES data (Optional) -------------
### without WGS/WES data --------------
# Generate segmentation table for each chromosome arm.
chrarm <- read.table("Clonalscope-main/data-raw/cytoarm_table_hg38.txt", stringsAsFactors = F, sep='\t', header=T)
chrarm <- chrarm[order(as.numeric(chrarm[,1]),as.numeric(chrarm[,3])),]
bin_bed <- chrarm[,-2]

seg_table_filtered <- data.frame("chr"=bin_bed[,1], 'start'=as.numeric(bin_bed[,2]),
                                 'end'=as.numeric(bin_bed[,3]), 'states'=1, 'length'=as.numeric(bin_bed[,3])-as.numeric(bin_bed[,2]),
                                 'mean'=0, 'var'=0, 'Var1'=1:nrow(bin_bed),'Freq'=50000,
                                 'chrr'=paste0(bin_bed[,1],":", bin_bed[,2]), stringsAsFactors = F)
seg_table_filtered %>% head()




## Step2. Subclone tracing and novel subclone detection in spatial transcriptomic data -----------
celltype <-
  celltype0 %>%
  mutate(celltype = case_when(
    celltype %in% c('Club cell', 'NE cell') ~ 'tumor',
    TRUE ~ 'normal'
  )) 
clustering_barcodes <- celltype$barcode[celltype$celltype == "tumor"]
celltype0$celltype %>% table()
celltype$celltype %>% table()


# Filter input files
Input_filtered <- FilterFeatures(mtx=mtx, barcodes=barcodes, features=features, cyclegenes=cyclegenes)

# Remove raw inputs
rm(mtx); rm(barcodes); rm(features)

# Subclone tracing in ST data
# Read in prior object from another associated dataset
Cov_obj0 <- readRDS("clonalscope_result/ADPC_subclone_detection/ADPC_Cov_obj.rds") # <------------
clustering0 <- Cov_obj0$result_final$clustering

Cov_obj0$result_final$result$Zest %>% table()
Zest <- readRDS('clonalscope_result/ADPC_subclone_detection/Zest_ADPC_subclone_detection.rds')
Cov_obj0$result_final$result$Zest <- Zest
Cov_obj0$result_final$result$Zest %>% table()
result0 <- Cov_obj0$result_final$result
result0$Zest %>% table()

# # Trace subclone in the new dataset
# undebug(RunCovCluster)
# 
# seg_table_filtered_validation <- Obj_filtered$seg_table
# seg_table_filtered_validation$mean <- 0
# seg_table_filtered_validation$var <- 0


set.seed(2022)

subset_cell <- intersect(Input_filtered$barcodes$V1, celltype$barcode)
subset_idx  <- Input_filtered$barcodes$V1 %in% subset_cell
subset_idx2 <- celltype$barcode %in% subset_cell
subset_idx3 <- clustering_barcodes %in% subset_cell


# undebug(RunCovCluster)
Cov_obj <- 
  RunCovCluster(
    mtx       = Input_filtered$mtx[, subset_idx],  # 
    barcodes  = Input_filtered$barcodes[subset_idx, ,drop = FALSE], # 
    features  = Input_filtered$features,  
    bed       = bed,
    seg_table_filtered  = seg_table_filtered, 
    celltype0           = celltype[subset_idx2,], # 
    clustering_barcodes = clustering_barcodes[subset_idx3], # 
    var_pt              = 0.99, 
    var_pt_ctrl         = 0.99,
    include             = 'all',
    alpha_source        = 'all',
    ctrl_region         = NULL,
    size                = size,
    dir_path            = dir_path,
    breaks              = 50,
    clustering0=clustering0, 
    result0=result0, 
    # est_cap             = 2, 
    # clust_mode          = 'cna_only',
    prep_mode           = 'intersect'
  )


# save the object
saveRDS(Cov_obj,paste0(dir_path,"/NEPC_Cov_obj.rds"))




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




## Step4. Integrate with Spatial information (Optional) ------------------
library(Seurat)
vs <- readRDS('/mnt/NexData_file/home/szt/server203/20250318_FWSP20242434/result_8um/NEPC4_annotation.rds')

vs$orig.ident %>% table()

# Zest2 <- 
  Zest %>% 
  as.data.frame() %>% 
  rownames_to_column('cell_id')
colnames(Zest2) <- c('cell_id', 'Zest')

vs_meta <- 
  vs@meta.data %>% 
  rownames_to_column('cell_id') %>% 
  left_join(Zest2, by = 'cell_id') %>% 
  mutate(
    subclone = ifelse(
      annotation_full %in% c('Club cell', 'NE cell'), 
      Zest, 
      'Other cell type'
      )
  ) %>% 
  column_to_rownames('cell_id')

vs <- AddMetaData(vs, metadata = vs_meta[c('Zest', 'subclone')])

Idents(vs) <- 'subclone'
library(gplots)
cols <- c(colors()[c(536, 652, 611, 463, 498, 71, 258, 84, 56, 26, 154, 59, 134, 78, 116, 85, 20, 259)])[1:length(Idents(vs) %>% levels())]
names(cols) <- Idents(vs) %>% levels()
pdf(file.path(dir_path, 'NEPC_subclone_detection_ct.pdf'), width = 8, height = 6)
SpatialDimPlot(vs, label.size = 3, pt.size.factor = 10, image.scale = "hires", cols=cols)
dev.off()    

# vs$subclone2 <- case_when(
#   vs$subclone == '1' ~ 'subclone_1', 
#   vs$subclone == '3' ~ 'subclone_3', 
#   vs$subclone == '' ~ 'subclone_', 
#   vs$subclone == '12' ~ 'subclone_12', 
#   vs$subclone == '16' ~ 'subclone_16', 
#   vs$subclone == '17' ~ 'subclone_17', 
#   vs$subclone == '28' ~ 'subclone_28', 
#   vs$subclone == '40' ~ 'subclone_40', 
#   TRUE ~ 'Others'
# )
vs$subclone2 <- str_replace_all(vs$subclone, '^', 'subclone_')
vs$subclone2 <- str_replace_all(vs$subclone2, 'subclone_Other cell type', 'Other cell type')

vs$subclone2 %>% table()
vs$subclone %>% table()

Idents(vs) <- 'subclone2'
library(gplots)
cols <- c("#FFFF00","#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#8B4513","#F08080","#1E90FF","#7CFC00",
          "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB",
          "#40E0D0","#5F9EA0","#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4",
          "#32CD32","#F0E68C","#FFFFE0","#EE82EE","#FF6347","#6A5ACD","#9932CC","#8B008B","#DEB887")
names(cols) <- Idents(vs) %>% levels()
cols
cols <- str_replace_all(cols, c("#20B2AA" = 'lightgrey'))
names(cols) <- Idents(vs) %>% levels()

SpatialDimPlot(vs, label.size = 3, pt.size.factor = 10, image.scale = "hires", cols=cols, group.by = 'subclone2', image.alpha = 0.1)
ggsave(file.path(dir_path, 'NEPC_subclone_detection_3clusters.pdf'), width = 8, height = 6)
