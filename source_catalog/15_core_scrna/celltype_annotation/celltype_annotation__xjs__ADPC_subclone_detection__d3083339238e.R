library(Clonalscope)
library(tidyverse)

size=read.table("Clonalscope-main/data-raw/sizes.cellranger-GRCh38-1.0.0.txt", stringsAsFactors = F)
size

# List of cyclegenes retrieved from the "CopyKAT"package (https://github.com/navinlabcode/copykat)
cyclegenes=readRDS("Clonalscope-main/data-raw/cyclegenes.rds")
cyclegenes %>% head()

# bed file indicating gene positions (hg19 and GRCh38 are provided.)
bed=read.table("Clonalscope-main/data-raw/hg38.genes.bed", sep='\t', header = T)
bed %>% head()




## Step0. Load the input files ------
dir_path <- "clonalscope_result/ADPC_subclone_detection"
dir.create(dir_path)

# files for scRNA-seq data
mtx=readMM("/data/outputdata/outputdata-sp/FSWP20202434/ADPC/outs/binned_outputs/square_008um/filtered_feature_bc_matrix/matrix.mtx.gz")
barcodes=read.table("/data/outputdata/outputdata-sp/FSWP20202434/ADPC/outs/binned_outputs/square_008um/filtered_feature_bc_matrix/barcodes.tsv.gz", stringsAsFactors = F, sep='\t', header=F)
features=read.table("/data/outputdata/outputdata-sp/FSWP20202434/ADPC/outs/binned_outputs/square_008um/filtered_feature_bc_matrix/features.tsv.gz", stringsAsFactors = F, sep='\t', header=F)
celltype0=readRDS("ADPC_celltype.rds")
celltype0 %>% head()
barcodes$V1 <- str_c('ADPC_', barcodes$V1)

### without WGS/WES data --------------
# Generate segmentation table for each chromosome arm.
chrarm=read.table("Clonalscope-main/data-raw/cytoarm_table_hg38.txt", stringsAsFactors = F, sep='\t', header=T)
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




## Step2. Subclone detection in spatial transcriptomic data ---------------
# celltype <- celltype0
# celltype$celltype %>% table()
# celltype %>% head()
# celltype[which(!grepl('Luminal cell', celltype[,2])),2] <- 'normal' # 非上皮细胞定为normal
# 
# clustering_barcodes <- celltype[which(grepl("Luminal cell",celltype[,2])),1] # 只保留上皮细胞做后续分析

celltype <-
  celltype0 %>%
  mutate(celltype = case_when(
    celltype == 'Luminal cell' ~ 'tumor',
    TRUE ~ 'normal'
  )) 
clustering_barcodes <- celltype$barcode[celltype$celltype == "tumor"]


### Filter input files -------
Input_filtered=FilterFeatures(mtx=mtx, barcodes=barcodes, features=features, cyclegenes=cyclegenes)

# Remove raw inputs
rm(mtx); rm(barcodes); rm(features)


### Subclone detection ------------
set.seed(2022)

# with WGS/WES
# set.seed(1234)
# subset_cell <- Input_filtered$barcodes %>% .[[1]] %>% sample(size = 10000, replace = F)

# subset_cell <- Input_filtered$barcodes$V1[Input_filtered$barcodes$V1 %in% celltype$barcode]
# subset_idx  <- Input_filtered$barcodes$V1 %in% subset_cell
# subset_idx2 <- clustering_barcodes %in% subset_cell


# inter_cells <- 

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
    # clust_mode          = 'cna_only',
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
    # est_cap             = 2, 
    prep_mode           = 'intersect'
  )

saveRDS(Cov_obj,file.path(dir_path,"ADPC_Cov_obj.rds"))




## Step3. visiualization --------------
clustering= Cov_obj$result_final$clustering
clustering2= Cov_obj$result_final$clustering2
result=Cov_obj$result_final$result
Zest=result$Zest
table(result$Zest)

# 调整resolution，mincell越大cluster数越少
result=AssignCluster(clustering2, mincell = 800)
Zest=result$Zest
table(result$Zest)
saveRDS(Zest, file.path(dir_path, 'Zest_ADPC_subclone_detection_mincell_500_cluster19.rds'))

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
vs <- readRDS('/mnt/NexData_file/home/szt/server203/20250318_FWSP20242434/result_8um/ADPC7_annotation.rds')

vs$orig.ident %>% table()

Zest2 <- 
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
      annotation_full != 'Luminal cell', 
      'Other cell type', 
      Zest)
  ) %>% 
  column_to_rownames('cell_id')

vs <- AddMetaData(vs, metadata = vs_meta[c('Zest', 'subclone')])

Idents(vs) <- 'subclone'
library(gplots)
cols <- c(colors()[c(536, 652, 611, 463, 498, 71, 258, 84, 56, 26, 154, 59, 134, 78, 116, 85, 20, 259)])[1:length(Idents(vs) %>% levels())]
names(cols) <- Idents(vs) %>% levels()
pdf(file.path(dir_path, 'ADPC_subclone_detection_ct.pdf'), width = 8, height = 6)
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
cols <- str_replace_all(cols, c("#0000FF" = 'lightgrey'))
names(cols) <- Idents(vs) %>% levels()

SpatialDimPlot(vs, label.size = 3, pt.size.factor = 10, image.scale = "hires", cols=cols, group.by = 'subclone2', image.alpha = 0.1)
ggsave(file.path(dir_path, 'ADPC_subclone_detection_19clusters.pdf'), width = 8, height = 6)
