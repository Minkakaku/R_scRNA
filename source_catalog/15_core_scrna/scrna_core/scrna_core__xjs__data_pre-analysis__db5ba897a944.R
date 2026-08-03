library(tidyverse)
library(Seurat)

# setting variable   --------------------
dr_outdir <- "2. dimension reduction and cluster"  
dir.create(dr_outdir)
sp_output <- "3. spatial plot" 
dir.create(sp_output)

dir.create("R_OUTPUT")
dir.create("R_OUTPUT/data_for_analysis")


all_data <- readRDS("/mnt/NexData_file/home/xjs/case/VH2601040012-01/output/data/adata_filtered_add_celltype.rds")
all_data

all_data[['RNA']] <- as(all_data[['RNA']], "Assay5")
all_data@assays$RNA$data <- NULL
all_data[['RNA']] <- split(all_data[['RNA']], all_data$project, layers = "counts")


# normalize ----------------------------------------------------------------------
all_data <- NormalizeData(all_data, scale.factor = median(all_data$nCount_RNA))
all_data <- FindVariableFeatures(all_data)

# PCA ----------------------------------------------------------------------------
## PCA --------
#DefaultAssay(all_data) <- "sketch"
#all_data <- FindVariableFeatures(all_data)
VariableFeaturePlot(all_data)

all_data <- ScaleData(all_data, vars.to.regress = 'nCount_RNA')
all_data <- RunPCA(all_data)


### pca inspector ------
pdf(file.path(dr_outdir, 'pca plot.pdf'), height = 6, width = 8)
ElbowPlot(all_data, ndims = 50)
dev.off()

all_data@meta.data %>% head()

all_data$cell_id2 <- NULL

DimPlot(all_data, group.by = "sample", label = TRUE, reduction = "umap") + NoLegend()

# export 
saveRDS(all_data, file = file.path("R_OUTPUT", "all_data.rds"))

DimPlot(all_data, group.by = "sample", label = TRUE, reduction = "umap")
ggsave(file.path(dr_outdir, 'Dimplot_sample.pdf'), width = 8, height = 6)
