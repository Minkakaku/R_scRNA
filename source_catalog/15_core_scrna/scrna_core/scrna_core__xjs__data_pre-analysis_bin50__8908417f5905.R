library(tidyverse)
library(Seurat)

# setting variable   --------------------
dr_outdir <- "2. dimension reduction and cluster"  
dir.create(dr_outdir)
sp_output <- "3. spatial plot" 
dir.create(sp_output)

dir.create("R_OUTPUT")
dir.create("R_OUTPUT/data_for_analysis")


all_data <- readRDS("/mnt/NexData_file/home/xjs/case/VH2601040012-01/output_bin50/data/adata_filtered.rds")
all_data

all_data[['RNA']] <- as(all_data[['RNA']], "Assay5")
all_data@assays$RNA$data <- NULL
all_data[['RNA']] <- split(all_data[['RNA']], all_data$sample, layers = "counts")


# normalize ----------------------------------------------------------------------
all_data <- NormalizeData(all_data, scale.factor = median(all_data$total_counts))
all_data <- FindVariableFeatures(all_data)

# PCA ----------------------------------------------------------------------------
## PCA --------
#DefaultAssay(all_data) <- "sketch"
#all_data <- FindVariableFeatures(all_data)
VariableFeaturePlot(all_data)

all_data <- ScaleData(all_data, vars.to.regress = 'total_counts')
all_data <- RunPCA(all_data)


### pca inspector ------
pdf(file.path(dr_outdir, 'pca plot_bin50.pdf'), height = 6, width = 8)
ElbowPlot(all_data, ndims = 50)
dev.off()

all_data@meta.data %>% head()

all_data <- RenameCells(all_data, new.names = all_data$cell_id_sample)

DimPlot(all_data, group.by = "sample", label = TRUE, reduction = "umap") + NoLegend()

# export 
saveRDS(all_data, file = file.path("R_OUTPUT", "all_data_bin50.rds"))

# DimPlot(all_data, group.by = "sample", label = TRUE, reduction = "umap")
# ggsave(file.path(dr_outdir, 'Dimplot_sample_bin50.pdf'), width = 8, height = 6)
