library(SeuratDisk)
Luminal2 <- readRDS('result_8um/Luminal_reannotation/Luminal2.rds')

Luminal2@assays$RNA <- as(Luminal2@assays$Spatial.008um, "Assay")
DefaultAssay(Luminal2) <- 'RNA'
saveRDS(Luminal2, 'result_8um/Luminal_reannotation/Luminal2_v4.rds')

## sceasy转换格式，conda env tools
# convertFormat(Luminal2, from = "seurat",
#               to = "anndata",
#               outFile = "Luminal2_direct.h5ad"
# )








### reticulate #####
library(reticulate)
library(Seurat)

use_python('/home/szt/miniforge3/envs/rpy/bin/python')
ad <- import('anndata')


## ADPC -------------------------------------------------------------------------
Luminal4 <- readRDS('result_8um/Luminal_reannotation3/Luminal4_final.rds')
ADPC7 <- readRDS('result_8um/ADPC7_annotation.rds')

Luminal4$clusters %>% table()
ADPC7$annotation_full %>% table()

Luminal_meta <- 
  Luminal4@meta.data %>% 
  rownames_to_column('cell_id') %>% 
  dplyr::select(cell_id, clusters)
ADPC_meta <-
  ADPC7@meta.data %>% 
  rownames_to_column('cell_id') %>% 
  dplyr::select(cell_id, annotation_full) %>% 
  left_join(Luminal_meta, by = 'cell_id') %>% 
  mutate(annotation_inferCNV = ifelse(annotation_full == 'Luminal cell', 
                                      clusters, 
                                      annotation_full) %>% 
           str_replace_all(c('1' = 'Luminal_0', 
                             '2' = 'Luminal_1', 
                             '3' = 'Luminal_2', 
                             '4' = 'Luminal_3', 
                             '5' = 'Luminal_4', 
                             '6' = 'Luminal_5'))
  ) %>% 
  column_to_rownames('cell_id')

ADPC7 <- AddMetaData(ADPC7, ADPC_meta['annotation_inferCNV'])

adata = ad$AnnData(
  t(LayerData(JoinLayers(ADPC7), "counts")),
  obs=ADPC7[[]], 
  var=data.frame(row.names = rownames(ADPC7), gene = rownames(ADPC7))
)

adata$write('result_8um/Luminal_reannotation3/ADPC7.manualy.h5ad')


ADPC7@assays$RNA <- as(ADPC7@assays$Spatial.008um, "Assay")
DefaultAssay(ADPC7) <- 'RNA'
saveRDS(ADPC7, 'result_8um/Luminal_reannotation3/ADPC7_direct.rds')

## sceasy转换格式，conda env tools
convertFormat(ADPC7, from = "seurat",
              to = "anndata",
              outFile = "/home/szt/20250318_FWSP20242434/result_8um/Luminal_reannotation3/ADPC7_direct.h5ad"
)
VariableFeatures(ADPC7) %>% write.csv('result_8um/Luminal_reannotation3/ADPC7_HVG.csv')
VariableFeatures(Luminal4) %>% write.csv('result_8um/Luminal_reannotation3/Luminal_HVG.csv')




#NEPC ---------------------------------------------------------------------------
NE <- readRDS('result_8um/NE_subtype/filtered/NE_filterd.rds')
NEPC4 <- readRDS('result_8um/NEPC4_annotation.rds')
NE$Spatial.008um_snn_res.0.25 %>% table()
NE$clusters <- as.character(NE$Spatial.008um_snn_res.0.25)
NE$clusters <- factor(NE$clusters, levels = c(0:4))
NE$clusters %>% table()
NEPC4$annotation_full %>% table()

NE_meta <- 
  NE@meta.data %>% 
  rownames_to_column('cell_id') %>% 
  dplyr::select(cell_id, clusters)
NEPC_meta <- 
  NEPC4@meta.data %>% 
  rownames_to_column('cell_id') %>% 
  dplyr::select(cell_id, annotation_full) %>% 
  left_join(NE_meta, by = 'cell_id') %>% 
  mutate(annotation_inferCNV = ifelse(annotation_full == 'NE cell', 
                                      clusters, 
                                      annotation_full) %>% 
           str_replace_all(c('1' = 'NE_0', 
                             '2' = 'NE_1', 
                             '3' = 'NE_2', 
                             '4' = 'NE_3', 
                             '5' = 'NE_4'))
  ) %>% 
  column_to_rownames('cell_id')
NEPC_meta$annotation_inferCNV %>% table()
NEPC4 <- AddMetaData(NEPC4, NEPC_meta['annotation_inferCNV'])
NEPC5 <- subset(NEPC4, annotation_inferCNV != 'NA')
saveRDS(NEPC5, 'result_8um/NEPC5.rds')

adata = ad$AnnData(
  t(LayerData(NEPC5, layers = 'counts')), 
  obs = NEPC5[[]], 
  var = data.frame(row.names = rownames(NEPC5), gene = rownames(NEPC5))
)
adata$write('result_8um/NE_subtype/filtered/NEPC5.manual.h5ad')

NEPC5@assays$RNA <- as(NEPC5@assays$Spatial.008um, "Assay")
DefaultAssay(NEPC5) <- 'RNA'
saveRDS(NEPC5, 'result_8um/NE_subtype/filtered/NEPC5_direct.rds')

## sceasy转换格式，conda env tools
convertFormat(NEPC5, from = "seurat",
              to = "anndata",
              outFile = "/home/szt/20250318_FWSP20242434/result_8um/NE_subtype/filtered/NEPC5_direct.h5ad"
)

VariableFeatures(NEPC5) %>% write.csv('result_8um/NE_subtype/filtered/NEPC5_HVG.csv')
VariableFeatures(NE) %>% write.csv('result_8um/NE_subtype/filtered/NE_HVG.csv')
