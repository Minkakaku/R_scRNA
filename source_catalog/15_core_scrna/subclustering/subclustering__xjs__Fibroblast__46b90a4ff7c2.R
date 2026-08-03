library(Seurat)
library(tidyverse)

Fibro_out_dir <- 'Fibroblast_subtype'
dir.create(Fibro_out_dir)
all_data <- readRDS(file.path("R_OUTPUT", "all_data.rds"))
all_data$sample %>% table()
DimPlot(all_data, group.by = 'annotation', label = T) + NoLegend()
FeaturePlot(all_data, features = 'Acta2')
VlnPlot(subset(all_data, Acta2 > 0), features = 'Acta2', group.by = 'annotation', pt.size = 0)

HG_would_id <- 
  read.csv('ROI/area1.csv') %>% 
  filter(area1 == 'area1') %>%
  select(Barcode) %>% 
  unlist()
HKPG_would_id <- 
  read.csv('ROI/area2.csv') %>% 
  filter(area2 == 'area2') %>%
  select(Barcode) %>% 
  unlist()
HG_wound <- subset(all_data, sample == 'HG' & cell_id_raw %in% HG_would_id)
HKPG_wound <- subset(all_data, sample == 'HKPG' & cell_id_raw %in% HKPG_would_id)

all_wound <- merge(HG_wound, HKPG_wound)
SpatialFeaturePlot(all_wound, features = "Phb2")

all_wound@images$HG <- all_wound@images$FWSP20251583.HKPG.HG.008um
all_wound@images$HKPG <- all_wound@images$FWSP20251583.HKPG.HG.008um.2
all_wound@images$FWSP20251583.HKPG.HG.008um <- NULL
all_wound@images$FWSP20251583.HKPG.HG.008um.2 <- NULL

all_wound$annotation %>% table()
all_wound_fibroblast <- subset(all_wound, annotation == 'Fibroblasts')


## PCA --------
all_wound_fibroblast <- FindVariableFeatures(all_wound_fibroblast)
all_wound_fibroblast <- ScaleData(all_wound_fibroblast, vars.to.regress = 'nCount_Spatial.008um')
all_wound_fibroblast <- RunPCA(all_wound_fibroblast)

pdf(file.path(Fibro_out_dir, 'pca plot.pdf'), height = 6, width = 8)
ElbowPlot(all_wound_fibroblast, ndims = 50)
dev.off()

## dimension reduction --------
dims <- 10

library(future)
plan(multisession, workers = 10)
options(future.globals.maxSize = 10 * 1024 ^ 3) # 10Gb

all_wound_fibroblast <- RunUMAP(all_wound_fibroblast, dims = 1:dims,return.model = TRUE)

pdf(file.path(Fibro_out_dir, str_glue("Dimplot dim{dims} sample.pdf")), height = 6, width = 8)
DimPlot(all_wound_fibroblast, group.by = "sample", label = TRUE) + NoLegend()
dev.off()

all_wound_fibroblast <- all_wound_fibroblast %>% 
  harmony::RunHarmony(
    dims.use = 1:dims,
    group.by.vars = "sample",
    reduction.save = "harmony"
    # theta = 3, # default 2, bigger value -- more diversity per cluster
    # lambda = 0.5, # default, smaller value -- more integrated
  )

all_wound_fibroblast <- RunUMAP(
  all_wound_fibroblast, 
  dims = 1:dims, 
  reduction = "harmony", 
  reduction.name = "harmonyumap", 
  reduction.key  = "harmonyUMAP_",
  return.model = TRUE
)

pdf(file.path(Fibro_out_dir, str_glue("Dimplot dim{dims} sample - harmony.pdf")), height = 6, width = 8)
DimPlot(all_wound_fibroblast, group.by = "sample", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()

reduction.use <- 'harmony' # 'harmony' 'pca'

### cluster -------
all_wound_fibroblast <- 
  FindNeighbors(
    all_wound_fibroblast ,
    reduction = reduction.use,
    dims = 1:dims,
    k.param = 20
  )

all_wound_fibroblast <- 
  FindClusters(
    all_wound_fibroblast,
    resolution = c(0.2, 0.4, 0.6, 0.8, 1.0)
  )

# dimplot
pdf(file.path(Fibro_out_dir, str_glue("Dimplot dim{dims} cluster.pdf")), height = 6, width = 8)
DimPlot(all_wound_fibroblast, group.by = "Spatial.008um_snn_res.0.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(all_wound_fibroblast, group.by = "Spatial.008um_snn_res.0.4", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(all_wound_fibroblast, group.by = "Spatial.008um_snn_res.0.6", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(all_wound_fibroblast, group.by = "Spatial.008um_snn_res.0.8", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(all_wound_fibroblast, group.by = "Spatial.008um_snn_res.1", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(all_wound_fibroblast, group.by = "Spatial.008um_snn_res.1.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(all_wound_fibroblast, group.by = "Spatial.008um_snn_res.1.5", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()


# HALLMARK scoring ----------------------
Fibro_TGFB <- 
  read_tsv('hallmark_gene_set/HALLMARK_TGF_BETA_SIGNALING.v2025.1.Mm.tsv') %>% 
  filter(STANDARD_NAME == 'GENE_SYMBOLS') %>% 
  .[1,2] %>% 
  unnest() %>% 
  deframe() %>% 
  .[1] %>% 
  strsplit(split = ',') %>% 
  .[[1]]
Fibro_WNT <- 
  read_tsv('hallmark_gene_set/HALLMARK_WNT_BETA_CATENIN_SIGNALING.v2025.1.Mm.tsv') %>% 
  filter(STANDARD_NAME == 'GENE_SYMBOLS') %>% 
  .[1,2] %>% 
  unnest() %>% 
  deframe() %>% 
  .[1] %>% 
  strsplit(split = ',') %>% 
  .[[1]]
Fibro_JAKSTAT <- 
  read_tsv('hallmark_gene_set/HALLMARK_IL6_JAK_STAT3_SIGNALING.v2025.1.Mm.tsv') %>% 
  filter(STANDARD_NAME == 'GENE_SYMBOLS') %>% 
  .[1,2] %>% 
  unnest() %>% 
  deframe() %>% 
  .[1] %>% 
  strsplit(split = ',') %>% 
  .[[1]]
GLYCOLYSIS <- 
  read_tsv('hallmark_gene_set/HALLMARK_GLYCOLYSIS.v2025.1.Mm.tsv') %>% 
  filter(STANDARD_NAME == 'GENE_SYMBOLS') %>% 
  .[1,2] %>% 
  unnest() %>% 
  deframe() %>% 
  .[1] %>% 
  strsplit(split = ',') %>% 
  .[[1]]
OXIDATIVE <- 
  read_tsv('hallmark_gene_set/HALLMARK_OXIDATIVE_PHOSPHORYLATION.v2025.1.Mm.tsv') %>% 
  filter(STANDARD_NAME == 'GENE_SYMBOLS') %>% 
  .[1,2] %>% 
  unnest() %>% 
  deframe() %>% 
  .[1] %>% 
  strsplit(split = ',') %>% 
  .[[1]]

pathway <- list(Fibro_TGFB, Fibro_WNT, Fibro_JAKSTAT, GLYCOLYSIS, OXIDATIVE)
names(pathway) <- c('Fibro_TGFB', 'Fibro_WNT', 'Fibro_JAKSTAT', 'GLYCOLYSIS', 'OXIDATIVE')

n_na_idx <- pathway %>% map_lgl(function(x) all(is.na(x))) %>% {! .} %>%  which
pathway <- pathway[n_na_idx] %>% map(function(x) x[!is.na(x)])

for ( i in seq_along(names(pathway))){
  all_wound_fibroblast <- AddModuleScore(JoinLayers(all_wound_fibroblast), pathway[i], name = names(pathway)[i], nbin=5, ctrl=10)
}
for ( i in seq_along(names(pathway))){
  all_wound <- AddModuleScore(JoinLayers(all_wound), pathway[i], name = names(pathway)[i], nbin=5, ctrl=10)
}

p <- VlnPlot(all_wound_fibroblast,
             features = c('Fibro_TGFB1', 'Fibro_WNT1', 'Fibro_JAKSTAT1', 'GLYCOLYSIS1', 'OXIDATIVE1'),
             group.by='Spatial.008um_snn_res.0.6',
             pt.size = 0)
p
FeaturePlot(all_wound_fibroblast, 
            reduction = 'harmonyumap', 
            features = c('Fibro_TGFB1', 'Fibro_WNT1', 'Fibro_JAKSTAT1', 'GLYCOLYSIS1', 'OXIDATIVE1'))
ggsave(
  p,
  filename = file.path(Fibro_out_dir, 'Hallmark_score.pdf'),
  height = 10,
  width  = 15
) 

VlnPlot(all_wound,
        features = c('Fibro_TGFB1', 'Fibro_WNT1', 'Fibro_JAKSTAT1', 'GLYCOLYSIS1', 'OXIDATIVE1'),
        group.by='annotation',
        pt.size = 0)
