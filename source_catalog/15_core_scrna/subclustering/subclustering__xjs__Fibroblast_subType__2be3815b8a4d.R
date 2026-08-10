library(Seurat)
library(tidyverse)

major_marker_ls <- c(
  'epithelials'  = 'EPCAM',
  'epithelials'  = 'KRT5',
  'epithelials'  = 'KRT14',
  
  "Foveolar cells" = "GKN1", #ACT
  "Foveolar cells" = "MUC5AC",
  "Foveolar cells" = "TFF1",
  
  'endothelials' = 'PECAM1', 
  'endothelials' = 'VWF', 
  
  "Fibroblasts" = "PDGFRA", 
  "Fibroblasts" = "CCDC80",
  "Fibroblasts" = "SFRP4",
  "Fibroblasts" = "COL1A1",
  "Fibroblasts" = "COL1A2",
  "Fibroblasts" = "DCN",
  "Fibroblasts" = "ACTA2",
  
  "B cells" = "CD19",
  "B cells" = "CD79A",
  "B cells" = "MS4A1",
  "Plasma" = "TNFRSF17",
  "Plasma" = "JCHAIN",
  "Plasma" = "MZB1",
  "Plasma" = "XBP1",
  
  'T cells'           = 'CD3D', 
  'T cells'           = 'CD3E', 
  'T cells'           = 'CD3G', 
  'T cells'           = 'TRAC', 
  
  "NK" = "GNLY",
  "NK" = "NKG7",
  "NK" = "NCAM1",
  
  'Proliferating' = 'MKI67',
  'Proliferating' = 'TOP2A',
  
  "Pericytes" = "PDGFRB",
  "Pericytes" = "RGS5",
  "Pericytes" = "MCAM",
  
  "Monocytes" = "CD14", 
  "Monocytes" = "FCGR3A",
  
  "Macrophage" = "CD68",
  "Macrophage" = "CD163",
  "Macrophage" = "C1QA",
  "Macrophage" = "C1QC", 
  'Macrophage'  = 'CSF1R',
  
  "Dendritic cell" = "CD83",
  'Dendritic cell' = "LILRA4", 
  'pDCs' = "IRF4", 
  'pDCs' = "IRF7", 
  'pDCs' = "IRF8", 
  'Dendritic cell' = "CLEC9A",
  'Dendritic cell' = "IDO1",
  'Dendritic cell' = "CD1C",
  'Dendritic cell' = "FCER1A",
  'Lamp3_DCs' = "LAMP3",
  
  'netrophil' = "CSF3R",
  'netrophil' = "G0S2",
  'netrophil' = "RETNLG",
  
  "Myocytes" = "TNNC1", 
  "Myocytes" = "MYH7", 
  
  "Mast cell" = "CPA3",
  "Mast cell" = "CTSG",
  "Mast cell" = "TPSAB1",
  "Mast cell" = "KIT"
)

theme_90angle <-
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
    strip.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )

all_data <- readRDS(file.path("R_OUTPUT", "all_data.rds"))
all_data$sample %>% table()
all_data$group %>% table()
all_data$annotation %>% table()
all_data$group2 %>% table()


# Fibroblasts --------------------------------------------------------------
Fibro_out_dir<-'Fibroblast_subtype_remove_D2' 
dir.create(Fibro_out_dir)
Fibroblast <- subset(all_data, annotation == 'Fibroblast' & sample != 'D2')

## PCA --------
DefaultAssay(Fibroblast) <- "RNA"
Fibroblast <- FindVariableFeatures(Fibroblast)
Fibroblast <- ScaleData(Fibroblast, vars.to.regress = 'nCount_RNA')
Fibroblast <- RunPCA(Fibroblast)


### pca inspector ------
pdf(file.path(Fibro_out_dir, 'pca plot.pdf'), height = 6, width = 8)
ElbowPlot(Fibroblast, ndims = 50)
dev.off()


## dimension reduction --------
dims <- 15

library(future)
plan(multisession, workers = 10)
options(future.globals.maxSize = 10 * 1024 ^ 3) # 10Gb

Fibroblast <- RunUMAP(Fibroblast, dims = 1:dims,return.model = TRUE)

pdf(file.path(Fibro_out_dir, str_glue("Dimplot dim{dims} sample.pdf")), height = 6, width = 8)
DimPlot(Fibroblast, group.by = "sample", label = TRUE) + NoLegend()
dev.off()

Fibroblast <- Fibroblast %>% 
  harmony::RunHarmony(
    dims.use = 1:dims,
    group.by.vars = "sample",
    reduction.save = "harmony"
    # theta = 3, # default 2, bigger value -- more diversity per cluster
    # lambda = 0.5, # default, smaller value -- more integrated
  )

Fibroblast <- RunUMAP(
  Fibroblast, 
  dims = 1:dims, 
  reduction = "harmony", 
  reduction.name = "harmonyumap", 
  reduction.key  = "harmonyUMAP_",
  return.model = TRUE
)

pdf(file.path(Fibro_out_dir, str_glue("Dimplot dim{dims} sample - harmony.pdf")), height = 6, width = 8)
DimPlot(Fibroblast, group.by = "sample", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()

reduction.use <- 'harmony' # 'harmony' 'pca'

### cluster -------
Fibroblast  <-
  FindNeighbors(
    Fibroblast ,
    reduction = reduction.use,
    dims = 1:dims,
    k.param = 15
  )

Fibroblast <-
  FindClusters(
    Fibroblast,
    resolution = c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.5)
  )

# dimplot
pdf(file.path(Fibro_out_dir, str_glue("Dimplot dim{dims} cluster.pdf")), height = 6, width = 8)
DimPlot(Fibroblast, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Fibroblast, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Fibroblast, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Fibroblast, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Fibroblast, group.by = "RNA_snn_res.1", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Fibroblast, group.by = "RNA_snn_res.1.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Fibroblast, group.by = "RNA_snn_res.1.5", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()


### dotplot -------
Fibro_major_marker_ls1 <- c(
  'mCAF'='MMP11',
  'mCAF'='MYL9',
  'mCAF'='POSTN',
  'mCAF'='TPM2',
  'mCAF'='TAGLN',
  'mCAF'='HOPX',
  # 'mCAF'='DCN',
  # 'mCAF'='CDH11',
  # 'mCAF'='LUM',
  'iCAF'='PLA2G2A',
  'iCAF'='CFD',
  'iCAF'='CD34',
  'iCAF'='CXCL14',
  'iCAF'='PI16',
  'iCAF'='PRG4',
  'iCAF'='APOD',
  'iCAF'='DPT',
  'iCAF'='CXCL12',
  'iCAF'='CXCL1',
  'iCAF'='IL6',
  'iCAF'='CCL2',
  'apCAF'='HLA-DRA',
  'apCAF'='HLA-DRB1',
  'apCAF'='CD74',
  'apCAF'='IGLC2',
  'apCAF'='IGKC',
  'apCAF'='CXCR4',
  'vCAF'='ACTA2',
  'vCAF'='MCAM',
  'vCAF'='RGS5',
  'vCAF'='MYH11',
  'vCAF'='ADIRF',
  'tCAF'='ENO1',
  'tCAF'='GAPDH',
  'tCAF'='NDRG1',
  # 'tCAF'='HSPH1',
  # 'tCAF'='VEGFA',
  'tCAF'='MME',
  'dCAF'='MKI67',
  'ifnCAF'='IDO1',
  'ifnCAF'='CXCL9',
  'ifnCAF'='CXCL10',
  'ifnCAF'='CXCL11',
  'rCAF'='CCL19',
  'rCAF'='CCL21'
  # 'hsptcaf'='HSPH1',
  # 'hsptcaf'='VEGFA'
)

Fibro_major_marker_ls2 <- c(
  NF = 'IGFBP6', 
  NF = 'PI16', 
  # NF = 'DPT', 
  # NF = 'CXCL12', 
  NF = 'GPX3', 
  NF = 'TNXB', 
  CAF = 'FAP',
  CAF = 'PDGFRA',
  CAF = 'PDGFRB',
  CAF = 'HSF1',
  CAF = 'STAT3',
  CAF = 'MYC',
  CAF = 'CD10',
  CAF = 'GPR77',
  iCAF = 'DPT',
  iCAF = 'IGF1',
  iCAF = 'ITM2A',
  iCAF = 'CXCL12',
  iCAF = 'CCL11',
  iCAF = 'CD34',
  iCAF = 'IL6',
  myCAF = 'MMP11',
  myCAF = 'CTHRC1',
  myCAF = 'COL11A1',
  myCAF = 'POSTN',
  myCAF = 'COL10A1',
  myCAF = 'TNC',
  meCAF = 'SLC2A1',
  meCAF = 'BNIP3',
  meCAF = 'ADM',
  meCAF = 'ANGPTL4',
  apCAF = 'FSP1',
  pCAF = 'NUSAP1',
  pCAF = 'TOP2A',
  pCAF = 'CENPF',
  pCAF = 'PTTG1',
  pCAF = 'SCG2',
  pCAF = 'STMN1',
  pCAF = 'PLAU',
  pCAF = 'IGFBP2',
  vCAF = 'NOTCH3',
  vCAF = 'COL18A1',
  vCAF = 'MCAM',
  vCAF = 'RGS5',
  vCAF = 'EPAS1',
  tCAF = 'PDPN',
  tCAF = 'MME',
  tCAF = 'NDRG1',
  tCAF = 'VEGFA',
  ifnCAF = 'IDO1',
  ifnCAF = 'CXCL9',
  ifnCAF = 'CXCL10',
  ifnCAF ='CXCL11',
  rCAF = 'CCL19',
  dCAF = 'MKI67',
  dCAF = 'SOX9',
  dCAF = 'SOX10'
)

Fibro_major_marker_ls3 <- c(
  'iCAF' = "IL11", 
  'iCAF' = "IL24", 
  'iCAF' = "MME", 
  'iCAF' = "CD10", 
  'iCAF' = "ANPEP", 
  'iCAF' = "CD13", 
  
  "POSTN",
  "CST1", 
  "SFRP1", 
  'mCAF' = "LUM",
  'mCAF' = "LDC",
  'mCAF' = "SLRR2D",
  
  "MYLK",
  "RGS5",
  'vCAF' = "MCAM"
)

Idents(Fibroblast)<-'RNA_snn_res.1'
Fibroblast_markers <-
  FindAllMarkers(
    JoinLayers(Fibroblast),
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    only.pos = TRUE
  )

Fibroblast_top5_markers <- 
  Fibroblast_markers %>% 
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 30) %>%
  arrange(desc(pct.1)) %>%
  slice_head(n = 5) %>%
  nest() %>%
  deframe() %>% 
  map(pull, 'gene') %>%
  imap(function(genes, name){
    re <- genes %>%
      map(function(gene){
        # Fibroblastrowser()
        names(gene) = str_glue("cluster_{name}")
        gene
      })
    re
  }) %>%
  unlist() %>% 
  .[!duplicated(.)]

names(Fibroblast_top5_markers) <- sub('.*\\.','',names(Fibroblast_top5_markers))
Fibroblast_markers %>% write_csv(file.path(Fibro_out_dir, 'Fibroblast_DEGs_RNA_snn_res.1.csv'))


# dotplot
pdf(file.path(Fibro_out_dir, str_glue("dotplot dims{dims} major_marker.pdf")), width = 25, height = 6)
DotPlot(Fibroblast, group.by = "RNA_snn_res.1", features = major_marker_ls)+RotatedAxis()
DotPlot(Fibroblast, group.by = "RNA_snn_res.1", features = Fibro_major_marker_ls1)+RotatedAxis()
DotPlot(Fibroblast, group.by = "RNA_snn_res.1", features = Fibro_major_marker_ls2)+RotatedAxis()
DotPlot(Fibroblast, group.by = "RNA_snn_res.1", features = Fibro_major_marker_ls3)+RotatedAxis()
DotPlot(Fibroblast, group.by = "RNA_snn_res.1", features = Fibroblast_top5_markers)+RotatedAxis()
# DotPlot(Fibroblast, group.by = "RNA_snn_res.0.8", features = macrophage_major_marker_ls)+RotatedAxis()
dev.off()

# Idents(Fibroblast) <- "RNA_snn_res.0.6"
# Fibroblast <- FindSubCluster(Fibroblast, cluster = 1, graph.name = "RNA_snn", algorithm = 1,resolution = 0.3)
# Idents(Fibroblast) <- "sub.cluster"
# DimPlot(Fibroblast, label = TRUE, reduction = 'harmonyumap') + NoLegend()
# # DimPlot(Fibroblast, label = TRUE, reduction = 'harmonyumap', split.by = "sub.cluster") + NoLegend()
# # DimPlot(Fibroblast, label = TRUE, reduction = 'harmonyumap', split.by = "RNA_snn_res.0.6") + NoLegend()
# DotPlot(Fibroblast, group.by = "sub.cluster", features = Fibro_major_marker_ls2)+RotatedAxis()

###score--------------------
iCAF <- c('DPT','IGF1','ITM2A','CXCL12','CCL11','CD34','IL6')
mCAF <-c('MMP11','COL11A1','POSTN','COL10A1')
myCAF <- c('ACTA2','TAGLN','MMP11', 'MYL9','HOPX','POSTN','TPM1','TPM2')
meCAF <- c('SLC2A1','BNIP3','ADM','ANGPTL4')
pCAF <- c('NUSAP1','TOP2A','CENPF','PTTG1','SCG2','STMN1','PLAU','IGFBP2')
vCAF <- c('NOTCH3','COL18A1','MCAM','RGS5','EPAS1')
tCAF <- c('PDPN','MME','NDRG1','VEGFA')
infCAF <- c('IDO1','CXCL9','CXCL10','CXCL11')
rCAF <- c('CCL19')
dCAF <- c('MKI67','SOX9','SOX10')

Macrophages <- c("CD68","CD163","C1QA","C1QC",'CSF1R')
Neutrophils <- c("FCGR3B","IFITM2","CSF3R")
Plasma <- c("TNFRSF17","JCHAIN","MZB1","XBP1","IGHG1","IGHA1","IGHG3","IGKC")
Tcells <- c("CD3D","CD3E","CD3G","LCK","TRAC","CD4","CD8A","GZMK")
Fibroblasts <- c("SFRP4","PDGFRA","CCDC80","COL1A1","ACTA2",'DCN',"FAP","TNC","PDPN","PDGFRB","COL4A2","COL5A2","STC2","LOXL2","LUM","VIM","COL1A2")
Bcells <- c("CD19","CD79A","MS4A1")
DC <- c("CD83","LILRA4","CLEC9A","IDO1","CD1C","FCER1A")

pathway<-list(iCAF,pCAF,mCAF,myCAF,meCAF,vCAF,tCAF,infCAF,rCAF,dCAF,Macrophages,Neutrophils,Plasma,Tcells,Fibroblasts,Bcells,DC)
names(pathway)<-c('iCAF','pCAF','mCAF','myCAF','meCAF','vCAF','tCAF','infCAF','rCAF','dCAF','Macrophages','Neutrophils','Plasma','Tcells','Fibroblasts','Bcells','DC')

# M1_M2<-read.table("/home/lr/ATAC_thirteen/Immune/MG/M1_M2.txt",header=TRUE,fill=TRUE)
n_na_idx <- pathway %>% map_lgl(function(x) all(is.na(x))) %>% {! .}%>%  which
pathway <- pathway[n_na_idx] %>% map(function(x) x[!is.na(x)])

for ( i in seq_along(names(pathway))){
  Fibroblast <- AddModuleScore(JoinLayers(Fibroblast), pathway[i], name =  names(pathway)[i],nbin=5,ctrl=10)
}

p<-VlnPlot(Fibroblast,
           features = c('CAF1','iCAF1','pCAF1','mCAF1','myCAF1','meCAF1','vCAF1','tCAF1','infCAF1','rCAF1','dCAF1',
                        'Macrophage1','Neutrophils1','Plasma1','Tcells1','Fibroblasts1','Bcells1','DC1'),
           group.by='RNA_snn_res.1',
           combine = T,
           pt.size = 0)

FeaturePlot(Fibroblast,
            features = c('CAF1','iCAF1','pCAF1','mCAF1','myCAF1','meCAF1','vCAF1','tCAF1','infCAF1','rCAF1','dCAF1',
                         'Macrophage','Neutrophils1','Plasma1','Tcells1','Fibroblasts1','Bcells1','DC1'))

ggsave(
  p,
  filename = file.path(Fibro_out_dir, 'CAF score.pdf'),
  height = 20,
  width  = 15
) 


## annotation ----
VlnPlot(Fibroblast, features = 'nCount_RNA', group.by = 'RNA_snn_res.1', pt.size = 0)

cell_cluster <- c(
  '0'  = "LENG8-mCAF", 
  '1'  = "mix", # Epi 
  '2'  = "low-nCounts", 
  '3'  = "mix", # Macro
  '4'  = "THBS1-mCAF", 
  '5'  = "iCAF", 
  '6'  = "tCAF", 
  '7'  = "ifnCAF", 
  '8'  = "vCAF", 
  '9'  = "mix", # Macro
  '10' = "mix", # Epi 
  '11' = "mix", # B, T
  '12' = "mix" # Epi
)

Idents(Fibroblast) <- "RNA_snn_res.1"

Fibroblast <- RenameIdents(Fibroblast, cell_cluster)
Fibroblast$subtype <- Idents(Fibroblast)

pdf(file.path(Fibro_out_dir,str_glue("Dimplot dim{dims} annotation.pdf")), height = 6, width = 8)
# DimPlot(Fibroblast, group.by = "subtype", label = TRUE, reduction = 'umap',repel = TRUE) + NoLegend()
DimPlot(Fibroblast, group.by = "subtype", label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

pdf(file.path(Fibro_out_dir, str_glue("Dimplot dim{dims} annotation split.pdf")), height = 6, width = 14)
DimPlot(Fibroblast, group.by = "subtype",split.by = 'group', label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

saveRDS(Fibroblast, file.path(Fibro_out_dir, 'Fibroblast.rds'))


Fibroblast_pure <- subset(JoinLayers(Fibroblast), subtype %in% c("LENG8-mCAF", "THBS1-mCAF", "vCAF", "ifnCAF", "tCAF", "iCAF"))
Fibroblast_pure <- RunUMAP(Fibroblast_pure, dims = 1:dims, return.model = TRUE)
DimPlot(Fibroblast_pure, group.by = "subtype", label = TRUE, reduction = 'harmonyumap',repel = TRUE)

Fibroblast_pure <- FindVariableFeatures(Fibroblast_pure)
Fibroblast_pure <- ScaleData(Fibroblast_pure, vars.to.regress = 'nCount_RNA')
Fibroblast_pure <- RunPCA(Fibroblast_pure)

pdf(file.path(Fibro_out_dir, 'pca plot - pure.pdf'), height = 6, width = 8)
ElbowPlot(Fibroblast_pure, ndims = 50)
dev.off()

Fibroblast_pure <- RunUMAP(Fibroblast_pure, dims = 1:dims, return.model = TRUE)
DimPlot(Fibroblast_pure, group.by = "subtype", label = TRUE, reduction = 'harmonyumap',repel = TRUE)

Fibroblast_pure <- Fibroblast_pure %>% 
  harmony::RunHarmony(
    dims.use = 1:dims,
    group.by.vars = "sample",
    reduction.save = "harmony"
    # theta = 3, # default 2, bigger value -- more diversity per cluster
    # lambda = 0.5, # default, smaller value -- more integrated
  )

Fibroblast_pure <- RunUMAP(
  Fibroblast_pure, 
  dims = 1:dims, 
  reduction = "harmony", 
  reduction.name = "harmonyumap", 
  reduction.key  = "harmonyUMAP_",
  return.model = TRUE
)

DimPlot(Fibroblast_pure, group.by = "RNA_snn_res.1", label = TRUE, reduction = 'harmonyumap',repel = TRUE)
DimPlot(Fibroblast_pure, group.by = "subtype", label = TRUE, reduction = 'harmonyumap',repel = TRUE)

pdf(file.path(Fibro_out_dir,str_glue("Dimplot dim{dims} annotation - pure.pdf")), height = 6, width = 9)
DimPlot(Fibroblast_pure, group.by = "subtype", label = TRUE, reduction = 'harmonyumap',repel = TRUE)
dev.off()

pdf(file.path(Fibro_out_dir, str_glue("Dimplot dim{dims} annotation split - pure.pdf")), height = 7, width = 16)
DimPlot(Fibroblast_pure, group.by = "subtype",split.by = 'group', label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()





## cell props barplot  --------
Fibroblast_meta <- Fibroblast_pure[[]]

p_col_proj <-
  Fibroblast_meta %>%
  group_by(sample, subtype) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = sample, y = freq, fill = subtype)) + 
  geom_col() +
  ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )


# export
ggsave(
  p_col_proj,
  filename = file.path(Fibro_out_dir, "cell_freq_sample.pdf"),
  height = 6,
  width  = 6
)

p_col_proj <-
  Fibroblast_meta %>%
  group_by(group, subtype) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = group, y = freq, fill = subtype)) + 
  geom_col() +
  ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )

ggsave(
  p_col_proj,
  filename = file.path(Fibro_out_dir, "cell_freq_group.pdf"),
  height = 6,
  width  = 4
)

p_col_proj <-
  Fibroblast_meta %>%
  group_by(group2, subtype) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = group2, y = freq, fill = subtype)) + 
  geom_col() +
  ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )

ggsave(
  p_col_proj,
  filename = file.path(Fibro_out_dir, "cell_freq_group2.pdf"),
  height = 6,
  width  = 4
)


cellType_marker_ls <- c(
  'mCAF' = 'MMP11', 
  'mCAF' = 'POSTN', 
  "LENG8-mCAF" = 'LENG8', 
  "THBS1-mCAF" = 'THBS1', 
  
  "iCAF" = 'IGF1', 
  "iCAF" = 'ITM2A', 
  "iCAF" = 'CXCL12', 
  "iCAF" = 'CD34', 
  
  "tCAF" = 'NDRG1', 
  "tCAF" = 'MME', 
  "tCAF" = 'SLC2A1', 
  "tCAF" = 'BNIP3', 
  "tCAF" = 'ADM', 
  "tCAF" = 'ANGPTL4', 
  
  "ifnCAF" = 'IDO1', 
  "ifnCAF" = 'CXCL9', 
  "ifnCAF" = 'CXCL10', 
  "ifnCAF" = 'CXCL11', 
  
  "vCAF" = 'NOTCH3', 
  "vCAF" = 'COL18A1', 
  "vCAF" = 'MCAM', 
  "vCAF" = 'RGS5'
)

pdf(file.path(Fibro_out_dir, "dotplot celltype main_marker.pdf"), width = 12, height = 6)
DotPlot(Fibroblast, group.by = "subtype", features = cellType_marker_ls, dot.scale = 10) + theme_90angle
dev.off()

pdf(file.path(Fibro_out_dir, "dotplot celltype main_marker - pure.pdf"), width = 12, height = 5)
DotPlot(Fibroblast_pure, group.by = "subtype", features = cellType_marker_ls, dot.scale = 10) + theme_90angle
dev.off()

saveRDS(Fibroblast_pure, file.path(Fibro_out_dir, 'Fibroblast_pure.rds'))




# enrichment -------------
Idents(Fibroblast_pure) <- "subtype"
Fibroblast_pure_markers <-
  FindAllMarkers(
    JoinLayers(Fibroblast_pure),
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    only.pos = TRUE
  )
write.csv(Fibroblast_pure_markers,file.path(Fibro_out_dir, 'Fibroblast_pure_DEGs_subtype.csv'))

library(EnhancedVolcano)
type<-unique(Fibroblast_pure_markers$cluster)
volcano_ls <- type %>% 
  setNames(.,.) %>%
  imap(function(type,name){
    markers<- Fibroblast_pure_markers[Fibroblast_pure_markers$cluster == type, ]
    p<-EnhancedVolcano(markers,
                       lab = rownames(markers),
                       x = 'avg_log2FC',
                       y = 'p_val_adj',
                       title = str_glue('{name} volcano plot'),
                       pCutoff = 10e-2,
                       FCcutoff = 0,
                       pointSize = 3.0,
                       labSize = 6.0)
    p
  })
volcano_ls %>%
  iwalk(function(x, name){
    ggsave(
      x,
      filename = file.path(Fibro_out_dir, str_glue("{name} diffgene valcano.pdf")),
      height = 6,
      width = 8
    )
  })


all_markers_major <- list(Fibroblast_pure = Fibroblast_pure_markers)
# scatter
p_volcanoplot_major_ls <- 
  all_markers_major %>%
  imap(function(dat, name){
    
    set.seed(1234)
    rect_height <- 1
    dodge_width <- 0.2
    
    df_scatter <-
      dat %>%
      mutate(x_base = as.numeric(factor(cluster))) %>%
      slice_sample(prop = 1) %>%
      group_by(cluster) %>%
      
      mutate(
        x_delta = 1:n(),
        x_delta = scales::rescale(x_delta, to = c(dodge_width, 1 - dodge_width)),
        xx = x_base - 1 + x_delta,
        y = ifelse(avg_log2FC > 0, avg_log2FC + rect_height, avg_log2FC - rect_height)
      ) %>%
      mutate(
        group = ifelse(p_val_adj > 0.01, "p > 0.01", "p < 0.01"),
        group = factor(group, levels = c("p > 0.01", "p < 0.01"))
      ) %>%
      dplyr::arrange(group)
    
    # df_label 
    df_label_up <- df_scatter %>%
      dplyr::filter(as.character(group) == "p < 0.01") %>%
      group_by(cluster) %>% 
      slice_max(avg_log2FC, n = 10)
    
    df_label_down <- df_scatter %>%
      dplyr::filter(as.character(group) == "p < 0.01") %>%
      group_by(cluster) %>% 
      slice_min(avg_log2FC, n = 10)
    
    
    
    # grp rect  
    df_grp_rect <-
      dat %>%
      group_by(cluster) %>%
      summarise() %>%
      mutate(x_base = as.numeric(factor(cluster))) %>%
      mutate(xmin = x_base - 1, xmax = x_base) %>%
      mutate(ymin = -rect_height, ymax = rect_height)
    
    
    # background rect
    df_bg_rect <-
      dat %>%
      group_by(cluster) %>%
      summarise(ymax = max(avg_log2FC) + rect_height, ymin = min(avg_log2FC) - rect_height) %>%
      mutate(x_base = as.numeric(factor(cluster))) %>%
      mutate(xmin = x_base - 1 + dodge_width, xmax = x_base - dodge_width)
    
    
    
    df_scatter %>%
      ggplot(aes(x = xx, y = y)) +
      geom_rect(
        data = df_bg_rect, 
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
        fill = "gray90",
        inherit.aes = FALSE
      ) +
      geom_rect(
        data = df_grp_rect,
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = cluster),
        inherit.aes = FALSE
      ) +
      geom_point(aes(color = group), size = 1) +
      
      ggrepel::geom_text_repel(data = df_label_up, aes(label = gene)) +
      ggrepel::geom_text_repel(data = df_label_down, aes(label = gene)) +
      scale_color_manual(values = c(
        'p < 0.01' = "red",
        'p > 0.01' = "gray70"
      )) +
      # scale_fill_manual(values = grp_cols_sc) +
      theme_minimal() +
      theme(
        panel.grid = element_blank(),
        axis.line.y = element_line(color = "black"),
        axis.text.x = element_blank()
      ) +
      labs(
        x = NULL,
        y = "Log2FC",
        title = name
      )
    
  })

p_volcanoplot_major_ls %>%
  iwalk(function(x, name){
    ggsave(
      x,
      filename = file.path(Fibro_out_dir, str_glue("{name} diffgene dot.pdf")),
      height = 8,
      width = 15
    )
  })


Fibroblast_pure_de_genes_for_enrichmenet <- 
  Fibroblast_pure_markers %>% 
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  dplyr::filter(pct.1 > 0.1) %>%
  dplyr::filter(abs(avg_log2FC) > 0.5) %>% 
  group_by(cluster) %>%
  nest() %>%
  deframe() %>%
  map("gene")


enrich_ls <-
  Fibroblast_pure_de_genes_for_enrichmenet %>%
  map(function(x){
    clusterProfiler::enrichGO(
      x,
      OrgDb = org.Hs.eg.db::org.Hs.eg.db, 
      ont='ALL',
      pAdjustMethod = 'BH',
      pvalueCutoff = 0.05, 
      qvalueCutoff = 0.05,
      keyType = 'SYMBOL'
    )
  })


Fibroblast_pure_er_output <- file.path(Fibro_out_dir, 'enrichment') # dir.create(Fibroblast_pure_er_output)
# expot enrichment table
enrich_ls %>% 
  # discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    dir.create(file.path(Fibroblast_pure_er_output, name))
    
    write.csv(
      x@result,
      file = file.path(Fibroblast_pure_er_output, name, "Go enriment table.csv"),
      row.names = FALSE
    )
  })

# dotplot
enrich_dot_ls <-
  enrich_ls %>%
  map(function(x){
    # if(nrow(as.data.frame(x)) < 1){
    #   return(NULL)
    # }
    enrichplot::dotplot(
      x,
      showCategory = 10,
      label_format = 50
    )
  })

enrich_bar_ls <-
  enrich_ls %>%
  map(function(x){
    # if(nrow(as.data.frame(x)) < 1){
    #   return(NULL)
    # }
    barplot(
      x,
      showCategory = 10,
      label_format = 50
    )
  })

# export dotplot
enrich_dot_ls %>%
  # discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    # cluster_name <- paste0("cluster_", name)
    # dir.create(file.path(er_output, cluster_name))
    
    ggsave(
      x,
      filename = file.path(Fibroblast_pure_er_output, name, "go dotplot.pdf"),
      height = 6,
      width = 8
    )
  })

enrich_bar_ls %>%
  # discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    # cluster_name <- paste0("cluster_", name)
    # dir.create(file.path(er_output, cluster_name))
    
    ggsave(
      x,
      filename = file.path(Fibroblast_pure_er_output, name, "go barplot.pdf"),
      height = 6,
      width = 8
    )
  })

Fibroblast_pure_de_genes_id_for_enrichmenet <-
  Fibroblast_pure_de_genes_for_enrichmenet %>%
  map(function(x){
    clusterProfiler::bitr(
      x,
      OrgDb = "org.Hs.eg.db",
      fromType = "SYMBOL",
      toType = "ENTREZID"
    ) %>%
      .$ENTREZID
  }) %>%
  purrr::discard(function(x) length(x) < 1)


enrich_ls <-
  Fibroblast_pure_de_genes_id_for_enrichmenet %>%
  map(function(x){
    clusterProfiler::enrichKEGG(
      x,
      organism = "hsa", 
      keyType = "ncbi-geneid", 
      pvalueCutoff = 0.05, 
      pAdjustMethod = "BH", 
      minGSSize = 10,
      maxGSSize = 500, 
      qvalueCutoff = 0.2,
      use_internal_data = FALSE
    )
  })

# saveRDS(enrich_ls,file.path(Fibroblast_pure_er_output, 'enrich_ls.rds'))

# expot enrichment table
enrich_ls %>%
  # discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    dir.create(file.path(Fibroblast_pure_er_output, name))
    
    write.csv(
      x@result,
      file = file.path(Fibroblast_pure_er_output, name, "KEGG enriment table.csv"),
      row.names = FALSE
    )
  })

# dotplot
enrich_dot_ls <-
  enrich_ls %>%
  map(function(x){
    if(nrow(as.data.frame(x)) < 1){
      return(NULL)
    }
    enrichplot::dotplot(
      x,
      showCategory = 10,
      label_format = 50
    )
  })

enrich_bar_ls <-
  enrich_ls %>%
  map(function(x){
    if(nrow(as.data.frame(x)) < 1){
      return(NULL)
    }
    barplot(
      x,
      showCategory = 10,
      label_format = 50
    )
  })

# export dotplot
enrich_dot_ls %>%
  purrr::discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    # cluster_name <- paste0("cluster_", name)
    # dir.create(file.path(er_output, cluster_name))
    
    ggsave(
      x,
      filename = file.path(Fibroblast_pure_er_output, name, "kegg dotplot.pdf"),
      height = 6,
      width = 8
    )
  })

enrich_bar_ls %>%
  purrr::discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    # cluster_name <- paste0("cluster_", name)
    # dir.create(file.path(er_output, cluster_name))
    
    ggsave(
      x,
      filename = file.path(Fibroblast_pure_er_output, name, "kegg barplot.pdf"),
      height = 6,
      width = 8
    )
  })



