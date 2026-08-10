library(Seurat)
library(tidyverse)
library(BPCells)
library(clusterProfiler)
library(org.Hs.eg.db)

all_data <- readRDS(file.path("R_OUTPUT", "all_data_remove_Tumor3.rds"))

# markers: https://cloud.tencent.com/developer/article/2115987
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




all_data$annotation %>% table()
all_data$sample %>% table()
Fibro_out_dir<-'Fibroblast_subtype' # dir.create(Fibro_out_dir)
Fibroblast<-subset(all_data, annotation == 'Fibroblast')
Fibroblast$sample %>% table()

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
dims <- 20

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
    k.param = 20
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
  NF = 'IGFBP6', 
  # NF = 'PI16', 
  # NF = 'DPT', 
  # NF = 'CXCL12', 
  NF = 'GPX3', 
  NF = 'TNXB', 
  CAF = 'DCN',
  CAF = 'COL1A1',
  CAF = 'COL3A1',
  CAF = 'COL6A1',
  # CAF = 'ACTA2',
  CAF = 'FAP',
  CAF = 'CD10',
  CAF = 'MMP1', 
  # CAF = 'MMP11',
  CAF = 'RUNX1',
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
  'infCAF'='IDO1',
  'infCAF'='CXCL9',
  'infCAF'='CXCL10',
  'infCAF'='CXCL11',
  'rCAF'='CCL19',
  'rCAF'='CCL21', 
  # 'hsptcaf'='HSPH1',
  # 'hsptcaf'='VEGFA', 
  # 'myCAF' = 'ACTA2', 
  # 'myCAF' = 'TAGLN', 
  # 'myCAF' = 'MYL9', 
  # 'myCAF' = 'HOPX', 
  'myCAF' = 'TPM1'
  # 'myCAF' = 'TPM2'
)

Fibro_major_marker_ls2 <- c(
  iCAF = 'DPT',
  iCAF = 'IGF1',
  iCAF = 'ITM2A',
  iCAF = 'CXCL12',
  iCAF = 'CCL11',
  iCAF = 'CD34',
  iCAF = 'IL6',
  mCAF = 'MMP11',
  mCAF = 'CTHRC1',
  mCAF = 'COL11A1',
  mCAF = 'POSTN',
  mCAF = 'COL10A1',
  mCAF = 'TNC',
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
  infCAF = 'IDO1',
  infCAF = 'CXCL9',
  infCAF = 'CXCL10',
  infCAF ='CXCL11',
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

pdf(file.path(Fibro_out_dir, str_glue('Dotplot dim{dims} major marker_multi res.pdf')), width = 15, height = 6)
DotPlot(Fibroblast, features = Fibro_major_marker_ls2, group.by = "RNA_snn_res.0.2") + theme_90angle
DotPlot(Fibroblast, features = Fibro_major_marker_ls2, group.by = "RNA_snn_res.0.4") + theme_90angle
DotPlot(Fibroblast, features = Fibro_major_marker_ls2, group.by = "RNA_snn_res.0.6") + theme_90angle
DotPlot(Fibroblast, features = Fibro_major_marker_ls2, group.by = "RNA_snn_res.0.8") + theme_90angle
DotPlot(Fibroblast, features = Fibro_major_marker_ls2, group.by = "RNA_snn_res.1") + theme_90angle
DotPlot(Fibroblast, features = Fibro_major_marker_ls2, group.by = "RNA_snn_res.1.2") + theme_90angle
DotPlot(Fibroblast, features = Fibro_major_marker_ls2, group.by = "RNA_snn_res.1.5") + theme_90angle
dev.off()

Idents(Fibroblast)<-'RNA_snn_res.0.8'
Fibroblast_markers <-
  FindAllMarkers(
    JoinLayers(Fibroblast),
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    only.pos = TRUE
  )
Fibroblast_markers %>% openxlsx::write.xlsx(file.path(Fibro_out_dir, 'Fibroblast_markers_res0.8.xlsx'))

Fibroblast_top5_markers <- 
  Fibroblast_markers %>% 
  group_by(cluster) %>%
  dplyr::filter(
    p_val_adj < 0.05, 
    avg_log2FC > 1, 
    pct.1-pct.2 > 0.1
    ) %>%
  # arrange(desc(avg_log2FC)) %>%
  # slice_head(n = 30) %>%
  arrange(desc(pct.1 - pct.2)) %>%
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

# dotplot
pdf(file.path(Fibro_out_dir, str_glue("dotplot dims{dims} major_marker.pdf")), width = 25, height = 6)
DotPlot(Fibroblast, group.by = "RNA_snn_res.0.8", features = major_marker_ls)+RotatedAxis()
DotPlot(Fibroblast, group.by = "RNA_snn_res.0.8", features = Fibro_major_marker_ls1)+RotatedAxis()
DotPlot(Fibroblast, group.by = "RNA_snn_res.0.8", features = Fibro_major_marker_ls2)+RotatedAxis()
DotPlot(Fibroblast, group.by = "RNA_snn_res.0.8", features = Fibro_major_marker_ls3)+RotatedAxis()
DotPlot(Fibroblast, group.by = "RNA_snn_res.0.8", features = Fibroblast_top5_markers)+RotatedAxis()
dev.off()

# Idents(Fibroblast) <- "RNA_snn_res.0.6"
# Fibroblast <- FindSubCluster(Fibroblast, cluster = 1, graph.name = "RNA_snn", algorithm = 1,resolution = 0.3)
# Idents(Fibroblast) <- "sub.cluster"
# DimPlot(Fibroblast, label = TRUE, reduction = 'harmonyumap') + NoLegend()
# # DimPlot(Fibroblast, label = TRUE, reduction = 'harmonyumap', split.by = "sub.cluster") + NoLegend()
# # DimPlot(Fibroblast, label = TRUE, reduction = 'harmonyumap', split.by = "RNA_snn_res.0.6") + NoLegend()
# DotPlot(Fibroblast, group.by = "sub.cluster", features = Fibro_major_marker_ls2)+RotatedAxis()

###score--------------------
CAF <- c('FAP','PDGFRA','PDGFRB','HSF1','STAT3','MYC','CD10','GPR77')
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

pathway<-list(CAF,iCAF,pCAF,mCAF,myCAF,meCAF,vCAF,tCAF,infCAF,rCAF,dCAF,Macrophages,Neutrophils,Plasma,Tcells,Fibroblasts,Bcells,DC)
names(pathway)<-c('CAF','iCAF','pCAF','mCAF','myCAF','meCAF','vCAF','tCAF','infCAF','rCAF','dCAF','Macrophages','Neutrophils','Plasma','Tcells','Fibroblasts','Bcells','DC')

# M1_M2<-read.table("/home/lr/ATAC_thirteen/Immune/MG/M1_M2.txt",header=TRUE,fill=TRUE)
n_na_idx <- pathway %>% map_lgl(function(x) all(is.na(x))) %>% {! .}%>%  which
pathway <- pathway[n_na_idx] %>% map(function(x) x[!is.na(x)])

for ( i in seq_along(names(pathway))){
  Fibroblast <- AddModuleScore(JoinLayers(Fibroblast), pathway[i], name =  names(pathway)[i],nbin=5,ctrl=10)
}

p<-VlnPlot(Fibroblast,
           features = c('CAF1','iCAF1','pCAF1','mCAF1','myCAF1','meCAF1','vCAF1','tCAF1','infCAF1','rCAF1','dCAF1',
                        'Macrophage1','Neutrophils1','Plasma1','Tcells1','Fibroblasts1','Bcells1','DC1'),
           group.by='RNA_snn_res.0.6',
           combine = T,
           pt.size = 0)

FeaturePlot(Fibroblast,
            features = c('CAF1','iCAF1','pCAF1','mCAF1','myCAF1','meCAF1','vCAF1','tCAF1','infCAF1','rCAF1','dCAF1',
                         'Macrophage','Neutrophils1','Plasma1','Tcells1','Fibroblasts1','Bcells1','DC1'))

ggsave(
  p,
  filename = file.path(Fibro_out_dir, str_glue('CAF score.pdf')),
  height = 20,
  width  = 15
) 


## annotation ----
cell_cluster <- c(
  '0'  = "mCAF",
  '1' = "mCAF",
  '2'  = "low-counts",
  '3'  = "Normal fibroblast", # LTBP2_CAF, https://doi.org/10.1042/CS20180435, https://doi.org/10.1158/0008-5472.CAN-24-1816
  '4'  = "mCAF",
  '5'  = "mCAF",
  '6'  = "Exocrine glandular cell", # ERP27
  '7'  = "M2 macrophage", # CD163, MRC1（CD206）
  '8'  = "Normal fibroblast", # TNXB
  '9'  = "tCAF", # MME, NDRG1
  '10' = "Exocrine glandular cell", # ERP27
  '11'  = "mCAF", 
  '12'  = "Ductal cell", # CFTR
  '13'  = "dCAF", # TOP2A, BIRC5
  '14'  = "mCAF",
  '15'  = "Normal fibroblast",
  '16'  = "mCAF",
  '17'  = "low-counts"
)

# VlnPlot(Fibroblast, features = 'nCount_RNA', group.by = "RNA_snn_res.0.8", pt.size = 0)
Idents(Fibroblast) <- "RNA_snn_res.0.8"
#Idents(Fibroblast) <- "sub.cluster"

Fibroblast <- RenameIdents(Fibroblast, cell_cluster)
Fibroblast$subtype <- Idents(Fibroblast)

pdf(file.path(Fibro_out_dir,str_glue("Dimplot dim{dims} annotation.pdf")), height = 6, width = 8)
# DimPlot(Fibroblast, group.by = "subtype", label = TRUE, reduction = 'umap',repel = TRUE) + NoLegend()
DimPlot(Fibroblast, group.by = "subtype", label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

Fibroblast$group <- factor(Fibroblast$group, levels = c('N', 'CP', 'T'))
pdf(file.path(Fibro_out_dir, str_glue("Dimplot dim{dims} annotation split.pdf")), height = 6, width = 14)
DimPlot(Fibroblast, group.by = "subtype",split.by = 'group', label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()
saveRDS(Fibroblast, file.path(Fibro_out_dir, 'Fibroblast.rds'))




# Fibroblast subtype dotplot ----------------------------------------------------
Fibroblast$subtype %>% table()
all_Fib <- subset(Fibroblast, subtype %in% c("mCAF", "tCAF", "dCAF", "Normal fibroblast"))
Fibroblast_cellTye_markers <- c(  
  NF = 'IGFBP6', 
  NF = 'PI16', 
  NF = 'DPT',
  NF = 'CXCL12',
  NF = 'GPX3', 
  NF = 'TNXB', 
  CAF = 'DCN',
  CAF = 'COL1A1',
  CAF = 'COL3A1',
  CAF = 'COL6A1',
  CAF = 'ACTA2',
  CAF = 'FAP',
  CAF = 'CD10',
  CAF = 'MMP1', 
  # CAF = 'MMP11',
  CAF = 'RUNX1',
  
  'mCAF'='MMP11',
  'mCAF'='MYL9',
  'mCAF'='POSTN',
  'mCAF'='TPM2',
  'mCAF'='TAGLN',
  'mCAF'='HOPX', 
  # 'mCAF'='DCN',
  'mCAF'='CDH11',
  'mCAF'='LUM', 
  
  'tCAF' = 'PDPN',
  'tCAF' = 'MME',
  'tCAF' = 'NDRG1',
  'tCAF' = 'VEGFA',
  
  'dCAF' = 'MKI67',
  'dCAF' = 'NUSAP1',
  'dCAF' = 'TOP2A',
  'dCAF' = 'CENPF',
  'dCAF' = 'PTTG1',
  'dCAF' = 'SCG2',
  'dCAF' = 'STMN1',
  'dCAF' = 'PLAU',
  'dCAF' = 'IGFBP2'
)

pdf(file.path(Fibro_out_dir, "dotplot celltype main_marker.pdf"), width = 28, height = 6)
DotPlot(all_Fib, group.by = "subtype", features = Fibroblast_cellTye_markers, dot.scale = 10) + theme_90angle
dev.off()
