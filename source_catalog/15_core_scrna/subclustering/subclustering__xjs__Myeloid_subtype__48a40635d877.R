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

Myeloid_out_dir<-'Myeloid_subtype' # dir.create(Myeloid_out_dir)
Myeloid<-subset(all_data, annotation == 'Myeloid cell')
Myeloid$sample %>% table()

## PCA --------
DefaultAssay(Myeloid) <- "RNA"
Myeloid <- FindVariableFeatures(Myeloid)
Myeloid <- ScaleData(Myeloid, vars.to.regress = 'nCount_RNA')
Myeloid <- RunPCA(Myeloid)


### pca inspector ------
pdf(file.path(Myeloid_out_dir, 'pca plot.pdf'), height = 6, width = 8)
ElbowPlot(Myeloid, ndims = 50)
dev.off()


## dimension reduction --------
dims <- 20

library(future)
plan(multisession, workers = 10)
options(future.globals.maxSize = 10 * 1024 ^ 3) # 10Gb

Myeloid <- RunUMAP(Myeloid, dims = 1:dims,return.model = TRUE)

pdf(file.path(Myeloid_out_dir, str_glue("Dimplot dim{dims} sample.pdf")), height = 6, width = 8)
DimPlot(Myeloid, group.by = "sample", label = TRUE) + NoLegend()
dev.off()

Myeloid <- Myeloid %>% 
  harmony::RunHarmony(
    dims.use = 1:dims,
    group.by.vars = "sample",
    reduction.save = "harmony"
    # theta = 3, # default 2, bigger value -- more diversity per cluster
    # lambda = 0.5, # default, smaller value -- more integrated
  )

Myeloid <- RunUMAP(
  Myeloid, 
  dims = 1:dims, 
  reduction = "harmony", 
  reduction.name = "harmonyumap", 
  reduction.key  = "harmonyUMAP_",
  return.model = TRUE
)

pdf(file.path(Myeloid_out_dir, str_glue("Dimplot dim{dims} sample - harmony.pdf")), height = 6, width = 8)
DimPlot(Myeloid, group.by = "sample", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()

reduction.use <- 'harmony' # 'harmony' 'pca'

### cluster -------
Myeloid  <-
  FindNeighbors(
    Myeloid ,
    reduction = reduction.use,
    dims = 1:dims,
    k.param = 10
  )

Myeloid <-
  FindClusters(
    Myeloid,
    resolution = c(
      0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.5
    )
  )

# dimplot
pdf(file.path(Myeloid_out_dir, str_glue("Dimplot dim{dims} cluster.pdf")), height = 6, width = 8)
DimPlot(Myeloid, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Myeloid, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Myeloid, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Myeloid, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Myeloid, group.by = "RNA_snn_res.1", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Myeloid, group.by = "RNA_snn_res.1.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Myeloid, group.by = "RNA_snn_res.1.5", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()


### dotplot -------
Myeloid_major_marker_ls1 <- c(
  "Monocytes" = "CD14", 
  "Monocytes" = "FCGR3A",
  
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
  
  'Netrophil' = "CSF3R",
  'Netrophil' = "G0S2",
  'Netrophil' = "RETNLG",
  
  "Myocytes" = "TNNC1", 
  "Myocytes" = "MYH7", 
  
  "Mast cell" = "CPA3",
  "Mast cell" = "CTSG",
  "Mast cell" = "TPSAB1",
  "Mast cell" = "KIT",
  
  M1 = 'NOS2',
  M1 = 'IL12',
  M1 = 'FCGR1A',
  M1 = 'FCGR1B',
  M1 = 'FCGR1C',
  M1 = 'CD80',
  M1 = 'CXCR10',
  M1 = 'IL23',
  M1 = 'CXCL9',
  M1 = 'CXCL10',
  M1 = 'CXCL11',
  M1 = 'CD86',
  M1 = 'IL1A',
  M1 = 'IL1B',
  M1 =  'IL6',
  M1 =  'TNF',
  M1 = 'MHCII',
  M1 = 'CCL5',
  M1 = 'IRF5',
  M1 = 'IRF1',
  M1 = 'CD40',
  M1 = 'KYNU',
  M1 = 'CCR7',
  M1 =  'VCAN',
  M1 =  'CD300E',
  M1 =  'INHBA',
  M1 =  'TNIP3',
  M2 = 'ARG1',
  M2 = 'IL10',
  M2 = 'CD32',
  M2 = 'CD163',
  M2 = 'FCER2',
  M2 = 'CD200R1',
  M2 = 'PDCD1LG2',
  M2 = 'CD274',
  M2 = 'MARCO',
  M2 = 'CSF1R',
  M2 = 'MRC1',
  M2 = 'IL1RN',
  M2 = 'Il1R2',
  M2 = 'IL4R',
  M2 = 'CCL4',
  M2 = 'CCL13',
  M2 = 'CCL20',
  M2 = 'CCL17',
  M2 = 'CCL18',
  M2 = 'CCL22',
  M2 = 'CCL24',
  M2 = 'LYVE1',
  M2 = 'VEGFA',
  M2 = 'VEGFB',
  M2 = 'VEGFC',
  M2 = 'VEGFD',
  M2 = 'EGF',
  M2 = 'CTSA',
  M2 = 'CTSB',
  M2 = 'CSTC',
  M2 = 'CTSD',
  M2 = 'TGFB1',
  M2 = 'CD209',
  M2 = 'CD163L1',
  M2 = 'FOLR2',
  M0 = 'CD68',
  M0 = 'C1QA',
  M0 = 'C1QB',
  M0 = 'SELENOP',
  M0 = 'LYZ',
  M0 = 'HLA-DPB1',
  M0 = 'HLA-DPA1',
  M0 = 'AIF1', 
  TAM = 'SPP1', 
  TAM = 'CQC'
)

Myeloid_major_marker_ls2 <- c(
  'Mast_KIT' = 'KIT', 
  'Mast_KIT' = 'TPSAB1', 
  'Mast_KIT' = 'CPA3', 
  
  'pDC_LILRA4' = 'LILRA4', 
  'pDC_LILRA4' = 'GZMB', 
  'pDC_LILRA4' = 'IL3RA',
  
  'cDC1_CLEC9A' = 'CLEC9A', 
  'cDC1_CLEC9A' = 'FLT3', 
  'cDC1_CLEC9A' = 'IDO1', 
  
  'cDC2_CD1C' = 'CD1C', 
  'cDC2_CD1C' = 'FCER1A', 
  'cDC2_CD1C' = 'HLA-DQA1', 
  
  'cDC3_LAMP3' = 'LAMP3', 
  'cDC3_LAMP3' = 'CCR7', 
  'cDC3_LAMP3' = 'FSCN1', 
  
  'Mono_CD14' = 'FCN1', 
  'Mono_CD14' = 'S100A9', 
  'Mono_CD14' = 'S100A8', 
  
  'Mono_CD16' = 'FCGR3A', 
  'Mono_CD16' = 'LST1', 
  'Mono_CD16' = 'LILRB2', 
  
  'Macro_INHBA' = 'INHBA', 
  'Macro_INHBA' = 'IL1RN', 
  'Macro_INHBA' = 'CCL4', 
  
  'Macro_NLRP3' = 'NLRP3', 
  'Macro_NLRP3' = 'EREG', 
  'Macro_NLRP3' = 'IL1B', 
  
  'Macro_LYVE1' = 'LYVE1', 
  'Macro_LYVE1' = 'PLTP', 
  'Macro_LYVE1' = 'SEPP1', 
  
  'Macro_C1QC' = 'C1QC', 
  'Macro_C1QC' = 'C1QA', 
  'Macro_C1QC' = 'APOE'
)

pdf(file.path(Myeloid_out_dir, str_glue('Dotplot dim{dims} major marker_multi res.pdf')), width = 15, height = 6)
DotPlot(Myeloid, features = Myeloid_major_marker_ls1, group.by = "RNA_snn_res.0.2") + theme_90angle
DotPlot(Myeloid, features = Myeloid_major_marker_ls1, group.by = "RNA_snn_res.0.4") + theme_90angle
DotPlot(Myeloid, features = Myeloid_major_marker_ls1, group.by = "RNA_snn_res.0.6") + theme_90angle
DotPlot(Myeloid, features = Myeloid_major_marker_ls1, group.by = "RNA_snn_res.0.8") + theme_90angle
DotPlot(Myeloid, features = Myeloid_major_marker_ls1, group.by = "RNA_snn_res.1") + theme_90angle
DotPlot(Myeloid, features = Myeloid_major_marker_ls1, group.by = "RNA_snn_res.1.2") + theme_90angle
DotPlot(Myeloid, features = Myeloid_major_marker_ls1, group.by = "RNA_snn_res.1.5") + theme_90angle
dev.off()

Idents(Myeloid)<-'RNA_snn_res.0.6'
Myeloid_markers <-
  FindAllMarkers(
    JoinLayers(Myeloid),
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    only.pos = TRUE
  )
Myeloid_markers %>% openxlsx::write.xlsx(file.path(Myeloid_out_dir, 'Myeloid_markers_res0.6.xlsx'))

Myeloid_top5_markers <- 
  Myeloid_markers %>% 
  group_by(cluster) %>%
  dplyr::filter(
    p_val_adj < 0.05, 
    avg_log2FC > 1, 
    pct.1-pct.2 > 0.2
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
        # Myeloidrowser()
        names(gene) = str_glue("cluster_{name}")
        gene
      })
    re
  }) %>%
  unlist() %>% 
  .[!duplicated(.)]

names(Myeloid_top5_markers) <- sub('.*\\.','',names(Myeloid_top5_markers))

# dotplot
pdf(file.path(Myeloid_out_dir, str_glue("dotplot dims{dims} major_marker.pdf")), width = 25, height = 6)
DotPlot(Myeloid, group.by = "RNA_snn_res.0.6", features = major_marker_ls)+RotatedAxis()
DotPlot(Myeloid, group.by = "RNA_snn_res.0.6", features = Myeloid_major_marker_ls)+RotatedAxis()
DotPlot(Myeloid, group.by = "RNA_snn_res.0.6", features = Myeloid_top5_markers)+RotatedAxis()
dev.off()

# Idents(Myeloid) <- "RNA_snn_res.0.6"
# Myeloid <- FindSubCluster(Myeloid, cluster = 1, graph.name = "RNA_snn", algorithm = 1,resolution = 0.3)
# Idents(Myeloid) <- "sub.cluster"
# DimPlot(Myeloid, label = TRUE, reduction = 'harmonyumap') + NoLegend()
# # DimPlot(Myeloid, label = TRUE, reduction = 'harmonyumap', split.by = "sub.cluster") + NoLegend()
# # DimPlot(Myeloid, label = TRUE, reduction = 'harmonyumap', split.by = "RNA_snn_res.0.6") + NoLegend()
# DotPlot(Myeloid, group.by = "sub.cluster", features = Myeloid_major_marker_ls2)+RotatedAxis()

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
Myeloids <- c("SFRP4","PDGFRA","CCDC80","COL1A1","ACTA2",'DCN',"FAP","TNC","PDPN","PDGFRB","COL4A2","COL5A2","STC2","LOXL2","LUM","VIM","COL1A2")
Bcells <- c("CD19","CD79A","MS4A1")
DC <- c("CD83","LILRA4","CLEC9A","IDO1","CD1C","FCER1A")

pathway<-list(CAF,iCAF,pCAF,mCAF,myCAF,meCAF,vCAF,tCAF,infCAF,rCAF,dCAF,Macrophages,Neutrophils,Plasma,Tcells,Myeloids,Bcells,DC)
names(pathway)<-c('CAF','iCAF','pCAF','mCAF','myCAF','meCAF','vCAF','tCAF','infCAF','rCAF','dCAF','Macrophages','Neutrophils','Plasma','Tcells','Myeloids','Bcells','DC')

# M1_M2<-read.table("/home/lr/ATAC_thirteen/Immune/MG/M1_M2.txt",header=TRUE,fill=TRUE)
n_na_idx <- pathway %>% map_lgl(function(x) all(is.na(x))) %>% {! .}%>%  which
pathway <- pathway[n_na_idx] %>% map(function(x) x[!is.na(x)])

for ( i in seq_along(names(pathway))){
  Myeloid <- AddModuleScore(JoinLayers(Myeloid), pathway[i], name =  names(pathway)[i],nbin=5,ctrl=10)
}

p<-VlnPlot(Myeloid,
           features = c('CAF1','iCAF1','pCAF1','mCAF1','myCAF1','meCAF1','vCAF1','tCAF1','infCAF1','rCAF1','dCAF1',
                        'Macrophage1','Neutrophils1','Plasma1','Tcells1','Myeloids1','Bcells1','DC1'),
           group.by='RNA_snn_res.0.6',
           combine = T,
           pt.size = 0)

FeaturePlot(Myeloid,
            features = c('CAF1','iCAF1','pCAF1','mCAF1','myCAF1','meCAF1','vCAF1','tCAF1','infCAF1','rCAF1','dCAF1',
                         'Macrophage','Neutrophils1','Plasma1','Tcells1','Myeloids1','Bcells1','DC1'))

ggsave(
  p,
  filename = file.path(Myeloid_out_dir, str_glue('CAF score.pdf')),
  height = 20,
  width  = 15
) 


## annotation ----
# cell_cluster <- c(
#   '0'  = "M2 macrophage", # CD163, MRC1（CD206）
#   '1' = "mCAF", # CTHRC1, POSTN
#   '2'  = "Neutropil",
#   '3'  = "M2 macrophage",
#   '4'  = "Exocrine glandular cell", # ERP27
#   '5'  = "cDC2", # FCER1A, CD1C
#   '6'  = "low-counts",
#   '7'  = "iCAF", # CXCL12, DPT1
#   '8'  = "M2 macrophage", # CCL18
#   '9'  = "Ductal cell", # CFTR
#   '10' = "M1 macrophage",# CD300E
#   '11'  = "IFN-responsive macrophage",
#   '12'  = "Proliferating macrophage",
#   '13'  = "LAMP3 DC",
#   '14'  = "Endothelial cell"
# )

cell_cluster <- c(
  '0'  = "M2 macrophage", # CD163, MRC1（CD206）
  '1' = "M2 macrophage", # CTSL
  '2'  = "mCAF", # CTHRC1, POSTN
  '3'  = "Neutropil",
  '4'  = "Exocrine glandular cell", # ERP27
  '5'  = "low-counts",
  '6'  = "cDC2", # FCER1A, CD1C
  '7'  = "iCAF", # CXCL12, DPT1 
  '8'  = "M1 macrophage",# CD300E CCR2,https://zhuanlan.zhihu.com/p/634193461 CLEC5A,10.1155/2022/9926305
  '9'  = "Ductal cell", # CFTR
  '10' = "Proliferating macrophage",
  '11'  = "IFN-responsive macrophage",
  '12' = "LAMP3 DC",
  '13' = "Endothelial cell"
)

# VlnPlot(Myeloid, features = 'nCount_RNA', group.by = "RNA_snn_res.0.6", pt.size = 0)
Idents(Myeloid) <- "RNA_snn_res.0.6"
#Idents(Myeloid) <- "sub.cluster"

Myeloid <- RenameIdents(Myeloid, cell_cluster)
Myeloid$subtype <- Idents(Myeloid)

pdf(file.path(Myeloid_out_dir,str_glue("Dimplot dim{dims} annotation.pdf")), height = 6, width = 8)
# DimPlot(Myeloid, group.by = "subtype", label = TRUE, reduction = 'umap',repel = TRUE) + NoLegend()
DimPlot(Myeloid, group.by = "subtype", label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

Myeloid$group <- factor(Myeloid$group, levels = c('N', 'CP', 'T'))
pdf(file.path(Myeloid_out_dir, str_glue("Dimplot dim{dims} annotation split.pdf")), height = 6, width = 14)
DimPlot(Myeloid, group.by = "subtype",split.by = 'group', label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()
saveRDS(Myeloid, file.path(Myeloid_out_dir, 'Myeloid.rds'))




# Myeloid subtype dotplot ----------------------------------------------------
Myeloid$subtype %>% table()
all_Mye <- subset(Myeloid, subtype %in% c('M2 macrophage', 'Neutropil', 'cDC2', 'M1 macrophage', 'Proliferating macrophage', 'IFN-responsive macrophage', 'LAMP3 DC'))
Myeloid_cellTye_markers <- c(  
  'M2 macrophage' = "CD163", 
  'M2 macrophage' = "MRC1", 
  'M2 macrophage' = "CTSL", 
  
  'Neutropil' = "CSF3R", 
  
  'DC' = "CD83", 
  
  'cDC2' = "FCER1A", 
  'cDC2' = "CD1C", 
  
  'M1 macrophage' = "CD300E", 
  'M1 macrophage' = "CCR2", 
  'M1 macrophage' = "CLEC5A", 
  
  'Proliferating macrophage' = "TOP2A", 
  'Proliferating macrophage' = "MKI67", 

  'IFN-responsive macrophage' = "IFIT3", 
  'IFN-responsive macrophage' = "MX1", 
  'IFN-responsive macrophage' = "IFIT1", 
  'IFN-responsive macrophage' = "IFIT2", 
  
  'LAMP3 DC' = "IDO1", 
  'LAMP3 DC' = "LAMP3"
)

pdf(file.path(Myeloid_out_dir, "dotplot celltype main_marker.pdf"), width = 28, height = 6)
DotPlot(all_Mye, group.by = "subtype", features = Myeloid_cellTye_markers, dot.scale = 10) + theme_90angle
dev.off()




