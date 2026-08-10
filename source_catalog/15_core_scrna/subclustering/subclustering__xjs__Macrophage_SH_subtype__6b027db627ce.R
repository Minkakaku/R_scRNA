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

CD8_major_marker_ls <- c(  
  'CD8'='CD8A',
  'CD8'='CD8B',
  'CD8_Tex'='HAVCR2',
  'CD8_Tex'='CTLA4',
  'CD8_Tex'='LAYN',
  'CD8_Tex'='LAG3',
  'CD8_Tex'='TIGIT',
  'CD8_Tex'='ENPD1',
  'CD8_Tem'='GZMK',
  'CD8_Tem'='GZMA',
  'CD8_Tem'='IFNG',
  'CD8_Tem'='KLRG1',
  'CD8_Tem'='GZMB',
  'CD8_Tem'='PRF1',
  'CD8_Tem'='NKG7',
  'CD8_Tnaive'='SELL',
  'CD8_Tnaive'='LEF1',
  'CD8_Tnaive'='TCF7',
  'CD8_Tnaive'='CCR7',
  'CD8_MAIT'='NCR3',
  'CD8_MAIT'='SLC4A10',
  'CD8_MAIT'='KLRB1',
  'CD8_Trm'='MYADM',
  'CD8_Trm'='CD69',
  'CD8_Trm'='ITGAE',
  'CD8_Trm'='CXCR6',
  'CD8_Trm'='ITGA1',
  'CD8_Tcm'='IL7R',
  'CD8_Tcm'='CD28',
  'CD8_Tcm'='GPR183'
)

Epi_major_marker_ls <- c(
  'Tumor cell' = 'MIAT', 
  'Tumor cell' = 'DSG3', 
  'Tumor cell' = 'MUC4', 
  
  "Anti-tumor cell" = 'FAM107A', 
  "Anti-tumor cell" = 'HSD17B2', 
  "Anti-tumor cell" = 'SERPINA3', 
  
  "Mucous cell" = 'MUC5B', 
  "Basal cell" = 'TP63', 
  "AT1 cell" = 'AGER', 
  
  # 'Mix' = 'POSTN', 
  # 'Mix' = 'EPCAM', 
  
  "Ciliated cell" = 'DLEC1', 
  "Ciliated cell" = 'DNAI1'
)

theme_90angle <-
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
    strip.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )

all_data <- readRDS(file.path("R_OUTPUT", "all_data.rds"))
all_data$sample %>% table()
all_data$annotation %>% table()

# Macrophages -------------------------------------------------------------------
Macrophage_out_dir <- "Macrophage_SH_subtype" 
dir.create(Macrophage_out_dir)
table(all_data$group, all_data$annotation)
Macrophage <- subset(all_data, annotation == 'Macrophage' & group == 'SH')

## PCA --------
Macrophage <- FindVariableFeatures(Macrophage)
Macrophage <- ScaleData(Macrophage, vars.to.regress = 'nCount_RNA')
Macrophage <- RunPCA(Macrophage)

### pca inspector ------
pdf(file.path(Macrophage_out_dir, 'pca plot.pdf'), height = 6, width = 8)
ElbowPlot(Macrophage, ndims = 50)
dev.off()

## dimension reduction --------
dims <- 20

library(future)
options(future.globals.maxSize = 10 * 1024 ^ 3) # 10Gb
plan(multisession, workers = 10)

Macrophage <- RunUMAP(Macrophage, dims = 1:dims,return.model = TRUE)

pdf(file.path(Macrophage_out_dir, str_glue("Dimplot dim{dims} sample.pdf")), height = 6, width = 8)
DimPlot(Macrophage, group.by = "sample", label = TRUE) + NoLegend()
dev.off()

Macrophage <- Macrophage %>% 
  harmony::RunHarmony(
    dims.use = 1:dims,
    group.by.vars = "sample",
    # theta = 4, # default 2, bigger value -- more diversity per cluster
    # sigma = 0.5,
    # lambda = 0.5, # default, smaller value -- more integrated
    reduction.save = "harmony"
  )

Macrophage <- RunUMAP(
  Macrophage, 
  dims = 1:dims, 
  reduction = "harmony", 
  reduction.name = "harmonyumap", 
  reduction.key  = "harmonyUMAP_",
  return.model = TRUE
)

pdf(file.path(Macrophage_out_dir, str_glue("Dimplot dim{dims} sample - harmony.pdf")), height = 6, width = 8)
DimPlot(Macrophage, group.by = "sample", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()

reduction.use <- 'harmony' # 'harmony' 'pca'


### cluster -------
Macrophage  <-
  FindNeighbors(
    Macrophage ,
    reduction = reduction.use,
    dims = 1:dims,
    k.param = 20
  )


Macrophage <- 
  FindClusters(
    Macrophage,
    resolution = c(0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.5)
  )



# dimplot
pdf(file.path(Macrophage_out_dir, str_glue("Dimplot dim{dims} cluster.pdf")), height = 6, width = 8)
DimPlot(Macrophage, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Macrophage, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Macrophage, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Macrophage, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Macrophage, group.by = "RNA_snn_res.1", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Macrophage, group.by = "RNA_snn_res.1.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Macrophage, group.by = "RNA_snn_res.1.5", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()

Idents(Macrophage) <- "RNA_snn_res.1.2"


### dotplot -------
macrophage_major_marker_ls <- c(
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


Idents(Macrophage)<-'RNA_snn_res.1.2'
macrophage_markers <-
  FindAllMarkers(
    JoinLayers(Macrophage),
    test.use = "wilcox",
    logfc.threshold = 0,
    min.pct = 0,
    only.pos = TRUE
  )

macrophage_top5_markers <- 
  macrophage_markers %>% 
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
        # browser()
        names(gene) = str_glue("cluster_{name}")
        gene
      })
    re
  }) %>%
  unlist() %>% 
  .[!duplicated(.)]

names(macrophage_top5_markers) <- sub('.*\\.','',names(macrophage_top5_markers))

# dotplot
pdf(file.path(Macrophage_out_dir, str_glue("dotplot dims{dims} major_marker.pdf")), width = 25, height = 6)
DotPlot(Macrophage, group.by = "RNA_snn_res.1.2", features = major_marker_ls)+RotatedAxis()
DotPlot(Macrophage, group.by = "RNA_snn_res.1.2", features = macrophage_major_marker_ls)+RotatedAxis()
DotPlot(Macrophage, group.by = "RNA_snn_res.1.2", features = CD8_major_marker_ls)+RotatedAxis()
DotPlot(Macrophage, group.by = "RNA_snn_res.1.2", features = Epi_major_marker_ls)+RotatedAxis()
DotPlot(Macrophage, group.by = "RNA_snn_res.1.2", features = macrophage_top5_markers)+RotatedAxis()
dev.off()


###score--------------------
M1<-c( 'NOS2',
       'IL12',
       'FCGR1A',
       'FCGR1B',
       'FCGR1C',
       'CD80',
       'CXCR10',
       'IL23',
       'CXCL9',
       'CXCL10',
       'CXCL11',
       'CD86',
       'IL1A',
       'IL1B',
       'IL6',
       'TNF',
       'MHCII',
       'CCL5',
       'IRF5',
       'IRF1',
       'CD40',
       'KYNU',
       'CCR7',
       'VCAN',
       'CD300E',
       'INHBA',
       'TNIP3')
M2<-c( 'ARG1',
       'IL10',
       'CD32',
       'CD163',
       'FCER2',
       'CD200R1',
       'PDCD1LG2',
       'CD274',
       'MARCO',
       'CSF1R',
       'MRC1',
       'IL1RN',
       'Il1R2',
       'IL4R',
       'CCL4',
       'CCL13',
       'CCL20',
       'CCL17',
       'CCL18',
       'CCL22',
       'CCL24',
       'LYVE1',
       'VEGFA',
       'VEGFB',
       'VEGFC',
       'VEGFD',
       'EGF',
       'CTSA',
       'CTSB',
       'CSTC',
       'CTSD',
       'TGFB1',
       'CD209',
       'CD163L1',
       'FOLR2')
Neutrophils <- c("FCGR3B","IFITM2","CSF3R")
Plasma <- c("TNFRSF17","JCHAIN","MZB1","XBP1","IGHG1","IGHA1","IGHG3","IGKC")
Tcells <- c("CD3D","CD3E","CD3G","LCK","TRAC","CD4","CD8A","GZMK")
Fibroblasts <- c("SFRP4","PDGFRA","CCDC80","COL1A1","ACTA2",'DCN',"FAP","TNC","PDPN","PDGFRB","COL4A2","COL5A2","STC2","LOXL2","LUM","VIM","COL1A2")
Bcells <- c("CD19","CD79A","MS4A1")
DC <- c("CD83","LILRA4","CLEC9A","IDO1","CD1C","FCER1A")

pathway<-list(M1,M2,Neutrophils,Plasma,Tcells,Fibroblasts,Bcells,DC)
names(pathway)<-c('M1','M2','Neutrophils','Plasma','Tcells','Fibroblasts','Bcells','DC')

# M1_M2<-read.table("/home/lr/ATAC_thirteen/Immune/MG/M1_M2.txt",header=TRUE,fill=TRUE)
n_na_idx <- pathway %>% map_lgl(function(x) all(is.na(x))) %>% {! .}%>%  which
pathway <- pathway[n_na_idx] %>% map(function(x) x[!is.na(x)])

for ( i in seq_along(names(pathway))){
  Macrophage <- AddModuleScore(JoinLayers(Macrophage), pathway[i], name =  names(pathway)[i],nbin=5,ctrl=10)
}

p<-VlnPlot(Macrophage,
           features = c('M11','M21','Neutrophils1','Plasma1','Tcells1','Fibroblasts1','Bcells1','DC1'),
           group.by='RNA_snn_res.0.6',
           pt.size = 0)

FeaturePlot(Macrophage,
            features = c('M11','M21','Neutrophils1','Plasma1','Tcells1','Fibroblasts1','Bcells1','DC1'))

ggsave(
  p,
  filename = file.path(Macrophage_out_dir, '/M1_M2 score.pdf'),
  height = 10,
  width  = 15
) 




## annotation ----------------------------------------------------------------
VlnPlot(Macrophage, features = 'nCount_RNA', pt.size = 0)
# cell_cluster <- c(
#   '0'  = "M2 macrophage", 
#   '1'  = "M2 macrophage",  
#   '2'  = "mix", # T
#   '3'  = "low-nCounts", 
#   '4'  = "mix", # DC
#   '5'  = "Tumor cell", # MIAT
#   '6'  = "Tumor cell", # MIAT
#   '7'  = "mix", # Fibro
#   '8'  = "M2 macrophage", 
#   '9'  = "M1 macrophage", 
#   '10' = "low-nCounts", 
#   '11' = "M2 macrophage", 
#   '12'  = "mix", # Tumor
#   '13'  = "mix", # Fibro
#   '14'  = "mix", # Endo
#   '15'  = "LAMP3 DC", # CD83, LAMP3 
#   '16'  = "Mast cell", # KIT 
#   '17'  = "mix", # T
#   '18'  = "Tumor cell" # EEF1A2, https://www.nature.com/articles/s41416-020-01250-4
# )

cell_cluster <- c(
  '0' = "",
  '1' = "",
  '2' = "",
  '3' = "",
  '4' = "",
  '5' = "",
  '6' = "",
  '7' = "",
  '8' = "",
  '9' = "",
  '10' = "",
  '11' = "",
  '12' = "",
  '13' = "",
  '14' = "",
  '15' = "",
  '16' = "",
  '17' = "", 
  '18' = "",
  '19' = "", 
  '20' = "",
  '21' = "",
  '22' = "",
  '23' = ""
)
 
# cell_cluster <- c(
#   '0' = "",
#   '1' = "",
#   '2' = "",
#   '3' = "",
#   '4' = "",
#   '5' = "",
#   '6' = "",
#   '7' = "",
#   '8' = "",
#   '9' = "",
#   '10' = "",
#   '11' = "",
#   '12' = "",
#   '13' = "",
#   '14' = "",
#   '15' = "",
#   '16' = "",
#   '17' = ""
# )
Idents(Macrophage) <- "RNA_snn_res.1.2"

Macrophage <- RenameIdents(Macrophage, cell_cluster)
Macrophage$subtype <- Idents(Macrophage)


pdf(file.path(Macrophage_out_dir, str_glue("Dimplot dim{dims} annotation.pdf")), height = 6, width = 8)
DimPlot(Macrophage, group.by = "subtype", label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

pdf(file.path(Macrophage_out_dir, str_glue("Dimplot dim{dims} annotation split by treat_result.pdf")), height = 6, width = 15)
DimPlot(Macrophage, group.by = "subtype",split.by = 'treat_result', label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

pdf(file.path(Macrophage_out_dir, str_glue("Dimplot dim{dims} annotation - pure.pdf")), height = 6, width = 8)
DimPlot(subset(Macrophage, subtype %in% c("M2 macrophage", "M1 macrophage", "LAMP3 DC", "Mast cell", "Tumor cell")), group.by = "subtype", label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

pdf(file.path(Macrophage_out_dir, str_glue("Dimplot dim{dims} annotation split by treat_result - pure.pdf")), height = 6, width = 15)
DimPlot(subset(Macrophage, subtype %in% c("M2 macrophage", "M1 macrophage", "LAMP3 DC", "Mast cell", "Tumor cell")), group.by = "subtype",split.by = 'treat_result', label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()


Macrophage_cellTye_markers <- c( 
  "Macrophage" = "CD68",
 
  "M2 macrophage" = "MRC1",
  "M2 macrophage" = "MARCO",

  # "mix" = "POSTN",
  # "mix" = "CD3E",
  # "mix" = "CD83",
  # "mix" = "PECAM1",
  
  "Tumor cell" = "MIAT",
  
  "M1 macrophage" = "CXCL9",
  "M1 macrophage" = "CXCL10", 

  "Mast cell" = "KIT"
)

pdf(file.path(Macrophage_out_dir, "dotplot celltype main_marker.pdf"), width = 16, height = 6)
DotPlot(Macrophage, group.by = "subtype", features = Macrophage_cellTye_markers, dot.scale = 10) + theme_90angle
dev.off()

pdf(file.path(Macrophage_out_dir, "dotplot celltype main_marker - pure.pdf"), width = 10, height = 5)
DotPlot(subset(Macrophage, subtype %in% c("M2 macrophage", "M1 macrophage", "LAMP3 DC", "Mast cell", "Tumor cell")), group.by = "subtype", features = Macrophage_cellTye_markers, dot.scale = 10) + theme_90angle
dev.off()

saveRDS(Macrophage, file.path(Macrophage_out_dir, 'Macrophage.rds'))
