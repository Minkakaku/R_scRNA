library(Seurat)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(BPCells)
library(clusterProfiler)
library(org.Hs.eg.db)


all_data <- readRDS('/mnt/NexData_file/home/lxy/case/XE2603262474-01/stRNA_2/output/data/all_data.rds')

# Endothelials -------------------------------------------------------------------
#all_data$annotation <- str_replace(all_data$annotation, 'Proliferating T cell', 'T cell')
all_data$annotation %>% table()
Endothelial_out_dir <- "Endothelial_subtype" 
dir.create(Endothelial_out_dir)
Endothelial <- subset(all_data, annotation == 'Endothelial')

all_data$annotation %>% table()
all_data$celltype %>% table()
## PCA --------
Endothelial <- FindVariableFeatures(Endothelial)

Endothelial <- ScaleData(Endothelial, vars.to.regress = 'nCount_RNA')
Endothelial <- RunPCA(Endothelial)

### pca inspector ------
pdf(file.path(Endothelial_out_dir, 'pca plot.pdf'), height = 6, width = 8)
ElbowPlot(Endothelial, ndims = 50)
dev.off()


## dimension reduction --------
dims <- 17 # 18

Endothelial <- RunUMAP(Endothelial, dims = 1:dims,return.model = TRUE)

pdf(file.path(Endothelial_out_dir, str_glue("Dimplot dim{dims} sample.pdf")), height = 6, width = 8)
DimPlot(Endothelial, group.by = "sample", label = TRUE) + NoLegend()
dev.off()

Endothelial <- Endothelial %>% 
  harmony::RunHarmony(
    dims.use = 1:dims,
    group.by.vars = "project",
    # theta = 4, # default 2, bigger value -- more diversity per cluster
    # sigma = 0.5,
    # lambda = 0.5, # default, smaller value -- more integrated
    reduction.save = "harmony"
  )

Endothelial <- RunUMAP(
  Endothelial, 
  dims = 1:dims, 
  reduction = "harmony", 
  reduction.name = "harmonyumap", 
  reduction.key  = "harmonyUMAP_",
  return.model = TRUE
)

pdf(file.path(Endothelial_out_dir, str_glue("Dimplot dim{dims} sample - harmony.pdf")), height = 6, width = 8)
DimPlot(Endothelial, group.by = "sample", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()

reduction.use <- 'harmony' # 'harmony' 'pca'


### cluster -------
Endothelial  <-
  FindNeighbors(
    Endothelial ,
    reduction = reduction.use,
    dims = 1:dims,
    k.param = 20
  )


Endothelial <-
  FindClusters(
    Endothelial,
    resolution = c(0.2, 0.4, 0.6, 0.8)
  )



# dimplot
pdf(file.path(Endothelial_out_dir, str_glue("Dimplot dim{dims} cluster.pdf")), height = 6, width = 8)
DimPlot(Endothelial, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Endothelial, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Endothelial, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Endothelial, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(Endothelial, group.by = "RNA_snn_res.1", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(Endothelial, group.by = "RNA_snn_res.1.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(Endothelial, group.by = "RNA_snn_res.1.5", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()

Idents(Endothelial) <- "RNA_snn_res.0.4"

EC_subtype_marker_ls <- c(
  # 通用内皮 Endothelial Cell (EC)
  "EC" = "CD34",
  "EC" = "PECAM1",
  "EC" = "PTPRB",
  "EC" = "MEIS2",
  "EC" = "FLT1",
  #"EC" = "EMCN",
  
  # 肾小球毛细血管内皮 Glomerular Capillary Endothelial Cell (EC-GC)
  "EC-GC" = "EMCN",
  "EC-GC" = "HECW2",
  "EC-GC" = "PLAT",
  "EC-GC" = "ITGA8",
  "EC-GC" = "EHD3",
  "EC-GC" = "KDR",
  "EC-GC" = "SOST",
  
  # 入/出球小动脉内皮 Afferent / Efferent Arteriole Endothelial Cell (EC-AEA)
  
  "EC-AEA" = "VEGFC",
  "EC-AEA" = "CCDC3",
  "EC-AEA" = "CDH5",
  "EC-AEA" = "SERPINE2",
  "EC-AEA" = "FBLN5",
  "EC-AEA" = "CXCL12",
  "EC-AEA" = "SOX17",
  
  "EC-AEA-DVR" = "BTNL9",
  "EC-AEA-DVR" = "ADAMTS6",
  "EC-AEA-DVR" = "PALMD",
  "EC-AEA-DVR" = "AQP1",
  "EC-AEA-DVR" = "TM4SF1",
  # 直小血管降支内皮 Descending Vasa Recta Endothelial Cell (EC-DVR)
  "EC-DVR" = "MCTP1",
  "EC-DVR" = "SLC14A1",
  "EC-DVR" = "ENPP2",
  "EC-DVR" = "LYPD6B",
  
  # 管周毛细血管内皮 Peritubular Capilary Endothelial Cell (EC-PTC)
  "EC-PTC" = "PITPNC1",
  "EC-PTC" = "GRB10",
  "EC-PTC" = "SLCO2A1",
  "EC-PTC" = "RAPGEF4",
  
  "EC-PTC-AVR" = "CEACAM1",
  "EC-PTC-AVR" = "DNASE1L3",
  "EC-PTC-AVR" = "PLVAP",
  # 直小血管升支内皮 Ascending Vasa Recta Endothelial Cell (EC-AVR)
  "EC-AVR" = "GPM6A",
  "EC-AVR" = "EDIL3",
  "EC-AVR" = "TLL1",
  "EC-AVR" = "ZNF385D",
  "EC-AVR" = "NR2F2",
  
  # 淋巴管内皮 Lymphatic Cell (EC-LYM)
  "EC-LYM" = "MMRN1",
  "EC-LYM" = "CD36",
  "EC-LYM" = "TBX1",
  "EC-LYM" = "PKHD1L1",
  "EC-LYM" = "PROX1"
)

EC_subtype_marker_ls_1 <- c(
  # 通用内皮 EC
  "EC" = "CD34",
  "EC" = "CLDN5",
  "EC" = "CLEC14A",
  "EC" = "EGFL7",
  "EC" = "FLT1",
  "EC" = "HSPG2",
  "EC" = "MCAM",
  "EC" = "PECAM1",
  "EC" = "PLVAP",
  "EC" = "PODXL",
  "EC" = "CALCRL",
  "EC" = "CDH5",
  "EC" = "CTLA2A",
  "EC" = "EMCN",
  "EC" = "ENG",
  "EC" = "ESAM",
  "EC" = "LEF1",
  "EC" = "LY6C1",
  "EC" = "PTPRB",
  
  # 动脉内皮 Arterial EC
  "Arterial EC" = "GJA5",
  "Arterial EC" = "FBLN5",
  "Arterial EC" = "HEY1",
  "Arterial EC" = "GJA4",
  "Arterial EC" = "MECOM",
  "Arterial EC" = "SOX17",
  
  # 静脉内皮 Venous EC
  "Venous EC" = "VEGFC",
  "Venous EC" = "PRSS23",
  "Venous EC" = "VWF",
  "Venous EC" = "VCAM1",
  "Venous EC" = "NR2F2",
  "Venous EC" = "CFH",
  
  # 淋巴管内皮 Lymphatic EC
  "Lymphatic EC" = "FLT4",
  "Lymphatic EC" = "CCL21",
  "Lymphatic EC" = "PDPN",
  "Lymphatic EC" = "PROX1",
  "Lymphatic EC" = "LYVE1",
  
  # 毛细血管内皮 Capillary EC
  "Capillary EC" = "PLPP3",
  "Capillary EC" = "KDR",
  "Capillary EC" = "NRP1",
  "Capillary EC" = "EDNRB",
  "Capillary EC" = "CA4",
  
  # 普通毛细血管 general capillary EC
  "general capillary EC" = "GPIHBP1",
  "general capillary EC" = "KIT",
  
  # 增殖型普通毛细血管 proliferating general capillary EC
  "proliferating general capillary EC" = "TOP2A",
  "proliferating general capillary EC" = "MKI67",
  
  
  # 毛细血管后内皮 Postcapillary
  "Postcapillary" = "ACKR1",
  "Postcapillary" = "CPE",
  
  # 未成熟内皮 Immature
  "Immature" = "APLNR",
  "Immature" = "VWA1",
  
  # 活化内皮 Activated
  "Activated" = "THY1",
  "Activated" = "ACTA2"
)

EC_subtype_marker_ls_2 <- c(
  # 肾小球毛细血管内皮 EC-GC
  "EC-GC" = "BTNL9",
  "EC-GC" = "HECW2",
  "EC-GC" = "EHD3",
  "EC-GC" = "CRHBP",
  "EC-GC" = "KDR",
  "EC-GC" = "CMTM8",
  "EC-GC" = "SYNE1",
  "EC-GC" = "SOST",
  "EC-GC" = "PLPP3",
  #"EC-GC" = "LDB2",
  "EC-GC" = "TGFBR2",
  "EC-GC" = "FCN3",
  
  # 入/出球小动脉内皮 EC-AEA
  "EC-AEA" = "CHRM3",
  "EC-AEA" = "ADAMTS6",
  "EC-AEA" = "VEGFC",
  "EC-AEA" = "ADGRL4",
  "EC-AEA" = "PLPP1",
  "EC-AEA" = "SERPINE2",
  "EC-AEA" = "ABI3BP",
  "EC-AEA" = "JAG2",
  "EC-AEA" = "SOX17",
  "EC-AEA" = "ADAM15",
  "EC-AEA" = "GLUL",
  "EC-AEA" = "LTBP4",
  "EC-AEA" = "KCTD12",
  
  # 直小血管降支内皮 EC-DVR
  "EC-DVR" = "PALMD",
  "EC-DVR" = "SLC14A1",
  #"EC-DVR" = "ADAMTS6",
  "EC-DVR" = "MCTP1",
  "EC-DVR" = "CYYR1",
  "EC-DVR" = "NOTCH4",
  "EC-DVR" = "CD34",
  "EC-DVR" = "LYPD6B",
  "EC-DVR" = "THSD4",
  "EC-DVR" = "ENPP2",
  "EC-DVR" = "RBP7",
  "EC-DVR" = "LEPR",
  "EC-DVR" = "SEMA3G",
  "EC-DVR" = "TIMP3",
  
  # 管周毛细血管内皮 EC-PTC
  "EC-PTC" = "SLCO2A1",
  "EC-PTC" = "RAPGEF4",
  "EC-PTC" = "RAMP3",
  "EC-PTC" = "PLVAP",
  "EC-PTC" = "FLT1",
  "EC-PTC" = "SH3RF3",
  "EC-PTC" = "HSPG2",
  "EC-PTC" = "PITPNC1",
  "EC-PTC" = "RGCC",
  "EC-PTC" = "CD93",
  
  # 分化型管周毛细血管内皮 dEC-PTC
  "dEC-PTC" = "EGFL7",
  #"dEC-PTC" = "FLT1",
  "dEC-PTC" = "LDB2",
  "dEC-PTC" = "FTL",
  "dEC-PTC" = "DHFR",
  
  # 直小血管升支内皮 EC-AVR
  "EC-AVR" = "GPM6A",
  #"EC-AVR" = "PLVAP",
  "EC-AVR" = "TLL1",
  "EC-AVR" = "SEMA3D",
  "EC-AVR" = "RANBP3L",
  "EC-AVR" = "MIR99AHG",
  "EC-AVR" = "PMP22",
  "EC-AVR" = "CALCRL",
  "EC-AVR" = "IGFBP4",
  "EC-AVR" = "MARCKS",
  
  # 分化内皮 dEC
  "dEC" = "IGFBP5",
  "dEC" = "B2M",
  "dEC" = "TMSB4X",
  "dEC" = "TMSB10",
  "dEC" = "PECAM1",
  "dEC" = "FAU",
  "dEC" = "RPL34",
  "dEC" = "RPL26",
  "dEC" = "RPS2",
  "dEC" = "SRP14",
  
  # 增殖内皮 cycEC
  "cycEC" = "TEK",
  "cycEC" = "APOLD1",
  "cycEC" = "PTPRB",
  #"cycEC" = "FLT1",
  "cycEC" = "PRC1",
  "cycEC" = "MXD3",
  "cycEC" = "SMC4",
  "cycEC" = "CDKN2C",
  "cycEC" = "NCAPG",
  "cycEC" = "LMNB2",
  "cycEC" = "CDH5",
  #"cycEC" = "EGFL7",
  
  # 淋巴管内皮 EC-LYM
  "EC-LYM" = "PKHD1L1",
  "EC-LYM" = "MMRN1",
  "EC-LYM" = "TBX1",
  "EC-LYM" = "SEMA3A",
  "EC-LYM" = "CD36",
  "EC-LYM" = "ZNF521",
  "EC-LYM" = "GRAPL",
  "EC-LYM" = "SMAD1",
  "EC-LYM" = "PTPRE",
  "EC-LYM" = "PGM5",
  "EC-LYM" = "DOCK5",
  "EC-LYM" = "KALRN",
  "EC-LYM" = "CCL21",
  "EC-LYM" = "TFF3",
  "EC-LYM" = "PROX1",
  "EC-LYM" = "LAMA4",
  "EC-LYM" = "NRP2",
  "EC-LYM" = "GJA1",
  "EC-LYM" = "FLT4",
  "EC-LYM" = "FABP4",
  "EC-LYM" = "TFPI",
  "EC-LYM" = "FABP5",
  "EC-LYM" = "LAPTM5",
  
  # 系膜细胞 MC（Mesangial cell）
  "MC" = "PIEZO2",
  "MC" = "ROBO1",
  "MC" = "PIP5K1B",
  "MC" = "DAPK2",
  "MC" = "ITGA11",
  "MC" = "DAAM2",
  "MC" = "PDE3A",
  "MC" = "PHTF2",
  "MC" = "PLCL1",
  "MC" = "PRKCA",
  "MC" = "IL1RL1",
  "MC" = "AGTR1",
  "MC" = "POSTN",
  "MC" = "STK17B",
  "MC" = "LMO7",
  "MC" = "SLCO3A1",
  "MC" = "ITGA8"
)


Idents(Endothelial)<-'RNA_snn_res.0.4'
Endothelial_markers <-
  FindAllMarkers(
    JoinLayers(Endothelial),
    test.use = "wilcox",
    logfc.threshold = 0,
    min.pct = 0,
    only.pos = TRUE
  )
write.csv(Endothelial_markers, file = "Endothelial_subtype/Endothelial_all_markers.csv", row.names = TRUE)

Endothelial_top5_markers <- 
  Endothelial_markers %>% 
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

names(Endothelial_top5_markers) <- sub('.*\\.','',names(Endothelial_top5_markers))


# dotplot
# 如果已做其他细胞类型的亚群分析，这里就可以顺便画一下其他类型的细胞亚群基因，比如EC_subtype_marker_ls2，否则只用画T细胞
pdf(file.path(Endothelial_out_dir, str_glue("dotplot dims{dims} major_marker.pdf")), width = 25, height = 6)
DotPlot(Endothelial, group.by = "RNA_snn_res.0.4", features = major_marker_ls)+RotatedAxis()+theme_90angle
DotPlot(Endothelial, group.by = "RNA_snn_res.0.4", features = EC_subtype_marker_ls)+RotatedAxis()+theme_90angle
DotPlot(Endothelial, group.by = "RNA_snn_res.0.4", features = EC_subtype_marker_ls_1)+RotatedAxis()+theme_90angle
DotPlot(Endothelial, group.by = "RNA_snn_res.0.4", features = EC_subtype_marker_ls_2)+RotatedAxis()+theme_90angle
#DotPlot(Endothelial, group.by = "RNA_snn_res.0.6", features = Endothelial_markers)+RotatedAxis()
DotPlot(Endothelial, group.by = "RNA_snn_res.0.4", features = Endothelial_top5_markers)+RotatedAxis()+theme_90angle
#DotPlot(Endothelial, group.by = "RNA_snn_res.0.6", features = EC_subtype_marker_ls2)+RotatedAxis()
dev.off()

pdf(file.path(Endothelial_out_dir, str_glue("dotplot dims{dims} major_marker 0.8.pdf")), width = 25, height = 6)
DotPlot(Endothelial, group.by = "RNA_snn_res.0.8", features = major_marker_ls)+RotatedAxis()+theme_90angle
DotPlot(Endothelial, group.by = "RNA_snn_res.0.8", features = EC_subtype_marker_ls)+RotatedAxis()+theme_90angle
DotPlot(Endothelial, group.by = "RNA_snn_res.0.8", features = EC_subtype_marker_ls_1)+RotatedAxis()+theme_90angle
DotPlot(Endothelial, group.by = "RNA_snn_res.0.8", features = EC_subtype_marker_ls_2)+RotatedAxis()+theme_90angle
#DotPlot(Endothelial, group.by = "RNA_snn_res.0.6", features = Endothelial_markers)+RotatedAxis()
#DotPlot(Endothelial, group.by = "RNA_snn_res.0.8", features = Endothelial_top5_markers)+RotatedAxis()+theme_90angle
#DotPlot(Endothelial, group.by = "RNA_snn_res.0.6", features = EC_subtype_marker_ls2)+RotatedAxis()
dev.off()

saveRDS(Endothelial ,file = "Endothelial_subtype/Endothelial_cluster.rds")
Endothelial <- readRDS("Endothelial_subtype/Endothelial_cluster.rds")
## annotation ----------------------------------------------------------------
cell_cluster <- c(
  '0'  = "PTC_EC", # 
  '1'  = "PTC_EC", #
  '2'  = "PTC_EC", # 
  '3'  = "GC_EC", # 
  '4'  = "Podocyte", # 
  '5'  = "Mesangial cell",#  
  '6'  = "PTS2-like_EC", #  
  '7'  = "AEA_EC", #  
  '8'  = "LOH_DCT_DC-like_EC", #  
  '9'  = "LYM_EC"
)

Idents(Endothelial) <- "RNA_snn_res.0.4"




Endothelial <- RenameIdents(Endothelial, cell_cluster)
Endothelial$subtype <- Idents(Endothelial)


pdf(file.path(Endothelial_out_dir, str_glue("Dimplot dim{dims} annotation.pdf")), height = 6, width = 8)
DimPlot(Endothelial, group.by = "subtype", label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

pdf(file.path(Endothelial_out_dir, str_glue("Dimplot dim{dims} annotation split.pdf")), height = 6, width = 15)
DimPlot(Endothelial, group.by = "subtype",split.by = 'group', label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()


keep <- c(
  "PTC_EC",
  "AEA_EC",
  "GC_EC",
  "LYM_EC"
)

Endothelial_clean <- subset(Endothelial, subtype %in% keep)

Endothelial_clean$subtype <- droplevels(Endothelial_clean$subtype)
Idents(Endothelial_clean) <- Endothelial_clean$subtype

table(Endothelial_clean$subtype)
saveRDS(Endothelial_clean, "Endothelial_subtype/Endothelial_clean.rds")

Endothelial <- Endothelial_clean

Endothelial <- RunUMAP(Endothelial, dims = 1:17,n.neighbors = 50,return.model = TRUE)

pdf(file.path(Endothelial_out_dir, str_glue("Dimplot dim{dims} Endothelial annotation.pdf")), height = 6, width = 8)
DimPlot(Endothelial, group.by = "subtype", label = TRUE, reduction = 'umap',repel = TRUE) + NoLegend()
dev.off()

pdf(file.path(Endothelial_out_dir, str_glue("Dimplot dim{dims} Endothelial annotation split.pdf")), height = 6, width = 15)
DimPlot(Endothelial, group.by = "subtype",split.by = 'group', label = TRUE, reduction = 'umap',repel = TRUE) + NoLegend()
dev.off()




# 保留：0,1,7（下面的不跑，跳到700行）
Endothelial <- subset(Endothelial, RNA_snn_res.0.4 %in% c(0,1,2,3,7,9))
Endothelial$RNA_snn_res.0.4 %>% table()


## PCA --------
Endothelial <- FindVariableFeatures(Endothelial)

Endothelial <- ScaleData(Endothelial, vars.to.regress = 'nCount_RNA')
Endothelial <- RunPCA(Endothelial)

### pca inspector ------
pdf(file.path(Endothelial_out_dir, 'pca plot_1.pdf'), height = 6, width = 8)
ElbowPlot(Endothelial, ndims = 50)
dev.off()

## dimension reduction --------
dims <- 17 # 18

Endothelial <- RunUMAP(Endothelial, dims = 1:dims,return.model = TRUE)


Endothelial <- Endothelial %>% 
  harmony::RunHarmony(
    dims.use = 1:dims,
    group.by.vars = "project",
    # theta = 4, # default 2, bigger value -- more diversity per cluster
    # sigma = 0.5,
    # lambda = 0.5, # default, smaller value -- more integrated
    reduction.save = "harmony"
  )

Endothelial <- RunUMAP(
  Endothelial, 
  dims = 1:dims, 
  reduction = "harmony", 
  reduction.name = "harmonyumap", 
  reduction.key  = "harmonyUMAP_",
  return.model = TRUE
)



reduction.use <- 'harmony' # 'harmony' 'pca'


### cluster -------
Endothelial  <-
  FindNeighbors(
    Endothelial ,
    reduction = reduction.use,
    dims = 1:dims,
    k.param = 20
  )


Endothelial <-
  FindClusters(
    Endothelial,
    resolution = c(0.2, 0.4, 0.6, 0.8)
  )



# dimplot
pdf(file.path(Endothelial_out_dir, str_glue("Dimplot dim{dims} cluster_1.pdf")), height = 6, width = 8)
DimPlot(Endothelial, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Endothelial, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Endothelial, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Endothelial, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(Endothelial, group.by = "RNA_snn_res.1", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(Endothelial, group.by = "RNA_snn_res.1.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(Endothelial, group.by = "RNA_snn_res.1.5", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()


Idents(Endothelial)<-'RNA_snn_res.0.6'
Endothelial_markers <-
  FindAllMarkers(
    JoinLayers(Endothelial),
    test.use = "wilcox",
    logfc.threshold = 0,
    min.pct = 0,
    only.pos = TRUE
  )
write.csv(Endothelial_markers, file = "Endothelial_subtype/Endothelial_all_markers_1.csv", row.names = TRUE)

Endothelial_top5_markers <- 
  Endothelial_markers %>% 
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

names(Endothelial_top5_markers) <- sub('.*\\.','',names(Endothelial_top5_markers))

#Idents(Endothelial) <- "RNA_snn_res.0.4"
pdf(file.path(Endothelial_out_dir, str_glue("dotplot dims{dims} major_marker_1.pdf")), width = 25, height = 6)
DotPlot(Endothelial, group.by = "RNA_snn_res.0.6", features = major_marker_ls)+RotatedAxis()+theme_90angle
DotPlot(Endothelial, group.by = "RNA_snn_res.0.6", features = EC_subtype_marker_ls)+RotatedAxis()+theme_90angle
DotPlot(Endothelial, group.by = "RNA_snn_res.0.6", features = EC_subtype_marker_ls_1)+RotatedAxis()+theme_90angle
#DotPlot(Endothelial, group.by = "RNA_snn_res.0.6", features = Endothelial_markers)+RotatedAxis()
DotPlot(Endothelial, group.by = "RNA_snn_res.0.6", features = EC_subtype_marker_ls_2)+RotatedAxis()+theme_90angle
DotPlot(Endothelial, group.by = "RNA_snn_res.0.6", features = Endothelial_top5_markers)+RotatedAxis()+theme_90angle
dev.off()


saveRDS(Endothelial ,file = "Endothelial_subtype/Endothelial_cluster_1.rds")
Endothelial <- readRDS("Endothelial_subtype/Endothelial_cluster_1.rds")
## annotation ----------------------------------------------------------------

cell_cluster <- c(
  '0'  = "PTC_EC", # 
  '1'  = "PTC_EC", #
  '2'  = "GC_EC", # 
  '3'  = "PTC_EC", # 
  '4'  = "AEA_EC", # 
  '5'  = "Immune-like_EC",#  
  '6'  = "AEA_EC", #  
  '7'  = "Pericyte", #  
  '8'  = "LYM_EC", #  
  '9'  = "aPTC_EC"
)

Idents(Endothelial) <- "RNA_snn_res.0.6"




Endothelial <- RenameIdents(Endothelial, cell_cluster)
Endothelial$subtype <- Idents(Endothelial)


pdf(file.path(Endothelial_out_dir, str_glue("Dimplot dim{dims} annotation_1.pdf")), height = 6, width = 8)
DimPlot(Endothelial, group.by = "subtype", label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

pdf(file.path(Endothelial_out_dir, str_glue("Dimplot dim{dims} annotation split_1.pdf")), height = 6, width = 15)
DimPlot(Endothelial, group.by = "subtype",split.by = 'group', label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

table(Endothelial$subtype)
#saveRDS(Endothelial, "Endothelial_subtype/Endothelial_clean_1.rds")

keep <- c(
  "PTC_EC",
  "aPTC_EC",
  "AEA_EC",
  "GC_EC",
  "LYM_EC"
)

Endothelial_clean <- subset(Endothelial, subtype %in% keep)

Endothelial_clean$subtype <- droplevels(Endothelial_clean$subtype)
Idents(Endothelial_clean) <- Endothelial_clean$subtype

table(Endothelial_clean$subtype)
saveRDS(Endothelial_clean, "Endothelial_subtype/Endothelial_clean_1.rds")

Endothelial <- Endothelial_clean

Endothelial <- RunUMAP(Endothelial, dims = 1:17,n.neighbors = 50,return.model = TRUE)

pdf(file.path(Endothelial_out_dir, str_glue("Dimplot dim{dims} Endothelial annotation.pdf")), height = 6, width = 8)
DimPlot(Endothelial, group.by = "subtype", label = TRUE, reduction = 'umap',repel = TRUE) + NoLegend()
dev.off()

pdf(file.path(Endothelial_out_dir, str_glue("Dimplot dim{dims} Endothelial annotation split.pdf")), height = 6, width = 15)
DimPlot(Endothelial, group.by = "subtype",split.by = 'group', label = TRUE, reduction = 'umap',repel = TRUE) + NoLegend()
dev.off()




Endothelial_1 <- readRDS("Endothelial_subtype/Endothelial_cluster.rds")
cell_cluster <- c(
  '0'  = "PTC_EC", # 
  '1'  = "PTC_EC", #
  '2'  = "PTC_EC", # 
  '3'  = "GC_EC", # 
  '4'  = "Podocyte", # 
  '5'  = "Mesangial cell",#  
  '6'  = "PTS2-like_EC", #  
  '7'  = "AEA_EC", #  
  '8'  = "LOH_DCT_DC-like_EC", #  
  '9'  = "LYM_EC"
)


Idents(Endothelial_1) <- "RNA_snn_res.0.4"
Endothelial_1 <- RenameIdents(Endothelial_1, cell_cluster)
Endothelial_1$subtype <- Idents(Endothelial_1)

#Endothelial_1$subtype2 <- as.character(Endothelial_1$subtype)
#cell_subtype <- Endothelial$subtype
#names(cell_subtype) <- colnames(Endothelial)
#Endothelial_1$subtype2[colnames(Endothelial)] <- as.character(cell_subtype)
#table(Endothelial_1$subtype2)
table(Endothelial_1$subtype)
#table(Endothelial$subtype)
saveRDS(Endothelial_1, "Endothelial_subtype/Endothelial_allannotation.rds")
######下面合并到alldata的注释（先不运行）
all_data$annotation %>% table()
#all_data <- readRDS("output/data/all_data.rds")
# 1. 先把 all_data$annotation 复制到新列 T_annotation
all_data$annotation1 <- as.character(all_data$annotation)

# 2. 创建映射：细胞名 → 精细T亚型
cell_subtype <- Endothelial_1$subtype2
names(cell_subtype) <- colnames(Endothelial_1)

# 3. 只把 "T cell" 替换成精细亚型，其他细胞保持原样
all_data$annotation1[colnames(Endothelial_1)] <- as.character(cell_subtype)

# 查看结果（正确）
table(all_data$annotation1)
table(all_data$annotation)
table(Endothelial_1$subtype2)
saveRDS(all_data, "Endothelial_subtype/all_data_T_annotation.rds")
saveRDS(Endothelial_1 ,file = "Endothelial_subtype/Endothelial_annotation.rds")



Endothelial <- readRDS("Endothelial_subtype/Endothelial_clean.rds")
table(Endothelial$subtype)
#画一个注释完的点图
EC_subtype_marker_ls_3 <- c(
  "EC" = "PECAM1",
  "EC" = "CD34",
  "EC" = "CDH5",
  
  "PTC_EC" = "CEACAM1",
  "PTC_EC" = "PLVAP",
  "PTC_EC" = "HSPG2",
  "PTC_EC" = "CD93",
  
  "GC_EC" = "KDR",
  "GC_EC" = "SOST",
  
  "AEA_EC" = "VEGFC",
  "AEA_EC" = "ADGRL4",
  "AEA_EC" = "SOX17",
  "AEA_EC" = "GJA5", #https://pubmed.ncbi.nlm.nih.gov/20530546/
  
  
  "LYM_EC" = "PROX1",
  "LYM_EC" = "FLT4"
)

pdf(file.path(Endothelial_out_dir, str_glue("dotplot Endothelial subtype major_marker.pdf")), width = 12, height = 6)
DotPlot(Endothelial, group.by = "subtype", features = EC_subtype_marker_ls_3)+RotatedAxis()+theme_90angle
dev.off()


Endothelial$subtype %>% table()
Endothelial$sample %>% table()
Endothelial$group %>% table()
## cell props barplot  ----------------------------------------------------------
Endothelial_meta <- Endothelial[[]]

Endothelial_meta$sample <- factor(Endothelial_meta$sample,levels=c(Endothelial_meta$sample %>% unique()))

p_col_proj <-
  Endothelial_meta %>%
  group_by(sample, subtype) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = sample, y = freq, fill = subtype)) + 
  geom_col() +
  scale_fill_brewer(palette = "Set2", name = "Endothelial Subtype") +
  #ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )


# export
ggsave(
  p_col_proj,
  filename = file.path("Endothelial_subtype", "cell_freq_sample.pdf"),
  height = 6,
  width  = 8
)

p_col_proj <-
  Endothelial_meta %>%
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
  filename = file.path(Endothelial_out_dir, "cell_freq_group.pdf"),
  height = 6,
  width  = 4
)

#改图
# 提取meta
Endothelial_meta <- Endothelial[[]]

group_map <- c(
  "N"       = "Normal",
  "M_igAN"  = "Mild",
  "S_igAN"  = "Severe"
)
# 新分组显示顺序（左→右：Normal, Mild, Severe）
group_levels <- c("Normal", "Mild", "Severe")
# 对应配色（保持你指定色值不变）
group_col <- c(
  "Normal"  = "#2E86AB",
  "Mild"    = "#A23B72",
  "Severe"  = "#F18F01"
)


# 确保meta里的group是因子，顺序和颜色对应
Endothelial_meta$Group <- recode(Endothelial_meta$group, !!!group_map)
Endothelial_meta$Group <- factor(Endothelial_meta$Group, levels = group_levels)


# ========== 2. 计算样本水平百分比（用于箱线图 & 均值堆叠图） ==========
# 每个样本内各细胞亚群百分比
data_long <- Endothelial_meta %>%
  group_by(sample, Group, subtype) %>%
  summarise(Count = n(), .groups = "drop_last") %>%
  group_by(sample) %>%
  mutate(Percentage = Count / sum(Count) * 100) %>%
  ungroup()

# 保存原始百分比数据
write.csv(data_long, 
          file.path("Endothelial_subtype", "Endothelial_sample_percentage.csv"), 
          row.names = FALSE)

# 按Group + CellType 计算均值（用于堆叠柱状图）
category_avg <- data_long %>%
  group_by(Group, subtype) %>%
  summarise(MeanPercentage = mean(Percentage), .groups = "drop")

# 保存均值统计数据
write.csv(category_avg, 
          file.path("Endothelial_subtype", "Endothelial_group_mean_percentage.csv"), 
          row.names = FALSE)

# ========== 堆叠柱状图 ==========
p_category_stack <- ggplot(category_avg, 
                           aes(x = Group, y = MeanPercentage, fill = subtype)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  scale_fill_brewer(palette = "Set2", name = "Endothelial Subtype") +
  labs(title = "Endothelial Subtype Composition by Group",
       x = "Group", y = "Percentage (%)") +
  coord_cartesian(ylim = c(0, 100)) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 10, face = "bold"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
  )

ggsave(file.path("Endothelial_subtype", "Endothelial_StackedBar.pdf"), 
       p_category_stack, width = 8, height = 6, dpi = 300)
print(p_category_stack)


# ==================== 按细胞类型分组的箱线图 ====================
comp_list <- list(c("Normal","Mild"), c("Normal","Severe"), c("Mild","Severe"))
p <- ggplot(data_long, aes(x = Group, y = Percentage, fill = Group)) +
  geom_boxplot(width = 0.7, alpha = 0.8, outlier.shape = 19, outlier.size = 1.5) +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.7, aes(color = Group)) +
  facet_wrap(~subtype, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = group_col) +
  scale_color_manual(values = group_col) +
  labs(title = "Endothelial - Cell Type Percentages",
       x = "Group", y = "Percentage (%)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 10, face = "bold"),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5))+
  stat_compare_means(
    comparisons = comp_list,
    method = "wilcox.test",
    label = "p.format",
    size = 3
  )

# 保存+展示
ggsave(file.path("Endothelial_subtype","Endothelial_Boxplot.pdf"), p, width = 12, height = 8, dpi = 300)
print(p)