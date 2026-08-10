library(Seurat)
library(tidyverse)
library(BPCells)
library(clusterProfiler)
library(org.Hs.eg.db)
library(future)
#all_data <- readRDS(file.path("R_OUTPUT", "all_data.rds"))

all_data <- readRDS('/mnt/NexData_file/home/lzl/case/XE2602040886-02/output/data/adata_filtered_add_celltype.rds')

# change assay
all_data[['RNA']] <- as(all_data[['RNA']], "Assay5")
all_data@assays$RNA$data <- NULL
all_data[['RNA']] <- split(all_data[['RNA']], all_data$project, layers = "counts")

# normalize
all_data <- all_data %>% NormalizeData()

# markers: https://cloud.tencent.com/developer/article/2115987
major_marker_ls <- c(
    "Epithelial" = "EPCAM", 
    "Epithelial" = "KRT5", 
    "Epithelial" = "KRT14", 
    "Epithelial" = "KRT20", 

    "Endothelial" = "PECAM1", 
    "Endothelial" = "VWF", 
    "Endothelial" = "AQP1", 
    "Endothelial" = "CDH5", 
    "Endothelial" = "PLVAP", 

    "Fibroblast" = "PDGFRA", 
    "Fibroblast" = "CCDC80", 
    "Fibroblast" = "COL1A1", 
    "Fibroblast" = "COL1A2", 
    "Fibroblast" = "DCN", 
    "Fibroblast" = "ACTA2", 

    "B cell" = "CD19", 
    "B cell" = "CD79A", 
    "B cell" = "MS4A1", 
    "Plasma" = "TNFRSF17", 
    "Plasma" = "JCHAIN", 
    "Plasma" = "IGJ", 
    "Plasma" = "MZB1", 
    "Plasma" = "XBP1", 

    "T cell" = "CD3D", 
    "T cell" = "CD3E", 
    "T cell" = "CD3G", 
    "T cell" = "TRAC", 
    "T cell" = "BCL11B", 

    "NK" = "GNLY", 
    "NK" = "NKG7", 
    "NK" = "NCAM1", 
    "NK" = "UNC13D", 

    "Proliferating" = "MKI67", 
    "Proliferating" = "TOP2A", 

    "Pericyte" = "PDGFRB", 
    "Pericyte" = "RGS5", 
    "Pericyte" = "MCAM", 
    "Pericyte" = "PLVAP", 

    "Monocyte" = "CD14", 
    "Monocyte" = "FCGR3A", 
    "Monocyte" = "ITGAM", 

    "Macrophage" = "CD68", 
    "Macrophage" = "CD163", 
    "Macrophage" = "C1QA", 
    "Macrophage" = "C1QC", 
    "Macrophage" = "CSF1R", 
    "Macrophage" = "IL1B", 

    "Dendritic cell" = "CD83", 
    "Dendritic cell" = "LILRA4", 
    "Dendritic cell" = "CLEC9A", 
    "Dendritic cell" = "IDO1", 
    "Dendritic cell" = "CD1C", 
    "Dendritic cell" = "FCER1A", 
    "Dendritic cell" = "CST3", 
    "pDC" = "IRF4", 
    "pDC" = "IRF7", 
    "pDC" = "IRF8", 
    "Lamp3_DC" = "LAMP3", 

    "Netrophil" = "CSF3R", 
    "Netrophil" = "FCGR3B", 
    "Netrophil" = "S100A8", 
    "Netrophil" = "S100A9", 
    "Netrophil" = "G0S2", 
    "Netrophil" = "RETNLG", 

    "Mast cell" = "MS4A2", 
    "Mast cell" = "KIT", 
    "Mast cell" = "HDC", 
    "Mast cell" = "TPSAB1", 
    "Mast cell" = "CTSG", 
    "Mast cell" = "CPA3", 

    "Smooth muscle cell" = "PLN", 
    "Smooth muscle cell" = "MYLK", 
    "Smooth muscle cell" = "MYH11", 
    "Smooth muscle cell" = "ACTG2", 
    "Smooth muscle cell" = "TAGLN", 

    "Schwann cell" = "MPZ", 
    "Schwann cell" = "SOX10", 

    "Neuron" = "TUBB3", 
    "Neuron" = "SNAP25", 

    "Mesothelial" = "MSLN", 

    "Erythrocytes" = "HBA-A1", 
    "Erythrocytes" = "HBA-A2"
)

theme_90angle <-
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
    strip.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )




# Macrophages -------------------------------------------------------------------
all_data$annotation <- str_replace(all_data$annotation, 'Monocytes', 'Macrophages')

VEC_out_dir <- "VEC_subtype"  
dir.create(VEC_out_dir)
VEC_data <- subset(all_data, annotation == 'Vascular Endothelial Cell(VEC)')

VEC_data <- JoinLayers(VEC_data)
## PCA --------
VEC_data <- FindVariableFeatures(VEC_data)
length(VariableFeatures(VEC_data)) #2000
VEC_data <- ScaleData(VEC_data)
VEC_data <- RunPCA(VEC_data)

### pca inspector ------
pdf(file.path(VEC_out_dir, 'pca plot.pdf'), height = 6, width = 8)
ElbowPlot(VEC_data, ndims = 50)
dev.off()

pdf(file.path(VEC_out_dir, 'dimheatmap 1-5.pdf'), height = 6, width = 10)
DimHeatmap(VEC_data, dims = 1:5)
dev.off()
pdf(file.path(VEC_out_dir, 'dimheatmap 6-10.pdf'), height = 6, width = 10)
DimHeatmap(VEC_data, dims = 6:10)
dev.off()
pdf(file.path(VEC_out_dir, 'dimheatmap 11-15.pdf'), height = 6, width = 10)
DimHeatmap(VEC_data, dims = 11:15)
dev.off()
pdf(file.path(VEC_out_dir, 'dimheatmap 16-20.pdf'), height = 6, width = 10)
DimHeatmap(VEC_data, dims = 16:20)
dev.off()

## dimension reduction --------
dims <- 20 # 18

library(future)
options(future.globals.maxSize = 10 * 1024 ^ 3) # 10Gb
plan(multisession, workers = 10)

VEC_data <- RunUMAP(VEC_data, dims = 1:dims,return.model = TRUE)

pdf(file.path(VEC_out_dir, str_glue("Dimplot dim{dims} sample.pdf")), height = 6, width = 8)
DimPlot(VEC_data, group.by = "sample", label = TRUE) + NoLegend()
dev.off()

VEC_data <- VEC_data %>% 
  harmony::RunHarmony(
    dims.use = 1:dims,
    group.by.vars = "project",
    # theta = 4, # default 2, bigger value -- more diversity per cluster
    # sigma = 0.5,
    # lambda = 0.5, # default, smaller value -- more integrated
    reduction.save = "harmony"
  )

VEC_data <- RunUMAP(
  VEC_data, 
  dims = 1:dims, 
  reduction = "harmony", 
  reduction.name = "harmonyumap", 
  reduction.key  = "harmonyUMAP_",
  return.model = TRUE
)

pdf(file.path(VEC_out_dir, str_glue("Dimplot dim{dims} sample - harmony.pdf")), height = 6, width = 8)
DimPlot(VEC_data, group.by = "sample", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()

reduction.use <- 'harmony' # 'harmony' 'pca'


### cluster -------
VEC_data  <-
  FindNeighbors(
    VEC_data ,
    reduction = reduction.use,
    dims = 1:dims,
    k.param = 20
  )


VEC_data <-
  FindClusters(
    VEC_data,
    resolution = c(0.2, 0.4, 0.6, 0.8)
  )



# dimplot
pdf(file.path(VEC_out_dir, str_glue("Dimplot dim{dims} cluster.pdf")), height = 6, width = 8)
DimPlot(VEC_data, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(VEC_data, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(VEC_data, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(VEC_data, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(Macrophage, group.by = "RNA_snn_res.1", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(Macrophage, group.by = "RNA_snn_res.1.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(Macrophage, group.by = "RNA_snn_res.1.5", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()

Idents(VEC_data) <- "RNA_snn_res.0.4"


### dotplot -------
VEC_major_marker_ls <- c(
  Aorta_coronary_ec = 'BGN',
  Aorta_coronary_ec = 'ELN',
  Aorta_coronary_ec = 'SULF1',
  Aorta_coronary_ec = 'OGN',
  Aorta_coronary_ec = 'OMD',
  Art_ec_1 = 'BGN',
  Art_ec_1 = 'ELN',
  Art_ec_1 = 'SULF1',
  Art_ec_1 = 'HEY1',
  Art_ec_1 = 'DLL4',
  Art_ec_2 = 'HEY1',
  Art_ec_2 = 'DLL4',
  Art_ec_2 = 'NEBL',
  Art_ec_2 = 'GPIHBP1',
  Ven_ec_1 =  'ACKR1',
  Ven_ec_1 =  'POSTN',
  Ven_ec_1 = 'PLVAP',
  Ven_ec_2 = 'ICAM4',
  Ven_ec_2 = 'SELE',
  Lymphatic_ECs = 'PROX1',
  Lymphatic_ECs = 'LYVE1',
  Lymphatic_ECs = 'PDPN'
)

VEC_major_marker_ls <- intersect(VEC_major_marker_ls,rownames(VEC_data))

VEC_major_marker_ls <- c(
  Aorta_coronary_ec = 'ELN',
  Aorta_coronary_ec = 'SULF1',
  Art_ec_1 = 'HEY1',
  Art_ec_1 = 'DLL4',
  Art_ec_2 = 'NEBL',
  Ven_ec_1 =  'POSTN',
  Ven_ec_1 = 'PLVAP',
  Ven_ec_2 = 'ICAM4',
  Ven_ec_2 = 'NR2F2',
  Lymphatic_ECs = 'PROX1',
  Lymphatic_ECs = 'PDPN'
)



Idents(VEC_data)<-'RNA_snn_res.0.4'
VEC_data_markers <-
  FindAllMarkers(
    JoinLayers(VEC_data),
    test.use = "wilcox",
    logfc.threshold = 0,
    min.pct = 0,
    only.pos = TRUE
  )

VEC_data_top5_markers <- 
  VEC_data_markers %>% 
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

names(VEC_data_top5_markers) <- sub('.*\\.','',names(VEC_data_top5_markers))

# dotplot
pdf(file.path(Macrophage_out_dir, str_glue("dotplot dims{dims} major_marker.pdf")), width = 25, height = 6)

table(VEC_data$RNA_snn_res.0.4)
levels(factor(VEC_data$RNA_snn_res.0.4))
DotPlot(VEC_data, group.by = "RNA_snn_res.0.4", features = VEC_major_marker_ls)+RotatedAxis()
DotPlot(VEC_data, group.by = "RNA_snn_res.0.6", features = macrophage_major_marker_ls)+RotatedAxis()
DotPlot(VEC_data, group.by = "RNA_snn_res.0.6", features = CD4_major_marker_ls)+RotatedAxis()
DotPlot(VEC_data, group.by = "RNA_snn_res.0.6", features = macrophage_top5_markers)+RotatedAxis()
DotPlot(VEC_data, group.by = "RNA_snn_res.0.6", features = Fibro_major_marker_ls2)+RotatedAxis()
dev.off()


###score--------------------
Aorta_coronary_ec<-c("BGN","ELN","SULF1","OGN","OMD")
Art_ec_1<-c("BGN","ELN","SULF1","HEY1","DLL4")
Art_ec_2<-c("HEY1","DLL4","NEBL","GPIHBP1")
Ven_ec_1 <- c("ACKR1","POSTN","PLVAP")
Ven_ec_2 <- c("ICAM4","SELE")
Lymphatic_ECs <- c("PROX1","LYVE1","PDPN")


pathway<-list(Aorta_coronary_ec,Art_ec_1,Art_ec_2,Ven_ec_1,Ven_ec_2,Lymphatic_ECs)
names(pathway)<-c('Aorta_coronary_ec','Art_ec_1','Art_ec_2','Ven_ec_1','Ven_ec_2','Lymphatic_ECs')

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
  filename = 'Macrophage_subtype/M1_M2 score.pdf',
  height = 10,
  width  = 15
) 


## annotation ----------------------------------------------------------------
cell_cluster <- c(
  '0'  = "vCAFs", # POSTN, PDGFRA, COL18A1
  '1'  = "M2 macrophages", #
  '2'  = "M2 macrophages",
  '3'  = "Epithelials", # https://zhuanlan.zhihu.com/p/549337063
  '4'  = "Epithelial-like macrophages", # PGA3
  '5'  = "M2 macrophages", # MMP9,10.1007/s11033-022-07184-9, MMP12,https://doi.org/10.1016/j.jhazmat.2023.132733
  '6'  = "cDC2s", # 
  '7'  = "Proliferating macrophages",
  '8'  = "LAMP3 DCs",
  '9'  = "IFN-responsive macrophages", # 10.1002/jcsm.12584
  '10' = "Epithelials", # 
  '11' = "cDC1s"
)

Idents(Macrophage) <- "RNA_snn_res.0.6"

Macrophage <- RenameIdents(Macrophage, cell_cluster)
Macrophage$subtype <- Idents(Macrophage)


pdf(file.path(Macrophage_out_dir, str_glue("Dimplot dim{dims} annotation.pdf")), height = 6, width = 8)
DimPlot(Macrophage, group.by = "subtype", label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

pdf(file.path(Macrophage_out_dir, str_glue("Dimplot dim{dims} annotation split.pdf")), height = 6, width = 15)
DimPlot(Macrophage, group.by = "subtype",split.by = 'group', label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

Macrophage_cellTye_markers <- c( 
  "Macrophage" = "CD68",
  # "Macrophage" = "CD163",
  "Macrophage" = "C1QA",
  "Macrophage" = "C1QC", 
  # 'Macrophage'  = 'CSF1R',
  
  "Fibroblasts" = "PDGFRA", 
  "Fibroblasts" = "CCDC80",
  "Fibroblasts" = "POSTN",
  # "Fibroblasts" = "SFRP4",
  "Fibroblasts" = "COL1A1",
  "Fibroblasts" = "DCN",
  "Fibroblasts" = "ACTA2",
  
  'vCAFs' = "MCAM", 
  'vCAFs' = "NOTCH3", 
  'vCAFs' = "COL18A1",
  
  'M2 macrophages' = "CD163",
  'M2 macrophages' = "CSF1R",
  'M2 macrophages' = "MRC1",
  'M2 macrophages' = "CCL24",
  'M2 macrophages' = "VEGFB",
  'M2 macrophages' = "CTSA",
  'M2 macrophages' = "CD209",
  'M2 macrophages' = "FOLR2",
  
  'Epithelials' = "EPCAM",
  
  'Epithelial-like macrophages' = 'PGA3',

  'cDC2s' = 'CD1C',
  'cDC2s' = 'CLEC10A',
  
  'Proliferating' = 'MKI67',
  'Proliferating' = 'TOP2A',

  'LAMP3 DCs' = "LAMP3",
  'LAMP3 DCs' = 'CD83', 
  
  'IFN-responsive macrophages' = "IFIT3", 
  'IFN-responsive macrophages' = "IFIT2",
  'IFN-responsive macrophages' = "IFIT1",

  'cDC1s' = 'CLEC9A',
  'cDC1s' = 'IDO1'
)

pdf(file.path(Macrophage_out_dir, "dotplot celltype main_marker.pdf"), width = 20, height = 6)
DotPlot(Macrophage, group.by = "subtype", features = Macrophage_cellTye_markers, dot.scale = 10) + theme_90angle
dev.off()


## cell props barplot  ----------------------------------------------------------
Macrophage_meta <- Macrophage[[]]

Macrophage_meta$sample <- factor(Macrophage_meta$sample,levels=c(Macrophage_meta$sample %>% unique()))

p_col_proj <-
  Macrophage_meta %>%
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
  filename = file.path("Macrophage_subtype", "cell_freq_sample.pdf"),
  height = 6,
  width  = 8
)

p_col_proj <-
  Macrophage_meta %>%
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
  filename = file.path(Macrophage_out_dir, "cell_freq_group.pdf"),
  height = 6,
  width  = 4
)


## vocanol -------------------
Macrophage_markers <-
  FindAllMarkers(
    JoinLayers(Macrophage),
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    only.pos = TRUE
  )
write.csv(Macrophage_markers,'Macrophage_subtype/Macrophage_diff_gene.csv')

library(EnhancedVolcano)
type<-unique(Macrophage_markers$cluster)
volcano_ls <- type %>% 
  setNames(.,.) %>%
  imap(function(type,name){
    markers<- Macrophage_markers[Macrophage_markers$cluster == type, ]
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
      filename = file.path(Macrophage_out_dir, str_glue("{name} diffgene valcano.pdf")),
      height = 6,
      width = 8
    )
  })


Macrophage_markers_major<-list(Macrophage = Macrophage_markers)
# scatter
p_volcanoplot_major_ls <- 
  Macrophage_markers_major %>%
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
      filename = file.path(Macrophage_out_dir, str_glue("{name} diffgene dot.pdf")),
      height = 10,
      width = 25
    )
  })


## enrichment ---------------------------------------
Macrophage_de_genes_for_enrichmenet <- 
  Macrophage_markers %>% 
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  dplyr::filter(pct.1 > 0.15) %>%
  dplyr::filter(abs(avg_log2FC) > 0.75) %>% 
  group_by(cluster) %>%
  nest() %>%
  deframe() %>%
  map("gene")

Macrophage_er_output<-'Macrophage_subtype/enrichment' # dir.create(Macrophage_er_output)
enrich_ls <-
  Macrophage_de_genes_for_enrichmenet %>%
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


# expot enrichment table
enrich_ls %>%
  # discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    dir.create(file.path(Macrophage_er_output, name))
    
    write.csv(
      x@result,
      file = file.path(Macrophage_er_output, name, "Go enriment table.csv"),
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
      filename = file.path(Macrophage_er_output, name, "go dotplot.pdf"),
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
      filename = file.path(Macrophage_er_output, name, "go barplot.pdf"),
      height = 6,
      width = 8
    )
  })


Macrophage_de_genes_id_for_enrichmenet <-
  Macrophage_de_genes_for_enrichmenet %>%
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
  Macrophage_de_genes_id_for_enrichmenet %>%
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

# saveRDS(Macrophage_er_output, 'enrich_ls.rds')

# expot enrichment table
enrich_ls %>%
  # discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    dir.create(file.path(Macrophage_er_output, name))
    
    write.csv(
      x@result,
      file = file.path(Macrophage_er_output, name, "KEGG enriment table.csv"),
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
      filename = file.path(Macrophage_er_output, name, "kegg dotplot.pdf"),
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
      filename = file.path(Macrophage_er_output, name, "kegg barplot.pdf"),
      height = 6,
      width = 8
    )
  })

# saveRDS(Macrophage, 'Macrophage_subtype/all_Macrophage.rds')



