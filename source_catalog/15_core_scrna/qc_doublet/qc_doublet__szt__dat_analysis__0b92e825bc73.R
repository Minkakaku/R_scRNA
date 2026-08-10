library(Seurat)
library(tidyverse)

file_dir <- "data"

all_files_dir <- list.files(file_dir, pattern = "^[A-Za-z0-9_]*$")
all_files <- str_glue("data/{all_files_dir}")
names(all_files) <- all_files_dir

sample_info <-
  data.frame(
    sample = all_files_dir,
    id = seq_along(all_files_dir)
  )


# 6 samples
stopifnot(length(all_files) == 6)


all_data_ls <- 
  all_files %>%
  map(Read10X) %>%
  imap(~CreateSeuratObject(
    .x,
    project = .y,
    assay = "RNA",
    names.field = 1,
    names.delim = "_",
    min.cells = 10,
    min.features = 10
  ))

all_data_ls <- 
  all_data_ls %>%
  map(function(x){
    x$project <- x@project.name
    x
  })

### QC
all_data_ls <- 
  all_data_ls %>%
  map(function(x){
    x[["percent.mt"]]  <- PercentageFeatureSet(x, pattern = "^MT-")
    x[["percent.ribo"]]  <- PercentageFeatureSet(x, pattern = "^RP[SL]")
    x$log10GenesPerUMI <- log10(x$nFeature_RNA) / log10(x$nCount_RNA)
    x
  })


p_qc_ls <-
  all_data_ls %>%
  imap(function(x, name){
    VlnPlot(x, features = c("nFeature_RNA", "nCount_RNA", 
                            "percent.mt", "percent.ribo", "log10GenesPerUMI"), 
            ncol = 5, log = FALSE) +
      patchwork::plot_annotation(title = name)
  })

p_qc_log_ls <-
  all_data_ls %>%
  imap(function(x, name){
    VlnPlot(x, features = c("nFeature_RNA", "nCount_RNA", 
                            "percent.mt", "percent.ribo", "log10GenesPerUMI"), 
            ncol = 5, log = TRUE) +
      patchwork::plot_annotation(title = str_glue('{name} + log'))
  })

pdf(
  "qc plot.pdf",
  height = 6, 
  width = 20
)
p_qc_ls %>% walk(print)
p_qc_log_ls %>% walk(print)
dev.off()




### filter
qc_metric <- read.csv("qc_metric.csv") %>% map_dfc(str_trim)
qc_metric_variable <- c(
  'nCount_RNA', 'nFeature_RNA', 'percent.mt', 'log10GenesPerUMI'
)
qc_metric_ls <-
  qc_metric %>% 
  apply(1, name = names(.), simplify = FALSE, function(x, name = name){
    qc_ls <- map2(x, name, function(num, variable){
      if(!variable %in% qc_metric_variable) return()
      
      if(str_detect(num, ":")){
        num <- str_split(num, pattern = ":", n = 2) %>% .[[1]]
        qc <- str_glue("{variable} {num[1]} & {variable} {num[2]}")
      }else{
        qc <- str_glue("{variable} {num}")
      }
      qc
    })
    
    res <- qc_ls %>% discard(function(x) length(x) < 1) %>% str_c(collapse = " & ")
    res
  })

names(qc_metric_ls) <- qc_metric$sample

# filter
all_data_ls_pass <-
  names(all_data_ls) %>%
  set_names(.,.) %>%
  map(function(sample){ #browser()
    all_data <- all_data_ls[[sample]]
    qc <- qc_metric_ls[[sample]]
    
    e <- str_glue("subset(
      all_data,
      {qc}
    )")
    
    eval(parse(text = e))
  })


# qc plot
qc_plot_filter_ls <-
  names(all_data_ls) %>%
  map(function(sample){
    all_data <- all_data_ls[[sample]]
    all_data_pass <- all_data_ls_pass[[sample]]
    
    # all_data
    plot1 <- FeatureScatter(all_data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
      patchwork::plot_annotation(title = sample)
    plot2 <- FeatureScatter(all_data, feature1 = "nCount_RNA", feature2 = "percent.mt") +
      patchwork::plot_annotation(title = sample)
    plot3 <- FeatureScatter(all_data, feature1 = "nCount_RNA", feature2 = "log10GenesPerUMI") +
      patchwork::plot_annotation(title = sample)
    p1 <- plot1 + plot2 + plot3 + patchwork::plot_layout(ncol = 3)
    
    # alll_data_pass
    plot1_pass <- FeatureScatter(all_data_pass, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")+
      patchwork::plot_annotation(title = str_glue('{sample} + log'))
    plot2_pass <- FeatureScatter(all_data_pass, feature1 = "nCount_RNA", feature2 = "percent.mt")+
      patchwork::plot_annotation(title = str_glue('{sample} + log'))
    plot3_pass <- FeatureScatter(all_data_pass, feature1 = "nCount_RNA", feature2 = "log10GenesPerUMI")+
      patchwork::plot_annotation(title = str_glue('{sample} + log'))
    p2 <- plot1_pass + plot2_pass + plot3_pass  + patchwork::plot_layout(ncol = 3)
    
    p1 / p2
  })

# export 
pdf(
  "qc plot - pass.pdf",
  height = 12, 
  width = 18
)
qc_plot_filter_ls %>% walk(print)
dev.off()


# cells filter
cell_counts_df <-
  names(all_data_ls) %>%
  set_names(., .) %>%
  map(function(sample){
    all_data <- all_data_ls[[sample]]
    all_data_pass <- all_data_ls_pass[[sample]]
    
    structure(
      c(ncol(all_data), ncol(all_data_pass)),
      names = c("Before filter", "After filter")
    )
  }) %>%
  as.data.frame(check.names = FALSE)

#### add cell_id prefix
all_data_ls_pass <-
  all_data_ls_pass %>%
  names() %>%
  set_names(., .) %>%
  map(function(x){
    all_data <- all_data_ls_pass[[x]]
    ids <- sample_info$id[sample_info$sample == x]
    
    all_data <- RenameCells(all_data, add.cell.id = ids)
    all_data
  })

# all_data <- merge(all_data_ls_pass[[1]], all_data_ls_pass[-1])
all_data <- all_data_ls_pass[[1]]

### normalize
all_data <- NormalizeData(all_data)
all_data <- FindVariableFeatures(all_data)
all_data <- ScaleData(all_data, vars.to.regress = 'nCount_RNA')


### dimension reduction
all_data <- RunPCA(all_data)

pdf('pca plot.pdf', height = 6, width = 8)
ElbowPlot(all_data, ndims = 30)
dev.off()

pdf('dimheatmap 1-5.pdf', height = 6, width = 10)
DimHeatmap(all_data, dims = 1:5)
dev.off()
pdf('dimheatmap 6-10.pdf', height = 6, width = 10)
DimHeatmap(all_data, dims = 6:10)
dev.off()
pdf('dimheatmap 11-15.pdf', height = 6, width = 10)
DimHeatmap(all_data, dims = 11:15)
dev.off()
pdf('dimheatmap 16-20.pdf', height = 6, width = 10)
DimHeatmap(all_data, dims = 16:20)
dev.off()

saveRDS(all_data, "all_data_merge.rds")

# ndims 8 10   ---------------------
# ndims 8   ---------------------

all_data <- RunUMAP(all_data, dims = 1:5)

pdf("Dimplot dim8 project.pdf", height = 6, width = 8)
DimPlot(all_data, group.by = "project", label = TRUE) + NoLegend()
dev.off()


all_data <- all_data %>% 
  harmony::RunHarmony(
    dims.use = 1:8,
    group.by.vars = "project",
    reduction.save = "harmony"
  )

all_data <- RunUMAP(
  all_data, 
  dims = 1:8, 
  reduction = "harmony", 
  reduction.name = "harmonyumap", 
  reduction.key  = "harmonyUMAP_"
)


pdf("Dimplot dim8 project - harmony.pdf", height = 6, width = 8)
DimPlot(all_data, group.by = "project", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()

### cluster
all_data <-
  FindNeighbors(
    all_data,
    reduction = "harmony",
    dims = 1:8,
    k.param = 15
  )
all_data <-
  FindClusters(
    all_data,
    resolution = c(0.2, 0.4, 0.6, 0.8)
  )


# dimplot
pdf("Dimplot dim8 cluster - harmony.pdf", height = 6, width = 8)
DimPlot(all_data, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(all_data, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(all_data, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(all_data, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(all_data, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(all_data, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(all_data, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(all_data, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()


# dotplot
marker_ls <- c(
  "Cardiomyocytes" = "Myh6",
  "Cardiomyocytes" = "Ttn",
  "Cardiomyocytes" = "Tnnt2",
  "Cardiomyocytes" = "Actn2",
  "Cardiomyocytes" = "Gja1",
  "Smooth muscle cells" = "Myh11",
  "Smooth muscle cells" = "Mylk",
  "Smooth muscle cells" = "Cald1",
  "Smooth muscle cells" = "Pdgfrb",
  "Smooth muscle cells" = "Rgs5",
  "Fibroblasts" = "Col1a1",
  "Fibroblasts" = "Col1a2",
  "Fibroblasts" = "Col13a1",
  "Fibroblasts" = "Pdgfra",
  "Fibroblasts" = "Dcn",
  "Immune cells" = "Ctsc",
  "Immune cells" = "Mrc1",
  "Immune cells" = "Cd74",
  "Immune cells" = "F13a1",
  "Immune cells" = "Csf1r",
  "Endothelial cells" = "Fabp4",
  "Endothelial cells" = "Egfl7",
  "Endothelial cells" = "Flt1",
  "Endothelial cells" = "Cdh5",
  "Endothelial cells" = "Emcn"
)


pdf("dotplot dims8.pdf", width = 25, height = 6)
DotPlot(all_data, group.by = "RNA_snn_res.0.2", features = marker_ls)
DotPlot(all_data, group.by = "RNA_snn_res.0.4", features = marker_ls)
DotPlot(all_data, group.by = "RNA_snn_res.0.6", features = marker_ls)
DotPlot(all_data, group.by = "RNA_snn_res.0.8", features = marker_ls)
dev.off()

# saveRDS(all_data, 'all_data_dim8.rds')

# ndims 10 ---------------

pdf("Dimplot dim8 project.pdf", height = 6, width = 8)
DimPlot(all_data, group.by = "project", label = TRUE) + NoLegend()
dev.off()


all_data <- all_data %>% 
  harmony::RunHarmony(
    dims.use = 1:8,
    group.by.vars = "project",
    reduction.save = "harmony"
  )

all_data <- RunUMAP(
  all_data, 
  dims = 1:8, 
  reduction = "harmony", 
  reduction.name = "harmonyumap", 
  reduction.key  = "harmonyUMAP_"
)


pdf("Dimplot dim8 project - harmony.pdf", height = 6, width = 8)
DimPlot(all_data, group.by = "project", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()

### cluster
all_data <-
  FindNeighbors(
    all_data,
    reduction = "harmony",
    dims = 1:8,
    k.param = 15
  )
all_data <-
  FindClusters(
    all_data,
    resolution = c(0.2, 0.4, 0.6, 0.8)
  )


# dimplot
pdf("Dimplot dim8 cluster - harmony.pdf", height = 6, width = 8)
DimPlot(all_data, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(all_data, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(all_data, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(all_data, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(all_data, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(all_data, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(all_data, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(all_data, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()


# dotplot
marker_ls <- c(
  "Cardiomyocytes" = "Myh6",
  "Cardiomyocytes" = "Ttn",
  "Cardiomyocytes" = "Tnnt2",
  "Cardiomyocytes" = "Actn2",
  "Cardiomyocytes" = "Gja1",
  "Smooth muscle cells" = "Myh11",
  "Smooth muscle cells" = "Mylk",
  "Smooth muscle cells" = "Cald1",
  "Smooth muscle cells" = "Pdgfrb",
  "Smooth muscle cells" = "Rgs5",
  "Fibroblasts" = "Col1a1",
  "Fibroblasts" = "Col1a2",
  "Fibroblasts" = "Col13a1",
  "Fibroblasts" = "Pdgfra",
  "Fibroblasts" = "Dcn",
  "Immune cells" = "Ctsc",
  "Immune cells" = "Mrc1",
  "Immune cells" = "Cd74",
  "Immune cells" = "F13a1",
  "Immune cells" = "Csf1r",
  "Endothelial cells" = "Fabp4",
  "Endothelial cells" = "Egfl7",
  "Endothelial cells" = "Flt1",
  "Endothelial cells" = "Cdh5",
  "Endothelial cells" = "Emcn"
)


pdf("dotplot dims8.pdf", width = 25, height = 6)
DotPlot(all_data, group.by = "RNA_snn_res.0.2", features = marker_ls)
DotPlot(all_data, group.by = "RNA_snn_res.0.4", features = marker_ls)
DotPlot(all_data, group.by = "RNA_snn_res.0.6", features = marker_ls)
DotPlot(all_data, group.by = "RNA_snn_res.0.8", features = marker_ls)
dev.off()

# saveRDS(all_data, 'all_data_dim8.rds')





###### dim15 --------------------
# script sc_analysis.R --dims 15 --rds all_data_dim8.rds --save --k.param 15

library(Seurat)
library(tidyverse)

all_data <- readRDS('all_data_dim15.rds')

Idents(all_data) <- "RNA_snn_res.0.2"

cell_cluster <- c(
  '0'  = "Cardiomyocytes",
  '1'  = "Fibroblasts",
  '2'  = "Endothelial cells",
  '3'  = "Immune cells",
  '4'  = "Smooth muscle cells",
  '5'  = "Endothelial cells",
  '6'  = "Fibroblasts",
  '7'  = "cluster7",
  '8'  = "Cardiomyocytes",
  '9'  = "cluster9",
  '10' = "cluster10",
  '11' = "cluster11"
)

all_data <- RenameIdents(all_data, cell_cluster)
all_data$annotation <- Idents(all_data)



# cell props
all_data$group <- all_data$project %>% str_remove("[0-9]$")
plot_data <- 
  all_data[[]] %>%
  dplyr::select(
    group,
    annotation
  ) %>% 
  group_by(group, annotation) %>% 
  summarise(count = n())%>%
  group_by(group)%>%
  mutate(cell_props = count/sum(count)) %>%
  mutate(cell_props_label = cumsum(rev(cell_props)))
 

p_cell_props <-
  plot_data %>%
  ggplot(aes(x = group, y = cell_props, fill = annotation)) +
  geom_col() 
  # geom_text(aes(label = cell_props, y = cell_props_label), vjust = 1)
  

pdf("cell props - dims 15.pdf", height = 6, width = 8)
print(p_cell_props)
dev.off()

#### 
pdf(str_glue("dim15/Dimplot dim15 group - harmony.pdf"), height = 5, width = 15)
DimPlot(all_data, group.by = "group", label = TRUE, reduction = "harmonyumap") 
DimPlot(all_data, group.by = "annotation", split.by = "group", label = TRUE, reduction = "harmonyumap") 
FeaturePlot(all_data, features = "nCount_RNA",split.by = "group",  label = TRUE, reduction = "harmonyumap") 
FeaturePlot(all_data, features = "nFeature_RNA",split.by = "group",  label = TRUE, reduction = "harmonyumap") 
FeaturePlot(all_data, features = "percent.mt", split.by = "group", label = TRUE, reduction = "harmonyumap") 
FeaturePlot(all_data, features = "log10GenesPerUMI",split.by = "group",  label = TRUE, reduction = "harmonyumap")
dev.off()


pdf(str_glue("dim15/Dimplot dim15 project - harmony.pdf"), height = 5, width = 15)
DimPlot(all_data, group.by = "project", label = TRUE, reduction = "harmonyumap") 
DimPlot(all_data, group.by = "annotation", split.by = "project", label = TRUE, reduction = "harmonyumap") 
FeaturePlot(all_data, features = "nCount_RNA",split.by = "project",  label = TRUE, reduction = "harmonyumap") 
FeaturePlot(all_data, features = "nFeature_RNA",split.by = "project",  label = TRUE, reduction = "harmonyumap") 
FeaturePlot(all_data, features = "percent.mt", split.by = "project", label = TRUE, reduction = "harmonyumap") 
FeaturePlot(all_data, features = "percent.ribo", split.by = "project", label = TRUE, reduction = "harmonyumap") 
FeaturePlot(all_data, features = "log10GenesPerUMI",split.by = "project",  label = TRUE, reduction = "harmonyumap")
dev.off()

  
# all_data[["percent.ribo"]]  <- PercentageFeatureSet(all_data, pattern = "^Rp[sl]")



# cycle2  ----------------
## dim8 -----
## Rscript sc_analysis.R --dims 8 --rds all_data_merge.rds --save --k.param 20


library(Seurat)
library(tidyverse)

all_data <- readRDS('all_data_dim8.rds')

Idents(all_data) <- "RNA_snn_res.0.2"

cell_cluster <- c(
  '0'  = "Fibroblasts",
  '1'  = "Cardiomyocytes",
  '2'  = "Endothelial cells",
  '3'  = "Immune cells",
  '4'  = "Smooth muscle cells",
  '5'  = "Endothelial cells",
  '6'  = "Fibroblasts",
  '7'  = "Cardiomyocytes - Cluster7",
  '8'  = "Cardiomyocytes",
  '9'  = "Fibroblasts"
)

all_data <- RenameIdents(all_data, cell_cluster)
all_data$annotation <- Idents(all_data)



# cell props
all_data$group <- all_data$project %>% str_remove("[0-9]$")
plot_data <- 
  all_data[[]] %>%
  dplyr::select(
    group,
    annotation
  ) %>% 
  group_by(group, annotation) %>% 
  summarise(count = n())%>%
  group_by(group)%>%
  mutate(cell_props = count/sum(count)) %>%
  mutate(cell_props_label = cumsum(rev(cell_props)))


p_cell_props <-
  plot_data %>%
  ggplot(aes(x = group, y = cell_props, fill = annotation)) +
  geom_col() 
# geom_text(aes(label = cell_props, y = cell_props_label), vjust = 1)


pdf("cell props - dims 8.pdf", height = 6, width = 8)
print(p_cell_props)
dev.off()

#### 
pdf(str_glue("dim8/Dimplot dim8 group - harmony.pdf"), height = 5, width = 15)
DimPlot(all_data, group.by = "group", label = TRUE, reduction = "harmonyumap") 
DimPlot(all_data, group.by = "annotation", split.by = "group", label = TRUE, reduction = "harmonyumap") 
FeaturePlot(all_data, features = "nCount_RNA",split.by = "group",  label = TRUE, reduction = "harmonyumap") 
FeaturePlot(all_data, features = "nFeature_RNA",split.by = "group",  label = TRUE, reduction = "harmonyumap") 
FeaturePlot(all_data, features = "percent.mt", split.by = "group", label = TRUE, reduction = "harmonyumap") 
FeaturePlot(all_data, features = "log10GenesPerUMI",split.by = "group",  label = TRUE, reduction = "harmonyumap")
dev.off()


pdf(str_glue("dim8/Dimplot dim8 project - harmony.pdf"), height = 5, width = 15)
DimPlot(all_data, group.by = "project", label = TRUE, reduction = "harmonyumap") 
DimPlot(all_data, group.by = "annotation", split.by = "project", label = TRUE, reduction = "harmonyumap") 
FeaturePlot(all_data, features = "nCount_RNA",split.by = "project",  label = TRUE, reduction = "harmonyumap") 
FeaturePlot(all_data, features = "nFeature_RNA",split.by = "project",  label = TRUE, reduction = "harmonyumap") 
FeaturePlot(all_data, features = "percent.mt", split.by = "project", label = TRUE, reduction = "harmonyumap") 
FeaturePlot(all_data, features = "percent.ribo", split.by = "project", label = TRUE, reduction = "harmonyumap") 
FeaturePlot(all_data, features = "log10GenesPerUMI",split.by = "project",  label = TRUE, reduction = "harmonyumap")
dev.off()


# all_data[["percent.ribo"]]  <- PercentageFeatureSet(all_data, pattern = "^Rp[sl]")


plot_data %>%
  dplyr::filter(annotation == "Cardiomyocytes")


plot_data %>%
  dplyr::filter(annotation == "Cardiomyocytes")



pdf("test.pdf")
all_data[[]] %>% 
  # dplyr::select(nFeature_RNA, group) %>%
  dplyr::filter(annotation == "Cardiomyocytes") %>%
  ggplot(aes(x = group, y = nFeature_RNA, fill = group)) +
  geom_boxplot()
  
dev.off()


pdf("testw.pdf", height = 5, width = 12)
VlnPlot(all_data, features = "nCount_RNA",group.by = "annotation", split.by = "group")
dev.off()

saveRDS(all_data, 'all_data_dim8_annotation.rds')


## dim8 - not use harmony -----

### cluster
all_data <-
  FindNeighbors(
    all_data,
    reduction = "pca",
    dims = 1:8,
    k.param = 15
  )
all_data <-
  FindClusters(
    all_data,
    resolution = c(0.2, 0.4, 0.6, 0.8)
  )


# dimplot
pdf(str_glue("dim8_no_harmony/Dimplot dim8 cluster - harmony.pdf"), height = 6, width = 8)
DimPlot(all_data, group.by = "project", label = TRUE, reduction = "umap") + NoLegend()

DimPlot(all_data, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(all_data, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(all_data, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(all_data, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = "umap") + NoLegend()
dev.off()


# dotplot
marker_ls <- c(
  "Cardiomyocytes" = "Myh6",
  "Cardiomyocytes" = "Ttn",
  "Cardiomyocytes" = "Tnnt2",
  "Cardiomyocytes" = "Actn2",
  "Cardiomyocytes" = "Gja1",
  "Smooth muscle cells" = "Myh11",
  "Smooth muscle cells" = "Mylk",
  "Smooth muscle cells" = "Cald1",
  "Smooth muscle cells" = "Pdgfrb",
  "Smooth muscle cells" = "Rgs5",
  "Fibroblasts" = "Col1a1",
  "Fibroblasts" = "Col1a2",
  "Fibroblasts" = "Col13a1",
  "Fibroblasts" = "Pdgfra",
  "Fibroblasts" = "Dcn",
  "Immune cells" = "Ctsc",
  "Immune cells" = "Mrc1",
  "Immune cells" = "Cd74",
  "Immune cells" = "F13a1",
  "Immune cells" = "Csf1r",
  "Endothelial cells" = "Fabp4",
  "Endothelial cells" = "Egfl7",
  "Endothelial cells" = "Flt1",
  "Endothelial cells" = "Cdh5",
  "Endothelial cells" = "Emcn"
)


pdf(str_glue("dim8_no_harmony/dotplot dims8.pdf"), width = 25, height = 6)
DotPlot(all_data, group.by = "RNA_snn_res.0.2", features = marker_ls)
DotPlot(all_data, group.by = "RNA_snn_res.0.4", features = marker_ls)
DotPlot(all_data, group.by = "RNA_snn_res.0.6", features = marker_ls)
DotPlot(all_data, group.by = "RNA_snn_res.0.8", features = marker_ls)
dev.off()



