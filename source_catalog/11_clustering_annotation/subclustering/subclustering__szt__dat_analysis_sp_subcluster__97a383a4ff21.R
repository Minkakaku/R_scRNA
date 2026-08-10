library(tidyverse)
library(Seurat)
library(patchwork)

source("subcluster/sc_analysis_tools.R")

all_data <- readRDS('output/data/adata_filtered_add_celltype.rds')

out_dir <- "sp_subcluster"
dir.create(out_dir)

all_data$annotation_full <- all_data$annotation
all_data$annotation_full %>% table
# B              Endothelial Epithelial proliferative             Epithelial_1 
# 89515                   119845                   328204                    61475 
# Epithelial_2               Fibroblast             Fibroblast-T        Immune-Epithelial 
# 347381                   538602                    48274                    30859 
# Macrophage                 Pericyte                   Plasma                        T 
# 169641                   123042                    80867                   212926 


# Epi subcluster ----------------
sp_epi <- all_data[, all_data$annotation_full %in% c("Epithelial_1", "Epithelial_2", "Epithelial proliferative")]

DefaultAssay(sp_epi) <- 'RNA'

# saveRDS(sp_epi, 'sp_epi.rds')
# sp_epi <- readRDS('sp_epi.rds')

# sp_epi <- sc_normalize(sp_epi, scale.factor = "median")
sp_epi <- sc_normalize(sp_epi, scale.factor = "default")

## sketch data --------
sp_epi <-
  SketchData(
    sp_epi,
    ncells         = 50000L,
    sketched.assay = "sketch",
    method         = "LeverageScore",
    var.name       = "leverage.score",
    seed = 123L
  )
sp_epi <- sc_normalize(sp_epi, scale.factor = "default", nfeatures = 1500)


sp_epi <- sc_pca(
  sp_epi, 
  sc_out_dir = file.path(out_dir, "epi/sc.1 pca"),
  dir.create.force = TRUE
)

set.seed(1234)
sp_epi2 <- sp_epi[, sample(colnames(sp_epi), 1e4)]

pdf(file.path(out_dir, "epi/sc.1 pca", 'dimheatmap 1-5 - downsample.pdf'), height = 6, width = 10)
DimHeatmap(sp_epi2, dims = 1:5)
dev.off()
pdf(file.path(out_dir, "epi/sc.1 pca", 'dimheatmap 6-10 - downsample.pdf'), height = 6, width = 10)
DimHeatmap(sp_epi2, dims = 6:10)
dev.off()
pdf(file.path(out_dir, "epi/sc.1 pca", 'dimheatmap 11-15 - downsample.pdf'), height = 6, width = 10)
DimHeatmap(sp_epi2, dims = 11:15)
dev.off()
pdf(file.path(out_dir, "epi/sc.1 pca", 'dimheatmap 16-20 - downsample.pdf'), height = 6, width = 10)
DimHeatmap(sp_epi2, dims = 16:20)
dev.off()


ndims_epi <- 10

sp_epi <- sc_dr(
  sp_epi, 
  ndims = ndims_epi,
  return.model = TRUE,
  sc_out_dir = file.path(out_dir, "epi/sc.2 dr"),
  dir.create.force = TRUE
)


sp_epi <- sc_cluster(
  sp_epi, 
  ndims = ndims_epi,
  reduction =  "harmony", # "harmony", "pca"
  sc_out_dir = file.path(out_dir, "epi/sc.3 cluster"),
  dir.create.force = TRUE,
  
  harmony.reduction.name = "harmonyumap",
  resolution = c(0.2, 0.4, 0.6, 0.8)
)


sp_epi %>%
  VlnPlot(
    features = "total_counts",
    group.by = "sketch_snn_res.0.2",
    pt.size = 0
  ) 



pdf(
  file.path(out_dir, 'epi', 'clustree.pdf'),
  height = 8,
  width = 8
)
library(clustree)
clustree::clustree(
  sp_epi,
  prefix = "RNA_snn_res."
)
dev.off()



epi_markers <- c( 
  'EPCAM',
  "CD3D",
  "CD3E",
  "CD8A",
  
  "squamous cell carcinoma" = "KRT5",
  "squamous cell carcinoma" = "KRT6A",
  
  "squamocolumnar junction" = "KRT7",
  "squamocolumnar junction" = "KRT17",
  
  "columnar epithelium" = "MUC6",
  "columnar epithelium" = "MUC5B",
  
  "EMT" = "COL17A1",
  "EMT" = "POSTN",
  "EMT" = "SNAI2",
  
  "Prof." = "MKI67",
  "Prof." = "TOP2A",
  "Prof." = "UBE2C",
  "Prof." = "CDC20",
  
  "MHC I" = "HLA-A",
  "MHC I" = "HLA-B",
  "MHC I" = "HLA-C"
  
)

sp_epi <- sc_dotplot(
  sp_epi, 
  ndims = ndims_epi,
  marker = epi_markers,
  width = 12,
  sc_out_dir = file.path(out_dir, "epi/sc.4 dotplot"), 
  dir.create.force = TRUE
)

sp_epi <- sc_dotplot(
  sp_epi,
  ndims = ndims_epi,
  marker = major_marker_for_subcluster,
  width = 20,
  sc_out_dir = file.path(out_dir, "epi/sc.4 dotplot - major"),
  dir.create.force = TRUE
)

for (i in str_glue("RNA_snn_res.0.{c(2,4,6,8)}")){
  all_markers_epi <- 
    sc_findmarker(
      sp_epi,
      idents = i, # "RNA_snn_res.0.2",
      test.use = "wilcox",
      logfc.threshold = 0.25,
      min.pct = 0.1,
      only.pos = TRUE,
      # method = "conserved",
      grouping.var = "group"
    )
  write.csv(
    all_markers_epi,
    file = file.path(out_dir, 'epi', str_glue("de gene table  - {i}.csv") ),
    row.names = F
  )
}
  


sc_findmarker(sp_epi, all_markers = all_markers_epi, sort_by = "pct")


# spatial plot
# dat_sp_subcluster <-
#   Embeddings( all_data[['spatial']] ) %>%
#   merge(
#     sp_epi[[c(str_glue("RNA_snn_res.0.{c(2,4,6,8)}"), 'sample')]],
#     all.x = TRUE,
#     by = "row.names"
#   )
dat_sp_coords <- cbind(Embeddings( all_data[['spatial']] ), all_data[['sample']] )
dat_sp_epi    <- sp_epi[[c(str_glue("RNA_snn_res.0.{c(2,4,6,8)}"), 'sample')]]

p_sp_epi_sub_ls <-
  dat_sp_coords$sample %>% 
  unique() %>%
  as.character() %>%
  setNames(., .) %>% # .[1] %>%
  map(function(sam){
    dat1 <- dat_sp_coords[dat_sp_coords$sample ==  sam, ] 
    
    height <- 5
    width <- 5
    p <-
      dat1 %>%
      
      ggplot(aes(x = SPATIAL_1, SPATIAL_2)) +
      geom_tile(
        fill = "gray90",
        height = height, 
        width = width
      ) +
      theme_void() +
      coord_equal() +
      scale_y_reverse() +
      ggsci::scale_fill_igv(na.translate = T) +
      theme(panel.background = element_rect(fill = "black", color = NA)) 
    
    p_ls <-
      str_glue("RNA_snn_res.0.{c(2,4,6,8)}") %>%
      set_names(., .) %>%
      map(function(x){
        dat2 <- dat_sp_epi[dat_sp_epi$sample ==  sam, ] 
        dat2 <- dat2 %>% merge(dat1, by = "row.names")
        # print(head(dat2))
        p +
          geom_tile(
            data = dat2,
            aes(fill = .data[[x]]),
            height = height, 
            width = width,
            show.legend = c(fill = TRUE, alpha = FALSE)
          ) +
          ggtitle(x)
      })
    p_ls
    
  })

dir.create(file.path(out_dir, "epi/sc.7 spatial"), recursive = TRUE)
p_sp_epi_sub_ls %>%
  iwalk(function(p_ls, name){
    pdf(
      file.path(out_dir, "epi/sc.7 spatial", str_glue("spatial plot - {name}.pdf") ),
      height = 8,
      width = 8
    )
    p_ls %>% walk(print)
    dev.off()
  })




# #############33333
# dir.create(file.path(out_dir, "epi/sc.6 heatmap"), recursive = TRUE)
# pdf(
#   file.path(out_dir, "epi/sc.6 heatmap", 'heatmap.pdf'),
#   height = 4,
#   width  = 8
# )
# sc_mark_heatmap(
#   sp_epi,
#   markers = all_markers_epi,
#   idents = "RNA_snn_res.0.2",
#   p_val_adj = 0.05,
#   pct.1 = 0.15,
#   n_head = 5,
#   row_cluster = FALSE
# )
# dev.off()
# 
# 
# sc_mark_heatmap(
#   sp_epi,
#   markers = all_markers_epi,
#   idents = "RNA_snn_res.0.2",
#   p_val_adj = 0.05,
#   pct.1 = 0.15,
#   n_head = 5,
#   row_cluster = FALSE,
#   scale = "none",
#   log_trans = TRUE
# )
# 
# epi_cell_cluster <- c(
#   '0'  =  "Proximal mature enterocytes",
#   '1'  =  "Peneth cell",
#   '2'  =  "Transit amplifying cell",
#   '3'  =  "Others",
#   '4'  =  "Enterocytes",
#   '5'  =  "Others",
#   '6'  =  "Others",
#   '7'  =  "CD8T",
#   '8'  =  "Transit amplifying cell",
#   '9'  =  "Others",
#   '10'  = "Others",
#   '11'  = "Others",
#   '12'  = "Others",
#   '13'  = "Proximal mature enterocytes",
#   '14'  = "Others"
# 
# )
# 
# 
# 
# Idents(sp_epi) <- "RNA_snn_res.0.4"
# sp_epi <- RenameIdents(sp_epi, epi_cell_cluster)
# sp_epi$annotation_epi <- Idents(sp_epi)
# 
# dir.create(file.path(out_dir, "epi/sc.5 annotation"))
# pdf(file.path(out_dir, 'epi/sc.5 annotation', str_glue("Dimplot dim{ndims_epi} epi annotation .pdf")), height = 6, width = 8)
# DimPlot(sp_epi, group.by = "annotation_epi", label = TRUE, reduction = "umap") + NoLegend()
# DimPlot(sp_epi, group.by = "annotation_epi", label = TRUE, reduction = "harmonyumap") + NoLegend()
# dev.off()


##  project  data --------------------
DefaultAssay(sp_epi) <- "sketch"
# sp_epi[['full.harmonyumap']] <- NULL
sp_epi <-
  sp_epi %>%
  ProjectData(
    assay = "RNA",
    full.reduction = "harmony.full",

    sketched.assay = "sketch",
    sketched.reduction = "harmony", # "pca",

    umap.model = "harmonyumap",
    dims = 1:ndims_epi,

    refdata = list(
      # annotation_full = "annotation",
      RNA_snn_res.0.2 = 'sketch_snn_res.0.2',
      RNA_snn_res.0.4 = 'sketch_snn_res.0.4',
      RNA_snn_res.0.6 = 'sketch_snn_res.0.6',
      RNA_snn_res.0.8 = 'sketch_snn_res.0.8'
    )
  )

Idents(sp_epi) <- 'RNA_snn_res.0.2' # 'annotation_full'
DefaultAssay(sp_epi) <- "RNA"

# pdf(file.path(out_dir, 'epi', str_glue("Dimplot dim{ndims_epi} annotation - projected.pdf")), height = 6, width = 8)
# DimPlot(sp_epi, group.by = "annotation_full", label = TRUE, reduction = 'full.harmonyumap') + NoLegend()
# dev.off()

pdf(file.path(out_dir, 'epi', str_glue("Dimplot dim{ndims_epi} annotation - projected.pdf")), height = 6, width = 8)
DimPlot(sp_epi, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = 'full.harmonyumap') + NoLegend()
DimPlot(sp_epi, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = 'full.harmonyumap') + NoLegend()
DimPlot(sp_epi, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = 'full.harmonyumap') + NoLegend()
DimPlot(sp_epi, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = 'full.harmonyumap') + NoLegend()
dev.off()

# saveRDS(sp_epi, file.path(out_dir, 'sp_epi_projected.rds'))


# DimPlot(sp_epi, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = 'full.harmonyumap') + NoLegend()



pdf(
  file.path(out_dir, 'epi', 'clustree.pdf'),
  height = 8,
  width = 8
)
library(clustree)
# clustree::clustree(
#   sp_epi,
#   prefix = "RNA_snn_res."
# )
sp_epi[[]] %>%
  dplyr::select(
    RNA_snn_res.0.2,
    RNA_snn_res.0.4,
    RNA_snn_res.0.6,
    RNA_snn_res.0.8
  ) %>%
  clustree::clustree(
    prefix = "RNA_snn_res."
  )
dev.off()

## project data 2 --------------
# 见`dat_analysis_260518.R`
# remove 7 in sketch assay, and project to RNA assay again.



# de analysis -------
for (i in str_glue("RNA_snn_res.0.{c(2,4,6,8)}")){
  all_markers_epi <- 
    sc_findmarker(
      sp_epi,
      idents = i, # "RNA_snn_res.0.2",
      test.use = "wilcox",
      logfc.threshold = 0.25,
      min.pct = 0.1,
      only.pos = TRUE,
      # method = "conserved",
      grouping.var = "group"
    )
  write.csv(
    all_markers_epi,
    file = file.path(out_dir, 'epi', str_glue("de gene table  - {i}.csv") ),
    row.names = F
  )
}


# cluster中HPV gene dotplot  -------
hpv_gene <- sp_epi %>% rownames() %>% str_subset("^V[0-9]{2}[EL]")

# p_dot_hpv <-
dot_p_hpv_ls <-
  sp_epi %>%
  Assays() %>%
  setNames(., .) %>%
  map(function(x){
    "RNA_snn_res.0.{c(2,4,6,8)}" %>%
      str_glue() %>%
      map(function(y){
        DotPlot(
          sp_epi,
          assay = x,
          features = hpv_gene,
          group.by = y
        ) +
        # RotatedAxis()
        coord_flip() +
        ggtitle(y)
      })
    
  })

# export
pdf(
  file.path(out_dir, 'epi', 'dotplot of HPV gene - sketch assay.pdf'),
  height = 15,
  width  = 8
)
dot_p_hpv_ls$sketch %>% walk(print)
dev.off()

pdf(
  file.path(out_dir, 'epi', 'dotplot of HPV gene - RNA assay.pdf'),
  height = 15,
  width  = 8
)
dot_p_hpv_ls$RNA %>% walk(print)
dev.off()




# cluster中sample占比  -------
sample_color <- c(
  'NC' = '#87CEEB',
  'A'	 = '#F4D03F',
  'CH' = '#F5B7B1',
  'HL' = '#A9DFBF'
)
sample_levs <- c("NC", "A", "CH", "HL")

p_col_sketch_ls <-
  "sketch_snn_res.0.{c(2,4,6,8)}" %>%
  str_glue() %>%
  setNames(., .) %>%
  map(function(x){
    p1 <-
      sp_epi[[]] %>%
      group_by(.data[[x]], project) %>%
      
      dplyr::filter(!is.na(.data[[x]])) %>% 
      
      summarise(count = n()) %>%
      group_by(.data[[x]]) %>%
      mutate(
        project2 = project %>% str_remove("[0-9].*$") %>% factor(levels = sample_levs),
        freq = count / sum(count)
      ) %>%
      
      ggplot(aes(x = .data[[x]], y = freq, fill = project2)) +
      geom_col(color = "gray") +
      scale_fill_manual(values = sample_color) +
      theme_classic()  +
      labs(
        title = x,
        fill = "sample_group"
      ) +
      RotatedAxis()
    
    p2 <-
      sp_epi[[]] %>%
      group_by(.data[[x]], project) %>%
      
      dplyr::filter(!is.na(.data[[x]])) %>% 
      
      summarise(count = n()) %>%
      group_by(project) %>%
      mutate(
        project2 = project %>% str_remove("[0-9].*$") %>% factor(levels = sample_levs),
        freq = count / sum(count)
      ) %>%
      
      ggplot(aes(x = project, color = project2, y = freq, fill = .data[[x]])) +
      geom_col(width = 0.5, linewidth = 1.5) +
      scale_color_manual(values = sample_color) +
      ggsci::scale_fill_igv() +
      theme_classic()  +
      labs(
        title = x,
        color = "sample_group"
      )+
      RotatedAxis()
    
    # return
    list(
      p1 + theme(plot.margin = margin(2, 100, 2, 100)),
      p2
    )
  })


p_col_rna_ls <-
  "RNA_snn_res.0.{c(2,4,6,8)}" %>%
  str_glue() %>%
  setNames(., .) %>%
  map(function(x){
    p1 <-
      sp_epi[[]] %>%
      group_by(.data[[x]], project) %>%
      
      dplyr::filter(!is.na(.data[[x]])) %>% 
      
      summarise(count = n()) %>%
      group_by(.data[[x]]) %>%
      mutate(
        project2 = project %>% str_remove("[0-9].*$") %>% factor(levels = sample_levs),
        freq = count / sum(count)
      ) %>%
      
      ggplot(aes(x = .data[[x]], y = freq, fill = project2)) +
      geom_col(color = "gray") +
      scale_fill_manual(values = sample_color) +
      theme_classic() +
      labs(
        title = x,
        fill = "sample_group"
      )+
      RotatedAxis()
    
    p2 <-
      sp_epi[[]] %>%
      group_by(.data[[x]], project) %>%
      
      dplyr::filter(!is.na(.data[[x]])) %>% 
      
      summarise(count = n()) %>%
      group_by(project) %>%
      mutate(
        project2 = project %>% str_remove("[0-9].*$") %>% factor(levels = sample_levs),
        freq = count / sum(count)
      ) %>%
      
      ggplot(aes(x = project, color = project2, y = freq, fill = .data[[x]])) +
      geom_col(width = 0.5, linewidth = 1.5) +
      scale_color_manual(values = sample_color) +
      ggsci::scale_fill_igv() +
      theme_classic()  +
      labs(
        title = x,
        color = "sample_group"
      )+
      RotatedAxis()
    
    # return
    list(
      p1 + theme(plot.margin = margin(2, 100, 2, 100)),
      p2
    )
  })


# export
pdf(
  file.path(out_dir, 'epi', 'barplot of sample - sketch assay.pdf'),
  height = 8,
  width  = 10
)
p_col_sketch_ls %>% print
dev.off()

pdf(
  file.path(out_dir, 'epi', 'barplot of sample - RNA assay.pdf'),
  height = 8,
  width  = 10
)
p_col_rna_ls %>% print
dev.off()


# 20260423 export cell props -----
dat_props_sketch_ls <- p_col_sketch_ls  %>%  map(~ .[[1]]$data)
dat_props_rna_ls    <- p_col_rna_ls     %>%  map(~ .[[1]]$data)

dat_props_sketch_ls %>%
  openxlsx::write.xlsx(
    file = file.path(out_dir, "epi", "cell props table of sketch assay - sample.xlsx"),
    asTable = TRUE,
    withFilter =  FALSE
  )
dat_props_rna_ls %>%
  openxlsx::write.xlsx(
    file = file.path(out_dir, "epi", "cell props table of RNA assay - sample.xlsx"),
    asTable = TRUE,
    withFilter =  FALSE
  )




# cluster中ROI占比  -------
## ROI from `dat_analysis_260331.R`

roi_color <- c(
  NC =	'#FBC04C',
  CC =	'#648FEB',
  CH =	'#52C7B7',
  CL =	'#CCE4FE',
  HH =	'#6E6FB2',
  HL =	'#DDC4E2'
)
roi_levs <- c('NC', 'CC', 'CH', 'CL', 'HH', 'HL')

roi_info <- readRDS('roi_cellid_group.rds')
roi_group_new <- 
  readxl::read_excel("26-4-8 样本组别.xlsx") %>%
  dplyr::select(
    roi_name_tidy = `3.30 ROI`,
    roi_name  = `ROI映射名称`,
    roi_group = `ROI映射组别`
  ) %>%
  mutate(roi_group = fct_recode(roi_group, NC = 'HPV-NC')) %>%
  mutate(roi_group = factor(roi_group, levels = roi_levs))

stopifnot( all(roi_group_new$roi_name_tidy %in% roi_info$ROI) )

roi_info2 <- roi_info %>% left_join(roi_group_new, by = c(ROI = "roi_name_tidy"))

dat_epi_roi <-
  sp_epi[[]] %>% 
  dplyr::filter(cell_id_sample %in% roi_info2$cell_id_sample) %>%
  left_join(
    roi_info2[c("cell_id_sample", "roi_name", "roi_group")],
    by = "cell_id_sample"
  )
stopifnot(all(! is.na(dat_epi_roi$roi_group )))

# plot 
p_col_sketch_roi_ls <-
  "sketch_snn_res.0.{c(2,4,6,8)}" %>%
  str_glue() %>%
  setNames(., .) %>%
  map(function(x){
    p1 <-
      dat_epi_roi %>%
      group_by(.data[[x]], roi_name) %>%
      
      dplyr::filter(!is.na(.data[[x]])) %>% 
      
      summarise(count = n()) %>%
      group_by(.data[[x]]) %>%
      mutate(freq = count / sum(count) ) %>%
      left_join(roi_group_new, by = 'roi_name') %>%
      
      ggplot(aes(x = .data[[x]], y = freq, fill = roi_group)) +
      geom_col(color = "gray") +
      scale_fill_manual(values = roi_color) +
      theme_classic()  +
      labs(
        title = x,
        fill = "roi_group"
      ) +
      RotatedAxis()
    
    p2 <-
      dat_epi_roi %>%
      group_by(.data[[x]], roi_name) %>%
      
      dplyr::filter(!is.na(.data[[x]])) %>% 
      
      summarise(count = n()) %>%
      group_by(roi_name) %>%
      mutate(freq = count / sum(count) ) %>%
      left_join(roi_group_new, by = 'roi_name') %>% 
      
      ggplot(aes(x = reorder(roi_name, as.numeric(roi_group)), color = roi_group, y = freq, fill = .data[[x]])) +
      geom_col(width = 0.5, linewidth = 1.5) +
      scale_color_manual(values = roi_color) +
      ggsci::scale_fill_igv() +
      theme_classic()  +
      labs(
        x = NULL,
        title = x,
        color = "roi_group"
      )+
      RotatedAxis()
    
    # return
    list(
      p1 + theme(plot.margin = margin(2, 200, 2, 200)), 
      p2
    )
  })


p_col_rna_roi_ls <-
  "RNA_snn_res.0.{c(2,4,6,8)}" %>%
  str_glue() %>%
  setNames(., .) %>%
  map(function(x){
    p1 <-
      dat_epi_roi %>%
      group_by(.data[[x]], roi_name) %>%
      
      dplyr::filter(!is.na(.data[[x]])) %>% 
      
      summarise(count = n()) %>%
      group_by(.data[[x]]) %>%
      mutate(freq = count / sum(count) ) %>%
      left_join(roi_group_new, by = 'roi_name') %>%
      
      ggplot(aes(x = .data[[x]], y = freq, fill = roi_group)) +
      geom_col(color = "gray") +
      scale_fill_manual(values = roi_color) +
      theme_classic()  +
      labs(
        title = x,
        fill = "roi_group"
      ) +
      RotatedAxis()
    
    p2 <-
      dat_epi_roi %>%
      group_by(.data[[x]], roi_name) %>%
      
      dplyr::filter(!is.na(.data[[x]])) %>% 
      
      summarise(count = n()) %>%
      group_by(roi_name) %>%
      mutate(freq = count / sum(count) ) %>%
      left_join(roi_group_new, by = 'roi_name') %>% 
      
      ggplot(aes(x = reorder(roi_name, as.numeric(roi_group)), color = roi_group, y = freq, fill = .data[[x]])) +
      geom_col(width = 0.5, linewidth = 1.5) +
      scale_color_manual(values = roi_color) +
      ggsci::scale_fill_igv() +
      theme_classic()  +
      labs(
        x = NULL,
        title = x,
        color = "roi_group"
      )+
      RotatedAxis()
    
    # return
    list(
      p1 + theme(plot.margin = margin(2, 200, 2, 200)),
      p2
    )
  })


# export
pdf(
  file.path(out_dir, 'epi', 'barplot of ROI - sketch assay.pdf'),
  height = 8,
  width  = 15
)
p_col_sketch_roi_ls %>% print
dev.off()

pdf(
  file.path(out_dir, 'epi', 'barplot of ROI - RNA assay.pdf'),
  height = 8,
  width  = 15
)
p_col_rna_roi_ls %>% print
dev.off()

# 20260423 export cell props -----
dat_props_sketch_roi_ls <- p_col_sketch_roi_ls  %>%  map(~ .[[1]]$data)
dat_props_rna_roi_ls    <- p_col_rna_roi_ls     %>%  map(~ .[[1]]$data)

dat_props_sketch_roi_ls %>%
  openxlsx::write.xlsx(
    file = file.path(out_dir, "epi", "cell props table of sketch assay - roi.xlsx"),
    asTable = TRUE,
    withFilter =  FALSE
  )
dat_props_rna_roi_ls %>%
  openxlsx::write.xlsx(
    file = file.path(out_dir, "epi", "cell props table of RNA assay - roi.xlsx"),
    asTable = TRUE,
    withFilter =  FALSE
  )



# 20260423 new ROI group ------------------
# cluster中ROI占比  
## ROI from `dat_analysis_260331.R`

roi_color <- c(
  'HPV-NC' = '#FBC04C',
  'HPV+NC' = '#648FEB',
  'Ca'	   = '#52C7B7',
  'HSIL'	 = '#CCE4FE',
  'LSIL'	 = '#6E6FB2'
)
roi_levs <- c("HPV-NC", "HPV+NC", "Ca", "HSIL", "LSIL")

roi_info <- readRDS('roi_cellid_group.rds')
roi_group_new <- 
  readxl::read_excel("26-4-17 样本组别.xlsx") %>%
  dplyr::select(
    roi_name_tidy = `3.30 ROI`,
    roi_name  = `ROI映射名称`,
    roi_group = `group4.17roi实际映射`
  ) %>%
  mutate(roi_group = factor(roi_group, levels = roi_levs))

stopifnot( all(roi_group_new$roi_name_tidy %in% roi_info$ROI) )

roi_info2 <- roi_info %>% left_join(roi_group_new, by = c(ROI = "roi_name_tidy"))

dat_epi_roi <-
  sp_epi[[]] %>% 
  dplyr::filter(cell_id_sample %in% roi_info2$cell_id_sample) %>%
  left_join(
    roi_info2[c("cell_id_sample", "roi_name", "roi_group")],
    by = "cell_id_sample"
  )
stopifnot(all(! is.na(dat_epi_roi$roi_group )))

# plot 
p_col_sketch_roi_ls2 <-
  "sketch_snn_res.0.{c(2,4,6,8)}" %>%
  str_glue() %>%
  setNames(., .) %>%
  map(function(x){
    p1 <-
      dat_epi_roi %>%
      group_by(.data[[x]], roi_name) %>%
      
      dplyr::filter(!is.na(.data[[x]])) %>% 
      
      summarise(count = n()) %>%
      group_by(.data[[x]]) %>%
      mutate(freq = count / sum(count) ) %>%
      left_join(roi_group_new, by = 'roi_name') %>%
      
      ggplot(aes(x = .data[[x]], y = freq, fill = roi_group)) +
      geom_col(color = "gray") +
      scale_fill_manual(values = roi_color) +
      theme_classic()  +
      labs(
        title = x,
        fill = "roi_group"
      ) +
      RotatedAxis()
    
    p2 <-
      dat_epi_roi %>%
      group_by(.data[[x]], roi_name) %>%
      
      dplyr::filter(!is.na(.data[[x]])) %>% 
      
      summarise(count = n()) %>%
      group_by(roi_name) %>%
      mutate(freq = count / sum(count) ) %>%
      left_join(roi_group_new, by = 'roi_name') %>% 
      
      ggplot(aes(x = reorder(roi_name, as.numeric(roi_group)), color = roi_group, y = freq, fill = .data[[x]])) +
      geom_col(width = 0.5, linewidth = 1.5) +
      scale_color_manual(values = roi_color) +
      ggsci::scale_fill_igv() +
      theme_classic()  +
      labs(
        x = NULL,
        title = x,
        color = "roi_group"
      )+
      RotatedAxis()
    
    # return
    list(
      p1 + theme(plot.margin = margin(2, 200, 2, 200)), 
      p2
    )
  })


p_col_rna_roi_ls2 <-
  "RNA_snn_res.0.{c(2,4,6,8)}" %>%
  str_glue() %>%
  setNames(., .) %>%
  map(function(x){
    p1 <-
      dat_epi_roi %>%
      group_by(.data[[x]], roi_name) %>%
      
      dplyr::filter(!is.na(.data[[x]])) %>% 
      
      summarise(count = n()) %>%
      group_by(.data[[x]]) %>%
      mutate(freq = count / sum(count) ) %>%
      left_join(roi_group_new, by = 'roi_name') %>%
      
      ggplot(aes(x = .data[[x]], y = freq, fill = roi_group)) +
      geom_col(color = "gray") +
      scale_fill_manual(values = roi_color) +
      theme_classic()  +
      labs(
        title = x,
        fill = "roi_group"
      ) +
      RotatedAxis()
    
    p2 <-
      dat_epi_roi %>%
      group_by(.data[[x]], roi_name) %>%
      
      dplyr::filter(!is.na(.data[[x]])) %>% 
      
      summarise(count = n()) %>%
      group_by(roi_name) %>%
      mutate(freq = count / sum(count) ) %>%
      left_join(roi_group_new, by = 'roi_name') %>% 
      
      ggplot(aes(x = reorder(roi_name, as.numeric(roi_group)), color = roi_group, y = freq, fill = .data[[x]])) +
      geom_col(width = 0.5, linewidth = 1.5) +
      scale_color_manual(values = roi_color) +
      ggsci::scale_fill_igv() +
      theme_classic()  +
      labs(
        x = NULL,
        title = x,
        color = "roi_group"
      )+
      RotatedAxis()
    
    # return
    list(
      p1 + theme(plot.margin = margin(2, 200, 2, 200)),
      p2
    )
  })


# export
pdf(
  file.path(out_dir, 'epi', 'barplot of ROI - sketch assay - new roi group.pdf'),
  height = 8,
  width  = 15
)
p_col_sketch_roi_ls2 %>% print
dev.off()

pdf(
  file.path(out_dir, 'epi', 'barplot of ROI - RNA assay - new roi group.pdf'),
  height = 8,
  width  = 15
)
p_col_rna_roi_ls2 %>% print
dev.off()

# export cell props
dat_props_sketch_roi_ls2 <- p_col_sketch_roi_ls2  %>%  map(~ .[[1]]$data)
dat_props_rna_roi_ls2    <- p_col_rna_roi_ls2     %>%  map(~ .[[1]]$data)

dat_props_sketch_roi_ls2 %>%
  openxlsx::write.xlsx(
    file = file.path(out_dir, "epi", "cell props table of sketch assay - roi - new roi group.xlsx"),
    asTable = TRUE,
    withFilter =  FALSE
  )
dat_props_rna_roi_ls2 %>%
  openxlsx::write.xlsx(
    file = file.path(out_dir, "epi", "cell props table of RNA assay - roi - new roi group.xlsx"),
    asTable = TRUE,
    withFilter =  FALSE
  )


# 20260423 new ROI group - add new cell props barplot ------------------
# x axis is group, fill is cell annotation  
p_col_sketch_roi_ls3 <-
  "sketch_snn_res.0.{c(2,4,6,8)}" %>%
  str_glue() %>%
  setNames(., .) %>%
  map(function(x){
    p <-
      dat_epi_roi %>%
      group_by(roi_group, .data[[x]]) %>%
      
      dplyr::filter(!is.na(.data[[x]])) %>% 
      
      summarise(count = n()) %>%
      group_by(roi_group) %>%
      mutate(freq = count / sum(count) ) %>%
      
      ggplot(aes(fill = .data[[x]], y = freq, x = roi_group)) +
      geom_col() +
      
      ggsci::scale_fill_igv() +
      theme_classic()  +
      labs(
        title = x,
        fill = "roi_group"
      ) +
      RotatedAxis()
    
    
    # return
    p
  })


p_col_rna_roi_ls3 <-
  "RNA_snn_res.0.{c(2,4,6,8)}" %>%
  str_glue() %>%
  setNames(., .) %>%
  map(function(x){
    p <-
      dat_epi_roi %>%
      group_by(roi_group, .data[[x]]) %>%
      
      dplyr::filter(!is.na(.data[[x]])) %>% 
      
      summarise(count = n()) %>%
      group_by(roi_group) %>%
      mutate(freq = count / sum(count) ) %>%
      
      ggplot(aes(fill = .data[[x]], y = freq, x = roi_group)) +
      geom_col() +
      # scale_fill_manual(values = roi_color) +
      ggsci::scale_fill_igv() +
      theme_classic()  +
      labs(
        title = x,
        fill = "roi_group"
      ) +
      RotatedAxis()
    
    # return
    p
  })


# export
pdf(
  file.path(out_dir, 'epi', 'barplot of ROI - sketch assay - new roi group - group is xaxis.pdf'),
  height = 8,
  width  = 6
)
p_col_sketch_roi_ls3 %>% print
dev.off()

pdf(
  file.path(out_dir, 'epi', 'barplot of ROI - RNA assay - new roi group - group is xaxis.pdf'),
  height = 8,
  width  = 6
)
p_col_rna_roi_ls3 %>% print
dev.off()


# export cell props
dat_props_sketch_roi_ls3 <- p_col_sketch_roi_ls3  %>%  map('data')
dat_props_rna_roi_ls3    <- p_col_rna_roi_ls3     %>%  map('data')

dat_props_sketch_roi_ls3 %>%
  openxlsx::write.xlsx(
    file = file.path(out_dir, "epi", "cell props table of sketch assay - roi - new roi group - group is xaxis.xlsx"),
    asTable = TRUE,
    withFilter =  FALSE
  )
dat_props_rna_roi_ls3 %>%
  openxlsx::write.xlsx(
    file = file.path(out_dir, "epi", "cell props table of RNA assay - roi - new roi group - group is xaxis.xlsx"),
    asTable = TRUE,
    withFilter =  FALSE
  )



# spatial plot of epi subcluster -----------
dat_sp <- cbind(
  Embeddings( all_data[['spatial']] ),
  all_data[[c('sample', 'fov', 'annotation', 'cell_id_sample')]]
)

dat_epi_cluster <- sp_epi[[c(str_glue("RNA_snn_res.0.{c(2,4,6,8)}"), "cell_id_sample")]]



p_sp_ct_epi_ls <-
  dat_sp$sample %>% 
  unique() %>%
  as.character() %>%
  setNames(., .) %>% # .[1] %>%
  map(function(sam){
    dat2 <- dat_sp[dat_sp$sample ==  sam, ] 
    dat2 <-  dat2 %>% left_join(dat_epi_cluster, by = "cell_id_sample")  
    
    
    height <- 5
    width <- 5
    
    str_glue("RNA_snn_res.0.{c(2,4,6,8)}") %>%
      setNames(., .) %>%
      map(function(x){
        cell_s <- dat2[[x]] %>% unique()
        cell_cols <- structure(
          names = cell_s,
          ggsci::pal_igv()(length(cell_s))
        )
        
        p <-
          dat2 %>%
          ggplot(aes(x = SPATIAL_1, SPATIAL_2)) +
          geom_tile(
            aes(fill = .data[[x]]),
            height = height, 
            width = width,
            show.legend = c(fill = TRUE, alpha = FALSE)
          ) +
          theme_void() +
          coord_equal() +
          scale_y_reverse() +
          # ggsci::scale_fill_igv() +
          scale_fill_manual(values = cell_cols) +
          theme(panel.background = element_rect(fill = "black", color = NA)) +
          ggtitle(sam)
        p
      })
    
  })


# export 
dir.create(file.path(out_dir, "spatial"))
p_sp_ct_epi_ls %>%
  iwalk(function(p_ls, name){
    pdf(
      file.path(out_dir, "spatial", str_glue("spatial plot of epi cluster - {name}.pdf")),
      height = 8,
      width = 8
    )
    p_ls %>% walk(print)
    dev.off()
  })


# spatial plot of epi subcluster - res0.2 - fixed color -----------
# dat_sp 
# dat_epi_cluster

cell_s_res_0.2 <- dat_epi_cluster$RNA_snn_res.0.2 %>% unique %>% sort
cell_cols_res_0.2 <- structure(
  names = cell_s_res_0.2,
  ggsci::pal_igv()(length(cell_s_res_0.2))
)

p_sp_ct_epi_res0.2_ls <-
  dat_sp$sample %>% 
  unique() %>%
  as.character() %>%
  setNames(., .) %>% # .[1] %>%
  map(function(sam){
    dat2 <- dat_sp[dat_sp$sample ==  sam, ] 
    dat2 <-  dat2 %>% left_join(dat_epi_cluster, by = "cell_id_sample")  
    
    height <- 5
    width <- 5
    
    
    p <-
      dat2 %>%
      ggplot(aes(x = SPATIAL_1, SPATIAL_2)) +
      geom_tile(
        aes(fill = RNA_snn_res.0.2),
        height = height, 
        width = width,
        show.legend = c(fill = TRUE, alpha = FALSE)
      ) +
      theme_void() +
      coord_equal() +
      scale_y_reverse() +
      # ggsci::scale_fill_igv() +
      scale_fill_manual(values = cell_cols_res_0.2) +
      theme(panel.background = element_rect(fill = "black", color = NA)) +
      ggtitle(sam)
    
    p
  })


# export 
pdf(
  file.path(out_dir, 'epi', str_glue("spatial plot of epi cluster - res_0.2.pdf")),
  height = 8,
  width = 8
)
p_sp_ct_epi_res0.2_ls %>% walk(print)
dev.off()



# Epi subcluster not use sketch ---------------------
# 20260610
out_dir <- "sp_subcluster_no_sketch"
dir.create(out_dir)

sp_epi <- all_data[, all_data$annotation_full %in% c("Epithelial_1", "Epithelial_2", "Epithelial proliferative")]
DefaultAssay(sp_epi) <- 'RNA'

sp_epi <- sc_normalize(sp_epi, scale.factor = "default")

sp_epi <- sc_pca(
  sp_epi, 
  sc_out_dir = file.path(out_dir, "epi/sc.1 pca"),
  dir.create.force = TRUE
)

set.seed(1234)
pdf(file.path(out_dir, "epi/sc.1 pca", 'dimheatmap 1-5 - downsample.pdf'), height = 6, width = 10)
DimHeatmap(sp_epi, dims = 1:5, cells = 5000)
dev.off()
pdf(file.path(out_dir, "epi/sc.1 pca", 'dimheatmap 6-10 - downsample.pdf'), height = 6, width = 10)
DimHeatmap(sp_epi, dims = 6:10, cells = 5000)
dev.off()
pdf(file.path(out_dir, "epi/sc.1 pca", 'dimheatmap 11-15 - downsample.pdf'), height = 6, width = 10)
DimHeatmap(sp_epi, dims = 11:15, cells = 5000)
dev.off()
pdf(file.path(out_dir, "epi/sc.1 pca", 'dimheatmap 16-20 - downsample.pdf'), height = 6, width = 10)
DimHeatmap(sp_epi, dims = 16:20, cells = 5000)
dev.off()


ndims_epi <- 10

sp_epi <- sc_dr(
  sp_epi, 
  ndims = ndims_epi,
  return.model = TRUE,
  sc_out_dir = file.path(out_dir, "epi/sc.2 dr"),
  dir.create.force = TRUE
)


sp_epi <- sc_cluster(
  sp_epi, 
  ndims = ndims_epi,
  reduction =  "harmony", # "harmony", "pca"
  sc_out_dir = file.path(out_dir, "epi/sc.3 cluster"),
  dir.create.force = TRUE,
  k.param = 50,
  
  harmony.reduction.name = "harmonyumap",
  resolution = c(0.2, 0.4, 0.6, 0.8)
)


pdf(
  file.path(out_dir, "epi/split dimplot - res0.2.pdf"),
  height = 20,
  width  = 12
)
DimPlot(
  sp_epi,
  group.by = "RNA_snn_res.0.2",
  split.by = "RNA_snn_res.0.2",
  reduction = "harmonyumap",
  ncol = 2
)
dev.off()


pdf(
  file.path(out_dir, "epi/vlnplot - res_0.2.pdf"),
  height = 4,
  width  = 8
)
sp_epi %>%
  VlnPlot(
    features = "total_counts",
    group.by = "RNA_snn_res.0.2",
    pt.size = 0
  ) 
dev.off()



pdf(
  file.path(out_dir, 'epi', 'clustree.pdf'),
  height = 8,
  width = 8
)
library(clustree)
clustree::clustree(
  sp_epi,
  prefix = "RNA_snn_res."
)
dev.off()

saveRDS(sp_epi, file.path(out_dir, 'epi', 'sp_epi.rds'))


epi_markers <- c( 
  'EPCAM',
  "CD3D",
  "CD3E",
  "CD8A",
  
  "squamous cell carcinoma" = "KRT5",
  "squamous cell carcinoma" = "KRT6A",
  
  "squamocolumnar junction" = "KRT7",
  "squamocolumnar junction" = "KRT17",
  
  "columnar epithelium" = "MUC6",
  "columnar epithelium" = "MUC5B",
  
  "EMT" = "COL17A1",
  "EMT" = "POSTN",
  "EMT" = "SNAI2",
  
  "Prof." = "MKI67",
  "Prof." = "TOP2A",
  "Prof." = "UBE2C",
  "Prof." = "CDC20",
  
  "MHC I" = "HLA-A",
  "MHC I" = "HLA-B",
  "MHC I" = "HLA-C"
  
)

sp_epi <- sc_dotplot(
  sp_epi, 
  ndims = ndims_epi,
  marker = epi_markers,
  width = 12,
  sc_out_dir = file.path(out_dir, "epi/sc.4 dotplot"), 
  dir.create.force = TRUE
)

sp_epi <- sc_dotplot(
  sp_epi,
  ndims = ndims_epi,
  marker = major_marker_for_subcluster,
  width = 20,
  sc_out_dir = file.path(out_dir, "epi/sc.4 dotplot - major"),
  dir.create.force = TRUE
)

for (i in str_glue("RNA_snn_res.0.{c(2,4,6,8)}")){
  all_markers_epi <- 
    sc_findmarker(
      sp_epi,
      idents = i, # "RNA_snn_res.0.2",
      test.use = "wilcox",
      logfc.threshold = 0.25,
      min.pct = 0.1,
      only.pos = TRUE,
      # method = "conserved",
      grouping.var = "group"
    )
  write.csv(
    all_markers_epi,
    file = file.path(out_dir, 'epi', str_glue("de gene table  - {i}.csv") ),
    row.names = F
  )
}



sc_findmarker(sp_epi, all_markers = all_markers_epi, sort_by = "pct")


# spatial plot
# dat_sp_subcluster <-
#   Embeddings( all_data[['spatial']] ) %>%
#   merge(
#     sp_epi[[c(str_glue("RNA_snn_res.0.{c(2,4,6,8)}"), 'sample')]],
#     all.x = TRUE,
#     by = "row.names"
#   )
dat_sp_coords <- cbind(Embeddings( all_data[['spatial']] ), all_data[['sample']] )
dat_sp_epi    <- sp_epi[[c(str_glue("RNA_snn_res.0.{c(2,4,6,8)}"), 'sample')]]

p_sp_epi_sub_ls <-
  dat_sp_coords$sample %>% 
  unique() %>%
  as.character() %>%
  setNames(., .) %>% # .[1] %>%
  map(function(sam){
    dat1 <- dat_sp_coords[dat_sp_coords$sample ==  sam, ] 
    
    height <- 5
    width <- 5
    p <-
      dat1 %>%
      
      ggplot(aes(x = SPATIAL_1, SPATIAL_2)) +
      geom_tile(
        fill = "gray90",
        height = height, 
        width = width
      ) +
      theme_void() +
      coord_equal() +
      scale_y_reverse() +
      ggsci::scale_fill_igv(na.translate = T) +
      theme(panel.background = element_rect(fill = "black", color = NA)) 
    
    p_ls <-
      str_glue("RNA_snn_res.0.{c(2,4,6,8)}") %>%
      set_names(., .) %>%
      map(function(x){
        dat2 <- dat_sp_epi[dat_sp_epi$sample ==  sam, ] 
        dat2 <- dat2 %>% merge(dat1, by = "row.names")
        # print(head(dat2))
        p +
          geom_tile(
            data = dat2,
            aes(fill = .data[[x]]),
            height = height, 
            width = width,
            show.legend = c(fill = TRUE, alpha = FALSE)
          ) +
          ggtitle(x)
      })
    p_ls
    
  })

dir.create(file.path(out_dir, "epi/sc.7 spatial"), recursive = TRUE)
p_sp_epi_sub_ls %>%
  iwalk(function(p_ls, name){
    pdf(
      file.path(out_dir, "epi/sc.7 spatial", str_glue("spatial plot - {name}.pdf") ),
      height = 8,
      width = 8
    )
    p_ls %>% walk(print)
    dev.off()
  })




# #############33333
# dir.create(file.path(out_dir, "epi/sc.6 heatmap"), recursive = TRUE)
# pdf(
#   file.path(out_dir, "epi/sc.6 heatmap", 'heatmap.pdf'),
#   height = 4,
#   width  = 8
# )
# sc_mark_heatmap(
#   sp_epi,
#   markers = all_markers_epi,
#   idents = "RNA_snn_res.0.2",
#   p_val_adj = 0.05,
#   pct.1 = 0.15,
#   n_head = 5,
#   row_cluster = FALSE
# )
# dev.off()
# 
# 
# sc_mark_heatmap(
#   sp_epi,
#   markers = all_markers_epi,
#   idents = "RNA_snn_res.0.2",
#   p_val_adj = 0.05,
#   pct.1 = 0.15,
#   n_head = 5,
#   row_cluster = FALSE,
#   scale = "none",
#   log_trans = TRUE
# )
# 
# epi_cell_cluster <- c(
#   '0'  =  "Proximal mature enterocytes",
#   '1'  =  "Peneth cell",
#   '2'  =  "Transit amplifying cell",
#   '3'  =  "Others",
#   '4'  =  "Enterocytes",
#   '5'  =  "Others",
#   '6'  =  "Others",
#   '7'  =  "CD8T",
#   '8'  =  "Transit amplifying cell",
#   '9'  =  "Others",
#   '10'  = "Others",
#   '11'  = "Others",
#   '12'  = "Others",
#   '13'  = "Proximal mature enterocytes",
#   '14'  = "Others"
# 
# )
# 
# 
# 
# Idents(sp_epi) <- "RNA_snn_res.0.4"
# sp_epi <- RenameIdents(sp_epi, epi_cell_cluster)
# sp_epi$annotation_epi <- Idents(sp_epi)
# 
# dir.create(file.path(out_dir, "epi/sc.5 annotation"))
# pdf(file.path(out_dir, 'epi/sc.5 annotation', str_glue("Dimplot dim{ndims_epi} epi annotation .pdf")), height = 6, width = 8)
# DimPlot(sp_epi, group.by = "annotation_epi", label = TRUE, reduction = "umap") + NoLegend()
# DimPlot(sp_epi, group.by = "annotation_epi", label = TRUE, reduction = "harmonyumap") + NoLegend()
# dev.off()

# 
# ##  project  data --------------------
# DefaultAssay(sp_epi) <- "sketch"
# # sp_epi[['full.harmonyumap']] <- NULL
# sp_epi <-
#   sp_epi %>%
#   ProjectData(
#     assay = "RNA",
#     full.reduction = "harmony.full",
#     
#     sketched.assay = "sketch",
#     sketched.reduction = "harmony", # "pca",
#     
#     umap.model = "harmonyumap",
#     dims = 1:ndims_epi,
#     
#     refdata = list(
#       # annotation_full = "annotation",
#       RNA_snn_res.0.2 = 'sketch_snn_res.0.2',
#       RNA_snn_res.0.4 = 'sketch_snn_res.0.4',
#       RNA_snn_res.0.6 = 'sketch_snn_res.0.6',
#       RNA_snn_res.0.8 = 'sketch_snn_res.0.8'
#     )
#   )
# 
# Idents(sp_epi) <- 'RNA_snn_res.0.2' # 'annotation_full'
# DefaultAssay(sp_epi) <- "RNA"

# pdf(file.path(out_dir, 'epi', str_glue("Dimplot dim{ndims_epi} annotation - projected.pdf")), height = 6, width = 8)
# DimPlot(sp_epi, group.by = "annotation_full", label = TRUE, reduction = 'full.harmonyumap') + NoLegend()
# dev.off()

pdf(file.path(out_dir, 'epi', str_glue("Dimplot dim{ndims_epi} annotation - projected.pdf")), height = 6, width = 8)
DimPlot(sp_epi, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = 'full.harmonyumap') + NoLegend()
DimPlot(sp_epi, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = 'full.harmonyumap') + NoLegend()
DimPlot(sp_epi, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = 'full.harmonyumap') + NoLegend()
DimPlot(sp_epi, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = 'full.harmonyumap') + NoLegend()
dev.off()

# saveRDS(sp_epi, file.path(out_dir, 'sp_epi_projected.rds'))


# DimPlot(sp_epi, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = 'full.harmonyumap') + NoLegend()


## spatial plot  ----------
dat_sp <- cbind(
  Embeddings( all_data[['spatial']] ),
  all_data[[c('sample', 'fov', 'annotation', 'cell_id_sample')]]
)

dat_epi_cluster <-
  sp_epi[[c("RNA_snn_res.0.2", "cell_id_sample")]] %>%
  mutate(RNA_snn_res.0.2 = as.character(RNA_snn_res.0.2)) %>% 
  dplyr::filter(! RNA_snn_res.0.2 %in% c("7", "8", "9"))


cell_s_res_0.2 <- dat_epi_cluster$RNA_snn_res.0.2 %>% unique %>% sort
cell_cols_res_0.2 <- structure(
  names = cell_s_res_0.2,
  ggsci::pal_igv()(length(cell_s_res_0.2))
)

p_sp_ct_epi_res0.2_ls <-
  dat_sp$sample %>% 
  unique() %>%
  as.character() %>%
  setNames(., .) %>% # .[1] %>%
  map(function(sam){
    dat2 <- dat_sp[dat_sp$sample ==  sam, ] 
    dat2 <-  dat2 %>% left_join(dat_epi_cluster, by = "cell_id_sample")  
    
    height <- 5
    width <- 5
    
    
    p <-
      dat2 %>%
      ggplot(aes(x = SPATIAL_1, SPATIAL_2)) +
      geom_tile(
        aes(fill = RNA_snn_res.0.2),
        height = height, 
        width = width,
        show.legend = c(fill = TRUE, alpha = FALSE)
      ) +
      theme_void() +
      coord_equal() +
      scale_y_reverse() +
      # ggsci::scale_fill_igv() +
      scale_fill_manual(values = cell_cols_res_0.2) +
      theme(panel.background = element_rect(fill = "black", color = NA)) +
      ggtitle(sam)
    
    p
  })


# export 
pdf(
  file.path(out_dir, 'epi', str_glue("spatial plot of epi cluster - res_0.2.pdf")),
  height = 8,
  width = 8
)
p_sp_ct_epi_res0.2_ls %>% walk(print)
dev.off()




