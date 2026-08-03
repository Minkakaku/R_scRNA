library(tidyverse)
library(Seurat)
library(patchwork)

source("subcluster/sc_analysis_tools.R")
source("subcluster/marker_for_subcluster.R")

all_data <- readRDS("all_data_full.rds")

out_dir <- "sp_subcluster"
dir.create(out_dir)


all_data$annotation %>% table
# .
# B Endothelial  Epithelial     Myeloid      Others     P_2_pos     Stromal           T 
# 103113       85375      454025      102258       93538       13482      136460       95523


# Epi subcluster ----------------
sp_epi <- all_data[, all_data$annotation == "Epithelial"]

DefaultAssay(sp_epi) <- 'RNA'

# sp_epi <- sc_normalize(sp_epi, scale.factor = "median")
sp_epi <- sc_normalize(sp_epi, scale.factor = "default")
####
sp_epi <- ScaleData(
  JoinLayers(sp_epi), 
  vars.to.regress = 'nCount_RNA',
  # features = de_marker_sc_ls %>% unlist() %>% .[. %in% rownames(sp_epi)]
)
####
sp_epi <- sc_pca(
  sp_epi, 
  sc_out_dir = file.path(out_dir, "epi/sc.1 pca"),
  dir.create.force = TRUE
)

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

epi_markers <- c(
  'EPCAM',
  "CD3D",
  "CD3E",
  "CD8A",
  
  'PECAM1',
  'FAP',
  'CCDC80',
  'CD68',
  'SPP1',
  'VCAN'
  
  # Enterocyte  = "SLC26A3",
  # 
  # Goblet = "MUC2",
  # Goblet = "TFF3",
  
  
  # Enteroendocrine = "CHGA",
  # Enteroendocrine = "TPH1",
  # Enteroendocrine = "NEUROD1",
  # 
  # 'Paneth cell' = 'DEFA6',
  # 'Paneth cell' = 'IL9',
  # 
  # 'Transit amplifying cell' = 'CCNA2',
  # 'Transit amplifying cell' = 'HMGB3',
  # 'Transit amplifying cell' = 'MCM5',
  # 'Transit amplifying cell' = 'PCNA',
  # 'Transit amplifying cell' = 'TOP2A',
  # 
  # 'Proximal Mature Enterocytes' = 'FABP2',
  # 'Proximal Mature Enterocytes' = 'ACE',
  # 'Proximal Mature Enterocytes' = 'FABP6',
  # 'Proximal Mature Enterocytes' = 'APOC3',
  # 
  # 'Proximal Enterocytes' = 'HEBP1',
  # 'Proximal Enterocytes' = 'CCL25'
)

sp_epi <- sc_dotplot(
  sp_epi, 
  ndims = ndims_epi,
  marker = epi_markers,
  sc_out_dir = file.path(out_dir, "epi/sc.4 dotplot"), 
  dir.create.force = TRUE
)
sp_epi <- sc_dotplot(
  sp_epi, 
  ndims = ndims_epi,
  marker = major_marker_for_subcluster,
  sc_out_dir = file.path(out_dir, "epi/sc.4 dotplot - major"), 
  dir.create.force = TRUE
)

all_markers_epi <- 
  sc_findmarker(
    sp_epi,
    idents = "RNA_snn_res.0.2",
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    only.pos = TRUE,
    # method = "conserved",
    grouping.var = "group"
  )

sc_findmarker(sp_epi, all_markers = all_markers_epi, sort_by = "pct")

# dir.create()
write.csv(
  all_markers_epi,
  file = file.path(out_dir, "epi/de markers - res 0.2.csv"),
  row.names = FALSE
)


dir.create(file.path(out_dir, "epi/sc.6 heatmap"), recursive = TRUE)
pdf(
  file.path(out_dir, "epi/sc.6 heatmap", 'heatmap.pdf'),
  height = 4,
  width  = 8
)
sc_mark_heatmap(
  sp_epi,
  markers = all_markers_epi,
  idents = "RNA_snn_res.0.4",
  p_val_adj = 0.05,
  pct.1 = 0.15,
  n_head = 5,
  row_cluster = FALSE
)
dev.off()


sc_mark_heatmap(
  sp_epi,
  markers = all_markers_epi,
  idents = "RNA_snn_res.0.4",
  p_val_adj = 0.05,
  pct.1 = 0.15,
  n_head = 5,
  row_cluster = FALSE,
  scale = "none",
  log_trans = TRUE
)

epi_cell_cluster <- c(
  '0'  =  "Proximal mature enterocytes",
  '1'  =  "Peneth cell",
  '2'  =  "Transit amplifying cell",
  '3'  =  "Others",
  '4'  =  "Enterocytes",
  '5'  =  "Others",
  '6'  =  "Others",
  '7'  =  "CD8T",
  '8'  =  "Transit amplifying cell",
  '9'  =  "Others",
  '10'  = "Others",
  '11'  = "Others",
  '12'  = "Others",
  '13'  = "Proximal mature enterocytes",
  '14'  = "Others"
  
)



Idents(sp_epi) <- "RNA_snn_res.0.4"
sp_epi <- RenameIdents(sp_epi, epi_cell_cluster)
sp_epi$annotation_epi <- Idents(sp_epi)

dir.create(file.path(out_dir, "epi/sc.5 annotation"))
pdf(file.path(out_dir, 'epi/sc.5 annotation', str_glue("Dimplot dim{ndims_epi} epi annotation .pdf")), height = 6, width = 8)
DimPlot(sp_epi, group.by = "annotation_epi", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(sp_epi, group.by = "annotation_epi", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()


# # project  data
# DefaultAssay(sp_epi) <- "sketch"
# sp_epi[['full.harmonyumap']] <- NULL
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
#       annotation_full = "annotation",
#       RNA_snn_res.0.2 = 'sketch_snn_res.0.2',
#       RNA_snn_res.0.4 = 'sketch_snn_res.0.4',
#       RNA_snn_res.0.6 = 'sketch_snn_res.0.6',
#       RNA_snn_res.0.8 = 'sketch_snn_res.0.8'
#     )
#   )
# 
# Idents(sp_epi) <- 'annotation_full'
# DefaultAssay(sp_epi) <- "RNA"
# 
# pdf(file.path(out_dir, 'epi', str_glue("Dimplot dim{ndims_epi} annotation - projected.pdf")), height = 6, width = 8)
# DimPlot(sp_epi, group.by = "annotation_full", label = TRUE, reduction = 'full.harmonyumap') + NoLegend()
# dev.off()

# saveRDS(sp_epi, file.path(out_dir, 'sp_epi_lognorm_hvg.rds'))

## annotation from single cell -------
sc_epi <- readRDS("output/singlecell/anno_epi.rds")
sc_epi <- NormalizeData(sc_epi)

Idents(sc_epi) <- "annotation"
all_markers_epi_sc <- FindAllMarkers(sc_epi)

de_marker_sc <- 
  all_markers_epi_sc %>%
  group_by(cluster) %>%
  dplyr::filter(
    avg_log2FC > 0.25,
    pct.1 > 0.1,
    pct.1 - pct.2 > 0.1,
    p_val_adj < 0.05
  ) 

# all de 
de_marker_sc_ls <- 
  de_marker_sc %>%
  arrange(desc(avg_log2FC), .by_group = TRUE) %>%
  with(split(gene, cluster))

de_marker_sc_ls %>%
  map(~ list(.x, .x[.x %in% rownames(sp_epi)]))  %>%
  # map(2)
  map(~ length(.x[[2]]) / length(.x[[1]]))
  
# top de
de_marker_sc_top <-
  de_marker_sc %>%
  # arrange(pct.1 - pct.2, .by_group = TRUE) %>%
  arrange(desc(avg_log2FC), .by_group = TRUE) %>%
  slice_head(n = 10) %>%
  with(structure(
    gene,
    names = as.character(cluster)
  ))
  
DotPlot(
  sc_epi,
  features = de_marker_sc_top,
  group.by = 'annotation'
)  +
  RotatedAxis()
DotPlot(
  sp_epi,
  features = de_marker_sc_top,
  group.by = 'RNA_snn_res.0.2'
)  +
  RotatedAxis()
  
DotPlot(
  sc_epi,
  features = epi_top5,
  group.by = 'annotation'
)  +
  RotatedAxis()


# add module score
sp_epi_score <- 
  sp_epi %>%
  JoinLayers() %>%
  AddModuleScore(features = de_marker_sc_ls %>% map(head, 100))
path_idx <- str_detect(colnames(sp_epi_score[[]]), "^Cluster")
colnames(sp_epi_score[[]])[path_idx] <- names(de_marker_sc_ls)

DotPlot(
  sp_epi_score, # [, ! sp_epi_score$RNA_snn_res.0.2 %in% c(4, 5)],
  features = names(de_marker_sc_ls),
  group.by = 'RNA_snn_res.0.2'
)  +
  RotatedAxis()

DotPlot(
  sp_epi_score,
  features = names(de_marker_sc_ls),
  group.by = 'RNA_snn_res.0.4'
)  +
  RotatedAxis()
DotPlot(
  sp_epi_score,
  features = names(de_marker_sc_ls),
  group.by = 'RNA_snn_res.0.6'
)  +
  RotatedAxis()


sp_epi_score[[]] %>%
  group_by(sample, RNA_snn_res.0.4) %>%
  summarise(count = n()) %>%
  group_by(sample) %>%
  mutate(freq = count / sum(count)) %>%
  
  ggplot(aes(x = sample, y = freq, fill = RNA_snn_res.0.4)) +
  geom_col() +
  ggsci::scale_fill_igv() 


sp_epi_score[[ ]] %>%
  group_by(group, RNA_snn_res.0.2) %>%
  summarise(count = n()) %>%
  mutate(freq = count / sum(count)) %>%
  ggplot(aes(x = group, y = freq, fill = RNA_snn_res.0.2)) +
  geom_col() +
  theme_bw() +
  ggsci::scale_fill_aaas()  

sp_epi[[]] %>%
  group_by(sample, RNA_snn_res.0.2) %>%
  summarise(count = n()) %>%
  mutate(freq = count / sum(count))  %>%
  
  ggplot(aes(x = sample, y = freq, fill = RNA_snn_res.0.2)) +
  geom_col()+
  ggsci::scale_fill_igv()
sp_epi[[]] %>%
  group_by(group, RNA_snn_res.0.2) %>%
  summarise(count = n()) %>%
  mutate(freq = count / sum(count))  %>%
  
  ggplot(aes(x = group, y = freq, fill = RNA_snn_res.0.2)) +
  geom_col() +
  ggsci::scale_fill_igv()


sc_epi[[]] %>%
  group_by(sample, annotation) %>%
  summarise(count = n()) %>%
  mutate(freq = count / sum(count))  %>%
  
  ggplot(aes(x = sample, y = freq, fill = annotation)) +
  geom_col()
  
  
sc_epi[[]] %>%
  group_by(group, annotation) %>%
  summarise(count = n()) %>%
  mutate(freq = count / sum(count))  %>%
  
  ggplot(aes(x = group, y = freq, fill = annotation)) +
  geom_col()

sc_epi[[c('sample', 'group')]] %>%  distinct()




sp_epi_score %>% head

FeaturePlot(
  sp_epi_score,
  feature = "tumor_1",
  reduction = "harmonyumap"
)  +
  scale_color

DimPlot(
  sp_epi_score,
  group.by = "RNA_snn_res.0.2",
  reduction = "harmonyumap"
)



## annotation to single cell -------
Idents(sp_epi) <- "RNA_snn_res.0.2"
all_markers_epi_sp <- FindAllMarkers(sp_epi)

de_marker_sp <- 
  all_markers_epi_sp %>%
  group_by(cluster) %>%
  dplyr::filter(
    avg_log2FC > 0.25,
    pct.1 > 0.1,
    # pct.1 - pct.2 > 0.1,
    p_val_adj < 0.05
  ) 



# all de 
de_marker_sp_ls <- 
  de_marker_sp %>%
  arrange(desc(avg_log2FC), .by_group = TRUE) %>%
  with(split(gene, cluster))


# add module score
sc_epi_score <- 
  sc_epi %>%
  # JoinLayers() %>%
  AddModuleScore(features = de_marker_sp_ls %>% map(head, 100))

sc_path_idx <- str_detect(colnames(sc_epi_score[[]]), "^Cluster")
colnames(sc_epi_score[[]])[sc_path_idx] <- names(de_marker_sp_ls)

DotPlot(
  sc_epi_score,
  features = names(de_marker_sp_ls),
  group.by = 'cell_type'
)  +
  RotatedAxis()





###########################################################
# #### FORMAL CODE  #### 
############################################################
# 20260409 using v4 and de_gene to match cell cluster 
## using epi - v4 default normalize and sc de gene ---------- 
## sp_epi <---- sp_subcluster/sp_epi_lognorm_hvg.rds (疑似是这个)
## sp_epi <---- sp_subcluster/sp_subcluster/sp_epi_used.rds (可定是这个)


### spatial plot of all celltypes - res0.2 --------------
cells_position_ls <-
  cbind(
    all_data[['sample']],
    Embeddings(all_data[['spatial']])
  ) %>%
  rownames_to_column('cell_id') %>%
  group_by(sample) %>%
  nest() %>%
  deframe()

annotation_col <- 'RNA_snn_res.0.2'
plot_dat_sub <- 
  FetchData(
    sp_epi, 
    vars = c(annotation_col)
  ) %>%
  rownames_to_column('cell_id')

p_sp_ct_ls <-
  cells_position_ls %>% # .[1] %>%
  imap(function(coords, sam){
    dat2 <- coords # %>% left_join(plot_dat_sub, by = "cell_id")
    dat3 <- coords %>% inner_join(plot_dat_sub, by = "cell_id")
    
    height <- 5
    width <- 5
    p <-
      dat2 %>%
      ggplot(aes(x = SPATIAL_1, SPATIAL_2)) +
      geom_tile(
        height = height, 
        width = width,
        fill = "gray90"
      ) +
      geom_tile(
        data = dat3,
        aes(fill = .data[[annotation_col]]),
        height = height, 
        width = width,
        show.legend = c(fill = TRUE, alpha = FALSE)
      ) +
      theme_void() +
      coord_equal() +
      scale_y_reverse() +
      ggsci::scale_fill_igv() +
      # scale_fill_manual(values = allcolour) +
      theme(panel.background = element_rect(fill = "black", color = NA)) +
      ggtitle(sam)
    p
    
  })


pdf(
  file.path(out_dir, "epi - v4 default normalize and sc de gene", "spatial plot of epi cluster in each sample.pdf"),
  height = 8,
  width = 8
)
p_sp_ct_ls %>% walk(print)
dev.off()

### spatial plot of all celltypes - res0.4 --------------
annotation_col <- 'RNA_snn_res.0.4'
plot_dat_sub <- 
  FetchData(
    sp_epi, 
    vars = c(annotation_col)
  ) %>%
  rownames_to_column('cell_id')

p_sp_ct_ls_res0.4 <-
  cells_position_ls %>% # .[1] %>%
  imap(function(coords, sam){
    dat2 <- coords # %>% left_join(plot_dat_sub, by = "cell_id")
    dat3 <- coords %>% inner_join(plot_dat_sub, by = "cell_id")
    
    height <- 5
    width <- 5
    p <-
      dat2 %>%
      ggplot(aes(x = SPATIAL_1, SPATIAL_2)) +
      geom_tile(
        height = height, 
        width = width,
        fill = "gray90"
      ) +
      geom_tile(
        data = dat3,
        aes(fill = .data[[annotation_col]]),
        height = height, 
        width = width,
        show.legend = c(fill = TRUE, alpha = FALSE)
      ) +
      theme_void() +
      coord_equal() +
      scale_y_reverse() +
      ggsci::scale_fill_futurama() +
      # scale_fill_manual(values = allcolour) +
      theme(panel.background = element_rect(fill = "black", color = NA)) +
      ggtitle(sam)
    p
    
  })


pdf(
  file.path(out_dir, "epi - v4 default normalize and sc de gene", "spatial plot of epi cluster in each sample - res0.4.pdf"),
  height = 8,
  width = 8
)
p_sp_ct_ls_res0.4 %>% walk(print)
dev.off()



table(
  sp_epi$RNA_snn_res.0.2,
  sp_epi$RNA_snn_res.0.4
)
### cell props ---------
p_col_grp <-
  sp_epi_score[[ ]] %>%
  group_by(group, RNA_snn_res.0.2) %>%
  summarise(count = n()) %>%
  mutate(freq = count / sum(count)) %>%
  ggplot(aes(x = group, y = freq, fill = RNA_snn_res.0.2)) +
  geom_col() +
  theme_bw() +
  ggsci::scale_fill_aaas()  

p_col_sample <-
  sp_epi[[]] %>%
  group_by(sample, RNA_snn_res.0.2) %>%
  summarise(count = n()) %>%
  mutate(freq = count / sum(count))  %>%
  
  ggplot(aes(x = sample, y = freq, fill = RNA_snn_res.0.2)) +
  geom_col()+
  ggsci::scale_fill_igv()

# export
pdf(
  file.path(out_dir, "epi - v4 default normalize and sc de gene", "cell prop of epi subcluster - group.pdf"),
  height = 6,
  width = 4
)
p_col_grp %>% print
dev.off()

pdf(
  file.path(out_dir, "epi - v4 default normalize and sc de gene", "cell prop of epi subcluster - sample.pdf"),
  height = 6,
  width = 8
)
p_col_sample %>% print
dev.off()

# res 0.4
p_col_grp_res0.4 <-
  sp_epi[[ ]] %>%
  group_by(group, RNA_snn_res.0.4) %>%
  summarise(count = n()) %>%
  mutate(freq = count / sum(count)) %>%
  ggplot(aes(x = group, y = freq, fill = RNA_snn_res.0.4)) +
  geom_col() +
  theme_bw() +
  ggsci::scale_fill_igv() 

p_col_sample_res0.4 <-
  sp_epi[[]] %>%
  group_by(sample, RNA_snn_res.0.4) %>%
  summarise(count = n()) %>%
  mutate(freq = count / sum(count))  %>%
  
  ggplot(aes(x = sample, y = freq, fill = RNA_snn_res.0.4)) +
  geom_col()+
  ggsci::scale_fill_igv()
# export
pdf(
  file.path(out_dir, "epi - v4 default normalize and sc de gene", "cell prop of epi subcluster - group - res0.4.pdf"),
  height = 6,
  width = 4
)
p_col_grp_res0.4 %>% print
dev.off()

pdf(
  file.path(out_dir, "epi - v4 default normalize and sc de gene", "cell prop of epi subcluster - sample - res0.4.pdf"),
  height = 6,
  width = 8
)
p_col_sample_res0.4 %>% print
dev.off()

### sc and sp link - tumor_1 ----------------------
sp_epi <- readRDS("~/sp_subcluster/sp_epi_lognorm.rds")

c(
  'AMTN',
  'KRT14',
  'SPINK6',
  'LGR5',
  'OLFM2',
  'SFRP1',
  'IRX1',
  'APCDD1',
  'SDC2',
  'BGN',
  'CSPG4',
  'RBP4',
  'SOSTDC1',
  'IL27RA',
  'LDLRAD4',
  'GPC4',
  'PCSK9',
  'IGFBP5',
  'FZD8',
  'CXCL14'
) %>% .[. %in% rownames(sp_epi)]


p <- DotPlot(
  sp_epi[, !sp_epi$RNA_snn_res.0.2 %in% c('3', '4', '5')],
  features = c("LGR5", "OLFM2", "SFRP1", "IRX1", "APCDD1", "SDC2", "CSPG4", 
               "SOSTDC1", "LDLRAD4", "PCSK9"),
  group.by = "RNA_snn_res.0.2"
)

pdf(
  file.path(out_dir, "epi - v4 default normalize and sc de gene", "tumor_1 marker gene in spatial epi cluster.pdf"),
  height = 6,
  width = 8
)
p %>% print
dev.off()



# plot_dat <-
#   sp_epi %>%
#   FetchData(
#     vars = c("LGR5", "OLFM2", "SFRP1", "IRX1", "APCDD1", "SDC2", "CSPG4", 
#              "SOSTDC1", "LDLRAD4", "PCSK9", 'RNA_snn_res.0.2')
#   )
# 
# plot_dat2 <-
#   plot_dat %>%
#   group_by(RNA_snn_res.0.2) %>%
#   dplyr::filter(!RNA_snn_res.0.2 %in% c('3', '4', '5')) %>%
#   
#   summarise(across(
#     everything(), 
#     list(
#       avg.exp = function(x) mean(expm1(x[x != 0])),
#       pct.exp = function(x) sum(x > 0)/length(x) * 100
#     )
#   )) %>%
#   pivot_longer(
#     -RNA_snn_res.0.2,
#     names_to = "name",
#     values_to = "value"
#   ) %>%
#   separate(
#     name,
#     sep = "_",
#     into = c("gene", "var")
#   ) %>%
#   pivot_wider(names_from = "var", values_from = "value")
# 
# 
# plot_dat2 %>%
#   ggplot(aes(x = gene, y = RNA_snn_res.0.2, size = pct.exp, color = avg.exp )) +
#   geom_point() +
#   theme_classic() +
#   scale_color_gradientn(colours = c("gray90", 'orangered', "red"))
# 
# 
# plot_dat3 <-
#   plot_dat %>%
#   group_by(RNA_snn_res.0.2) %>%
#   dplyr::filter(!RNA_snn_res.0.2 %in% c('3', '4', '5')) %>%
#   summarise(across(
#     everything(), 
#     list(
#       avg.exp = function(x) mean(expm1(x)),
#       pct.exp = function(x) sum(x > 0)/length(x) * 100
#     )
#   )) %>%
#   pivot_longer(
#     -RNA_snn_res.0.2,
#     names_to = "name",
#     values_to = "value"
#   ) %>%
#   separate(
#     name,
#     sep = "_",
#     into = c("gene", "var")
#   ) %>%
#   pivot_wider(names_from = "var", values_from = "value")
# 
# 
# plot_dat3 %>%
#   ggplot(aes(x = gene, y = RNA_snn_res.0.2, size = pct.exp, color = avg.exp )) +
#   geom_point() +
#   theme_classic() +
#   scale_color_gradientn(colours = c("gray90", 'orangered', "red"))



pdf(
  file.path("sp_subcluster/epi/marker gene of tumor_1.pdf"),
  height = 6,
  width  = 8
)
FeaturePlot(
  sp_epi,
  features = c("LGR5", "OLFM2", "SFRP1", "IRX1", "APCDD1", "SDC2", "CSPG4", 
               "SOSTDC1", "LDLRAD4", "PCSK9"), # [1:3],
  reduction = "harmonyumap",
  combine = FALSE,
  order = TRUE
)
dev.off()

### pathway validation of tumor 1 --------
# causing tumor_1 have metastatic invasion proliferation and cellcycle,
# so check these pathway scores in spatial data to validate the cluster0 is tumor_1

Metastatic <- c("ANXA1","ANXA2","CD44","CDH2","CXCL8","EZH2","FOXM1","HIF1A","SNAI2","SOX4")
Invasion <- c("LAMC1","UBE2V2","COL1A1","COL5A1","COL5A2","TGFBI","NNMT","CBX1","H2AFZ","YBX1")
# [DNA replication]
Proliferation <- c(
  "RNASEH2A", "POLD3", "DNA2", "FEN1", "POLA2", "RNASEH1", "RPA4", 
  "LIG1", "MCM2", "MCM3", "MCM4", "MCM5", "MCM6", "MCM7", "PCNA", 
  "POLE3", "POLA1", "POLD1", "POLD2", "POLE", "POLE2", "PRIM1", 
  "PRIM2", "POLE4", "POLD4", "RFC1", "RFC2", "RFC3", "RFC4", "RFC5", 
  "RPA1", "RPA2", "RPA3", "SSBP1", "RNASEH2B", "RNASEH2C"
)
CellCycle <- c(
  "CDK2", "CDK4", "CDK6", "CDK7", "CDKN1A", "CDKN1B", "STAG1", 
  "CDKN1C", "CDKN2A", "CDKN2B", "CDKN2C", "CDKN2D", "ANAPC10", 
  "NDC80", "MAD2L2", "STAG2", "PTTG2", "GADD45G", "DBF4", "YWHAQ", 
  "CHEK1", "CHEK2", "CDCA5", "ESCO1", "ANAPC16", "CREBBP", "SGO1", 
  "ESCO2", "GADD45A", "DDX11", "CDC14C", "E2F1", "E2F2", "E2F3", 
  "E2F4", "E2F5", "EP300", "PDS5B", "WAPL", "PDS5A", "MAU2", "ORC6", 
  "ORC3", "CDC26", "ABL1", "NIPBL", "ANAPC13", "ANAPC15", "FBXO5", 
  "MTBP", "SMC1B", "SFN", "GSK3B", "ANAPC2", "ANAPC4", "HDAC1", 
  "HDAC2", "MAD2L1", "SMAD2", "SMAD3", "SMAD4", "MCM2", "MCM3", 
  "MCM4", "MCM5", "MCM6", "MCM7", "MDM2", "MYC", "GADD45B", "ATM", 
  "WEE2", "ORC1", "ORC2", "ORC4", "ORC5", "PCNA", "FZR1", "ANAPC5", 
  "ANAPC7", "ANAPC11", "PLK1", "ATR", "ATRX", "PPP2CA", "PPP2CB", 
  "PPP2R1A", "PPP2R1B", "PPP2R5A", "PPP2R5B", "PPP2R5C", "PPP2R5D", 
  "PPP2R5E", "HDAC8", "PRKDC", "KNL1", "RAD21", "RB1", "RBL1", 
  "RBL2", "CCND1", "ANAPC1", "SKP1", "SKP2", "BUB1", "BUB1B", "TFDP1", 
  "TFDP2", "TGFB1", "TGFB2", "TGFB3", "TP53", "TTK", "WEE1", "YWHAB", 
  "YWHAE", "YWHAG", "YWHAH", "YWHAZ", "ZBTB17", "DBF4B", "CDT1", 
  "SMC1A", "CDC7", "CDC45", "MAD1L1", "CUL1", "CCNB3", "CDC14B", 
  "CDC14A", "CDC23", "CDC16", "CCNA2", "CCNA1", "CCNB1", "CCND2", 
  "CCND3", "CCNE1", "CCNH", "TICRR", "PKMYT1", "SMC3", "CCNB2", 
  "CCNE2", "BUB3", "AURKB", "PTTG1", "TRIP13", "MAD2L1BP", "ESPL1", 
  "CDK1", "CDC6", "CDC20", "CDC25A", "CDC25B", "CDC25C", "CDC27", 
  "RBX1"
)

path_to_validation_tumor1 <-
  list(
    Metastatic = Metastatic,
    Invasion = Invasion,
    Proliferation = Proliferation,
    CellCycle = CellCycle
  )


sp_epi_score2 <- 
  sp_epi %>%
  AddModuleScore(features = path_to_validation_tumor1)

path_idx_validation <- str_detect(colnames(sp_epi_score2[[]]), "^Cluster")
colnames(sp_epi_score2[[]])[path_idx_validation] <- names(path_to_validation_tumor1)

p_vln_path_validation <-
  sp_epi_score2 %>%
  .[, sp_epi_score2$RNA_snn_res.0.2 %in% as.character(0:2)] %>% 
  VlnPlot(
    group.by = "RNA_snn_res.0.2",
    features = names(path_to_validation_tumor1),
    pt.size = 0
  ) &
  ggsci::scale_fill_igv()

# export 
pdf(
  file.path("sp_subcluster/epi - v4 default normalize and sc de gene/vlnplot of some pahway of tumor_1.pdf"),
  height = 6,
  width  = 12
)
p_vln_path_validation %>% print()
dev.off()

### sc and sp link - normal ----------------------
all_markers_epi_by_scanpy <- read.table("all_cluster_marker.xls", header = TRUE)
sc_cluster_label <- c(
  'tumor_1' = '0', 
  'tumor_2' = '1', 
  'tumor_3' = '2', 
  'tumor_4' = '3', 
  'tumor_5' = '4', 
  'normal' = '5'
)

de_de_sc_scanpy <-
  all_markers_epi_by_scanpy %>%
  group_by(cluster) %>%
  dplyr::filter(logfoldchanges > 0.25, pvals_adj < 0.05, pct_nz_group - pct_nz_reference > 0.1) %>%
  dplyr::slice_max(logfoldchanges, n = 50) %>%
  mutate(cluster = fct_recode(as.character(cluster), !!! sc_cluster_label)) %>% 
  with(split(Symbol,  cluster))


p_dot_de_gene_sc <-
  de_de_sc_scanpy %>%
  imap(function(gene, name){
    DotPlot(
      sp_epi[, !sp_epi$RNA_snn_res.0.2 %in% c('3', '4', '5')],
      features = gene,
      group.by = "RNA_snn_res.0.2"
    ) +
    ggtitle(name) +
    RotatedAxis()
  })

pdf(
  file.path(out_dir, "epi - v4 default normalize and sc de gene", "single cell marker gene dotplot in spatial epi cluster.pdf"),
  height = 5,
  width = 8
)
p_dot_de_gene_sc %>% walk(print)
dev.off()

# marker gene plot of normal 
pdf(
  file.path("sp_subcluster/epi - v4 default normalize and sc de gene/marker gene of normal.pdf"),
  height = 6,
  width  = 8
)
FeaturePlot(
  sp_epi[, !sp_epi$RNA_snn_res.0.2 %in% c('3', '4', '5')],
  features = de_de_sc_scanpy$normal,
  reduction = "harmonyumap",
  combine = FALSE,
  order = TRUE
)
dev.off()


# res 0.4
p_dot_de_gene_sc_res0.4 <-
  de_de_sc_scanpy %>%
  imap(function(gene, name){
    DotPlot(
      sp_epi[, sp_epi$RNA_snn_res.0.4 %in% as.character(0:5)],
      features = gene,
      group.by = "RNA_snn_res.0.2"
    ) +
    ggtitle(name) +
    RotatedAxis()
  })
pdf(
  file.path(out_dir, "epi - v4 default normalize and sc de gene", "single cell marker gene dotplot in spatial epi cluster - res0.4.pdf"),
  height = 6,
  width = 8
)
p_dot_de_gene_sc_res0.4 %>% walk(print)
dev.off()

### cell annotation ---------
sp_epi2 <- sp_epi[, !sp_epi$RNA_snn_res.0.2 %in% c('3', '4', '5')]
cell_annotaiton <- c(
  '0' = "tumor_1",
  '1' = "tumor_others",
  '2' = "normal"
)
Idents(sp_epi2) <- 'RNA_snn_res.0.2'
sp_epi2 <- sp_epi2 %>% RenameIdents(cell_annotaiton) 
sp_epi2$annotation_epi <- Idents(sp_epi2)

saveRDS(sp_epi, "sp_subcluster/sp_epi_used.rds")
saveRDS(sp_epi2, "sp_subcluster/sp_epi2.rds")



### spatial plot of all celltypes - annotation --------------
cells_position_ls <-
  cbind(
    all_data[['sample']],
    Embeddings(all_data[['spatial']])
  ) %>%
  rownames_to_column('cell_id') %>%
  group_by(sample) %>%
  nest() %>%
  deframe()

annotation_col2 <- 'annotation_epi'
plot_dat_sub2 <- 
  FetchData(
    sp_epi2, 
    vars = c(annotation_col2)
  ) %>%
  rownames_to_column('cell_id')

p_sp_ct_ls2 <-
  cells_position_ls %>% # .[1] %>%
  imap(function(coords, sam){
    dat2 <- coords # %>% left_join(plot_dat_sub, by = "cell_id")
    dat3 <- coords %>% inner_join(plot_dat_sub2, by = "cell_id")
    
    height <- 5
    width  <- 5
    p <-
      dat2 %>%
      ggplot(aes(x = SPATIAL_1, SPATIAL_2)) +
      geom_tile(
        height = height, 
        width = width,
        fill = "gray90"
      ) +
      geom_tile(
        data = dat3,
        aes(fill = .data[[annotation_col2]]),
        height = height, 
        width = width,
        show.legend = c(fill = TRUE, alpha = FALSE)
      ) +
      theme_void() +
      coord_equal() +
      scale_y_reverse() +
      ggsci::scale_fill_igv() +
      # scale_fill_manual(values = allcolour) +
      theme(panel.background = element_rect(fill = "black", color = NA)) +
      ggtitle(sam)
    p
    
  })


pdf(
  file.path(out_dir, "epi - v4 default normalize and sc de gene", "spatial plot of epi annotation in each sample.pdf"),
  height = 8,
  width = 8
)
p_sp_ct_ls2 %>% walk(print)
dev.off()


### cell props - annotation ---------
p_col_grp2 <-
  sp_epi2[[ ]] %>%
  group_by(group, annotation_epi) %>%
  summarise(count = n()) %>%
  mutate(freq = count / sum(count)) %>%
  ggplot(aes(x = group, y = freq, fill = annotation_epi)) +
  geom_col() +
  theme_bw() +
  ggsci::scale_fill_aaas()  

p_col_sample2 <-
  sp_epi2[[]] %>%
  group_by(sample, annotation_epi) %>%
  summarise(count = n()) %>%
  mutate(freq = count / sum(count))  %>%
  
  ggplot(aes(x = sample, y = freq, fill = annotation_epi)) +
  geom_col()+
  ggsci::scale_fill_igv()

# export
pdf(
  file.path(out_dir, "epi - v4 default normalize and sc de gene", "cell prop of epi annotation - group.pdf"),
  height = 6,
  width = 4
)
p_col_grp2 %>% print
dev.off()

pdf(
  file.path(out_dir, "epi - v4 default normalize and sc de gene", "cell prop of epi annotation - sample.pdf"),
  height = 6,
  width = 8
)
p_col_sample2 %>% print
dev.off()


# 20260511 merge cluster 0 and 1 to tumor_others 
# '0' = "tumor_1",
# '1' = "tumor_others",
# '2' = "normal"

## v2 - using epi - v4 default normalize and sc de gene ---------- 
### cell annotation ---------
sp_epi3 <- sp_epi[, !sp_epi$RNA_snn_res.0.2 %in% c('3', '4', '5')]
cell_annotaiton2 <- c(
  '0' = "tumor_1",
  '1' = "tumor_others",
  '2' = "tumor_others"
)
Idents(sp_epi3) <- 'RNA_snn_res.0.2'
sp_epi3 <- sp_epi3 %>% RenameIdents(cell_annotaiton2) 
sp_epi3$annotation_epi <- Idents(sp_epi3)

saveRDS(sp_epi3, "sp_subcluster/sp_epi3.rds")

### spatial plot of all celltypes - annotation --------------
# cells_position_ls <-
#   cbind(
#     all_data[['sample']],
#     Embeddings(all_data[['spatial']])
#   ) %>%
#   rownames_to_column('cell_id') %>%
#   group_by(sample) %>%
#   nest() %>%
#   deframe()

annotation_col3 <- 'annotation_epi'
plot_dat_sub3 <- 
  FetchData(
    sp_epi3, 
    vars = c(annotation_col2)
  ) %>%
  rownames_to_column('cell_id')

p_sp_ct_ls3 <-
  cells_position_ls %>% # .[1] %>%
  imap(function(coords, sam){
    dat2 <- coords # %>% left_join(plot_dat_sub, by = "cell_id")
    dat3 <- coords %>% inner_join(plot_dat_sub3, by = "cell_id")
    
    height <- 5
    width <- 5
    p <-
      dat2 %>%
      ggplot(aes(x = SPATIAL_1, SPATIAL_2)) +
      geom_tile(
        height = height, 
        width = width,
        fill = "gray90"
      ) +
      geom_tile(
        data = dat3,
        aes(fill = .data[[annotation_col3]]),
        height = height, 
        width = width,
        show.legend = c(fill = TRUE, alpha = FALSE)
      ) +
      theme_void() +
      coord_equal() +
      scale_y_reverse() +
      ggsci::scale_fill_igv() +
      # scale_fill_manual(values = allcolour) +
      theme(panel.background = element_rect(fill = "black", color = NA)) +
      ggtitle(sam)
    p
    
  })


pdf(
  file.path(out_dir, "epi - v4 default normalize and sc de gene", "spatial plot of epi annotation in each sample - v2.pdf"),
  height = 8,
  width = 8
)
p_sp_ct_ls3 %>% walk(print)
dev.off()


### cell props - annotation ---------
p_col_grp3 <-
  sp_epi3[[ ]] %>%
  group_by(group, annotation_epi) %>%
  summarise(count = n()) %>%
  mutate(freq = count / sum(count)) %>%
  ggplot(aes(x = group, y = freq, fill = annotation_epi)) +
  geom_col() +
  theme_bw() +
  ggsci::scale_fill_aaas()  

p_col_sample3 <-
  sp_epi3[[]] %>%
  group_by(sample, annotation_epi) %>%
  summarise(count = n()) %>%
  mutate(freq = count / sum(count))  %>%
  
  ggplot(aes(x = sample, y = freq, fill = annotation_epi)) +
  geom_col()+
  ggsci::scale_fill_igv()

# export
pdf(
  file.path(out_dir, "epi - v4 default normalize and sc de gene", "cell prop of epi annotation - group - v2.pdf"),
  height = 6,
  width = 4
)
p_col_grp3 %>% print
dev.off()

pdf(
  file.path(out_dir, "epi - v4 default normalize and sc de gene", "cell prop of epi annotation - sample - v2.pdf"),
  height = 6,
  width = 8
)
p_col_sample3 %>% print
dev.off()


### validation of pathway score ----------------
# `sp_epi_score2` is created in `pathway validation of tumor 1` 

path_name <- c('Metastatic', 'Invasion', 'Proliferation', 'CellCycle')
dat_epi_score <- sp_epi_score2[[path_name]]  

p_vln_path_validation2 <-
  sp_epi3 %>%
  AddMetaData(dat_epi_score) %>% 
  VlnPlot(
    group.by = "annotation_epi",
    features = path_name,
    pt.size = 0
  ) &
  ggsci::scale_fill_igv()

p_vln_path_validation2_2 <-
  p_vln_path_validation2 &
  geom_boxplot(width = 0.15, outliers = FALSE) &
  ggpubr::stat_compare_means() & #comparisons = list(unique(sp_epi3$annotation_epi)))
  scale_y_continuous(limits = function(x) c(x[1], 1.1 * x[2]))

# export 
pdf(
  file.path("sp_subcluster/epi - v4 default normalize and sc de gene/vlnplot of some pahway of tumor_1 - v2.pdf"),
  height = 7,
  width  = 12
)
p_vln_path_validation2_2 %>% print()
dev.off()



# sc cluster --------------
sc_epi2 <- sc_epi
DefaultAssay(sc_epi2) <- 'RNA'

# sc_epi2 <- sc_normalize(sc_epi2, scale.factor = "median")
sc_epi2 <- sc_normalize(sc_epi2, scale.factor = "default")
# ####
# sc_epi2 <- ScaleData(
#   JoinLayers(sc_epi2), 
#   vars.to.regress = 'nCount_RNA',
#   features = de_marker_sc_ls %>% unlist() %>% .[. %in% rownames(sc_epi2)]
# )
# ####
sc_epi2 <- sc_pca(
  sc_epi2, 
  sc_out_dir = file.path(out_dir, "epi - sc/sc.1 pca"),
  dir.create.force = TRUE
)

ndims_epi2 <- 20

sc_epi2 <- sc_dr(
  sc_epi2, 
  ndims = ndims_epi2,
  return.model = TRUE,
  sc_out_dir = file.path(out_dir, "epi - sc/sc.2 dr"),
  dir.create.force = TRUE
)


sc_epi2 <- sc_cluster(
  sc_epi2, 
  ndims = ndims_epi,
  reduction =  "harmony", # "harmony", "pca"
  sc_out_dir = file.path(out_dir, "epi - sc/sc.3 cluster"),
  dir.create.force = TRUE,
  
  harmony.reduction.name = "harmonyumap",
  resolution = c(0.2, 0.4, 0.6, 0.8)
)


sc_epi2[[]] %>%
  dplyr::filter(sample != "H2502760") %>%
  group_by(sample, RNA_snn_res.0.2) %>%
  summarise(count = n()) %>%
  mutate(freq = count / sum(count)) %>%
  ggplot(aes(x = sample, y = freq, fill = RNA_snn_res.0.2)) +
  geom_col() +
  theme_bw() +
  ggsci::scale_fill_aaas() 

sc_epi2[[]] %>%
  # dplyr::filter(sample != "H2502760") %>%
  dplyr::filter(RNA_snn_res.0.2 != "1") %>%
  group_by(group, RNA_snn_res.0.2) %>%
  summarise(count = n()) %>%
  mutate(freq = count / sum(count)) %>%
  ggplot(aes(x = group, y = freq, fill = RNA_snn_res.0.2)) +
  geom_col() +
  theme_bw() +
  ggsci::scale_fill_aaas()


sc_epi2[[]] %>%
  group_by(sample, RNA_snn_res.0.4) %>%
  summarise(count = n()) %>%
  mutate(freq = count / sum(count)) %>%
  ggplot(aes(x = sample, y = freq, fill = RNA_snn_res.0.4)) +
  geom_col() +
  theme_bw() +
  ggsci::scale_fill_aaas() 



## tranfer to spatial  -----
Idents(sc_epi2) <- "RNA_snn_res.0.2"
all_markers_epi_sc2 <- FindAllMarkers(sc_epi2)

# all de 
de_marker_sc_ls2 <- 
  all_markers_epi_sc2 %>%
  dplyr::filter(p_val_adj < 0.05, avg_log2FC > 0.25, pct.1 - pct.2 > 0.1) %>%
  arrange(desc(avg_log2FC), .by_group = TRUE) %>%
  with(split(gene, cluster))

# add module score
sp_epi2_score <- 
  sp_epi %>%
  JoinLayers() %>%
  AddModuleScore(features = de_marker_sc_ls2 %>% map(head, 100))
path_idx2 <- str_detect(colnames(sp_epi2_score[[]]), "^Cluster")
colnames(sp_epi2_score[[]])[path_idx2] <- names(de_marker_sc_ls2)

DotPlot(
  sp_epi2_score,
  features = names(de_marker_sc_ls2),
  group.by = 'RNA_snn_res.0.2'
)  +
  RotatedAxis()

DotPlot(
  sp_epi2_score,
  features = names(de_marker_sc_ls2),
  group.by = 'RNA_snn_res.0.4'
)  +
  RotatedAxis()
DotPlot(
  sp_epi2_score,
  features = names(de_marker_sc_ls2),
  group.by = 'RNA_snn_res.0.6'
)  +
  RotatedAxis()




sp_epi2_score[[]] %>%
  group_by(sample, RNA_snn_res.0.2) %>%
  summarise(count = n()) %>%
  mutate(freq = count / sum(count))  %>%
  
  ggplot(aes(x = sample, y = freq, fill = RNA_snn_res.0.2)) +
  geom_col()


sp_epi2_score[[]] %>%
  group_by(group, RNA_snn_res.0.2) %>%
  summarise(count = n()) %>%
  mutate(freq = count / sum(count))  %>%
  
  ggplot(aes(x = group, y = freq, fill = RNA_snn_res.0.2)) +
  geom_col()


FeaturePlot(
  sp_epi2_score,
  features = "1",
  reduction = "harmonyumap",
  order = TRUE
)

FeaturePlot(
  sp_epi2_score,
  features = "0",
  reduction = "harmonyumap",
  order = TRUE
)
FeaturePlot(
  sp_epi2_score,
  features = "2",
  reduction = "harmonyumap",
  order = TRUE
)




DimPlot(
  sp_epi2_score, 
  group.by = "RNA_snn_res.0.2",
  reduction = "harmonyumap"
)





FeaturePlot(sp_epi, features = "COL7A1", reduction = "harmonyumap")
DimPlot(sp_epi, group.by = "RNA_snn_res.0.2", reduction = "harmonyumap")

FeaturePlot(sp_epi, features = "GJA1", reduction = "harmonyumap")
FeaturePlot(sp_epi, features = "GJA1", reduction = "harmonyumap")








