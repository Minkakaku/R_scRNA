library(tidyverse)
library(Seurat)


# only used fibroblast
all_data_rctd <- readRDS("/mnt/NexData_file/home/yrc/XE2511046228/sc_rds/counts_raw.RDS")

all_data_sp_rctd <- readRDS("/mnt/NexData_file/home/yrc/XE2511046228/sp_rds/ALL.RDS")

cells_position_ls <- readRDS('/mnt/NexData_file/home/yrc/XE2511046228/sp_rds/cells_position_ls.rds')


# RCTD ------------------
library(spacexr)
rctd_out_dir <- 'RCTD' 
dir.create(rctd_out_dir)

#1. sketch data  
# NOT RUN. ROI data, only 50k cells
# visium_roi_ls %>%
#   map(function(visium){
#     # sketch the cortical subset of the Visium HD dataset
#     DefaultAssay(visium) <- "Spatial.008um"
#     visium <- FindVariableFeatures(visium)
#     visium <- SketchData(
#       object = cortex,
#       ncells = 50000,
#       method = "LeverageScore",
#       sketched.assay = "sketch"
#     )
#     
#     DefaultAssay(visium) <- "sketch"
#     visium <- ScaleData(visium)
#     visium <- RunPCA(visium, assay = "sketch", reduction.name = "pca.sketch", verbose = T)
#     visium <- FindNeighbors(visium, reduction = "pca.sketch", dims = 1:50)
#     visium <- RunUMAP(visium, reduction = "pca.sketch", reduction.name = "umap.sketch", return.model = T, dims = 1:50, verbose = T)
#     
#     visium
#   })

#2. single data
# all_data_rctd <- readRDS("all_data_sc_dim15_anno.rds")
# all_data_rctd <- all_data


Idents(all_data_rctd) <- "celltype"
counts_sc  <- all_data_rctd@assays$RNA@counts

cluster_sc <- as.factor(all_data_rctd$celltype)
nUMI_sc    <- all_data_rctd$total_counts

sc_cell_filter_idx <- ! cluster_sc == "Unknown"

# create the RCTD reference object
reference <- Reference(
  counts_sc[, sc_cell_filter_idx], 
  droplevels(cluster_sc[sc_cell_filter_idx]), 
  nUMI_sc[sc_cell_filter_idx]
)

# 3. create the RCTD query object
all_data_sp_rctd_ls <- all_data_sp_rctd %>% JoinLayers() %>% SplitObject('sample')

rctd_query_ls <-
  all_data_sp_rctd_ls %>%
  imap(function(xenium, name){ # browser()
    counts_hd <- xenium %>% LayerData(layer = "counts")
    coords    <- 
      cells_position_ls[[name]] %>% 
      dplyr::select(
        cell_id,
        x = x_centroid, 
        y = y_centroid
      ) %>%
      column_to_rownames('cell_id')
    
    # create the RCTD query object
    query <- SpatialRNA(coords, counts_hd, colSums(counts_hd))
    query
  })

#4. run RCTD
rctd_result_ls <- 
  rctd_query_ls %>%
  imap(function(query, name){
    # run RCTD
    RCTD <- create.RCTD(query, reference, max_cores = 10, UMI_min = 10, UMI_min_sigma = 50)
    RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")
    
    RCTD
  })

saveRDS(rctd_result_ls, 'RCTD/rctd_result_ls.rds')

#5. addMetaData
all_data_sp_rctd_anno_ls <-
  all_data_sp_rctd_ls %>%
  imap(function(xenium, name){ # browser()
    rctd <- rctd_result_ls[[name]]
    xenium <- AddMetaData(xenium, metadata = rctd@results$results_df)
    
    xenium$first_type <- as.character(xenium$first_type)
    xenium$first_type[is.na(xenium$first_type)] <- "Unknown"
    xenium
  })




#6. dimplot 
# all_data_sp_anno <-
#   merge(
#     all_data_sp_rctd_anno_ls[[1]],
#     all_data_sp_rctd_anno_ls[-1]
#   )

meta_rctd_result <-
  rctd_result_ls %>%
  map(~ .@results$results_df) %>%
  bind_rows()

all_data_sp_rctd_anno <-
  all_data_sp_rctd %>% 
  AddMetaData(meta_rctd_result)



# #7. ProjectData
# local({
#   sp <- all_data_sp_rctd_anno_ls$X_202413076_32 
#   idx <- sp$first_type != 'Unknown'
#   cells <- colnames(sp)[idx]
#   
#   assay_sub <- sp[['RNA']] %>% subset(Cells = cells)
#   sp[['RCTD']] <- assay_sub
#   DefaultAssay(sp) <- "RCTD"
#   
#   sp@reductions <- list()
#   sp <- sp %>% RunPCA()
#   
#   ElbowPlot(sp, ndims = 30) %>% print()
#   # stop()
#   dims <- 10
#   
#   
#   sp <- RunUMAP(sp, dims = 1:dims, return.model = TRUE)
#   
#   # sp <- sp %>%
#   #   harmony::RunHarmony(
#   #     dims.use = 1:dims,
#   #     group.by.vars = "project",
#   #     reduction.save = "harmony"
#   #   )
#   # 
#   # sp <- RunUMAP(
#   #   sp,
#   #   dims = 1:dims,
#   #   reduction = "harmony",
#   #   reduction.name = "harmonyumap",
#   #   reduction.key  = "harmonyUMAP_", 
#   #   return.model = TRUE
#   # )
#   
#   
#   
#   sp_rctd <<- 
#     sp %>% 
#     ProjectData(
#       assay = "RNA",
#       full.reduction = "pca.full", # "harmony.full",
#       
#       sketched.assay = "RCTD",
#       sketched.reduction = "pca", # "harmony",# "pca",
#       
#       umap.model =  "umap", # "harmonyumap", 
#       dims = 1:dims, 
#       
#       refdata = list(
#         annotation_full = "first_type"
#       )
#     )
# })
# sp_rctd %>% DimPlot(reduction = "umap", group.by = "first_type")



#6. spatial plot 
# SpatialFeaturePlot(
#   visium_roi_ls[[1]], 
#   image.scale = 'hires', 
#   pt.size.factor = 5,
#   shape = 22,
#   alpha = c(0.5, 0.8),
#   features = "nFeature_Spatial.008um"
# )

# image_scale_ls <-
#   visium_roi_rctd_ls %>%
#   # visium_roi_ls %>% 
#   map(GetTissueCoordinates) %>%
#   map('x') %>%
#   map(range) %>%
#   map(diff) %>%
#   map(log) %>%
#   {map(., function(x, max_val)x/max_val, max_val = max(unlist(.)))}

# pdf(
#   file.path(rctd_out_dir, 'spatial plot of RCTD.pdf'), 
#   height = 15, 
#   width = 15
# )
# visium_roi_rctd_ls %>%
#   imap(function(visium, name){
#     SpatialDimPlot(
#       visium,
#       group.by       = 'first_type',
#       image.scale    = 'hires',
#       pt.size.factor = 5 * image_scale_ls[[name]],
#       shape          = 22,
#       alpha          = 0.5
#     ) +
#       ggtitle(name) +
#       ggsci::scale_fill_igv() +
#       guides(fill = guide_legend(override.aes = list(size = 5)))
#   })
# dev.off()



p_col_rctd <-
  all_data_sp_rctd_anno[[c('first_type', 'group')]] %>%
  rownames_to_column("cell_id") %>%
  dplyr::filter(! is.na(first_type)) %>%
  group_by(group, first_type) %>%
  summarise(count = n()) %>%
  group_by(group) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = group, y = freq)) +
  geom_col(aes(fill = first_type)) +
  scale_fill_manual(values = c(
    'S3 fibroblasts'                             = "#1F77B4FF",
    'Immediate early-response (IER) fibroblasts' = "#FF7F0EFF",
    'S1 fibroblasts'                             = "#2CA02CFF",
    'Peri-cryptal S2 fibroblasts'                = "#D62728FF"
  )) +
  scale_y_continuous(labels = scales::percent) + 
  theme_bw() +
  RotatedAxis() +
  labs(
    x = NULL,
    y = NULL
  )

ggsave(
  p_col_rctd,
  filename = file.path('sp_figure', 'RCTD cell proportion barplot.pdf'),
  height = 6,
  width  = 8
)


# export
saveRDS(all_data_sp_rctd_anno, 'RCTD/all_data_sp_rctd_anno.rds')
DimPlot(all_data_sp_rctd_anno,group.by = "celltype",label = T)
ggsave("/mnt/NexData_file/home/yrc/XE2511046228/RCTD/ALL_celltype.png", height = 5,width  = 7)
ggsave("/mnt/NexData_file/home/yrc/XE2511046228/RCTD/ALL_celltype.pdf", height = 5,width  = 7)

DimPlot(all_data_sp_rctd_anno,group.by = "first_type",label = T)
ggsave("/mnt/NexData_file/home/yrc/XE2511046228/RCTD/ALL_celltype_rctd.png", height = 5,width  = 7)
ggsave("/mnt/NexData_file/home/yrc/XE2511046228/RCTD/ALL_celltype_rctd.pdf", height = 5,width  = 7)

####空间定位####
library(tidyverse)
library(Seurat)

out_dir <- 'RCTD/spatial'
dir.create(out_dir)

all_data <- all_data_sp_rctd_anno
all_data$celltype <- all_data$first_type
all_data@meta.data %>% head()

position <- readRDS("./R_OUTPUT/data_for_analysis/cells_position_ls.rds")


all_coords <- bind_rows(position, .id = "sample")
all_coords <- all_coords %>% select(cell_id, x_centroid, y_centroid)

# 将 all_coords 转换为以 cell_id 为行名的矩阵
coords_matrix <- all_coords %>%
  column_to_rownames("cell_id") %>%
  select(x_centroid, y_centroid)

# 只保留 Seurat 对象中存在的细胞
coords_matrix <- coords_matrix[colnames(all_data), , drop = FALSE]

# 添加到 meta.data
all_data <- AddMetaData(all_data, coords_matrix)

plot_data <- all_data@meta.data %>%
  select(x_centroid, y_centroid, celltype, sample) %>%
  filter(!is.na(x_centroid))  # 去除可能缺失坐标的细胞


for (sm in unique(plot_data$sample)) {
  p <- ggplot(plot_data %>% filter(sample == sm), aes(x = x_centroid, y = y_centroid, color = celltype)) +
    geom_point(size = 0.2, alpha = 0.8) +
    theme_classic() +   # 白色背景，无网格线
    theme(
      axis.text.x = element_blank(),   # 去掉x轴数值
      axis.text.y = element_blank(),   # 去掉y轴数值
      axis.ticks = element_blank(),    # 去掉刻度线（可选）
      axis.line = element_blank(),     # 去掉默认的轴线（因为我们要加方框）
      panel.border = element_rect(color = "black", fill = NA, size = 1)) +
    labs(title = sm, x = "spatial_1", y = "spatial_2") +
    coord_fixed()+
    guides(color = guide_legend(override.aes = list(size = 3)))
  ggsave(file.path(out_dir, paste0("spatial_", sm, ".png")), p, width = 8, height = 6)
  ggsave(file.path(out_dir, paste0("spatial_", sm, ".pdf")), p, width = 8, height = 6)
}

####堆叠图####
p_col_proj <-
  all_data@meta.data %>%
  group_by(group1, celltype) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = group1, y = freq, fill = celltype)) + 
  geom_col() +
  ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )
ggsave(
  p_col_proj,
  filename = file.path("RCTD","prop_group1.png"),
  height = 6,
  width  = 4
)
ggsave(
  p_col_proj,
  filename = file.path("RCTD","prop_group1.pdf"),
  height = 6,
  width  = 4
)

p_col_proj <-
  all_data@meta.data %>%
  group_by(sample, celltype) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = sample, y = freq, fill = celltype)) + 
  geom_col() +
  ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )
ggsave(
  p_col_proj,
  filename = file.path("RCTD","prop_sample.png"),
  height = 6,
  width  = 18
)
ggsave(
  p_col_proj,
  filename = file.path("RCTD","prop_sample.pdf"),
  height = 6,
  width  = 18
)

####箱线图####
all_data$group <- all_data$group3

Ratio <- all_data@meta.data %>%
  group_by(sample,group,celltype,) %>%
  count() %>%
  group_by(sample) %>%
  mutate(Freq = n/sum(n)*100)

my_comparisons <- list(c("Con","NR"),c("Con","PR"),c("Con","CR"),c("NR","PR"),c("NR","CR"),c("PR","CR"))#差异比较组

library(ggpubr)
ggplot(Ratio,aes(group,Freq))+
  geom_boxplot(aes(fill = group))+
  geom_point(size=1)+
  scale_x_discrete(limits=levels(Ratio$group))+ #很重要！保证X轴的label和实际level顺序对应
  scale_fill_discrete()+
  stat_compare_means(method = "wilcox.test", comparisons = my_comparisons, aes(label = ..p.signif..)) +  
  scale_y_continuous("Proportion(%)")+
  facet_grid(~celltype)+
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())+
  xlab("")+
  theme(panel.grid = element_blank())

ggsave(filename = file.path("RCTD","prop_boxplot_group3.png"),height = 6,width  = 18)
ggsave(filename = file.path("RCTD","prop_boxplot_group3.pdf"),height = 6,width  = 18)



# save all cells -------------
all_data <- readRDS("all_data_grp.rds")
all_data_meta_rctd <-
  all_data[['annotation_full']] %>% 
  rownames_to_column("cell_id") %>%
  left_join(rownames_to_column(meta_rctd_result['first_type'], "cell_id"), by = "cell_id") %>%
  mutate(annotation3 = case_when(
    annotation_full == "Fibroblast" &  is.na(first_type) ~ "Fibroblast_others",
    annotation_full == "Fibroblast" & !is.na(first_type) ~ first_type,
    TRUE ~ annotation_full
  )) %>%
  dplyr::select(cell_id, annotation3) %>%
  column_to_rownames("cell_id")

all_data_rctd <- all_data %>% AddMetaData(all_data_meta_rctd)
saveRDS(all_data_rctd, 'all_data_rctd.rds')

# singler  ------------------
library(SingleR)

sc_ref_dat <- all_data_rctd %>% JoinLayers() %>% LayerData(layer = "data")


DefaultAssay(all_data_sp_rctd) <- "RNA"
Idents(all_data_sp_rctd) <- "RNA_snn_res.0.6"

all_markers <-
  all_data_sp_rctd %>% 
  JoinLayers() %>%
  FindAllConservedMarkers(
    grouping.var = 'project'
  )
# saveRDS(all_markers, "all_markers_res.0.6.rds")
# all_markers <- readRDS("all_markers_res.0.6.rds")

all_markers_df <- 
  all_markers %>%
  map(rownames_to_column, 'gene') %>%
  enframe("cluster") %>%
  unnest("value") %>%
  mutate(
    p_val_adj  = minimump_p_val,
    avg_log2FC = mean(c_across(ends_with('avg_log2FC'))),
    pct.1      = mean(c_across(ends_with('pct.1')))
  ) 

hvg_gene <- VariableFeatures(all_data_sp_rctd)

head200_marker_gene <-
  all_markers_df %>%
  group_by(cluster) %>%
  arrange(desc(pct.1)) %>%
  slice_head(n = 100) %>%
  dplyr::pull(gene) %>%
  unique() 
  
head200_marker_gene <- 
  head200_marker_gene %>%
  .[. %in% hvg_gene]

xenium_data_full <- all_data_sp_rctd %>% JoinLayers() %>% LayerData(layer = "data")

singler_main <-
  SingleR::SingleR(
    test     = xenium_data_full[head200_marker_gene, ] %>% as.matrix(),
    ref      = sc_ref_dat %>% as.matrix(),
    labels   = all_data_rctd$annotation, 
    clusters = all_data_sp_rctd$RNA_snn_res.0.6
  )

cluster_fibro <- 
  structure(
    singler_main$labels,
    names = rownames(singler_main)
  )

# main cluster
Idents(all_data_sp_rctd) <- 'RNA_snn_res.0.6'
all_data_sp_rctd <- all_data_sp_rctd %>% RenameIdents(cluster_fibro)
all_data_sp_rctd$annotation_fibro <- Idents(all_data_sp_rctd)

DimPlot(all_data_sp_rctd, group.by = "annotation_fibro")

saveRDS(all_data_sp_rctd, file = file.path('sp_subcluster', 'sp_fibro_singler.rds'))

# # subcluster 
# singler_sub_ls <-
#   all_data_sp_rctd$annotation_main %>%
#   unique() %>%
#   set_names(., .) %>%
#   map(function(cluster){
#     x_idx  <- all_data_sp_rctd$annotation_main == cluster
#     sc_idx <- all_data_sp_rctd$annotation   == cluster
#     
#     x_dat  <- xenium_data_full[head200_marker_gene, x_idx]
#     sc_dat <- sc_ref_dat[, sc_idx]
#     
#     singler_main <-
#       SingleR(
#         test     = x_dat  %>% as.matrix(),
#         ref      = sc_dat %>% as.matrix(),
#         labels   = all_data_sc$annotation2[sc_idx], 
#       )
#     singler_main
#   })







