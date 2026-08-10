#1 上皮细胞以nCounts <200，重新降维分析epi分亚群，dim15 r0.4的，本次亚群的marker基因列表
#2 本次亚群与之前亚群的对应关系（cluster consistency+桑葚图）
#3 本次亚群的hpv gene表达dotplot+HPV基因感染率（y轴是hpv gene count大于0的细胞数除以本组内的总细胞数，5.19那次更新的定义）
#4  roi中cluster占比柱状图（x轴是各个细胞亚群，y是各个roi）
#5 cluster中roi占比柱状图（x轴是ROI分组HPV-NE,HPV+NE,HPV-NA,HPV+NA,LSIL,HSIL,Ca，y是各个亚群占比）(ROI信息excel版本26-5-6)
#6 本次亚群的空间分布图
#7 各个亚群细胞数量以及在各个roi中的细胞数量
# 20260630 dim15 r0.4

library(tidyverse)
library(Seurat)
library(patchwork)

source("subcluster/sc_analysis_tools.R")

all_data <- readRDS('output/data/adata_filtered_add_celltype.rds')

out_dir <- "sp_subcluster_nCount200"
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

# QC
pdf(
  file.path(out_dir, 'epi', 'Vlnplot.pdf'),
  height = 6,
  width = 10
)
VlnPlot(
  sp_epi,
  features = c('nCount_RNA', 'nFeature_RNA', 'log10GenesPerUMI' ),
  pt.size = 0
)
VlnPlot(
  sp_epi,
  features = c('nCount_RNA', 'nFeature_RNA', 'log10GenesPerUMI' ),
  pt.size = 0,
  log = T
)
dev.off()

# filter
sp_epi <- sp_epi[, sp_epi$nCount_RNA > 200]

# norm
DefaultAssay(sp_epi) <- 'RNA'
sp_epi <- sc_normalize(sp_epi, scale.factor = "default")

## sketch data --------
sp_epi <-
  SketchData(
    sp_epi,
    ncells         = 120000L, # 1/5 - 1/3
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

# set.seed(1234)
# sp_epi2 <- sp_epi[, sample(colnames(sp_epi), 1e4)]

pdf(file.path(out_dir, "epi/sc.1 pca", 'dimheatmap 1-5 - downsample.pdf'), height = 6, width = 10)
DimHeatmap(sp_epi, dims = 1:5, cells = 1000)
dev.off()
pdf(file.path(out_dir, "epi/sc.1 pca", 'dimheatmap 6-10 - downsample.pdf'), height = 6, width = 10)
DimHeatmap(sp_epi, dims = 6:10, cells = 1000)
dev.off()
pdf(file.path(out_dir, "epi/sc.1 pca", 'dimheatmap 11-15 - downsample.pdf'), height = 6, width = 10)
DimHeatmap(sp_epi, dims = 11:15, cells = 1000)
dev.off()
pdf(file.path(out_dir, "epi/sc.1 pca", 'dimheatmap 16-20 - downsample.pdf'), height = 6, width = 10)
DimHeatmap(sp_epi, dims = 16:20, cells = 1000)
dev.off()


ndims_epi <- 15

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
  k.param = 20,
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
FeaturePlot(
  sp_epi,
  features = c("total_counts", "nFeature_RNA"),
  reduction = "harmonyumap"
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

# DefaultAssay(sp_epi) <- 'RNA'

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

Idents(sp_epi) <- 'RNA_snn_res.0.4' # 'annotation_full'
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

sp_epi$RNA_snn_res.0.4 <- fct_inseq(sp_epi$RNA_snn_res.0.4)

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



# 2. cluster consistency -------------
# sp_epi_origin <- readRDS('sp_epi_projected-v2.rds')  


plot_dat_sankey_ls <- 
  list(
    original = sp_epi_origin,
    ncount_gt_200 = sp_epi
  ) %>%
  map(`[[`) %>%
  imap(function(x, name){
    if(name == 'original')(
      select_col = 'RNA_snn_res.0.2'
    )else{
      select_col = 'RNA_snn_res.0.4'
    }
    res <-
      x[select_col] %>%
      rownames_to_column('cell_id') %>%
      mutate(RNA_snn_res.0.4 = fct_inseq(.data[[select_col]], ordere = NA))
    colnames(res) <- c('cell_id', str_glue("res.0.4_{name}"))
    res
  })
plot_dat_sankey <-
  plot_dat_sankey_ls$original %>%
  left_join(plot_dat_sankey_ls$ncount_gt_200) %>%
  map_dfc(str_replace_na, "NA") %>%
  mutate(res.0.4_original = fct_inseq(res.0.4_original)) %>%

  group_by(res.0.4_original, res.0.4_ncount_gt_200) %>%
  summarise(count = n())  

library(ggalluvial)
p_sankey <-
  ggplot(data = plot_dat_sankey,
       aes(axis1 = res.0.4_original,   # 第一个轴：舱位等级
           axis2 = res.0.4_ncount_gt_200,     # 第二个轴：性别
           y = count)) +     # y 轴数值：频数
  geom_alluvium(aes(fill = res.0.4_original)) + # 冲积层按生存状态着色
  geom_stratum(fill = NA) +                      # 绘制层
  geom_text(stat = "stratum",           # 添加层标签
            aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("original", "ncount_gt_200"), 
                   expand = c(.2, .05)) +
  labs(
  ) +
  theme_minimal() +
  ggsci::scale_fill_igv()

# export
ggsave(
  filename = file.path(out_dir, 'epi', '2. sankey plot of cluster changes - v2 - res0.2 to res0.4.pdf'),
  height = 6,
  width = 10
)

# 3. HPV gene dotplot and average count -------------
hpv_gene <- all_data %>% rownames() %>% str_subset("^V[0-9]{2}[EL]")
stopifnot(length(hpv_gene) == 100)


## dotplot 
p_dot_hpv <-
  DotPlot(
    sp_epi,
    features = hpv_gene,
    group.by = "RNA_snn_res.0.4"
  ) +
  RotatedAxis()

pdf(
  file = file.path(out_dir, 'epi', '3. Dotplot of HPV gene expression - res_0.4.pdf'),
  height = 4,
  width = 22
)
p_dot_hpv %>% print
dev.off()


## batplot of HPV positiva ratio
# x axis is group, HPV gene positive ratio 
dat_hpv_res_0.4 <- 
  sp_epi %>%
  FetchData(
    vars = c(hpv_gene, "RNA_snn_res.0.4"),
    layer = "data"
  )


p_col_roi_hpv_ratio_ls <-
  hpv_gene %>% 
  setNames(., .) %>%
  map(function(x){
    p <-
      dat_hpv_res_0.4 %>%
      group_by(RNA_snn_res.0.4) %>%
      # dplyr::filter(!is.na(.data[[x]])) %>% 
      summarise(AverageExpr = sum(.data[[x]])/n()) %>%
      
      ggplot(aes(fill = RNA_snn_res.0.4, y = AverageExpr, x = RNA_snn_res.0.4)) +
      geom_col() +
      ggsci::scale_fill_igv() +
      theme_classic()  +
      labs(
        title = x,
        fill = "RNA_snn_res.0.4"
      ) +
      RotatedAxis()
    
    # return
    p
  })


# export
pdf(
  file.path(out_dir, 'epi', '3 barplot of HPV gene positive ratio - group by RNA_snn_res.0.4.pdf'),
  height = 8,
  width  = 6
)
p_col_roi_hpv_ratio_ls %>% print
dev.off()

# 4. ROI: ROI props in cluster ------------- 
# x is cluster, y is roi
roi_color <- c(
  'HPV-NE' = '#FBC04C',
  'HPV+NE' = '#648FEB',
  'HPV-NA' = '#52C7B7',
  'HPV+NA' = '#CCE4FE',
  'LSIL'   = '#6E6FB2',
  'HSIL'   = '#DDC4E2',
  'Ca'     = '#F27C7C'
)

cluster_res_0.4  <- sp_epi$RNA_snn_res.0.4 %>% unique() %>% str_sort(numeric = T)
color_cluster_res_0.4 <- structure(
  ggsci::pal_igv()(length(cluster_res_0.4)),
  names = cluster_res_0.4
)


df_roi_cell_id2_uniq <- readRDS('roi_cellid_group_260518.rds')

df_roi_cell_epi <-
  sp_epi[['RNA_snn_res.0.4']] %>% 
  rownames_to_column("cell_id_sample") %>%
  inner_join(df_roi_cell_id2_uniq, by = "cell_id_sample")

dim(df_roi_cell_epi)
# [1] 582359      7



p_bar_roi_1 <-
  df_roi_cell_epi %>%
  group_by(RNA_snn_res.0.4, ROI) %>%
  
  dplyr::filter(!is.na(RNA_snn_res.0.4)) %>% 
  
  summarise(count = n()) %>%
  group_by(RNA_snn_res.0.4) %>%
  mutate(freq = count / sum(count) ) %>%
  left_join(distinct(df_roi_cell_epi[c("ROI", "group")]), by = 'ROI') %>%

  ggplot(aes(x = RNA_snn_res.0.4, y = freq, fill = group)) +
  geom_col(color = "gray") +
  scale_fill_manual(values = roi_color) +
  theme_classic()  +
  labs(
    title = "RNA_snn_res.0.4",
    fill = "roi_group"
  ) +
  RotatedAxis()
    
p_bar_roi_2 <-
  df_roi_cell_epi %>%
  group_by(RNA_snn_res.0.4, ROI) %>%
  
  dplyr::filter(!is.na(RNA_snn_res.0.4)) %>% 
  
  summarise(count = n()) %>%
  group_by(ROI) %>%
  mutate(freq = count / sum(count) ) %>%
  left_join(distinct(df_roi_cell_epi[c("ROI", "group")]), by = 'ROI') %>% 
  
  ggplot(aes(x = reorder(ROI, as.numeric(group)), y = freq, fill = RNA_snn_res.0.4)) +
  geom_col(width = 0.5) +
  geom_col(
    data = ~ dplyr::distinct(., ROI, group, .keep_all = TRUE),
    aes(x = ROI, color = group, y = 1),
    inherit.aes = FALSE,
    fill = NA,
    linewidth = 1.5
  ) + 
  scale_color_manual(values = roi_color, guide = guide_legend(override.aes = list(width = 5, fill = NA))) +
  ggsci::scale_fill_igv() +
  theme_classic()  +
  labs(
    x = NULL,
    y = "RNA_snn_res.0.4",
    color = "roi_group"
  )+
  RotatedAxis()
    
# export
pdf(
  file.path(out_dir, 'epi', '4 barplot of ROI - x is cluster y is roi.pdf'),
  height = 8,
  width  = 20
)
  list(
    p_bar_roi_1 + theme(plot.margin = margin(2, 200, 2, 200)),
    p_bar_roi_2
  ) %>%
  walk(print)
dev.off()

# 5. ROI: cluster props in ROI ------------- 
# x is roi, y is cluster
p_bar_roi_3 <-
  df_roi_cell_epi %>%
  group_by(group, RNA_snn_res.0.4) %>%
  
  summarise(count = n()) %>%
  group_by(group) %>%
  mutate(freq = count / sum(count) ) %>%
  
  ggplot(aes(fill = RNA_snn_res.0.4, y = freq, x = group)) +
  geom_col() +
  
  # ggsci::scale_fill_igv() +
  scale_fill_manual(values = color_cluster_res_0.4) +
  theme_classic()  +
  labs(
    # title = x,
    # fill = "roi_group"
  ) +
  RotatedAxis()

# export
pdf(
  file.path(out_dir,  'epi', '5 barplot of cluster - x is roi y is cluster.pdf'),
  height = 8,
  width  = 5
)
p_bar_roi_3 %>% print
dev.off()

# 6. spatial plot of epi cluster ------------
dat_sp <- cbind(
  Embeddings( all_data[['spatial']] ),
  all_data[[c('sample', 'fov', 'annotation', 'cell_id_sample')]]
)

dat_epi_cluster <- sp_epi[[c("RNA_snn_res.0.4", "cell_id_sample")]]

p_sp_ct_epi_res0.4_ls <-
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
        aes(fill = RNA_snn_res.0.4),
        height = height, 
        width = width,
        show.legend = c(fill = TRUE, alpha = FALSE)
      ) +
      theme_void() +
      coord_equal() +
      scale_y_reverse() +
      # ggsci::scale_fill_igv() +
      scale_fill_manual(values = color_cluster_res_0.4, na.value = "gray90") +
      theme(panel.background = element_rect(fill = "black", color = NA)) +
      ggtitle(sam)
    
    p
  })


# export 
pdf(
  file.path(out_dir, 'epi', str_glue("6 spatial plot of epi cluster - res_0.4.pdf")),
  height = 8,
  width = 8
)
p_sp_ct_epi_res0.4_ls %>% walk(print)
dev.off()


# 7. cluster count and cluster clunt in each ROI ----------------
list(
  count = sp_epi[['RNA_snn_res.0.4']] %>% 
    group_by(RNA_snn_res.0.4) %>%
    summarise(count = n()),
  count_roi = df_roi_cell_epi %>%
    group_by(RNA_snn_res.0.4, ROI) %>%
    summarise(count = n()),
  count_roi.group = df_roi_cell_epi %>%
    group_by(RNA_snn_res.0.4, group) %>%
    summarise(count = n()) 
) %>%
  openxlsx::write.xlsx(
    file = file.path(out_dir, "epi", "7. cluster count and cluster count in each ROI or roi_group - 3 sheets.xlsx"),
    asTable = TRUE
  )



# supp. export for xenium explorer -----------
sp_epi[[]] %>% 
  dplyr::select(
    fov, 
    cell_id = cell_id2, 
    group = 'RNA_snn_res.0.4'
  ) %>%
  group_by(fov) %>%
  nest() %>%
  deframe() %>%
  iwalk(function(dat, name){
    write.csv(
      dat,
      file = file.path(out_dir, str_glue("supp.epi.RNA_snn_res.0.4 - {name}.csv")),
      row.names = FALSE
    )
  })

 