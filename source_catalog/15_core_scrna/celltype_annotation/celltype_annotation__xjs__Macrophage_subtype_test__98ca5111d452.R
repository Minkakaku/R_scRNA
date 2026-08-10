Macrophage_out_dir <- "Macrophage_subtype/redction_pca" 
dir.create(Macrophage_out_dir)

reduction.use <- 'pca' # 'harmony' 'pca'


# cluster -------
Macrophage  <-
  FindNeighbors(
    Macrophage ,
    reduction = reduction.use,
    dims = 1:dims,
    k.param = 15
  )


Macrophage <- 
  FindClusters(
    Macrophage,
    resolution = c(0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.5)
  )



# dimplot
pdf(file.path(Macrophage_out_dir, str_glue("Dimplot dim{dims} cluster.pdf")), height = 6, width = 8)
DimPlot(Macrophage, group.by = "RNA_snn_res.0.2", label = TRUE) + NoLegend()
DimPlot(Macrophage, group.by = "RNA_snn_res.0.4", label = TRUE) + NoLegend()
DimPlot(Macrophage, group.by = "RNA_snn_res.0.6", label = TRUE) + NoLegend()
DimPlot(Macrophage, group.by = "RNA_snn_res.0.8", label = TRUE) + NoLegend()
DimPlot(Macrophage, group.by = "RNA_snn_res.1", label = TRUE) + NoLegend()
DimPlot(Macrophage, group.by = "RNA_snn_res.1.2", label = TRUE) + NoLegend()
DimPlot(Macrophage, group.by = "RNA_snn_res.1.5", label = TRUE) + NoLegend()
dev.off()


Idents(Macrophage) <- "RNA_snn_res.0.8"


## dotplot -------
macrophage_major_marker_ls <- c(
  # Interstitial macrophages
  "Interstitial macrophages" = "Mrc1",
  "Interstitial macrophages" = "Lyve1",
  "Interstitial macrophages" = "F13a1",
  "Interstitial macrophages" = "Cd163", 
  
  # Spp1 macrophages
  # "Spp1 macrophages" = "Spp1",
  "Spp1 macrophages" = "Apoe",
  # "Spp1 macrophages" = "Trem2",
  "Spp1 macrophages" = "Lpl",
  "Spp1 macrophages" = "Ctsd",
  # "Spp1 macrophages" = "Gpnmb",
  
  # Alveolar macrophages
  "Alveolar macrophages" = "Fabp4",
  "Alveolar macrophages" = "Siglecf",
  #"Alveolar macrophages" = "Fcgr1",   # 替代 "Cd64"
  "Alveolar macrophages" = "Mafb",
  "Alveolar macrophages" = "Pparg",
  "Alveolar macrophages" = "Marco", 
  
  # Activated macrophages
  "Activated macrophages" = "Nos2",
  #"Activated macrophages" = "Il1b",
  "Activated macrophages" = "Tnf",
  "Activated macrophages" = "Cxcl9",
  "Activated macrophages" = "Cxcl10",
  
  # Inflammatory macrophages
  "Inflammatory macrophages" = "Il1b",
  "Inflammatory macrophages" = "Ccl2",
  "Inflammatory macrophages" = "Ccl3",
  "Inflammatory macrophages" = "Ccl4",
  
  "Lactate export macrophages" = 'Ldha', 
  "Lactate export macrophages" = 'Slc16a3', 
  "Lactate export macrophages" = 'Mct4', 
  "Lactate export macrophages" = 'Hk2', 
  "Lactate export macrophages" = 'Pfkfb3', 
  "Lactate export macrophages" = 'Pkm', 
  "Lactate export macrophages" = 'Hif1a', 
  
  "Profibrotic macrophages" = 'Spp1', 
  "Profibrotic macrophages" = 'Lgals3', 
  "Profibrotic macrophages" = 'Trem2', 
  "Profibrotic macrophages" = 'Gpnmb', 
  "Profibrotic macrophages" = 'Mmp12', 
  "Profibrotic macrophages" = 'Tgfb1', 
  "Profibrotic macrophages" = 'Arg1', 
  
  "Monocyte-derived macrophages" = 'Ccr2',
  "Monocyte-derived macrophages" = 'Ccr5',
  # "Monocyte-derived macrophages" = 'Ccl2',
  "Monocyte-derived macrophages" = 'Ccl12',
  "Monocyte-derived macrophages" = 'Ccl7', 
  "Monocyte-derived macrophages" = 'Fcnb', 
  "Monocyte-derived macrophages" = 'Vcan', 
  "Monocyte-derived macrophages" = 'Ly6c2'
)


Idents(Macrophage)<-'RNA_snn_res.0.8'
macrophage_markers <-
  FindAllMarkers(
    JoinLayers(Macrophage),
    test.use = "wilcox",
    logfc.threshold = 0,
    min.pct = 0,
    only.pos = TRUE
  )
write.csv(macrophage_markers, file.path(Macrophage_out_dir, 'macrophage_markers_res.0.8.csv'))


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
DotPlot(Macrophage, group.by = "RNA_snn_res.0.8", features = major_marker_ls)+RotatedAxis()
DotPlot(Macrophage, group.by = "RNA_snn_res.0.8", features = macrophage_major_marker_ls)+RotatedAxis()
# DotPlot(Macrophage, group.by = "RNA_snn_res.0.8", features = CD8_major_marker_ls)+RotatedAxis()
# DotPlot(Macrophage, group.by = "RNA_snn_res.0.8", features = Epi_major_marker_ls)+RotatedAxis()
DotPlot(Macrophage, group.by = "RNA_snn_res.0.8", features = macrophage_top5_markers)+RotatedAxis()
dev.off()


##score--------------------
Interstitial_macrophage <- c( 
  'Mrc1',
  'Lyve1', 
  'F13a1', 
  'Cd163'
)
Spp1_macrophage <- c( 
  'Spp1',
  'Apoe', 
  'Trem2', 
  'Lpl', 
  'Ctsd', 
  'Gpnmb'
)
Alveolar_macrophage <- c(
  'Fabp4', 
  'Siglecf', 
  'Fcgr1', 
  'Mafb', 
  'Pparg', 
  'Marco'
)
Activated_macrophage <- c(
  'Nos2', 
  'Il1b', 
  'Tnf', 
  'Cxcl9', 
  'Cxcl10'
)
Inflammatory_macrophage <- c(
  'Il1b', 
  'Ccl2', 
  'Ccl3', 
  'Ccl4'
)
Lactate_export_macrophage <- c(
  'Ldha', 
  'Slc16a3', 
  'Mct4', 
  'Hk2', 
  'Pfkfb3', 
  'Pkm', 
  'Hif1a'
)
Profibrotic_macrophage <- c(
  'Spp1', 
  'Lgals3', 
  'Trem2', 
  'Gpnmb', 
  'Mmp12', 
  'Tgfb1', 
  'Arg1'
)
Monocyte_derived_macrophage <- c(
  'Ccr2',
  'Ccr5',
  'Ccl2',
  'Ccl12',
  'Ccl7'
)

pathway<-list(Interstitial_macrophage, Spp1_macrophage, Alveolar_macrophage, Activated_macrophage, Inflammatory_macrophage, Lactate_export_macrophage, Profibrotic_macrophage, Monocyte_derived_macrophage)
names(pathway)<-c('Interstitial_macrophage', 'Spp1_macrophage', 'Alveolar_macrophage', 'Activated_macrophage', 'Inflammatory_macrophage', 'Lactate_export_macrophage', 'Profibrotic_macrophage', 'Monocyte_derived_macrophage')

n_na_idx <- pathway %>% map_lgl(function(x) all(is.na(x))) %>% {! .}%>%  which
pathway <- pathway[n_na_idx] %>% map(function(x) x[!is.na(x)])

for ( i in seq_along(names(pathway))){
  Macrophage <- AddModuleScore(JoinLayers(Macrophage), pathway[i], name =  names(pathway)[i],nbin=5,ctrl=10)
}

p<-VlnPlot(Macrophage,
           features = str_c(names(pathway), '1'), 
           group.by='RNA_snn_res.0.8',
           pt.size = 0)

p2<-FeaturePlot(Macrophage,
                reduction = 'harmonyumap', 
                features = str_c(names(pathway), '1'))

ggsave(
  p,
  filename = file.path(Macrophage_out_dir, 'macrophage score - vlnplot.pdf'),
  height = 10,
  width  = 15
) 

ggsave(
  p2,
  filename = file.path(Macrophage_out_dir, 'macrophage score - featureplot.pdf'),
  height = 12,
  width  = 15
) 




# annotation ----------------------------------------------------------------
VlnPlot(Macrophage, features = 'nCount_RNA', pt.size = 0)
cell_cluster <- c(
  '0'  = "NK cell", # S1pr5, Cd244a 
  '1'  = "Endothelial", 
  '2'  = "Fibroblast", 
  '3'  = "AT2 cell", # Abca3
  '4'  = "Alveolar macrophage", # 
  '5'  = "Profibrotic macrophage", #
  '6'  = "Inflammatory macrophage", # Acod1, Il1rn
  '7'  = "cDC1", # Xcr1, Itgae, Batf3, Irf8
  '8'  = "AT1 cell", # Ager
  '9'  = "Profibrotic macrophage", # 
  '10' = "Interstitial macrophage", # 
  '11' = "Proliferating cDC1", # Mki67, Top2a
  '12' = "Monocyte-derived macrophage", # CCL*,CCR*
  '13'  = "Mature DC", # Cd40, Ccr7
  '14'  = "Megakaryocyte", # Itga2b, Gata1
  '15'  = "Alveolar macrophage"
)

Idents(Macrophage) <- "RNA_snn_res.0.8"

Macrophage <- RenameIdents(Macrophage, cell_cluster)
Macrophage$subtype <- Idents(Macrophage)


pdf(file.path(Macrophage_out_dir, str_glue("Dimplot dim{dims} annotation.pdf")), height = 6, width = 8)
DimPlot(Macrophage, group.by = "subtype", label = TRUE, repel = TRUE) + NoLegend()
dev.off()

pdf(file.path(Macrophage_out_dir, str_glue("Dimplot dim{dims} annotation split by group.pdf")), height = 6, width = 18)
DimPlot(Macrophage, group.by = "subtype", split.by = 'group', label = TRUE, repel = TRUE) + NoLegend()
dev.off()




## pure ----------
Macrophage_pure <- Macrophage[, Macrophage$subtype %in% c("Profibrotic macrophage", "Alveolar macrophage", "Monocyte-derived macrophage", "Interstitial macrophage", "Inflammatory macrophage")]
Macrophage_pure <- RunUMAP(
  Macrophage_pure, 
  dims = 1:dims, 
  return.model = TRUE
)
pdf(file.path(Macrophage_out_dir, str_glue("Dimplot dim{dims} annotation - pure.pdf")), height = 6, width = 8)
DimPlot(Macrophage_pure, group.by = "subtype", label = TRUE, repel = TRUE) + NoLegend()
dev.off()


p<-VlnPlot(Macrophage_pure,
           features = str_c(names(pathway), '1'), 
           group.by='subtype',
           pt.size = 0) &  
  theme(axis.text.x = element_text(angle = 90))

p2<-FeaturePlot(Macrophage_pure,
                features = str_c(names(pathway), '1'))

ggsave(
  p,
  filename = file.path(Macrophage_out_dir, 'macrophage score - Macrophage_pure - vlnplot.pdf'),
  height = 20,
  width  = 15
) 

ggsave(
  p2,
  filename = file.path(Macrophage_out_dir, 'macrophage score - Macrophage_pure - featureplot.pdf'),
  height = 12,
  width  = 15
) 




Macrophage_cellTye_markers <- c( 
  # Alveolar macrophages
  "Alveolar macrophages" = "Pparg",
  "Alveolar macrophages" = "Marco", 
  
  "Profibrotic macrophages" = 'Trem2', 
  "Profibrotic macrophages" = 'Gpnmb', 
  "Profibrotic macrophages" = 'Mmp12', 
  
  "Inflammatory macrophage" = 'Acod1',
  "Inflammatory macrophage" = 'Il1rn',
  
  # Interstitial macrophages
  "Interstitial macrophages" = "Mrc1",
  "Interstitial macrophages" = "Lyve1",
  "Interstitial macrophages" = "F13a1",
  "Interstitial macrophages" = "Cd163",   

  "Monocyte-derived macrophages" = 'Ccr2',
  "Monocyte-derived macrophages" = 'Ccr5',
  # "Monocyte-derived macrophages" = 'Ccl2',
  "Monocyte-derived macrophages" = 'Ccl12',
  "Monocyte-derived macrophages" = 'Ccl7'
)

pdf(file.path(Macrophage_out_dir, "dotplot celltype main_marker.pdf"), width = 20, height = 6)
DotPlot(Macrophage_pure, group.by = "subtype", features = Macrophage_cellTye_markers, dot.scale = 10) + theme_90angle
dev.off()


saveRDS(Macrophage, file.path(Macrophage_out_dir, 'Macrophage.rds'))
saveRDS(Macrophage_pure, file.path(Macrophage_out_dir, 'Macrophage_pure.rds'))




# 映射 ------------
library(tibble)
library(data.table)
subtype <- rbind(Macrophage[['subtype']] %>% rownames_to_column('cells'))
all_data$annotation %>% table()

meta.data.t <-
  subtype%>%
  right_join(
    all_data[['annotation']] %>% rownames_to_column('cells'),
    by = "cells"
  ) %>% 
  mutate(
    annotation_full = ifelse((as.character(annotation) %in% c('Macrophage')),
                             as.character(subtype), 
                             as.character(annotation)))%>%
  dplyr::select(cells, annotation_full) %>%
  column_to_rownames('cells')


all_data<-AddMetaData(all_data,meta.data.t['annotation_full'])
all_data$annotation %>% table()
all_data$annotation_full %>% table()
all_data$annotation %>% table() %>% sum()
all_data$annotation_full %>% table() %>% sum()




# spatial plot ----------
dir.create(file.path(Macrophage_out_dir, 'spatial_plot'))
spatial_theme <- 
  theme(panel.grid.major = element_blank(), #主网格线
        panel.grid.minor = element_blank(), #次网格线
        panel.border = element_rect(color = "black", size = 1, fill = NA), #边框
        #axis.title = element_blank(),  #轴标题
        axis.text = element_blank(), # 文本
        axis.title = element_text(size = 15),
        axis.ticks = element_blank(),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 20),
        panel.background = element_rect(fill = 'white'), #背景色
        plot.title = element_text(hjust = 0.5, size = 20))


allcolour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#8B4513","#F08080","#1E90FF","#7CFC00","#FFFF00",
            "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB",
            "#40E0D0","#5F9EA0","#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4",
            "#32CD32","#F0E68C","#FFFFE0","#EE82EE","#FF6347","#6A5ACD","#9932CC","#8B008B","#DEB887")
cells_boundary_ls <- readRDS('R_OUTPUT/data_for_analysis/cells_boundary_ls.rds')
all_data_meta <- all_data[[]] %>% rownames_to_column('cell_id')


p_dimPlot_ls <- 
  cells_boundary_ls %>%
  imap(function(cells_boundary, sample) { # browser()
    # if (!x %in% colnames(all_data_meta)) return(NULL)
    
    dat <-
      cells_boundary %>%
      right_join(all_data_meta, by = "cell_id")
    
    p <- dat %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) + 
      geom_polygon(aes(fill = annotation_full, group = cell_id)) + 
      ggtitle(sample) + 
      # ggrastr::geom_point_rast(size = 0.05 , alpha = 1, raster.dpi = 300)  +  
      scale_fill_manual(values = allcolour) +
      guides(color = guide_legend(override.aes = list(size = 3, alpha = 0.7))) + 
      coord_equal() +
      scale_y_reverse() + 
      spatial_theme
    
    p
  })

p_dimPlot_ls %>% 
  iwalk(function(p, grp){
    pdf(file.path(Macrophage_out_dir, 'spatial_plot', 
                  str_glue("spatial dim plot {grp}_subType.pdf")), 
        height = 10, 
        width = 15)
    print(p)
    dev.off()
  })   




## highlight -----------
### Macrophage ----------
all_data$annotation_full %>% table()
all_data$annotation4 <- fct_collapse(
  all_data$annotation_full, 
  'Others' = unique(all_data$annotation_full)[!unique(all_data$annotation_full) %in% c("Profibrotic macrophage", "Alveolar macrophage", "Monocyte-derived macrophage", "Interstitial macrophage", "Inflammatory macrophage")]
)
all_data$annotation4 %>% table()
all_data_meta <- all_data[[]] %>% rownames_to_column('cell_id')

p_dimPlot_ls <- 
  cells_boundary_ls %>%
  imap(function(cells_boundary, sample) { # browser()
    # if (!x %in% colnames(all_data_meta)) return(NULL)
    
    dat <-
      cells_boundary %>%
      right_join(all_data_meta, by = "cell_id")
    
    p <- dat %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) + 
      geom_polygon(aes(fill = annotation4, group = cell_id)) + 
      ggtitle(sample) + 
      # ggrastr::geom_point_rast(size = 0.05 , alpha = 1, raster.dpi = 300)  +  
      scale_fill_manual(values = c("lightgrey","#0000FF","#7CFC00","#FFFF00","#7B68EE","#D2691E")) +
      guides(color = guide_legend(override.aes = list(size = 3, alpha = 0.7))) + 
      coord_equal() +
      scale_y_reverse() + 
      spatial_theme
    
    p
  })

p_dimPlot_ls %>% 
  iwalk(function(p, grp){
    pdf(file.path(Macrophage_out_dir, 'spatial_plot',   
                  str_glue("spatial dim plot {grp}_macrophage highlight.pdf")), 
        height = 10, 
        width = 10)
    print(p)
    dev.off()
  })  




# cell prop ------------
table(Macrophage_pure$sample, Macrophage_pure$subtype) / 263.22
p_col_proj <-
  Macrophage_pure[[]] %>%
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
  filename = file.path(Macrophage_out_dir, "cell_freq_sample.pdf"),
  height = 6,
  width  = 8
)


p_col_proj <-
  Macrophage_pure[[]] %>%
  group_by(group, subtype) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = group, y = freq, fill = subtype)) + 
  geom_col() +
  ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )

# export
ggsave(
  p_col_proj,
  filename = file.path(Macrophage_out_dir, "cell_freq_group.pdf"),
  height = 6,
  width  = 8
)




# cluster prop ----------
p_col_proj <-
  Macrophage_pure[[]] %>%
  group_by(sample, RNA_snn_res.0.8) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = sample, y = freq, fill = RNA_snn_res.0.8)) + 
  geom_col() +
  ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )


# export
ggsave(
  p_col_proj,
  filename = file.path(Macrophage_out_dir, "cluster_freq_sample.pdf"),
  height = 6,
  width  = 8
)


p_col_proj <-
  Macrophage[[]] %>%
  group_by(sample, RNA_snn_res.0.8) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = sample, y = freq, fill = RNA_snn_res.0.8)) + 
  geom_col() +
  ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )


# export
ggsave(
  p_col_proj,
  filename = file.path(Macrophage_out_dir, "cluster_freq_sample_raw.pdf"),
  height = 6,
  width  = 8
)




# Annotation info output -------------------------------------------------------
all_data$annotation %>% table()
cellType_anno_ls <- 
  all_data@meta.data[c('project', 'annotation')] %>% 
  rownames_to_column('Barcode') %>% 
  mutate(
    cell_id = str_replace(Barcode, '^.*_', ''), 
    group = annotation
  ) %>% 
  dplyr::select(cell_id, project, group) %>% 
  group_by(project) %>%
  nest() %>% 
  deframe()

cellType_anno_ls %>% 
  iwalk(function(cellType, project){
    write.csv(cellType, str_glue('anno_info/cell type {project}.csv'), row.names = F)
  })


all_data$annotation_full %>% table()
subType_anno_ls <- 
  all_data@meta.data[c('project', 'annotation_full')] %>% 
  rownames_to_column('Barcode') %>% 
  mutate(
    cell_id = str_replace(Barcode, '^.*_', ''), 
    group = annotation_full
  ) %>% 
  dplyr::select(cell_id, project, group) %>% 
  group_by(project) %>%
  nest() %>% 
  deframe()

subType_anno_ls %>% 
  iwalk(function(subType, project){
    write.csv(subType, str_glue('anno_info/sub cell type {project}.csv'), row.names = F)
  })
