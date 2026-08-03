library(Seurat)
library(tidyverse)


pdf("sp_subcluster/epi/Dimplot - project.pdf", height = 6 * 2, width = 12 * 2)
DimPlot(
  sp_epi,
  group.by = c("sample", "group", "group2", "group3", "paired_sample", "SCT_snn_res.0.2"),
  reduction = "harmonyumap"
)
dev.off()

pdf("sp_subcluster/epi/Dimplot-dim20.pdf", height = 6 * 2, width = 12 * 2)
DimPlot(
  sp_epi,
  group.by = c("sample", "group", "group2", "group3", "paired_sample", "RNA_snn_res.0.2"),
  reduction = "harmonyumap"
)
dev.off()

pdf("sp_subcluster/epi/Dimplot-dim30.pdf", height = 6 * 2, width = 12 * 2)
DimPlot(
  sp_epi,
  group.by = c("sample", "group", "group2", "group3", "paired_sample", "RNA_snn_res.0.4"),
  reduction = "harmonyumap"
)
dev.off()


pdf("sp_subcluster/epi - v1/sc.1 pca/pca plot.pdf", height = 6, width = 8)
ElbowPlot(sp_epi, ndims = 50)
dev.off()


sp_epi[[]] %>%
  group_by(SCT_snn_res.0.2, group) %>%
  summarise(count = n()) %>%
  mutate(freq = count / sum(count)) %>%
  # dplyr::filter(! SCT_snn_res.0.2 %in% c(11:17)) %>%
  ggplot(aes(x = SCT_snn_res.0.2, y = freq, fill = group)) +
  geom_col()



pdf("sp_subcluster/epi/Dimplot-dim10.pdf", height = 6 * 2, width = 12 * 2)
DimPlot(
  sp_epi,
  group.by = c("sample", "group", "group2", "group3", "paired_sample", "RNA_snn_res.0.4"),
  reduction = "harmonyumap"
)
dev.off()



sp_epi_dim20[[]] %>%
  group_by(RNA_snn_res.0.2, group) %>%
  summarise(count = n()) %>%
  mutate(freq = count / sum(count)) %>%
  dplyr::filter(! RNA_snn_res.0.2 %in% c(11:17)) %>%
  ggplot(aes(x = RNA_snn_res.0.2, y = freq, fill = group)) +
  geom_col()
  
  
sp_epi_dim20$annotation2 <-
  sp_epi_dim20$RNA_snn_res.0.6 %>%
  fct_recode(!!!c(
    Normal = '0',
    LR = '1',
    PH = '10',
    LR = '11',
    PH = '12',
    Tumor = '13',
    Tumor = '2',
    Tumor = '3',
    Tumor = '4',
    Tumor = '5',
    PH = '6',
    Tumor = '7',
    PH = '8',
    Tumor = '9'
  ))

DimPlot(
  sp_epi_dim20,
  reduction = "harmonyumap",
  group.by = "paired_sample"
)
DimPlot(
  sp_epi_dim20,
  reduction = "harmonyumap",
  group.by = "annotation2"
)
  
### 
sp_epi_dim20[[]] %>%
  group_by(RNA_snn_res.0.6, paired_sample) %>%
  summarise(count = n()) %>%
  mutate(freq = count / sum(count)) %>%
  dplyr::filter(! RNA_snn_res.0.6 %in% c(14:17)) %>%
  ggplot(aes(x = RNA_snn_res.0.6, y = freq, fill = paired_sample)) +
  geom_col() + 
  ggsci::scale_fill_igv()


sp_epi_dim20[[]] %>%
  group_by(RNA_snn_res.0.6, group) %>%
  summarise(count = n()) %>%
  mutate(freq = count / sum(count)) %>%
  dplyr::filter(! RNA_snn_res.0.6 %in% c(14:17)) %>%
  ggplot(aes(x = RNA_snn_res.0.6, y = freq, fill = group)) +
  geom_col() + 
  ggsci::scale_fill_igv()



sp_epi_dim20[[]] %>%
  group_by(group, RNA_snn_res.0.6) %>%
  summarise(count = n()) %>%
  mutate(freq = count / sum(count)) %>%
  dplyr::filter(! RNA_snn_res.0.6 %in% c(14:17)) %>%
  ggplot(aes(x = group, y = freq, fill = RNA_snn_res.0.6 )) +
  geom_col() + 
  ggsci::scale_fill_igv()


sp_epi_dim20[[]] %>%
  mutate(group4 = str_c(paired_sample, group, sep = "_")) %>%
  group_by(RNA_snn_res.0.6, group4) %>%
  summarise(count = n()) %>%
  mutate(freq = count / sum(count)) %>%
  dplyr::filter(! RNA_snn_res.0.6 %in% c(14:17)) %>%
  ggplot(aes(x = RNA_snn_res.0.6, y = freq, fill =  group4)) +
  geom_col() + 
  ggsci::scale_fill_igv()



# find markers
# all_markers_res_0.6 <- 
#   FindAllMarkers(
#     sp_epi_dim20,
#     group.by = "RNA_snn_res.0.6"
#   )

all_markers_res_0.6 <- 
  FindAllMarkers(
    sp_epi_dim20,
    group.by = "RNA_snn_res.0.6",
    logfc.threshold = 0,
    return.thresh = 1.1,
    min.pct = 0.0
  )

write.csv(
  all_markers_res_0.6,
  file = "all markers - RNA_snn_res.0.6.csv",
  row.names = FALSE
)



