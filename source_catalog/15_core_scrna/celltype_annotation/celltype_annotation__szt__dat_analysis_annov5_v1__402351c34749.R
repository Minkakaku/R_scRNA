library(tidyverse)
library(Seurat)

# single cell rds
all_data_sc <- readRDS("data/gastric_sc_st_datas_20251117/scRNA/adata_raw.rds")
all_data_sc[['RNA']] <- as(all_data_sc[['RNA']], "Assay5")
all_data_sc <- NormalizeData(all_data_sc)

all_data_sc$celltype <- all_data_sc$major_celltype
all_data_sc$project <- all_data_sc$sample_id

# spatial rds
all_data_full <- readRDS('all_data_full_v5.rds')

out_dir <- "20251126 anno v5"
sp_output <- file.path(out_dir, "spatial plot")

dir.create(out_dir)
dir.create(sp_output)

# 1. spatial plot of CN ---------------

p_sp_plot_anno <-
  Images(all_data_full) %>%
  setNames(., str_remove(., "\\.[0-9]{3}um$")) %>%
  imap(function(img, img_name){
    p <- SpatialDimPlot(
      all_data_full,
      images = img, 
      group.by = "celltype_final3",
      pt.size.factor	= 4
    ) +
      # ggsci::scale_fill_aaas() +
      ggsci::scale_fill_igv() +
      guides(fill = guide_legend(override.aes = list(size = 5))) +
      ggtitle(img_name)
    p
  })

p_sp_plot_cn <-
  Images(all_data_full) %>%
  setNames(., str_remove(., "\\.[0-9]{3}um$")) %>%
  imap(function(img, img_name){
    p <- SpatialDimPlot(
      all_data_full,
      images = img, 
      group.by = "CN_cellcharter",
      pt.size.factor	= 4
    ) +
      ggsci::scale_fill_aaas() +
      # ggsci::scale_fill_igv() +
      guides(fill = guide_legend(override.aes = list(size = 5))) +
      ggtitle(img_name)
    p
  })
# export
pdf(
  file.path(sp_output, "spatial plot annotation.pdf"),
  height = 8, 
  width  = 8
)
p_sp_plot_anno %>% walk(print)
dev.off()

pdf(
  file.path(sp_output, "spatial plot CN.pdf"),
  height = 8, 
  width  = 8
)
p_sp_plot_cn %>% walk(print)
dev.off()


# 关注一下OLFM4+tumor_cells和pmc_infla和pmc周围细胞的变化
all_data_full$celltype_final4 <- 
  ifelse(
    all_data_full$celltype_final3 == "tumor_cells",
    as.character(all_data_full$celltype_final),
    as.character(all_data_full$celltype_final3)
  )


# 2. proximity analysis ---------------
# Proximity analysis -------
# 参数介绍
# x,y 指定细胞的坐标
# group.by 一般是切片或者样本编号，代表分别计算最近邻关系
# classification 一般是细胞的注释
# k 计算的最近邻的细胞数
# searchtype 邻近度使用radius

# radius是只查看当前细胞的多少um内的细胞
# mpp 坐标与um的转换比例，xenium是1，celldive是0.325
# returntype 可以每个细胞多少radius内的细胞计数count，也可以返回细胞频率freq
# cell_enough_ratio_cutoff 为了降低计算时间，k值不建议设置太高，
# 如果k值未能完全覆盖全部的邻近细胞，则cell_enough_ratio_cutoff会自动加k值，
# 直到默认90%的细胞均已经计算了全部邻近细胞



## 1. input data ------------
cells_position_ls <- 
  all_data_full@images %>%
  map(GetTissueCoordinates) %>%
  map(dplyr::rename, cell_id = cell) %>% 
  setNames(., str_remove(names(.), "\\..*$"))  

all_data_meta <- 
  all_data_full[[]] %>%
  mutate(annotation_full = celltype_final4)

# all_data_cn_meta
all_data_meta_pos <-
  cells_position_ls %>%
  imap(function(cells_position, slide){ # browser()
    
    proj_meta_sub <-
      all_data_meta %>%  
      rownames_to_column('cell_id') %>% 
      # dplyr::filter(group_stage != 'late_cancer') %>% #############
    dplyr::select(annotation_full, cell_id, group, project) %>%
      # rownames_to_column("cell_id") %>%
      dplyr::filter(project == slide)
    
    cells_position[, 1:3] %>%
      inner_join(proj_meta_sub, by = "cell_id") %>%  
      mutate(annotation_full = factor(annotation_full)) %>%
      dplyr::rename(Classification1 = annotation_full,
                    X = x,
                    Y = y,
                    Slide_ID = project)
  }) %>%
  enframe() %>%
  unnest(cols = value) %>%
  dplyr::select(-name) %>%
  column_to_rownames("cell_id")

all_data_meta_pos$Slide_ID <- str_c(all_data_meta_pos$group, all_data_meta_pos$Slide_ID, sep = "_")

all_data_meta_pos_subset <- all_data_meta_pos %>% rownames_to_column("cell_id")

# saveRDS(all_data_meta_pos_subset, file.path(out_dir, 'all_data_meta_pos_subset.rds'))
# readRDS('all_data_meta_pos_subset.rds')

all_data_meta_pos_subset %>% head

# cell_id        X        Y Classification1 group       Slide_ID
# 1 1_aaaaaagg-1 1670.389 185.1358      Macrophage  CD-F X_202413076_32
# 2 1_aaadlghn-1 1681.172 176.4662          Plasma  CD-F X_202413076_32
# 3 1_aaadmplc-1 1677.561 182.6699               T  CD-F X_202413076_32
# 4 1_aaadoogj-1 1677.109 165.6637          Plasma  CD-F X_202413076_32
# 5 1_aaadpljf-1 1671.887 176.8148          Plasma  CD-F X_202413076_32
# 6 1_aaafaiko-1 1719.466 192.8071          Plasma  CD-F X_202413076_32
# ......


## 2. proximity analysis -------
# 参数介绍
# x,y 指定细胞的坐标
# group.by 一般是切片或者样本编号，代表分别计算最近邻关系
# classification 一般是细胞的注释
# k 计算的最近邻的细胞数
# searchtype 邻近度使用radius

# radius是只查看当前细胞的多少um内的细胞
# mpp 坐标与um的转换比例，xenium是1，celldive是0.325
# returntype 可以每个细胞多少radius内的细胞计数count，也可以返回细胞频率freq
# cell_enough_ratio_cutoff 为了降低计算时间，k值不建议设置太高，
# 如果k值未能完全覆盖全部的邻近细胞，则cell_enough_ratio_cutoff会自动加k值，
# 直到默认90%的细胞均已经计算了全部邻近细胞

### 2.1 proximity analysis --------
source('spatial analysis.R')
res_pa <-
  ProximutyAnalysis_NN_multi(
    all_data_meta_pos_subset,
    query = all_data_meta_pos_subset,
    x = "X",
    y = "Y",
    group.by = "Slide_ID",
    classification = "Classification1", # 若为NULL: 提前返回window的细胞近邻的list
    k = 300, # 'auto',
    searchtype = "radius",
    radius = 50, #  50um半径
    mpp = 8 / all_data_full@images[[1]]@scale.factors$spot, # !important, 1 for xenium
    returntype = "count",
    cell_enough_ratio_cutoff = 90
  )
# saveRDS(res_pa, file = file.path(out_dir, "res_pa.rds"))

# return res_pa
# rows: all celltypes
# cols: query cell
res_pa_colsum <- colSums(res_pa)
res_pa <- sweep(res_pa, 2, res_pa_colsum, FUN = "/")
if(any(res_pa_colsum == 0)){
  res_pa[, res_pa_colsum == 0] <- 0
}

# res_pa <- readRDS(file.path(cn_out_dir, 'res_pa.rds'))


# summarized by group and celltype
pa_groups <-
  all_data_meta_pos_subset %>%
  group_by(Slide_ID) %>%
  nest() %>%
  deframe() %>%
  imap(function(dat, name){
    with(
      dat,
      str_c(Classification1, group, name, sep = '..')
    )
  }) %>%
  unlist(use.names = FALSE)


# pa_groups2 <-
#   all_data_meta_pos_subset %>%
#   left_join(all_data_meta[c('cell_id', 'CN2')], "cell_id") %>%
#   with(
#     str_c(Classification1, group, Slide_ID, CN2, sep = '..')
#   )

# 分组内取细胞比例的中位数
# apply对每行操作时会转置结果
res_pa_df <-  apply(res_pa, 1, function(x){tapply(x, INDEX = pa_groups, FUN = "median")})
# res_pa_df2 <-  apply(res_pa, 1, function(x){tapply(x, INDEX = pa_groups2, FUN = "median")})

res_pa_df  <-
  res_pa_df %>%
  as.data.frame() %>%
  rownames_to_column('group') %>%
  pivot_longer(-group, names_to = "Others", values_to = "cell_props") %>%
  separate(group, into = c("Target", "Group", "project"), sep = "\\.\\.") 


# proximity props / slide props
all_props_df <-
  all_data_meta_pos_subset %>%
  group_by(Slide_ID, Classification1) %>%
  summarise(count = n()) %>%
  group_by(Slide_ID) %>%
  mutate(freq = count/sum(count))


res_pa_df2 <-
  res_pa_df %>%
  left_join(
    all_props_df[c('Slide_ID', 'Classification1', 'freq')],
    by = c(Others = 'Classification1', project = 'Slide_ID')
  ) %>%
  mutate(cell_props2 = cell_props / freq)



# res_pa_df2  <-
#   res_pa_df2 %>%
#   as.data.frame() %>%
#   rownames_to_column('group') %>%
#   pivot_longer(-group, names_to = "Others", values_to = "cell_props") %>%
#   separate(group, into = c("Target", "Group", "project", "CN2"), sep = "\\.\\.") 


# plot boxplot
# Target:中心细胞
# Others: 其他细胞占比

p_res_pa <-
  res_pa_df %>%
  mutate(Group = factor(Group, levels = c(
    "1NAG", "CAG",  "IM",  "LGD", "HGD",  "early_cancer", "late_cancer", "TLS", "other"        
  ))) %>%
  # dplyr::filter(Target == "Epithelial cells") %>%
  ggplot(aes(x = Others, y = cell_props, fill = Group)) +
  geom_boxplot(width = 0.45, outlier.colour = NA, linewidth = 0.1) +
  # geom_jitter(width = 0.25) +
  facet_wrap(~ Target, scales = "free_y", ncol = 2)+
  ggpubr::theme_pubclean() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

p_res_pa2 <-
  res_pa_df2 %>%
  mutate(Group = factor(Group, levels = c(
    "1NAG", "CAG", "IM",  "LGD", "HGD",  "early_cancer", "late_cancer", "TLS", "other"        
  ))) %>%
  # dplyr::filter(Target == "Epithelial cells") %>%
  ggplot(aes(x = Others, y = cell_props2, fill = Group)) +
  geom_boxplot(width = 0.45, outlier.colour = NA, linewidth = 0.1) +
  # geom_jitter(width = 0.25) +
  facet_wrap(~ Target, scales = "free_y", ncol = 2)+
  scale_y_log10() +
  ggpubr::theme_pubclean() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))


# p_res_pa_ls <-
#   res_pa_df2 %>%
#   with(split(res_pa_df2, CN2)) %>%
#   imap(function(dat, name){
#     dat %>%
#   # dplyr::filter(Target == "Epithelial cells") %>%
#   ggplot(aes(x = Others, y = cell_props, fill = Group)) +
#   geom_boxplot(width = 0.45, outlier.colour = NA) +
#   # geom_jitter(width = 0.25) +
#   facet_wrap(~ Target, scales = "free_y", ncol = 2)+
#   ggpubr::theme_pubclean() +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
#   ggtitle(name)
#   })

# export plot
pdf(
  file.path(out_dir, "proximity analysis.pdf"),
  height = 50,
  width  = 30
)
p_res_pa %>% print()
p_res_pa2 %>% print()
dev.off()

pdf(
  file.path(out_dir, "proximity analysis2.pdf"),
  height = 50,
  width  = 30
)
p_res_pa2 %>% print()
dev.off()
# 
# pdf(
#   file.path(cn_out_dir, "proximity analysis - CN.pdf"),
#   height = 50,
#   width  = 25
# )
# p_res_pa_ls %>% print()
# dev.off()




# 3. cellchat analysis ------------------


# 4. single cell cell props -------------
p_col_major <- 
  all_data_sc[[]] %>%
  group_by(group, celltype) %>%
  summarise(count = n()) %>%
  group_by(group) %>%
  mutate(freq = count / sum(count))  %>%
  
  ggplot(aes(x = group, y = freq, fill = celltype)) +
  geom_col() +
  ggsci::scale_fill_igv() +
  RotatedAxis()


p_col_macro_ct <-
  all_data_sc[[]] %>%
  dplyr::filter(major_celltype == "Macrophage") %>% 
  group_by(group, sub_celltype) %>%
  summarise(count = n()) %>%
  group_by(group) %>%
  mutate(freq = count / sum(count))  %>%
  
  ggplot(aes(x = group, y = freq, fill = sub_celltype)) +
  geom_col() +
  ggsci::scale_fill_igv() +
  RotatedAxis()

DotPlot(
  all_data_sc[, all_data_sc$major_celltype == "Macrophage"],
  features = c('CCL5', 'CCR1'),
  group.by = "sub_celltype"
)


p_col_t_ct <-
  all_data_sc[[]] %>%
  dplyr::filter(major_celltype == "T") %>% 
  group_by(group, sub_celltype) %>%
  summarise(count = n()) %>%
  group_by(group) %>%
  mutate(freq = count / sum(count))  %>%
  
  ggplot(aes(x = group, y = freq, fill = sub_celltype)) +
  geom_col() +
  ggsci::scale_fill_igv() +
  RotatedAxis()

DotPlot(
  all_data_sc[, all_data_sc$major_celltype == "T"],
  features = c('CCL5', 'CCR1'),
  group.by = "sub_celltype"
)




# 5. hallmark signature -------------------
h_geneset <- msigdbr::msigdbr(category = "H")

h_geneset$gs_name %>%
  unique() %>% 
  str_to_lower() %>%
  str_remove("^hallmark_")

h_geneset_ls <-
  h_geneset[c('gs_name', 'gene_symbol')] %>%
  with(split(gene_symbol, gs_name)) %>%
  {set_names(
    .,
    names(.) %>% 
      str_to_lower() %>%
      str_remove("^hallmark_")
  )}

all_data_sc_score <-
  all_data_sc %>%
  AddModuleScore(
    features = h_geneset_ls,
    # name = names(h_geneset_ls)
  )

all_data_sc_score %>% head

# msigdb anno
h_score_grp <- 
  all_data_sc_score[[]] %>%
  dplyr::select(group, Cluster1:Cluster50) %>%
  group_by(group) %>%
  summarise(across(everything(), mean))

p_hallmark_1 <-
  pheatmap::pheatmap(
    h_score_grp %>% column_to_rownames("group") %>% t,
    scale = 'row',
    labels_row = names(h_geneset_ls),
    color = colorRampPalette(rev(RColorBrewer::brewer.pal(Inf, "RdBu")))(50),
    border_color = NA,
    silent = TRUE
  )

# export 
pdf(
  file.path(out_dir, "Hallmark pathway score of group - single cell.pdf"), 
  height = 12, 
  width  = 5
)
p_hallmark_1$gtable %>% cowplot::plot_grid() %>% print()
dev.off()


# 6. de gene ------------

de_pmc_infla_vs_pmc <-
  FindMarkers(
    `Idents<-`(all_data_sc, value = "sub_celltype4"),
    ident.1 = 'pmc_infla',
    ident.2 = 'pmc'
  )
de_gene_pmc_infla <-
  rownames(de_pmc_infla_vs_pmc)[
    de_pmc_infla_vs_pmc$p_val_adj < 0.05 &
      de_pmc_infla_vs_pmc$avg_log2FC > 1
  ]

all_de_marker_group <-
  FindAllMarkers(
    `Idents<-`(all_data_sc, value = "group"),
    return.thresh = 1.1
  )
all_de_gene_group_ls <-
  all_de_marker_group %>%
  group_by(cluster) %>%
  nest() %>%
  deframe() %>%
  map(function(dat){
    dat$gene[
      dat$p_val_adj < 0.05 &
        dat$avg_log2FC > 1
    ]
  })

cellcaht_lr <-
  CellChatDB.human %>% 
  subsetDB(
    search = c("Secreted Signaling", "ECM-Receptor", "Cell-Cell Contact"), 
    key    = "annotation"
  )
cellchat_ligand <- cellcaht_lr$interaction$ligand

de_gene_pmc_infla %>%
  intersect(all_de_gene_group_ls$IM) %>%
  intersect(cellchat_ligand)


ggvenn::ggvenn(
  list(
    pmc_infla_vs_pmc = de_gene_pmc_infla,
    de_gene_im = all_de_gene_group_ls$IM,
    cellchat_ligand = cellchat_ligand
  ),
  show_stats = 'c'
)


 




