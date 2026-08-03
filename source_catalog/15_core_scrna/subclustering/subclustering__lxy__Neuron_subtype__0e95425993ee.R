all_data <- readRDS('R_OUTPUT/all_data.rds')
all_data$annotation1 %>% table()

all_data$annotation %>% table()

Neuron_out_dir <- "Neuron_subtype" 
dir.create(Neuron_out_dir)

# 1. 读取你的CSV文件（细胞ID列表）
Neuron_cells <- read.csv("/mnt/NexData_file/home/lxy/case/XE2601300820-01/stRNA_2/analysis/20260623_analysis/Neuron_subtype/ROI/Neuron.csv", stringsAsFactors = FALSE)

# 2. 提取CSV里的所有细胞ID
csv_cell_ids <- Neuron_cells$Cell.ID

# 3. 提取当前Seurat对象里的所有细胞ID
seurat_cell_ids <- rownames(all_data@meta.data)  # 或 mTEC$cell_id2 都可以


# 4. 取交集
cells_to_keep <- intersect(seurat_cell_ids, csv_cell_ids)

# 5. 只保留这些细胞
all_data_filtered <- subset(all_data, cells = cells_to_keep)

# 7. 保存为新的RDS文件
saveRDS(all_data_filtered, "Neuron_subtype/all_data_filtered_Neuron.rds")

all_data$annotation %>% table()

all_data <- all_data_filtered
all_data$annotation %>% table()

Neuron <- subset(all_data, annotation == 'Neuron')

all_data$annotation %>% table()
all_data$celltype %>% table()
## PCA --------
Neuron <- FindVariableFeatures(Neuron)

Neuron <- ScaleData(Neuron, vars.to.regress = 'nCount_RNA')
Neuron <- RunPCA(Neuron)

### pca inspector ------
pdf(file.path(Neuron_out_dir, 'pca plot.pdf'), height = 6, width = 8)
ElbowPlot(Neuron, ndims = 50)
dev.off()


## dimension reduction --------
dims <- 22 # 18

#library(future)
#options(future.globals.maxSize = 10 * 1024 ^ 3) # 10Gb
#plan(multisession, workers = 10)

Neuron <- RunUMAP(Neuron, dims = 1:dims,return.model = TRUE)

pdf(file.path(Neuron_out_dir, str_glue("Dimplot dim{dims} sample.pdf")), height = 6, width = 8)
DimPlot(Neuron, group.by = "sample", label = TRUE) + NoLegend()
dev.off()

Neuron <- Neuron %>% 
  harmony::RunHarmony(
    dims.use = 1:dims,
    group.by.vars = "project",
    # theta = 4, # default 2, bigger value -- more diversity per cluster
    # sigma = 0.5,
    # lambda = 0.5, # default, smaller value -- more integrated
    reduction.save = "harmony"
  )

Neuron <- RunUMAP(
  Neuron, 
  dims = 1:dims, 
  reduction = "harmony", 
  reduction.name = "harmonyumap", 
  reduction.key  = "harmonyUMAP_",
  return.model = TRUE
)

pdf(file.path(Neuron_out_dir, str_glue("Dimplot dim{dims} sample - harmony.pdf")), height = 6, width = 8)
DimPlot(Neuron, group.by = "sample", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()

reduction.use <- 'harmony' # 'harmony' 'pca'


### cluster -------
Neuron  <-
  FindNeighbors(
    Neuron ,
    reduction = reduction.use,
    dims = 1:dims,
    k.param = 20
  )


Neuron <-
  FindClusters(
    Neuron,
    resolution = c(0.2, 0.4, 0.6, 0.8)
  )



# dimplot
pdf(file.path(Neuron_out_dir, str_glue("Dimplot dim{dims} cluster.pdf")), height = 6, width = 8)
DimPlot(Neuron, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Neuron, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Neuron, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Neuron, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(Neuron, group.by = "RNA_snn_res.1", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(Neuron, group.by = "RNA_snn_res.1.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(Neuron, group.by = "RNA_snn_res.1.5", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()

Idents(Neuron) <- "RNA_snn_res.0.8"

### dotplot -------
Neuron_major_marker_ls <- c(
  "Neuron" = "Tubb3",
  "Neuron" = "Snap25",
  
  "Neuroblast" = "Cd24a",
  "Neuroblast" = "Dcx",
  "Neuroblast" = "Stmn2",
  
  
  # Excitatory neuron 兴奋性神经元 - Glutamatergic neuron 谷氨酸能神经元
  "Glutamatergic_neuron" = "Slc17a6",
  "Glutamatergic_neuron" = "Slc17a7",
  "Glutamatergic_neuron" = "Slc17a8",
  
  # Inhibitory neuron 抑制性神经元 - GABAergic neuron γ-氨基丁酸能神经元
  "GABAergic_neuron" = "Slc32a1",
  #"GABAergic_neuron" = "Slc18a2",
  "GABAergic_neuron" = "Gad1",
  "GABAergic_neuron" = "Gad2",
  "GABAergic_neuron" = "Aldh1a1",
  
  # Inhibitory neuron 抑制性神经元 - Glycinergic neuron 甘氨酸能神经元
  "Glycinergic_neuron" = "Slc6a5",
  
  # Modulatory neuron 调制型神经元 - Cholinergic neuron 胆碱能神经元
  "Cholinergic_neuron" = "Slc18a3",
  "Cholinergic_neuron" = "Chat",
  
  # Modulatory neuron 调制型神经元 - Dopaminergic neuron 多巴胺能神经元
  "Dopaminergic_neuron" = "Slc6a3",
  "Dopaminergic_neuron" = "Slc18a2",
  "Dopaminergic_neuron" = "Th",
  "Dopaminergic_neuron" = "Ddc",
  #"Dopaminergic_neuron" = "Syt1",
  
  # Modulatory neuron 调制型神经元 - Serotonergic neuron 5-羟色胺能神经元
  "Serotonergic_neuron" = "Slc6a4",
  #"Serotonergic_neuron" = "Slc18a2",
  "Serotonergic_neuron" = "Tph2",
  #"Serotonergic_neuron" = "Ddc",
  
  # Modulatory neuron 调制型神经元 - Noradrenergic neuron 去甲肾上腺素能神经元
  "Noradrenergic_neuron" = "Slc6a2",
  #"Noradrenergic_neuron" = "Slc18a2",
  "Noradrenergic_neuron" = "Dbh",
  
  # Modulatory neuron 调制型神经元 - Histaminergic neuron 组胺能神经元
  #"Histaminergic_neuron" = "Slc18a2",
  "Histaminergic_neuron" = "Hdc"
)


Idents(Neuron)<-'RNA_snn_res.0.8'
Neuron_markers <-
  FindAllMarkers(
    JoinLayers(Neuron),
    test.use = "wilcox",
    logfc.threshold = 0,
    min.pct = 0,
    only.pos = TRUE
  )
write.csv(Neuron_markers, file = "Neuron_subtype/Neuron_all_markers.csv", row.names = TRUE)

Neuron_top5_markers <- 
  Neuron_markers %>% 
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

names(Neuron_top5_markers) <- sub('.*\\.','',names(Neuron_top5_markers))


# dotplot
# 如果已做其他细胞类型的亚群分析，这里就可以顺便画一下其他类型的细胞亚群基因，比如Fibro_major_marker_ls2，否则只用画T细胞
pdf(file.path(Neuron_out_dir, str_glue("dotplot dims{dims} major_marker.pdf")), width = 25, height = 10)
DotPlot(Neuron, group.by = "RNA_snn_res.0.8", features = major_marker_ls)+RotatedAxis()+theme_90angle
DotPlot(Neuron, group.by = "RNA_snn_res.0.8", features = Neuron_major_marker_ls)+RotatedAxis()+theme_90angle
DotPlot(Neuron, group.by = "RNA_snn_res.0.8", features = Neuron_top5_markers)+RotatedAxis()+theme_90angle
#DotPlot(Neuron, group.by = "RNA_snn_res.0.6", features = Fibro_major_marker_ls2)+RotatedAxis()
dev.off()


saveRDS(Neuron ,file = "Neuron_subtype/Neuron_cluster.rds")
Neuron <- readRDS("Neuron_subtype/Neuron_cluster.rds")
## annotation ----------------------------------------------------------------
cell_cluster <- c(
  '0'  = "GABAergic_neuron", # 
  '1'  = "GABAergic_neuron", #
  '2'  = "GABAergic_neuron", # 
  '3'  = "GABAergic_neuron", # 
  '4'  = "GABAergic_neuron", # 
  '5'  = "GABAergic_neuron",#  
  '6'  = "GABAergic_neuron", #  
  '7'  = "GABAergic_neuron", #  
  '8'  = "Glutamatergic_neuron", #  
  '9'  = "Glutamatergic_neuron", #  
  '10'  = "GABAergic_neuron", #  
  '11'  = "Glutamatergic_neuron", #  
  '12'  = "Glutamatergic_neuron", #  
  '13'  = "Astrocyte-like_neuron", #  
  '14'  = "GABAergic_neuron", #  
  '15'  = "GABAergic_neuron", #  
  '16'  = "Glutamatergic_neuron", #  
  '17'  = "GABAergic_neuron", #  
  '18'  = "Glutamatergic_neuron", #  
  '19'  = "GABAergic_neuron", #  
  '20'  = "Glutamatergic_neuron", #  
  '21'  = "Glutamatergic_neuron", #  
  '22'  = "Oligodendrocyte-like_neuron", #  
  '23'  = "GABAergic_neuron"
)

Idents(Neuron) <- "RNA_snn_res.0.8"


Neuron <- RenameIdents(Neuron, cell_cluster)
Neuron$subtype <- Idents(Neuron)

group_map <- c(
  "CON1"  = "CON",
  "CON2"  = "CON",
  "CON3"  = "CON",
  "CRS4"  = "CRS",
  "CRS5"  = "CRS",
  "CRS6"  = "CRS"
)


Neuron@meta.data <- Neuron@meta.data %>%
  mutate(
    Group = recode(group, !!!group_map)
  )

table(Neuron$Group)
table(Neuron$group)
table(Neuron$subtype)

pdf(file.path(Neuron_out_dir, str_glue("Dimplot dim{dims} annotation.pdf")), height = 6, width = 8)
DimPlot(Neuron, group.by = "subtype", label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

pdf(file.path(Neuron_out_dir, str_glue("Dimplot dim{dims} annotation split.pdf")), height = 6, width = 30)
DimPlot(Neuron, group.by = "subtype",split.by = 'group', label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

pdf(file.path(Neuron_out_dir, str_glue("Dimplot dim{dims} annotation split by group.pdf")), height = 6, width = 15)
DimPlot(Neuron, group.by = "subtype",split.by = 'Group', label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()


keep <- c(
  "Glutamatergic_neuron",
  "GABAergic_neuron"
)

Neuron_clean <- subset(Neuron, subtype %in% keep)

Neuron_clean$subtype <- droplevels(Neuron_clean$subtype)
Idents(Neuron_clean) <- Neuron_clean$subtype

table(Neuron_clean$subtype)
saveRDS(Neuron_clean, "Neuron_subtype/Neuron_clean.rds")

Neuron <- Neuron_clean

Neuron <- RunUMAP(Neuron, dims = 1:{dims},n.neighbors = 50,return.model = TRUE)

pdf(file.path(Neuron_out_dir, str_glue("Dimplot dim{dims} Neuron annotation.pdf")), height = 6, width = 8)
DimPlot(Neuron, group.by = "subtype", label = TRUE, reduction = 'umap',repel = TRUE) + NoLegend()
dev.off()

pdf(file.path(Neuron_out_dir, str_glue("Dimplot dim{dims} Neuron annotation split.pdf")), height = 6, width = 30)
DimPlot(Neuron, group.by = "subtype",split.by = 'group', label = TRUE, reduction = 'umap',repel = TRUE) + NoLegend()
dev.off()

pdf(file.path(Neuron_out_dir, str_glue("Dimplot dim{dims} Neuron annotation split by group.pdf")), height = 6, width = 15)
DimPlot(Neuron, group.by = "subtype",split.by = 'Group', label = TRUE, reduction = 'umap',repel = TRUE) + NoLegend()
dev.off()





Neuron_1 <- readRDS("Neuron_subtype/Neuron_cluster.rds")
## annotation ----------------------------------------------------------------
cell_cluster <- c(
  '0'  = "GABAergic_neuron", # 
  '1'  = "GABAergic_neuron", #
  '2'  = "GABAergic_neuron", # 
  '3'  = "GABAergic_neuron", # 
  '4'  = "GABAergic_neuron", # 
  '5'  = "GABAergic_neuron",#  
  '6'  = "GABAergic_neuron", #  
  '7'  = "GABAergic_neuron", #  
  '8'  = "Glutamatergic_neuron", #  
  '9'  = "Glutamatergic_neuron", #  
  '10'  = "GABAergic_neuron", #  
  '11'  = "Glutamatergic_neuron", #  
  '12'  = "Glutamatergic_neuron", #  
  '13'  = "Astrocyte-like_neuron", #  
  '14'  = "GABAergic_neuron", #  
  '15'  = "GABAergic_neuron", #  
  '16'  = "Glutamatergic_neuron", #  
  '17'  = "GABAergic_neuron", #  
  '18'  = "Glutamatergic_neuron", #  
  '19'  = "GABAergic_neuron", #  
  '20'  = "Glutamatergic_neuron", #  
  '21'  = "Glutamatergic_neuron", #  
  '22'  = "Oligodendrocyte-like_neuron", #  
  '23'  = "GABAergic_neuron"
)

Idents(Neuron_1) <- "RNA_snn_res.0.8"
Neuron_1 <- RenameIdents(Neuron_1, cell_cluster)
Neuron_1$subtype <- Idents(Neuron_1)

table(Neuron_1$subtype)

saveRDS(Neuron_1, "Neuron_subtype/Neuron_allannotation.rds")

## annotation ----------------------------------------------------------------
######下面合并到alldata的注释
all_data$annotation %>% table()
#all_data <- readRDS("output/data/all_data.rds")
# 1. 先把 all_data$annotation 复制到新列 T_annotation
all_data$annotation1 <- as.character(all_data$annotation)

# 2. 创建映射：细胞名 → 精细T亚型
cell_subtype <- Neuron_1$subtype
names(cell_subtype) <- colnames(Neuron_1)

# 3. 只把 "T cell" 替换成精细亚型，其他细胞保持原样
all_data$annotation1[colnames(Neuron_1)] <- as.character(cell_subtype)

# 查看结果（正确）
table(all_data$annotation1)
table(all_data$annotation)
table(Neuron_1$subtype)
saveRDS(all_data, "Neuron_subtype/all_data_Neuron_annotation.rds")




Neuron <- readRDS("Neuron_subtype/Neuron_clean.rds")
#画一个注释完的点图
Neuron_major_marker_ls_2 <- c(
  "Neuron" = "Tubb3",
  "Neuron" = "Snap25",
  
  
  # Excitatory neuron 兴奋性神经元 - Glutamatergic neuron 谷氨酸能神经元
  "Glutamatergic_neuron" = "Slc17a6",
  "Glutamatergic_neuron" = "Slc17a7",
  #"Glutamatergic_neuron" = "Slc17a8",
  
  # Inhibitory neuron 抑制性神经元 - GABAergic neuron γ-氨基丁酸能神经元
  "GABAergic_neuron" = "Slc32a1",
  #"GABAergic_neuron" = "Slc18a2",
  "GABAergic_neuron" = "Gad1",
  "GABAergic_neuron" = "Gad2")

pdf(file.path(Neuron_out_dir, str_glue("dotplot Neuron subtype major_marker.pdf")), width = 12, height = 6)
DotPlot(Neuron, group.by = "subtype", features = Neuron_major_marker_ls_2)+RotatedAxis()+theme_90angle
dev.off()


Neuron$subtype %>% table()
Neuron$sample %>% table()
Neuron$group %>% table()
## cell props barplot  ----------------------------------------------------------
Neuron_meta <- Neuron[[]]

Neuron_meta$sample <- factor(Neuron_meta$sample,levels=c(Neuron_meta$sample %>% unique()))

p_col_proj <-
  Neuron_meta %>%
  group_by(sample, subtype) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = sample, y = freq, fill = subtype)) + 
  geom_col() +
  scale_fill_brewer(palette = "Set2", name = "Neuron Subtype") +
  #ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )


# export
ggsave(
  p_col_proj,
  filename = file.path(Neuron_out_dir, "cell_freq_sample.pdf"),
  height = 6,
  width  = 8
)

p_col_proj <-
  Neuron_meta %>%
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
  filename = file.path(Neuron_out_dir, "cell_freq_group_1.pdf"),
  height = 6,
  width  = 4
)

p_col_proj <-
  Neuron_meta %>%
  group_by(Group, subtype) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = Group, y = freq, fill = subtype)) + 
  geom_col() +
  ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )

ggsave(
  p_col_proj,
  filename = file.path(Neuron_out_dir, "cell_freq_group_2.pdf"),
  height = 6,
  width  = 4
)


#改图
# 提取meta
Neuron_meta <- Neuron[[]]

# 对应配色（保持你指定色值不变）
group_col <- c(
  "CON"  = "#2E86AB",
  "CRS"    = "#A23B72"
)

# ========== 2. 计算样本水平百分比（用于箱线图 & 均值堆叠图） ==========
# 每个样本内各细胞亚群百分比
data_long <- Neuron_meta %>%
  group_by(sample, Group, subtype) %>%
  summarise(Count = n(), .groups = "drop_last") %>%
  group_by(sample) %>%
  mutate(Percentage = Count / sum(Count) * 100) %>%
  ungroup()

# 保存原始百分比数据
write.csv(data_long, 
          file.path(Neuron_out_dir, "Neuron_sample_percentage.csv"), 
          row.names = FALSE)


# ==================== 按细胞类型分组的箱线图 ====================
comp_list <- list(c("CON", "CRS"))
p <- ggplot(data_long, aes(x = Group, y = Percentage, fill = Group)) +
  geom_boxplot(width = 0.7, alpha = 0.8, outlier.shape = 19, outlier.size = 1.5) +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.7, aes(color = Group)) +
  facet_wrap(~subtype, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = group_col) +
  scale_color_manual(values = group_col) +
  labs(title = "Neuron - Cell Type Percentages",
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
ggsave(file.path(Neuron_out_dir,"Neuron_Boxplot.pdf"), p, width = 8, height = 4.5, dpi = 300)
print(p)

names(Neuron)

# 映射到空间切片
Neuron_spatial_out_dir <- "/mnt/NexData_file/home/lxy/case/XE2601300820-01/stRNA_2/analysis/20260623_analysis/Neuron_subtype/spatial" 
dir.create(Neuron_spatial_out_dir)
#all_data <- readRDS("Neuron_subtype/all_data_Neuron_annotation.rds")
all_data$annotation1 %>% table()

cells_boundary_ls <- readRDS("R_OUTPUT/data_for_analysis/cells_boundary_ls.rds")

head(cells_boundary_ls)
tail(cells_boundary_ls)


all_data_meta <- 
  all_data[[]] %>% 
  rownames_to_column('cell_id') 
head(all_data_meta)
tail(all_data_meta)


cn_sp_p_ls <- 
  cells_boundary_ls %>% 
  imap(function(cells_boundary, sample){ # browser()
    anno_data   <- all_data_meta[c("sample", "group")] %>% distinct()
    group_info  <- anno_data$Group[anno_data$sample == sample]
    
    dat <-
      cells_boundary %>%
      right_join(all_data_meta, by = "cell_id")
    
    dat <- dat %>%
      mutate(plot_group = ifelse(annotation1 %in% keep, annotation1, "Others"))
    
    p <- dat %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) +
      geom_polygon(aes(fill = annotation1, group = cell_id))  +
      scale_y_reverse() +
      coord_equal() + 
      # 自定义颜色：为keep列表中的细胞分配颜色，"Others"设为灰色
      scale_fill_manual(
        values = c(
          "GABAergic_neuron" = "#1f77b4",    # 蓝色
          "Glutamatergic_neuron" = "#ff7f0e",         # 橙色
          "Others" = "#dcdcdc"          # 浅灰色
        ),
        limits = c(keep, "Others")) +
      theme_void() +
      labs(fill = NULL, title = sample )
    p
  })

cn_sp_p_ls %>% 
  iwalk(function(p, name){
    if(all(is.na(p$data$vertex_x))) return(NULL)
    if(any(is.na(p$data$vertex_x))) warning("vertex_x has NAs ...", .immediate = T)
    pdf(
      file.path(Neuron_spatial_out_dir, str_glue("Subtype spatial plot {name}.pdf")),
      width  = 20,
      height = 20
    )
    print(p)
    dev.off()
  })

library(patchwork)

# 方法1：将所有图拼接成一个大的网格
combined_plot <- wrap_plots(cn_sp_p_ls, ncol = 4)  # 每行4张图

# 保存
ggsave(
  file.path(Neuron_out_dir, "Subtype spatial plot.pdf"),
  plot = combined_plot,
  width = 40,   # 根据子图数量调整
  height = 40,
  limitsize = FALSE
)

# 组间差异分析和通路富集
#nohup Rscript /mnt/NexData_file/home/lxy/case/XE2603262474-01/stRNA_2/deg_gokegg/deg_enrich.R -i /mnt/NexData_file/home/lxy/case/XE2601300820-01/stRNA_2/analysis/20260623_analysis/Neuron_subtype/Neuron_clean.rds --species mouse -g Group -c subtype --ident1 CRS --ident2 CON -o CRS_vs_CON_DEG &
