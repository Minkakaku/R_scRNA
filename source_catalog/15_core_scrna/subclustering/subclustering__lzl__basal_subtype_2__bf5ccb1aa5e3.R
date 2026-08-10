Basal_data <- readRDS(file.path(Basal_out_dir, 'Basal_data_all.rds'))
Basal_out_dir <- "Basal_subtype_all_harmony" 
dir.create(Basal_out_dir)
dims=15
###
pdf(file.path(Basal_out_dir, str_glue("Dimplot dim{dims} sample - harmony.pdf")), height = 6, width = 8)
DimPlot(Basal_data, group.by = "sample", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()

pdf(file.path(Basal_out_dir, str_glue("Dimplot dim{dims} sample.pdf")), height = 6, width = 8)
DimPlot(Basal_data, group.by = "sample", label = TRUE, reduction = "umap") + NoLegend()
dev.off()

########
df <- Basal_data@meta.data
p_col_proj <-
  df %>%
  group_by(sample, RNA_snn_res.0.6) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = sample, y = freq, fill = RNA_snn_res.0.6)) + 
  geom_col() +
  #scale_fill_manual(values = basal_colors) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )

# export
ggsave(
  p_col_proj,
  filename = file.path("Basal_subtype_all_harmony/cluster_freq_sample.pdf"),
  height = 6,
  width  = 8
)
common_clusters <- Basal_data@meta.data %>%
  group_by(RNA_snn_res.0.6) %>%
  summarise(n_samples = n_distinct(sample)) %>%  # 计算每个cluster跨越了多少个样本
  filter(n_samples == n_distinct(Basal_data$sample)) # 筛选出样本数等于总样本数的cluster

# 查看结果
print(common_clusters)

cluster_sample_stats <- Basal_data@meta.data %>%
  group_by(RNA_snn_res.0.6) %>%                     # 按照 Cluster 分组
  summarise(n_samples = n_distinct(sample)) %>%     # 统计每个 Cluster 包含的不重复样本数
  arrange(desc(n_samples))                          # 按样本数降序排列，方便查看

# 查看统计结果
print(cluster_sample_stats)
# dimplot
pdf(file.path(Basal_out_dir, str_glue("Dimplot dim{dims} cluster.pdf")), height = 6, width = 8)
DimPlot(Basal_data, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Basal_data, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Basal_data, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Basal_data, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(Basal_data, group.by = "RNA_snn_res.1", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(Basal_data, group.by = "RNA_snn_res.1.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(Basal_data, group.by = "RNA_snn_res.1.5", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()


pdf(file.path(Basal_out_dir, str_glue("Dimplot dim{dims} cluster noharmony.pdf")), height = 6, width = 8)
DimPlot(Basal_data, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = 'umap') + NoLegend()
DimPlot(Basal_data, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = 'umap') + NoLegend()
DimPlot(Basal_data, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = 'umap') + NoLegend()
DimPlot(Basal_data, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = 'umap') + NoLegend()
# DimPlot(Basal_data, group.by = "RNA_snn_res.1", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(Basal_data, group.by = "RNA_snn_res.1.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(Basal_data, group.by = "RNA_snn_res.1.5", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()
###########
#####剔除两个复发的大部分细胞
Basal_data_use <- subset(Basal_data,`RNA_snn_res.0.6`%in%c("11","13","12","8"))
table(Basal_data_use$sample)

## PCA --------
Basal_data_use <- FindVariableFeatures(Basal_data_use)

options(future.globals.maxSize = 20 * 1024 ^ 3) # 20Gb
Basal_data_use <- ScaleData(Basal_data_use, vars.to.regress = 'nCount_RNA')
Basal_data_use <- RunPCA(Basal_data_use)

### pca inspector ------
pdf(file.path(Basal_out_dir, 'pca plot use.pdf'), height = 6, width = 8)
ElbowPlot(Basal_data_use, ndims = 50)
dev.off()


Basal_data_use$subtype %>% table()
#Basal_data$celltype %>% table()
## dimension reduction --------
dims <- 10 # 18

library(future)
plan(multisession, workers = 10)
options(future.globals.maxSize = 20 * 1024 ^ 3) # 20Gb

Basal_data_use <- RunUMAP(Basal_data_use, dims = 1:dims,return.model = TRUE)

Basal_data_use <- Basal_data_use %>% 
  harmony::RunHarmony(
    dims.use = 1:dims,
    group.by.vars = "project",
    # theta = 4, # default 2, bigger value -- more diversity per cluster
    # sigma = 0.5,
    # lambda = 0.5, # default, smaller value -- more integrated
    reduction.save = "harmony"
  )

Basal_data_use <- RunUMAP(
  Basal_data_use, 
  dims = 1:dims, 
  reduction = "harmony", 
  reduction.name = "harmonyumap", 
  reduction.key  = "harmonyUMAP_",
  return.model = TRUE
)

pdf(file.path(Basal_out_dir, str_glue("Dimplot dim{dims} sample  - harmony use.pdf")), height = 6, width = 8)
DimPlot(Basal_data_use, group.by = "sample", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()

pdf(file.path(Basal_out_dir, str_glue("Dimplot dim{dims} sample use.pdf")), height = 6, width = 8)
DimPlot(Basal_data_use, group.by = "sample", label = TRUE, reduction = "umap") + NoLegend()
dev.off()

### cluster -------
Basal_data_use  <-
  FindNeighbors(
    Basal_data_use ,
    reduction = reduction.use,
    dims = 1:dims,
    k.param = 20
  )


Basal_data_use <-
  FindClusters(
    Basal_data_use,
    resolution = c(0.2, 0.4, 0.6, 0.8)
  )



# dimplot
pdf(file.path(Basal_out_dir, str_glue("Dimplot dim{dims} cluster use.pdf")), height = 6, width = 8)
DimPlot(Basal_data_use, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Basal_data_use, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Basal_data_use, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Basal_data_use, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(Basal_data, group.by = "RNA_snn_res.1", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(Basal_data, group.by = "RNA_snn_res.1.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(Basal_data, group.by = "RNA_snn_res.1.5", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()
###########
Idents(Basal_data_use) <- "RNA_snn_res.0.4"

###marker_ls
# 2. 过滤 marker 列表，只保留存在的基因
Basal_markers <-
  FindAllMarkers(
    JoinLayers(Basal_data_use),
    test.use = "wilcox",
    logfc.threshold = 0,
    min.pct = 0,
    only.pos = TRUE
  )
write.csv(Basal_markers, file = "Basal_subtype_all_harmony/Basal_all_markers.csv", row.names = TRUE)
Basal_markers <- read.csv("Basal_subtype_all_harmony/Basal_all_markers.csv")

Basal_top5_markers <- 
  Basal_markers %>% 
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

names(Basal_top5_markers) <- sub('.*\\.','',names(Basal_top5_markers))


# dotplot
# 如果已做其他细胞类型的亚群分析，这里就可以顺便画一下其他类型的细胞亚群基因，比如Fibro_major_marker_ls2，否则只用画T细胞
pdf(file.path(Basal_out_dir, str_glue("dotplot dims{dims} major_marker use.pdf")), width = 25, height = 6)
DotPlot(Basal_data_use, group.by = "RNA_snn_res.0.4", features = major_marker_ls)+RotatedAxis()+theme_90angle
DotPlot(Basal_data_use, group.by = "RNA_snn_res.0.4", features = Epi_major_marker_ls)+RotatedAxis()+theme_90angle
#DotPlot(Macrophage, group.by = "RNA_snn_res.0.4", features = Macrophage_major_marker_ls_1)+RotatedAxis()+theme_90angle
#DotPlot(Macrophage, group.by = "RNA_snn_res.0.6", features = Macrophage_markers)+RotatedAxis()
DotPlot(Basal_data_use, group.by = "RNA_snn_res.0.4", features = Basal_top5_markers)+RotatedAxis()+theme_90angle
#DotPlot(Macrophage, group.by = "RNA_snn_res.0.6", features = Fibro_major_marker_ls2)+RotatedAxis()
dev.off()


pdf(file.path(Basal_out_dir, "ncount_vlnplot.pdf"), width = 8, height = 4)
VlnPlot(Basal_data_use, features = c("nCount_RNA", "nFeature_RNA"), group.by = "RNA_snn_res.0.2", pt.size = 0)

VlnPlot(Basal_data_use, features = c("TP63"), group.by = "subtype", pt.size = 0)
dev.off()
#####
cell_cluster <- c(
  '0'  = "Basal_0", # 
  '1'  = "Basal_1", #
  '2'  = "Basal_2", #
  '3'  = "Basal_3", # 
  '4'  = "AT1/AT2-like", # 
  '5'  = "Pulmonary neuroendocrine cell", # 
  '6'  = "Unknown"
)

Idents(Basal_data_use) <- "RNA_snn_res.0.4"



Basal_data_use <- RenameIdents(Basal_data_use, cell_cluster)
Basal_data_use$subtype <- Idents(Basal_data_use)


pdf(file.path(Basal_out_dir, str_glue("Dimplot dim{dims} annotation use.pdf")), height = 6, width = 8)
DimPlot(Basal_data_use, group.by = "subtype", label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

pdf(file.path(Basal_out_dir, str_glue("Dimplot dim{dims} annotation split use.pdf")), height = 6, width = 15)
DimPlot(Basal_data_use, group.by = "subtype",split.by = 'group', label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

saveRDS(Basal_data_use, file.path(Basal_out_dir, 'Basal_data_use.rds'))
############
####
Basal_clean <- Basal_data_use[, !Basal_data_use$subtype %in% c('Unknown',"AT1/AT2-like","Pulmonary neuroendocrine cell")]
Basal_clean <- RunUMAP(Basal_clean, dims = 1:dims, return.model = TRUE)


pdf(file.path(Basal_out_dir,str_glue("Dimplot dim{dims} annotation - clean use.pdf")), height = 6, width = 8)
DimPlot(Basal_clean, group.by = "subtype", label = TRUE, repel = TRUE) + NoLegend()
dev.off()


Basal_markers <-
  FindAllMarkers(
    JoinLayers(Basal_clean),
    test.use = "wilcox",
    logfc.threshold = 0,
    min.pct = 0,
    only.pos = TRUE
  )
write.csv(Basal_markers, file = "Basal_subtype_all_harmony/Basal_clean_markers.csv", row.names = TRUE)

Basal_top5_markers <- 
  Basal_markers %>% 
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  dplyr::filter(pct.1 > 0.3) %>%  
  dplyr::filter(avg_log2FC > 0.5) %>% # 核心修改：限制表达比例大于 30%
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
        names(gene) = str_glue("cluster_{name}")
        gene
      })
    re
  }) %>%
  unlist() %>% 
  .[!duplicated(.)]

names(Basal_top5_markers) <- sub('.*\\.','',names(Basal_top5_markers))

#names(Basal_top5_markers) <- sub('.*\\.','',names(Basal_top5_markers))

pdf(file.path(Basal_out_dir, str_glue("dotplot dims{dims} major_marker use.pdf")), width = 25, height = 6)
DotPlot(Basal_clean, group.by = "subtype", features = major_marker_ls)+RotatedAxis()+theme_90angle
DotPlot(Basal_clean, group.by = "subtype", features = Epi_major_marker_ls)+RotatedAxis()+theme_90angle
#DotPlot(Macrophage, group.by = "RNA_snn_res.0.4", features = Macrophage_major_marker_ls_1)+RotatedAxis()+theme_90angle
#DotPlot(Macrophage, group.by = "RNA_snn_res.0.6", features = Macrophage_markers)+RotatedAxis()
DotPlot(Basal_clean, group.by = "subtype", features = Basal_top5_markers)+RotatedAxis()+theme_90angle
#DotPlot(Macrophage, group.by = "RNA_snn_res.0.6", features = Fibro_major_marker_ls2)+RotatedAxis()
dev.off()


cell_cluster <- c(
  '0'  = "EIF4G1+ AKR1C3+ basal", # 
  '1'  = "H3F3B+ basal", #
  '2'  = "SLC2A1+ ALDOA+ basal", #
  '3'  = "Proliferating basal"
)

Idents(Basal_clean) <- "RNA_snn_res.0.4"
Idents(Basal_clean) <- "subtype"

#####
Basal_clean <- RenameIdents(Basal_clean, cell_cluster)
Basal_clean$subtype <- Idents(Basal_clean)


pdf(file.path(Basal_out_dir, str_glue("Dimplot dim{dims} annotation clean.pdf")), height = 6, width = 8)
DimPlot(Basal_clean, group.by = "subtype", label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

pdf(file.path(Basal_out_dir, str_glue("Dimplot dim{dims} annotation split clean.pdf")), height = 6, width = 15)
DimPlot(Basal_clean, group.by = "subtype",split.by = 'group', label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()
pdf(file.path(Basal_out_dir, str_glue("Dimplot dim{dims} annotation split clean group_3.pdf")), height = 6, width = 15)
DimPlot(Basal_clean, group.by = "subtype",split.by = 'group_3', label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()


saveRDS(Basal_clean, file.path(Basal_out_dir, 'Basal_data_clean.rds'))

#################
df = Basal_clean[[]]
df$subtype <- as.character(df$subtype)
df$subtype <- as.factor(df$subtype)

#df <- df[df$group_3%in%c("A_P","A_R","A_S"),]
summary <- table(df[,c('subtype','group_3')])
roe <- as.data.frame(ROIE(summary))
long_df <- tidyr::gather(roe, key = "group_3", value = "Ro/e", colnames(roe)[1:ncol(roe)-1])
write.csv(roe, file = "/mnt/NexData_file/home/lzl/20260515_XE2510316140-03/Basal_subtype_all_harmony/Roe_group_3.csv",quote=F)


p = pheatmap::pheatmap(roe, 
                       color = colorRampPalette(brewer.pal(n = 7, name = "Reds"))(10),
                       cluster_rows = F, 
                       cluster_cols = F, 
                       display_numbers = TRUE,
                       number_color = "black",
                       cellwidth = 30, cellheight = 16, fontsize = 10,                    # fontsize字号大小
                       border_color = '#ffffff',
                       angle_col='45')
p$gtable$heights[2] = unit(15,"bigpts")

p$gtable[[2]][p$gtable[[2]]$name =="legend",1:4] = c(4,5,4,6)
p$gtable <- gtable_add_grob(p$gtable, textGrob("Ro/e",x = unit(0.32,"npc"),y = unit(0.7,"npc")), t = 2, l = 5, r= 5, b=2)
#pdf("heatmap.pdf",width = ncol(roe)*1.5,height = nrow(roe)*0.3)
pdf("/mnt/NexData_file/home/lzl/20260515_XE2510316140-03/Basal_subtype_all_harmony/heatmap_group_3.pdf",width = ncol(roe)*1.5,height = nrow(roe)*0.6)
p
dev.off()


summary_t <- as.data.frame(summary)
summary_prop <- summary_t %>%
  group_by(group_3) %>%             # 按分组计算
  mutate(total = sum(Freq), prop = Freq/total) %>%  # 计算组内总数、占比
  ungroup()


basal_colors <- c(
  "H3F3B+ basal"         = "#4A7C59",   # 深绿 - 基础/正常
  "SLC2A1+ ALDOA+ basal"         = "#F77F00",   # 橙红 - 侵袭/危险
  "EIF4G1+ AKR1C3+ basal"       = "#9B5DE5",   # 紫色 - 淋巴/特殊微环境
  "Proliferating basal"         = "#E63946"    # 鲜红 - 强增殖/高活力 (推荐)
)
p2 <- ggplot(summary_prop, aes(x = group_3, y = prop, fill = subtype)) +
  geom_col(position = "stack", width = 0.7) +  # 堆叠柱状图
  scale_y_continuous(labels = scales::percent_format(), expand = c(0,0)) + # y轴转为百分比
  # 配色
  scale_fill_manual(values = basal_colors) +
  # 坐标轴、标签主题设置
  labs(
    x = "Group",
    y = "Cell proportion",
    fill = "Cell type"
  ) +
  # 论文级主题美化
  theme_classic() +
  theme(
    axis.title = element_text(size = 12, color = "black"),
    axis.text = element_text(size = 11, color = "black"),
    legend.title = element_text(size = 13, color = "black"),
    legend.text = element_text(size = 12, color = "black"),
    legend.position = "right",
    plot.margin = margin(10,20,10,10)
  )

pdf("/mnt/NexData_file/home/lzl/20260515_XE2510316140-03/Basal_subtype_all_harmony/subtypetype_prop.pdf",width = 8,height = 8)
p2
dev.off()

##
p_col_proj <-
  df %>%
  group_by(sample, subtype) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = sample, y = freq, fill = subtype)) + 
  geom_col() +
  scale_fill_manual(values = basal_colors) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )


# export
ggsave(
  p_col_proj,
  filename = file.path("Basal_subtype_all_harmony/cell_freq_sample.pdf"),
  height = 6,
  width  = 8
)

p_col_proj <-
  df %>%
  group_by(group_3, subtype) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = group_3, y = freq, fill = subtype)) + 
  geom_col() +
  scale_fill_manual(values = basal_colors) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )


ggsave(
  p_col_proj,
  filename = file.path("Basal_subtype_all_harmony/cell_freq_group3.pdf"),
  height = 6,
  width  = 4
)

###
out_dir <- "spatial plot"
dir.create(out_dir)

all_data <- readRDS("/mnt/NexData_file/home/lzl/20260515_XE2510316140-03/data/all_data.rds")
all_data$subtype <- as.character(all_data$annotation)  # 或 NA
# 提取 KC_data 的细胞名和对应注释
Basal_cells <- colnames(Basal_clean)
Basal_subtype <- Basal_clean$subtype  # 如 KC_data$cell_type 或 KC_data$RNA_snn_res.0.4

unique(Basal_clean$subtype)
# 创建命名向量
subtype_map <- setNames(as.character(Basal_subtype), Basal_cells)

# 赋值到 all_data
all_data$subtype[match(Basal_cells, colnames(all_data))] <- subtype_map[Basal_cells]

# 非 KC 细胞保留原注释
non_kc <- setdiff(colnames(all_data), Basal_cells)
all_data$subtype[non_kc] <- as.character(all_data$annotation[non_kc])
unique(all_data$subtype)



annotation_col <- 'subtype'

cells_position_ls <-
  cbind(
    all_data[['sample']],
    Embeddings(all_data[['spatial']])
  ) %>%
  rownames_to_column('cell_id') %>%
  group_by(sample) %>%
  nest() %>%
  deframe()


plot_dat <-
  FetchData(
    all_data, 
    vars = c(annotation_col)
  ) %>%
  rownames_to_column('cell_id')


all_cell_types <- names(basal_colors)

# 遍历每种细胞类型
for(target_cell_type in all_cell_types){
  
  message("正在绘制：", target_cell_type)
  
  # 强制锁定当前细胞类型和颜色
  current_cell <- target_cell_type
  current_color <- basal_colors[[current_cell]]
  
  p_sp_ct_ls <- cells_position_ls %>%
    imap(function(coords, sam){
      
      ct <- current_cell
      
      dat2 <- coords %>% left_join(plot_dat, by = "cell_id")
      
      # ========== 关键修复1：检查该样本是否包含目标细胞类型 ==========
      cell_types_vec <- dat2 %>% pull(!!sym(annotation_col))
      
      if (!ct %in% cell_types_vec) {
        message("  跳过样本 ", sam, "：无 ", ct, " 细胞")
        return(NULL)
      }
      
      # ========== 关键修复2：统一字符型比较，避免 factor 问题 ==========
      cell_types_char <- as.character(cell_types_vec)
      ct_char <- as.character(ct)
      
      highlight_vec <- ifelse(cell_types_char == ct_char, ct_char, "Other")
      highlight_vec[is.na(highlight_vec)] <- "Other"
      
      dat2$highlight <- highlight_vec
      
      dat2$highlight <- factor(
        dat2$highlight,
        levels = c(ct_char, "Other")
      )
      
      # ========== 绘图 ==========
      fill_values <- setNames(c(current_color, "gray60"), c(ct_char, "Other"))
      alpha_values <- setNames(c(1, 0.3), c(ct_char, "Other"))
      
      p <- ggplot(dat2, aes(x = SPATIAL_1, y = SPATIAL_2)) +
        # 灰色背景点：大小不变
        geom_tile(
          data = dat2 %>% filter(highlight == "Other"),
          fill = "gray60", alpha = 0.3,
          height = 10, width = 10
        ) +
        # 目标细胞类型：放大
        geom_tile(
          data = dat2 %>% filter(highlight != "Other"),
          aes(fill = highlight, alpha = highlight),
          height = 16, width = 16   # ← 调大这个数值
        ) +
        theme_void() +
        coord_equal() +
        scale_y_reverse() +
        scale_alpha_manual(values = setNames(c(1, 0.3), c(ct_char, "Other"))) +
        scale_fill_manual(
          values = setNames(c(current_color, "gray60"), c(ct_char, "Other"))
        ) +
        theme(
          panel.background = element_rect(fill = NA, color = NA),
          plot.title = element_text(hjust = 0.5)
        ) +
        ggtitle(paste0(sam, " (", ct, ")"))
      
      return(p)
    }) %>%
    # ========== 关键修复3：过滤掉 NULL（跳过的样本）==========
  compact()
  
  # 如果所有样本都被跳过，给出提示并继续下一个细胞类型
  if (length(p_sp_ct_ls) == 0) {
    message("  警告：", target_cell_type, " 在所有样本中均未找到，跳过输出")
    next
  }
  
  # 输出PDF（根据样本数调整宽度）
  pdf(
    file.path(out_dir, paste0("spatial plot of ", target_cell_type, ".pdf")),
    height = 8, width = 8 
  )
  walk(p_sp_ct_ls, print)
  dev.off()
  
  message("  完成：", target_cell_type, "，共 ", length(p_sp_ct_ls), " 个样本")
}

library(clusterProfiler)
library(org.Hs.eg.db)


# 按 cluster 分组做富集
marker_list <- Basal_markers %>%
  filter(p_val_adj < 0.05, avg_log2FC > 0.5) %>%
  split(.$cluster)
names(marker_list) <- c("EIF4G1+ AKR1C3+ basal","H3F3B+ basal","SLC2A1+ ALDOA+ basal","Proliferating basal")
enrichGO_results <- imap(marker_list, function(df, cluster_name) {
  message("正在富集：", cluster_name, " (", nrow(df), " genes)")
  
  entrez <- bitr(df$gene, "SYMBOL", "ENTREZID", org.Hs.eg.db)$ENTREZID
  
  if (length(entrez) < 10) return(NULL)
  
  enrichGO(
    gene = entrez,
    OrgDb = org.Hs.eg.db,
    ont = "BP",
    readable = TRUE
  )
}) %>% compact()

names(enrichGO_results)
#enrichGO_results
# 可视化
iwalk(enrichGO_results, function(ego, cluster_name) {
  
  # 检查是否有显著结果
  if (is.null(ego) || nrow(as.data.frame(ego)) == 0) {
    message(cluster_name, " 无显著富集结果")
    return(NULL)
  }
  
  # 点图
  p <- dotplot(ego, showCategory = 15) +
    ggtitle(paste0(cluster_name, " - GO BP Enrichment")) +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.text.y = element_text(size = 10)
    )
  
  # 保存
  ggsave(
    filename = file.path(out_dir, paste0("GO_dotplot_", cluster_name, ".pdf")),
    plot = p,
    width = 10,
    height = 8
  )
  
  message("已保存: ", cluster_name)
})


####
df = all_data[[]]
df$subtype <- as.character(df$subtype)
df$subtype <- as.factor(df$subtype)

p_col_proj <-
  df %>%
  group_by(sample, subtype) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = sample, y = freq, fill = subtype)) + 
  geom_col() +
  scale_fill_manual(values = basal_colors) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )


# export
ggsave(
  p_col_proj,
  filename = file.path("Basal_subtype_all_harmony/all_cell_freq_sample.pdf"),
  height = 6,
  width  = 8
)
###########
pdf(file.path(Basal_out_dir,str_glue("Feature marker group_3.pdf")), height = 6, width = 30)
FeaturePlot(Basal_clean, features = "EIF4G1",split.by = "group_3", label = FALSE, repel = TRUE,reduction = 'harmonyumap')
FeaturePlot(Basal_clean, features = "AKR1C3",split.by = "group_3", label = FALSE, repel = TRUE,reduction = 'harmonyumap')
FeaturePlot(Basal_clean, features = "H3F3B",split.by = "group_3", label = FALSE, repel = TRUE,reduction = 'harmonyumap')
FeaturePlot(Basal_clean, features = "SLC2A1",split.by = "group_3", label = FALSE, repel = TRUE,reduction = 'harmonyumap')
FeaturePlot(Basal_clean, features = "ALDOA",split.by = "group_3", label = FALSE, repel = TRUE,reduction = 'harmonyumap')
FeaturePlot(Basal_clean, features = "DSG3",split.by = "group_3", label = FALSE, repel = TRUE,reduction = 'harmonyumap')
dev.off()


pdf(file.path(Basal_out_dir,str_glue("Feature marker.pdf")), height = 6, width = 8)
FeaturePlot(Basal_clean, features = "EIF4G1", label = FALSE, repel = TRUE,reduction = 'harmonyumap')
FeaturePlot(Basal_clean, features = "AKR1C3", label = FALSE, repel = TRUE,reduction = 'harmonyumap')
FeaturePlot(Basal_clean, features = "H3F3B", label = FALSE, repel = TRUE,reduction = 'harmonyumap')
FeaturePlot(Basal_clean, features = "SLC2A1", label = FALSE, repel = TRUE,reduction = 'harmonyumap')
FeaturePlot(Basal_clean, features = "ALDOA", label = FALSE, repel = TRUE,reduction = 'harmonyumap')
FeaturePlot(Basal_clean, features = "DSG3", label = FALSE, repel = TRUE,reduction = 'harmonyumap')
dev.off()
###########
all_data <- JoinLayers(all_data)
Idents(all_data) <- "group_3"

# 按 group3 的每个组 vs 其余组找差异基因
markers <- FindAllMarkers(
  all_data,
  only.pos = TRUE,      # 只保留上调基因
  min.pct = 0,       # 基因必须在至少25%的细胞中表达
  logfc.threshold = 0,
  test.use = "wilcox"   # 默认 wilcox，也可换 "MAST" 等
)


marker_list <- markers %>%
  filter(p_val_adj < 0.05, avg_log2FC > 0.5) %>%
  split(.$cluster)
enrichGO_results <- imap(marker_list, function(df, cluster_name) {
  message("正在富集：", cluster_name, " (", nrow(df), " genes)")
  
  entrez <- bitr(df$gene, "SYMBOL", "ENTREZID", org.Hs.eg.db)$ENTREZID
  
  if (length(entrez) < 10) return(NULL)
  
  enrichGO(
    gene = entrez,
    OrgDb = org.Hs.eg.db,
    ont = "BP",
    readable = TRUE
  )
}) %>% compact()

names(enrichGO_results)
#enrichGO_results
# 可视化
iwalk(enrichGO_results, function(ego, cluster_name) {
  
  # 检查是否有显著结果
  if (is.null(ego) || nrow(as.data.frame(ego)) == 0) {
    message(cluster_name, " 无显著富集结果")
    return(NULL)
  }
  
  # 点图
  p <- dotplot(ego, showCategory = 15) +
    ggtitle(paste0(cluster_name, " - GO BP Enrichment")) +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.text.y = element_text(size = 10)
    )
  
  # 保存
  ggsave(
    filename = file.path(out_dir, paste0("GO_dotplot_", cluster_name, ".pdf")),
    plot = p,
    width = 10,
    height = 8
  )
  
  message("已保存: ", cluster_name)
})
