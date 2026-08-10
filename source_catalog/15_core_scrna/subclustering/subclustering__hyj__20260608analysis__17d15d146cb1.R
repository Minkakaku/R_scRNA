result_dir <- "result20260608/"
dir.create(result_dir)

all_data <- readRDS("R_OUTPUT/all_data.rds")
all_data$sample %>% table()
all_data$annotation %>% table()

# 注释信息导出：一个FOV ---------
annotation_metadata <- all_data[[]] %>%
  rownames_to_column("cell_id") %>%
  mutate(cell_id = str_remove(cell_id, paste0(project, "_"))) %>%
  dplyr::select(cell_id, group = annotation)

write.csv(annotation_metadata, file = file.path(result_dir, "cell_annotation.csv"), row.names = FALSE)

all_data <- readRDS("result20260528/all_data_area_subtype.rds")
# 1、数据拆分肿瘤区域和边界区域--------

# WT组肿瘤
# "WT_1_2"   1-肿瘤核心-2-701（右侧圆圈）-肿瘤组织
# "WT_4_1"   0317-1-2-肿瘤组织
# "WT_5_1"   YH-2-肿瘤组织

# MUT组肿瘤
# "MUT_2_2"  240722高杨年-右侧圆圈-肿瘤组织
# "MUT_3_1"  240722刘丽-上部分-肿瘤组织
# "MUT_1_1"  240522赵诗宇-红线左下-肿瘤组织

table(all_data$group, all_data$sample)
table(all_data$area, all_data$sample)

## 确定肿瘤区域-------
all_data$area <- "other"
all_data@meta.data <- all_data[[]] %>%
  mutate(area = case_when(
    sample %in% c("WT_1_2", "WT_4_1", "WT_5_1", "MUT_2_2", "MUT_3_1", "MUT_1_1") ~ "Tumor",
    TRUE ~ area
  ))

table(all_data$sample, all_data$area)

# #WT组边界
# #WT_1_1
# #WT_4_2
# #WT_5_2
#
# #MUT组边界
# #MUT_2_1
#
## 确定边界区域------
# all_data@meta.data <- all_data[[]] %>%
#   mutate(area = case_when(
#     sample%in%c("WT_1_1","WT_4_2","WT_5_2","MUT_2_1") ~ "Boundary",
#     TRUE ~ area
#   ))
#
# table(all_data$sample,all_data$area)
#
# all_data$area_group <- paste0(all_data$group,"_",all_data$area)
# all_data$area_group %>% table()

# 2、加入拆分区域------
all_data2 <- all_data
all_data2@meta.data$ROI_old <- all_data2@meta.data$ROI
all_data2@meta.data$ROI <- NULL
all_data2$cell_id <- gsub("TMA_", "", rownames(all_data2@meta.data))

## 读取文件夹内所有CSV文件到列表-----
library(tidyverse)
folder_path <- "./result20260608/20260608TumorROI/"

csv_files <- list.files(
  path = folder_path,
  pattern = "\\.csv$",
  full.names = TRUE
)

# 读取所有CSV文件到命名列表
dataframes <- map(csv_files, ~ read_csv(.x, skip = 2)) %>%
  set_names(tools::file_path_sans_ext(basename(csv_files)))

print(head(dataframes[["WT1_2_Tumor2"]]))

## 更新ROI信息-----
all_data2@meta.data$ROI <- "other"

obs <- all_data2@meta.data %>%
  mutate(
    ROI = as.character(ROI),
    cell_id = as.character(cell_id)
  )

walk2(names(dataframes), dataframes, ~ {
  cell_ids <- .y[["Cell ID"]] %>% as.character()
  obs[obs$cell_id %in% cell_ids, "ROI"] <<- .x
})

obs <- obs %>%
  mutate(
    ROI = factor(ROI),
    cell_id = factor(cell_id)
  )
obs$ROI %>% table()

all_data2@meta.data <- obs
all_data2$ROI %>% unique()

all_data2$ROI <- gsub("_1_", "-", all_data2$ROI)
all_data2$ROI <- gsub("_2_", "-", all_data2$ROI)
all_data2$ROI %>% table()

all_data2$area_group %>% table()

all_data2$area_group <- paste0(all_data2$group, "_", all_data2$area)
all_data2$area_group %>% table()

table(all_data2$area_group, all_data2$ROI)

## area spatial plot -----
rownames(all_data2@meta.data) %>% head()
all_data2@meta.data$cell_id %>% head()
all_data2@meta.data$cell_id_raw <- all_data2@meta.data$cell_id
all_data2@meta.data$cell_id <- NULL

all_data_meta <- all_data2@meta.data %>%
  rownames_to_column("cell_id")

#sample_ls <- unique(all_data2$sample_raw)
sample_ls <- c("WT_1_2", "WT_4_1", "WT_5_1", "MUT_2_2", "MUT_3_1", "MUT_1_1")

cells_boundary_ls2 <- readRDS("R_OUTPUT/data_for_analysis/cells_boundary_ls.rds")

sp_p_ls <-
  cells_boundary_ls2 %>%
  # .[1] %>%
  #  .[names(all_data_ls_susbet)] %>%
  .[sample_ls] %>%
  imap(function(cells_boundary, project) { # browser()
    # image_wh_r <- image_wh_r_ls[[project]]
    dat <-
      cells_boundary %>%
      right_join(all_data_meta, by = "cell_id")

    p <- dat %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) +
      geom_polygon(aes(fill = ROI, group = cell_id)) +
      scale_y_reverse() +
      coord_equal() +
      # scale_fill_manual(values = plan_rctd) +
      ggsci::scale_fill_igv() +
      guides(fill = guide_legend(override.aes = list(size = 5, color = "black"))) +
      theme_void()
    p
  })


spatial_dir <- file.path(result_dir, "spatial_plot/")
dir.create(spatial_dir)

# export
sp_p_ls %>%
  iwalk(function(p, grp) {
    pdf(file.path(spatial_dir, str_glue("spatial plot area {grp}.pdf")), height = 18, width = 18)
    print(p)
    dev.off()
  })


# 3、排除六个tumor样本未在ROI内的细胞------
all_data2_clean <- all_data2[
  ,
  !(all_data2$area_group == "MUT_Tumor" & all_data2$ROI == "other") &
    !(all_data2$area_group == "WT_Tumor" & all_data2$ROI == "other")
]

table(all_data2_clean$area_group, all_data2_clean$ROI)

all_data2 <- all_data2_clean
rm(all_data2_clean)

all_data2$sample_raw <- all_data2$sample
all_data2$sample <- all_data2$ROI

table(all_data2$area_group, all_data2$ROI)

saveRDS(all_data2, file.path(result_dir, "all_data_TumorROI.rds"))
rm(all_data)


# 4、各区域间细胞比例统计------
prop_dir <- file.path(result_dir, "cell_prop/")
dir.create(prop_dir)

library(ggplot2)
library(ggpubr)
library(dplyr)

## WT Tumor vs MUT Tumor-----
boxplot_data <- all_data2[[]] %>%
  filter(area == "Tumor") %>%
  group_by(sample, area_group, subtype) %>%
  summarise(count = n()) %>%
  mutate(freq = count / sum(count)) %>%
  mutate(area_group = factor(area_group, levels = c("WT_Tumor", "MUT_Tumor")))

boxplot_data$area_group %>% table()

group_colors <- c(
  "WT_Tumor" = "#1f77b4",
  # "Post-Acq-R" = "#ff7f0e",
  "MUT_Tumor" = "#2ca02c"
)

### 绘制箱线图------
p <- boxplot_data %>%
  ggplot(aes(x = area_group, y = freq, fill = area_group)) +
  geom_boxplot(width = 0.6) +
  facet_wrap(~subtype, scales = "free_y", ncol = 5) +
  stat_compare_means(
    comparisons = list(c("WT_Tumor", "MUT_Tumor")),
    method = "wilcox.test",
    label = "p.signif", # p.signif, "p.value"
    p.adjust.method = "none"
  ) +
  theme_bw() +
  labs(
    title = "",
    x = "area",
    y = "freq"
  ) +
  scale_fill_manual(values = group_colors)

print(p)

ggsave(file.path(prop_dir, "WT_Tumor vs MUT_Tumor boxplot.pdf"), p, width = 15, height = 18)
ggsave(file.path(prop_dir, "WT_Tumor vs MUT_Tumor boxplot.png"), p, width = 15, height = 18)


### 柱状图-分样本------
p <- boxplot_data %>%
  ggplot(aes(x = sample, y = freq, fill = subtype)) +
  geom_col() +
  # facet_wrap(~subtype, scales = "free_y", ncol = 5) +
  ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )
ggsave(file.path(prop_dir, "WT_Tumor vs MUT_Tumor sample_prop.pdf"), p, width = 10, height = 10)
ggsave(file.path(prop_dir, "WT_Tumor vs MUT_Tumor sample_prop.png"), p, width = 10, height = 10)



###柱状图-分组------
boxplot_data <- all_data2[[]] %>%
  filter(area == "Tumor") %>%
  group_by(area_group, subtype) %>%
  summarise(count = n()) %>%
  mutate(freq = count / sum(count)) %>%
  mutate(area_group = factor(area_group, levels = c("WT_Tumor", "MUT_Tumor")))


p <- boxplot_data %>%
  ggplot(aes(x = area_group, y = freq, fill = subtype)) +
  geom_col() +
  # facet_wrap(~subtype, scales = "free_y", ncol = 5) +
  ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )
ggsave(file.path(prop_dir, "WT_Tumor vs MUT_Tumor areagroup_prop.pdf"), p, width = 10, height = 10)
ggsave(file.path(prop_dir, "WT_Tumor vs MUT_Tumor areagroup_prop.png"), p, width = 10, height = 10)



## WT Tumor vs MUT Tumor Tcell_Mac-----
all_data2$subtype %>% table()

target_subtype <- c("M0_Microglia", "M1_Microglia", "M2_Microglia", "M0_Macrophage", "M1_Macrophage", "M2_Macrophage", "DC", "Neutrophil", "CD4_T", "CD8_T")

boxplot_data <- all_data2[[]] %>%
  filter(area == "Tumor") %>%
  group_by(sample, area_group, subtype) %>%
  summarise(count = n()) %>%
  mutate(freq = count / sum(count)) %>%
  filter(subtype %in% target_subtype) %>%
  mutate(area_group = factor(area_group, levels = c("WT_Tumor", "MUT_Tumor")))


group_colors <- c(
  "WT_Tumor" = "#1f77b4",
  # "Post-Acq-R" = "#ff7f0e",
  "MUT_Tumor" = "#2ca02c"
)

### 绘制箱线图--------
p <- boxplot_data %>%
  ggplot(aes(x = area_group, y = freq, fill = area_group)) +
  geom_boxplot(width = 0.6) +
  facet_wrap(~subtype, scales = "free_y", ncol = 2) +
  stat_compare_means(
    comparisons = list(c("WT_Tumor", "MUT_Tumor")),
    method = "wilcox.test",
    label = "p.signif", # p.signif, "p.value"
    p.adjust.method = "none"
  ) +
  theme_bw() +
  labs(
    title = "",
    x = "area",
    y = "freq"
  ) +
  scale_fill_manual(values = group_colors)

print(p)

ggsave(file.path(prop_dir, "WT_Tumor vs MUT_Tumor Tcell_Mac boxplot.pdf"), p, width = 7, height = 15)
ggsave(file.path(prop_dir, "WT_Tumor vs MUT_Tumor Tcell_Mac boxplot.png"), p, width = 7, height = 15)


### 柱状图-分样本------
target_subtype <- c("M0_Microglia", "M1_Microglia", "M2_Microglia", "M0_Macrophage", "M1_Macrophage", "M2_Macrophage", "DC", "Neutrophil", "CD4_T", "CD8_T")

boxplot_data <- all_data2[[]] %>%
  filter(area == "Tumor") %>%
  filter(subtype %in% target_subtype) %>%
  group_by(sample, area_group, subtype) %>%
  summarise(count = n()) %>%
  mutate(freq = count / sum(count)) %>%
  mutate(area_group = factor(area_group, levels = c("WT_Tumor", "MUT_Tumor")))


p <- boxplot_data %>%
  ggplot(aes(x = sample, y = freq, fill = subtype)) +
  geom_col() +
  # facet_wrap(~subtype, scales = "free_y", ncol = 5) +
  ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )
ggsave(file.path(prop_dir, "WT_Tumor vs MUT_Tumor Tcell_Mac sample_prop.pdf"), p, width = 10, height = 10)
ggsave(file.path(prop_dir, "WT_Tumor vs MUT_Tumor Tcell_Mac sample_prop.png"), p, width = 10, height = 10)


### 柱状图-分组------
boxplot_data <- all_data2[[]] %>%
  filter(area == "Tumor") %>%
  filter(subtype %in% target_subtype) %>%
  group_by(area_group, subtype) %>%
  summarise(count = n()) %>%
  mutate(freq = count / sum(count)) %>%
  mutate(area_group = factor(area_group, levels = c("WT_Tumor", "MUT_Tumor")))


p <- boxplot_data %>%
  ggplot(aes(x = area_group, y = freq, fill = subtype)) +
  geom_col() +
  # facet_wrap(~subtype, scales = "free_y", ncol = 5) +
  ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )
ggsave(file.path(prop_dir, "WT_Tumor vs MUT_Tumor Tcell_Mac areagroup_prop.pdf"), p, width = 10, height = 10)
ggsave(file.path(prop_dir, "WT_Tumor vs MUT_Tumor Tcell_Mac areagroup_prop.png"), p, width = 10, height = 10)



#5、ROE分析---------
ROE_result <- file.path(result_dir,"ROE_result/")
dir.create(ROE_result)

library(dplyr)
library(ggplot2)
library(ggsci)
library(Seurat)
library(pheatmap)
library(RColorBrewer)
library(gtable)
library(ggpubr)
source("result20260528/ROE_result/source.R")

##WT Tumor vs MUT Tumor Tcell_Mac-----
target_subtype <- c("M0_Microglia", "M1_Microglia", "M2_Microglia", "M0_Macrophage", "M1_Macrophage", "M2_Macrophage", "DC", "Neutrophil", "CD4_T", "CD8_T")
df = all_data2@meta.data %>% 
  filter(area == "Tumor") %>% 
  filter(subtype %in% target_subtype) %>% 
  mutate(area_group = factor(area_group, levels = c("WT_Tumor","MUT_Tumor")))

summary <- table(df[,c('subtype','area_group')])
roe <- as.data.frame(ROIE(summary))
long_df <- tidyr::gather(roe, key = "area_group", value = "Ro/e", colnames(roe)[1:ncol(roe)-1])

p = pheatmap::pheatmap(roe, 
                       color = colorRampPalette(brewer.pal(n = 7, name = "Reds"))(10),
                       cluster_rows = F, 
                       cluster_cols = F, 
                       display_numbers = TRUE,
                       number_color = "black",
                       cellwidth = 30, cellheight = 16, fontsize = 10,                    # fontsize字号大小
                       border_color = '#ffffff',
                       angle_col='45')
dev.off()
p$gtable$heights[2] = unit(15,"bigpts")

p$gtable[[2]][p$gtable[[2]]$name =="legend",1:4] = c(4,5,4,6)
p$gtable <- gtable_add_grob(p$gtable, text_grob("Ro/e",x = unit(0.32,"npc"),y = unit(0.7,"npc")), t = 2, l = 5, r= 5, b=2)
pdf(file.path(ROE_result,"WT Tumor vs MUT Tumor roe heatmap.pdf"),width = ncol(roe)*1.8,height = nrow(roe)*0.4)
p
dev.off()


##test V3-----
ROE_result <- file.path(result_dir,"ROE_result/")
dir.create(ROE_result)

# AT1 
# roe_out_dir <- "subcluster/AT1"
# 
# anno_column <- "annotation"
# grp_column  <- "project"
# 
# all_data_meta_roe <- 
#   sc_at1[[]] %>%
#   dplyr::filter(.data[[anno_column]]  != "doublet") %>%
#   mutate({{anno_column}} := fct_drop(  .data[[anno_column]] ))

# AT2
target_subtype <- c("M0_Microglia", "M1_Microglia", "M2_Microglia", "M0_Macrophage", "M1_Macrophage", "M2_Macrophage", "DC", "Neutrophil", "CD4_T", "CD8_T")

df = all_data2@meta.data %>% 
  filter(area == "Tumor") %>% 
  filter(subtype %in% target_subtype) %>% 
  mutate(area_group = factor(area_group, levels = c("WT_Tumor","MUT_Tumor")))


anno_column <- "subtype"
grp_column  <- "area_group"

all_data_meta_roe <- df  #sc_at2_2[[]]

# Cell Preference
library(ComplexHeatmap)
# 生成细胞数目分布频数统计表observe.data
observe.data <- unclass(table(all_data_meta_roe[[anno_column]], all_data_meta_roe[[grp_column]]))
stopifnot(is.matrix(observe.data))

expected <-
  function(x){
    if(!is.matrix(x)){stop("Must be a matrix")}
    rtot <- margin.table(x, 1)
    ctot <- margin.table(x, 2)
    tot  <- margin.table(x)
    outer(rtot, ctot, "*")/tot
  }
expected.data <- expected(observe.data) # 生成期望频数表
plot.data <- observe.data/expected.data # 计算实际细胞数目与期望细胞数目的比值
# 绘图数据的列顺序
col.order <- ggplot2:::`%||%`(
  levels(all_data_meta_roe[[grp_column]]),
  unique(all_data_meta_roe[[grp_column]])
) 
plot.data <- plot.data[, col.order] # 对绘图数据列进行排序
plot.data <- plot.data[order(apply(plot.data, 1, which.max), decreasing = T), ] # 对绘图数据行进行排序
# 如果想深色色块在主对角线使用which.min，在副对角线则使用which.max

# 热图颜色设定
roe_cols <- structure(
  c("#FEE6CE", "#FDC08C", "#F5904B", "#E6550D"),
  names = c("1", "1.5", "2", ">3")
)

cell_fun = function(j, i, x, y, width, height, fill) { # browser()
  
  # print(c(j, i, x, y, width, height, fill))
  
  if(plot.data[i,j] < as.numeric(names(roe_cols)[1])) {
    grid::grid.rect(x = x, y = y, width = width, height = height, gp = grid::gpar(col = roe_cols[1], fill = roe_cols[1]))
    grid::grid.text("\u00B1", x = x, y = y)
  } else if(plot.data[i,j] < as.numeric(names(roe_cols)[2])) {
    grid::grid.rect(x = x, y = y, width = width, height = height, gp = grid::gpar(col = roe_cols[2], fill = roe_cols[2]))
    grid::grid.text("+", x = x, y = y)
  } else if(plot.data[i,j] < as.numeric(names(roe_cols)[3])) {
    grid::grid.rect(x = x, y = y, width = width, height = height, gp = grid::gpar(col = roe_cols[3], fill = roe_cols[3]))
    grid::grid.text("++", x = x, y = y)
  } else {
    grid::grid.rect(x = x, y = y, width = width, height = height, gp = grid::gpar(col = roe_cols[4], fill = roe_cols[4]))
    grid::grid.text("+++", x = x, y = y)
  }
  
}

# 生成热图
roi_hm <-
  ComplexHeatmap::Heatmap(
    plot.data, 
    cell_fun = cell_fun,
    cluster_rows = F, 
    cluster_columns = F,
    width = unit(5, "cm"), 
    height = unit(11, "cm"), 
    show_heatmap_legend = F
  )

# stat comparasion
chi2.T <- sum((observe.data - expected.data)^2/expected.data)
p_chisq <- pchisq(q = chi2.T, df = (nrow(observe.data)-1)*(ncol(observe.data)-1), lower.tail = FALSE)
print(p_chisq)

# 生成图例
lgd <-
  ComplexHeatmap::Legend(
    labels = names(roe_cols), 
    title = expression(R[o/e]), 
    legend_gp = grid::gpar(fill = roe_cols)
  )

pdf(
  file.path(ROE_result, "WT Tumor vs MUT Tumor - Roe plot.pdf"), 
  width = 8, 
  height = 6
)
ComplexHeatmap::draw(roi_hm, annotation_legend_list = lgd)
dev.off()

#6、各ROI内注释空间映射-----
source("/mnt/NexData_file/home/hyj/h5ad2rds(1).R")
library(Seurat)
library(SeuratDisk)
library(sceasy)
library(reticulate)
library(dplyr)
library(ggplot2)
library(argparser)
library(patchwork)
library(RColorBrewer)
library(arrow)

all_data <-  JoinLayers(all_data2,layers = 'count')
rownames(all_data@assays$RNA@layers$count) <- rownames(all_data)
colnames(all_data@assays$RNA@layers$count) <- colnames(all_data)
xenium <- CreateSeuratObject(
  all_data@assays$RNA@layers$count,
  assay = "Xenium")

#saveRDS(xenium,"xenium_pro.rds")
#xenium <- readRDS("xenium_pro.rds")

xenium@meta.data <- all_data@meta.data
##筛选肿瘤区域
xenium$sample %>% table()
xenium <- xenium[,xenium$sample != "other"]

###改行名
rownames(xenium@meta.data) <- gsub("TMA_","",rownames(xenium@meta.data))
xenium@meta.data$cell_id <- rownames(xenium@meta.data)
rownames(xenium@meta.data) <- paste0(xenium@meta.data$sample,"_",xenium@meta.data$cell_id)

#umap <- all_data@reductions$umap@cell.embeddings
#rownames(umap) <- rownames(xenium@meta.data)
metadata <- xenium@meta.data
xenium <- Add_meta_umapData_Xenium(object = xenium,umapfile = umap,metedata = metadata)
colnames(xenium@assays$Xenium) <-  colnames(xenium)
row.names(xenium@reductions$umap@cell.embeddings) <- colnames(xenium)
xenium <- RenameCells(xenium,new.names = colnames(xenium))

#saveRDS(xenium,"xenium_renames.rds")
#xenium_renames <- readRDS("xenium_renames.rds")

##所有fov在一个文件夹内-----
has_dt <- requireNamespace("data.table", quietly = TRUE) && requireNamespace("R.utils", quietly = TRUE)
data.dir = "/data/xenium/FWSP20255039/output-XETG00305__0077158__TMA__20251017__072720"
xenium@meta.data$label = xenium@meta.data$sample
centroids = {
  if (has_dt) {
    centroids <- as.data.frame(data.table::fread(file.path(data.dir,
                                                           "cells.csv.gz")))
  } else {
    centroids <- read.csv(file.path(data.dir, "cells.csv.gz"))
  }
  centroids <- data.frame(x = centroids$x_centroid,
                          y = centroids$y_centroid, cell = centroids$cell_id,
                          stringsAsFactors = FALSE)}

segmentations = {
  if (has_dt) {
    segmentations <- as.data.frame(data.table::fread(file.path(data.dir,
                                                               "cell_boundaries.csv.gz")))
  } else {
    segmentations <- read.csv(file.path(data.dir,
                                        "cell_boundaries.csv.gz"), stringsAsFactors = FALSE)
  }
  segmentations <- data.frame(x = segmentations$vertex_x,y = segmentations$vertex_y,
                              cell = segmentations$cell_id)}
mols.qv.threshold = 20
# microns = {
#   if (has_dt) {
#     microns <- data.table::fread(file.path(data.dir, "transcripts.csv"))
#     microns <- subset(microns, qv >= mols.qv.threshold)
#   } else {
#     microns <- read.csv(file.path(data.dir, "transcripts.csv"))
#     microns <- subset(microns, qv >= mols.qv.threshold)
#   }
#   microns <- data.frame(x = microns$x_location, y = microns$y_location,
#                         cell = microns$cell_id,
#                         gene = microns$feature_name, stringsAsFactors = FALSE)
# }

microns = {
  microns <- read_parquet(file.path(data.dir, "transcripts.parquet"))
  microns <- subset(microns, qv >= mols.qv.threshold)
  microns <- data.frame(x = microns$x_location, y = microns$y_location,
                        cell = microns$cell_id,
                        gene = microns$feature_name, stringsAsFactors = FALSE)
}

object <- xenium
assay = "Xenium"
xenium_fov <- function(object,label,assay){
  cells_data <- data.frame(cell_id = object@meta.data$cell_id,sample = object@meta.data$sample)
  cells <- subset(cells_data,sample == label)$cell_id
  features <- rownames(object)
  
  segmentations_sub <- subset(segmentations,cell %in% cells)
  centroids_sub <- subset(centroids,cell %in% cells)
  microns_sub <- subset(microns,gene %in% features & cell %in% cells)
  
  segmentations_sub$cell <- paste0(label,"_",segmentations_sub$cell)
  centroids_sub$cell <- paste0(label,"_",centroids_sub$cell)
  microns_sub$cell <- paste0(label,"_",microns_sub$cell)
  
  segmentations.data <- list(centroids = CreateCentroids(centroids_sub),
                             segmentation = CreateSegmentation(segmentations_sub))
  
  coords <- CreateFOV(coords = segmentations.data, type = c("segmentation",
                                                            "centroids"), molecules = microns_sub, assay = assay)
  return(coords)
  
}

sample_label <- unique(xenium@meta.data$sample)
for (i in sample_label) {
  fov <- xenium_fov(object = xenium,label = i,assay = "Xenium")
  xenium@images[[i]] <- fov
}

future::plan("sequential")
xenium <- NormalizeData(xenium)

names(xenium@images)
###测试一个样本----
DefaultBoundary(xenium@images[["WT1-Tumor1"]]) <- "segmentation"
p <- ImageFeaturePlot(
  object = xenium,
  fov = "WT1-Tumor1",
  nmols = 20000,       # 每个图像中显示的最大分子数（随机子抽样）
  border.color = "white",  # 细胞边界的颜色（深蓝色）
  border.size = 0.1,    # 细胞边界的线宽
  features = 'INHBA',
  max.cutoff = "q90",
  axes = FALSE,
  cols = c('#003366', "#E7475E"),
  size = 1.5,
  dark.background = T
)
ggsave(p,filename = paste0("feature_spatial.png"),width = 12 ,height = 12,dpi = 300)


###Imagedimplot_annotation----
all_colors <- c(
  "#A6CEE3", "#5B9EC9", "#2D82AF", "#7EBA98", "#98D277",
  "#52AF43", "#6F9E4C", "#DD9A88", "#FDBB69", "#FE982C",
  "#F78620", "#D9A295", "#feee7d", "#02c9c9", "#600473",
  "#B294C7", "#7D54A5", "#9E8099", "#F0EB99", "#DBB466",
  "#B15928", # 原始21色
  "#FF6B6B", "#4ECDC4", "#FF9F1C", "#7209B7", "#4361EE", 
  "#F15BB5", "#00BBF9", "#00F5D4", "#9B5DE5", "#00F5A7",
  "#FEE440", "#9C89B8", "#06D6A0"  # 新增13色
)

all_colors8 <- c(
  "#FF6B6B","#A6CEE3","#00F5A7", "#FDBB69","#feee7d", "#9E8099",
  "#00BBF9","#7209B7"
)

xenium$annotation %>% unique()

sample_label <- unique(xenium@meta.data$sample)
for( i in sample_label){
  #j = gsub("_","-",i)
  DefaultBoundary(xenium@images[[i]]) <- "segmentation"
  p <- ImageDimPlot(xenium,
                    fov = i,
                    nmols = 20000,
                    dark.background = F,
                    border.color = "black", 
                    border.size = 0.1,
                    cols = all_colors8,
                    shuffle.cols = F, #随机打乱颜色顺序
                    #size = 1.2,
                    group.by = "annotation"
                    #molecules = c("INHBA"),
                    #mols.size = 1,
                    #mols.cols = "red",
                    #mols.alpha = 0.5, #分子的透明度值
                    
  )
  ggsave(p,filename = file.path(result_dir,paste0("Imagedimplot/",i,"_anno_spatial.png")),width = 12 ,height = 12,dpi = 300)
  ggsave(p,filename = file.path(result_dir,paste0("Imagedimplot/",i,"_anno_spatial.pdf")),width = 12 ,height = 12)
}

xenium$subtype %>% unique()
xenium$subtype <- factor(xenium$subtype)

sample_label <- unique(xenium@meta.data$sample)
for( i in sample_label){
  #j = gsub("_","-",i)
  DefaultBoundary(xenium@images[[i]]) <- "segmentation"
  p <- ImageDimPlot(xenium,
                    fov = i,
                    nmols = 20000,
                    dark.background = F,
                    border.color = "black", 
                    border.size = 0.1,
                    cols = all_colors,
                    #shuffle.cols = F, #随机打乱颜色顺序
                    #size = 1.2,
                    group.by = "subtype"
                    #molecules = c("INHBA"),
                    #mols.size = 1,
                    #mols.cols = "red",
                    #mols.alpha = 0.5, #分子的透明度值
                    
  )
  ggsave(p,filename = file.path(result_dir,paste0("Imagedimplot/",i,"_sub_spatial.png")),width = 12 ,height = 12,dpi = 300)
  ggsave(p,filename = file.path(result_dir,paste0("Imagedimplot/",i,"_sub_spatial.pdf")),width = 12 ,height = 12)
}

saveRDS(xenium, file.path(result_dir,"xenium_TumorROI.rds"))


#7、邻域分析--------
library(Seurat)
library(tidyverse)
source("spatial analysis.R")

CN_out_dir <- file.path(result_dir,"CN analysis/")
dir.create(CN_out_dir)

cells_position_ls <- readRDS("R_OUTPUT/data_for_analysis/cells_position_ls.rds")
cells_position_df <- cells_position_ls %>% bind_rows()

## CN analysis ---------
# options(warn = 0)
# Classification1 must be a factor
# add a column: obj_id

all_data2$subtype %>% unique()

target_subtype <- c('M0_Microglia', 'M1_Microglia', 'M2_Microglia' , 'M0_Macrophage', 'M1_Macrophage', 'M2_Macrophage',"CD4_T", "CD8_T", "DC", "Neutrophil",
                    "Endothelial", "Fibroblast", "Neuron", "Oligodendrocyte", "Pericyte", "Tumor"
) 

all_data_meta_for_cn_analysis <- 
  all_data2[[c("sample", "group", "area", "area_group","annotation","subtype")]] %>% 
  rownames_to_column("cell_id") %>%
  filter(area != c("other")) %>% 
  filter(subtype %in% target_subtype) %>% 
  left_join(cells_position_df[c('cell_id', 'x_centroid', 'y_centroid')], by = "cell_id") %>%
  #rename(X = x_centroid, Y = y_centroid) %>%
  mutate(
    subtype = factor(subtype), # import!
    obj_id = cell_id # import!
  )

all_data_meta_for_cn_analysis$X <- all_data_meta_for_cn_analysis$x_centroid 
all_data_meta_for_cn_analysis$Y <- all_data_meta_for_cn_analysis$y_centroid 
all_data_meta_for_cn_analysis$x_centroid <- NULL
all_data_meta_for_cn_analysis$y_centroid <- NULL

saveRDS(all_data_meta_for_cn_analysis,file.path(CN_out_dir,"all_data_meta_for_cn_analysis.rds"))

### CN k10 ---------
res_cn_10 <- CN_compute_CN(
  all_data_meta_for_cn_analysis,
  x = "X",
  y = "Y",
  group.by       = "sample",
  classification = "subtype",
  k.nn           = 10,
  k.cluster      = 10
)

saveRDS(res_cn_10, file.path(CN_out_dir, "res_cn_10.rds"))

### heatmap ----------------------
# FoldChange cluster/allCell
cn10_cluster_mean <- res_cn_10$cluster_centers
cn10_baseMean     <- res_cn_10$base_mean

# log2(FC)
cn10_cluster_FC <- log2(1 + sweep(cn10_cluster_mean, MARGIN = 2, FUN = "/", STATS = cn10_baseMean))

# heatmap
p_cn10_heat <-
  pheatmap::pheatmap(
    cn10_cluster_FC,
    color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
    border_color = NA,
    treeheight_row = 0,
    labels_row = str_c("CN", rownames(cn10_cluster_FC), sep = ""),
    # clustering_method = 'ward.D2',
    cluster_rows = F,
    scale = "column",
    silent = TRUE
  )

pdf(
  file.path(CN_out_dir, "CN_10_heatmap.pdf"),
  height = 6,
  width  = 10
)
cowplot::plot_grid(p_cn10_heat$gtable)
dev.off()


### CN annotaton --------------------
# # CN k10
# cn_annotation <- c(
#   "CN1"  = "Macrophage",
#   "CN2"  = "T",
#   "CN3"  = "B",
#   "CN4"  = "Fibro_Stroma",
#   "CN5"  = "Myo_Stroma",
#   "CN6"  = "Pan immune",
#   "CN7"  = "NK_Endothelial_Stroma",
#   "CN8"  = "Epithelial",
#   "CN9"  = "T_Epithelial",
#   "CN10" = "T_mast_Epithelial"
# )
# 
# p_cn_heat2 <-
#   cn_cluster_FC %>%
#   pheatmap::pheatmap(
#     color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
#     border_color = NA,
#     treeheight_row = 0,
#     labels_row = cn_annotation,
#     # clustering_method = 'ward.D2',
#     cluster_rows = F,
#     scale = "column",
#     silent = TRUE
#   )
# 
# pdf(
#   file.path(CN_out_dir, "CN_heatmap_annotation.pdf"),
#   height = 6,
#   width  = 11
# )
# cowplot::plot_grid(p_cn_heat2$gtable)
# dev.off()


### CN spatial plot -----------------
cn10_sp_out_dir <- file.path(CN_out_dir, 'CN10 spatial plot/') 
dir.create(cn10_sp_out_dir)

cells_boundary_ls <- readRDS("R_OUTPUT/data_for_analysis/cells_boundary_ls.rds")
# image_wh_r_ls     <- readRDS('data_for_analysis/image_wh_r_ls.rds')

All_type_meta_cn10 <- 
  all_data2[[]] %>%
  rownames_to_column('cell_id') %>%
  # left_join(rownames_to_column(res_cn_10$cn_cluster, "cell_id"), by = "cell_id") %>%
  # mutate(CN2 = fct_recode(CN, !!!structure(names(cn_annotation), names = cn_annotation)))
  left_join(rownames_to_column(res_cn_10$cn_cluster, "cell_id"), by = "cell_id")

saveRDS(All_type_meta_cn10, file.path(CN_out_dir,'all_data_meta_cn10.rds'))

cn10_sp_p_ls <-
  cells_boundary_ls %>% 
  imap(function(cells_boundary, project){ # browser()
    # image_wh_r <- image_wh_r_ls[[project]]
    
    anno_data   <- All_type_meta_cn10[c("sample", "group")] %>% distinct()
    group_info  <- anno_data$group[anno_data$sample == project]
    
    dat <-
      cells_boundary %>%
      right_join(All_type_meta_cn10, by = "cell_id")
    
    p <- dat %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) +
      geom_polygon(aes(fill = CN, group = cell_id))  +
      scale_y_reverse() +
      # coord_fixed(ratio = 1/image_wh_r) +
      coord_equal() + 
      ggsci::scale_fill_igv() +
      theme_void() +
      labs(fill = NULL, title = paste(project, group_info, sep = "_") )
    p
  })


# export image
cn10_sp_p_ls %>% 
  iwalk(function(p, name){
    pdf(
      file.path(cn10_sp_out_dir, str_glue("CN10 spatial plot {name}.pdf")),
      width  = 20,
      height = 20
    )
    print(p)
    dev.off()
  })

###CN加入带images的xenium中-----
xenium$cell_id_raw %>% head()
xenium$cell_id <- paste0(xenium$project,"_",xenium$cell_id_raw)
xenium$cell_id %>% head()

xenium_meta_cn <- xenium[[]] %>%
  rownames_to_column('cell_id_rownames') %>%
  # left_join(rownames_to_column(res_cn_10$cn_cluster, "cell_id"), by = "cell_id") %>%
  # mutate(CN2 = fct_recode(CN, !!!structure(names(cn_annotation), names = cn_annotation)))
  left_join(rownames_to_column(res_cn_10$cn_cluster, "cell_id"), by = "cell_id")

xenium_meta_cn %>% head()

rownames(xenium_meta_cn) <- xenium_meta_cn$cell_id_rownames

xenium@meta.data <- xenium_meta_cn

###Imagedimplot_annotation----
all_colors <- c(
  "#A6CEE3", "#5B9EC9", "#2D82AF", "#7EBA98", "#98D277",
  "#52AF43", "#6F9E4C", "#DD9A88", "#FDBB69", "#FE982C",
  "#F78620", "#D9A295", "#feee7d", "#02c9c9", "#600473",
  "#B294C7", "#7D54A5", "#9E8099", "#F0EB99", "#DBB466",
  "#B15928", # 原始21色
  "#FF6B6B", "#4ECDC4", "#FF9F1C", "#7209B7", "#4361EE", 
  "#F15BB5", "#00BBF9", "#00F5D4", "#9B5DE5", "#00F5A7",
  "#FEE440", "#9C89B8", "#06D6A0"  # 新增13色
)

all_colors10 <- c(
  "#FF6B6B","#A6CEE3","#00F5A7", "#FDBB69","#feee7d", "#9E8099",
  "#00BBF9","#7209B7","#2D82AF", "#F15BB5"
)

xenium$CN %>% unique()
table(is.na(xenium$CN))
xenium_cn <- xenium[, !is.na(xenium$CN)]
table(is.na(xenium_cn$CN))

sample_label <- unique(xenium_cn@meta.data$sample)
for( i in sample_label){
  #j = gsub("_","-",i)
  DefaultBoundary(xenium_cn@images[[i]]) <- "segmentation"
  p <- ImageDimPlot(xenium_cn,
                    fov = i,
                    nmols = 20000,
                    dark.background = F,
                    border.color = "black", 
                    border.size = 0.1,
                    cols = all_colors10,
                    shuffle.cols = F, #随机打乱颜色顺序
                    #size = 1.2,
                    group.by = "CN"
                    #molecules = c("INHBA"),
                    #mols.size = 1,
                    #mols.cols = "red",
                    #mols.alpha = 0.5, #分子的透明度值
                    
  )
  ggsave(p,filename = file.path(cn10_sp_out_dir,paste0(i,"_CN_spatial.png")),width = 12 ,height = 12,dpi = 300)
  ggsave(p,filename = file.path(cn10_sp_out_dir,paste0(i,"_CN_spatial.pdf")),width = 12 ,height = 12)
}





### cell props ------
####WT vs MUT肿瘤区域-------
All_type_cn10_count_df <-
  All_type_meta_cn10 %>% 
  filter(area == c("Tumor")) %>% 
  filter(subtype %in% target_subtype) %>%
  group_by(sample, area_group, CN) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>% 
  mutate(area_group = factor(area_group, levels = c("WT_Tumor","MUT_Tumor")))

p_cn10_props_group <-
  All_type_cn10_count_df %>%
  
  ggplot(aes(x = CN, y = freq, fill = area_group)) +
  geom_boxplot() +
  ggsci::scale_fill_aaas() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) +
  labs(
    x = NULL
  )

pdf(
  file.path(CN_out_dir, "CN10_boxplot_Tumor.pdf"),
  height = 6,
  width  = 11
)
p_cn10_props_group %>% print
dev.off()


#8、邻近度--------
proximity_out_dir <- file.path(result_dir,"proximity analysis")
dir.create(proximity_out_dir)

all_data_meta_for_cn_analysis <- readRDS(file.path(CN_out_dir,"all_data_meta_for_cn_analysis.rds"))
all_data_meta_for_cn_analysis %>% head()

all_data_meta_pos <- all_data_meta_for_cn_analysis %>% 
  mutate(Slide_ID = sample)

## proximity analysis -------
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

res_pa <-
  ProximutyAnalysis_NN_multi(
    all_data_meta_pos,
    query = all_data_meta_pos,
    x = "X",
    y = "Y",
    group.by = "Slide_ID",
    classification = "subtype", # 若为NULL: 提前返回window的细胞近邻的list
    k = 750, # 'auto',
    searchtype = "radius",
    radius = 50, #  50um半径
    mpp = 1, # !important, 1 for xenium
    returntype = "count",
    cell_enough_ratio_cutoff = 90
  )

# return res_pa
# rows: all celltypes
# cols: query cell
res_pa_colsum <- colSums(res_pa)
res_pa <- sweep(res_pa, 2, res_pa_colsum, FUN = "/")
if(any(res_pa_colsum == 0)){
  res_pa[, res_pa_colsum == 0] <- 0
}

# res_pa <- readRDS('data/res_pa.rds')


# summarized by group and celltype
#pa_groups <-
#  with(
#    all_data_meta_pos_subset,
#    str_c(Classification1, Group, Slide_ID, sep = '..')
#  )


###WT vs MUT肿瘤区域-------
pa_groups <-
  all_data_meta_pos %>%
  group_by(Slide_ID) %>% # compute pa_groups per Slide_ID
  nest() %>%
  deframe() %>%
  imap(function(dat, Slide_ID){
    with(
      dat,
      str_c(subtype, area_group, Slide_ID, sep = '..')
    )
  }) %>%
  unlist(use.names = FALSE)


# 分组内取细胞比例的中位数
# apply对每行操作时会转置结果
res_pa_df <-  apply(res_pa, 1, function(x){tapply(x, INDEX = pa_groups, FUN = "median")})

res_pa_df  <-
  res_pa_df %>%
  as.data.frame() %>%
  rownames_to_column('group') %>%
  pivot_longer(-group, names_to = "Others", values_to = "cell_props") %>%
  separate(group, into = c("Target", "Group", "project"), sep = "\\.\\.") 

# plot boxplot
# Target:中心细胞
# Others: 其他细胞占比

p_res_pa <-
  res_pa_df %>%
  dplyr::filter(Group %in% c("WT_Tumor","MUT_Tumor")) %>% 
  # dplyr::filter(Target == "Epithelial cells") %>%
  ggplot(aes(x = Others, y = cell_props, fill = Group)) +
  geom_boxplot(width = 0.35, outlier.colour = NA) +
  geom_jitter(width = 0.25) +
  facet_wrap(~ Target, scales = "free_y")+
  ggpubr::theme_pubclean() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

# export plot
pdf(
  file.path(proximity_out_dir , "area_Tumor proximity analysis.pdf"),
  height = 30,
  width  = 30
)
p_res_pa %>% print()
dev.off()


p_res_pa <-
  res_pa_df %>%
  dplyr::filter(Group %in% c("WT_Tumor","MUT_Tumor")) %>% 
  dplyr::filter(Target == "Tumor") %>%
  ggplot(aes(x = Others, y = cell_props, fill = Group)) +
  geom_boxplot(width = 0.35, outlier.colour = NA) +
  geom_jitter(width = 0.25) +
  facet_wrap(~ Target, scales = "free_y")+
  ggpubr::theme_pubclean() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

# export plot
pdf(
  file.path(proximity_out_dir , "area_Tumor proximity analysis Target_Tumor.pdf"),
  height = 6,
  width  = 12
)
p_res_pa %>% print()
dev.off()



#9、蛋白基因集打分------
score_dir <- file.path(result_dir,"score_result")
dir.create(score_dir)

xenium$ROI %>% table()

genesets <- c(
  "ADGRL3","ANGPT2","APOD","ART3","C3","C4B","CADM4","COL3A1","CRYL1","CSPG4",
  "CTSF","DCT","DDAH2","DPEP2","EGFLAM","ENPP2","FBN1","FGFR3","GSDMD","GSTM2",
  "HTRA1","LGI4","MMP11","MOXD1","NTM","PCSK1N","PDLIM4","PGM5","PLXNB1","PTN",
  "RTN4RL2","S100B","SCRG1","SELENOP","SEMA4G","SMPDL3A","TGFB2"
)

valid_gene <- genesets[genesets %in% rownames(xenium)]
cat("有效基因:", length(valid_gene), "/", length(genesets))

missing_gene <- setdiff(genesets, valid_gene)

library(Seurat)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(ggsignif)
library(ggpubr)

xenium <- AddModuleScore(
  object = xenium,
  features = list(valid_gene),
  name = "portein_Score"
)

xenium$portein_Score <- xenium$portein_Score1

xenium$group %>% unique()
xenium$group <- factor(xenium$group,levels = c("WT","MUT"))
##小提琴-----
V <- VlnPlot(
  xenium,
  features = "portein_Score",
  group.by = "group",
  pt.size = 0,
  fill.by = "group"
) +
  geom_boxplot(width=0.1, fill="white") +
  theme_bw() +
  ggtitle("") +
  stat_compare_means(
    method = "wilcox.test",
    label = "p.signif",
    ref.group = "WT", # 以Con为对照组，自动和所有其他组比较
    hide.ns = F,
    size = 3,
    vjust = 1.1,
    hjust = 2
  )

ggsave(file.path(score_dir,"Score_Violin.pdf"),V, width=7, height=6)
ggsave(file.path(score_dir,"Score_Violin.png"),V, width=7, height=6)

xenium$area_group %>% unique()
xenium$area_group <- factor(xenium$area_group,levels = c("WT_Tumor","MUT_Tumor"))
##小提琴-----
V <- VlnPlot(
  xenium,
  features = "portein_Score",
  group.by = "area_group",
  pt.size = 0,
  fill.by = "area_group"
) +
  geom_boxplot(width=0.1, fill="white") +
  theme_bw() +
  ggtitle("") +
  stat_compare_means(
    method = "wilcox.test",
    label = "p.signif",
    ref.group = "WT_Tumor", # 以Con为对照组，自动和所有其他组比较
    hide.ns = F,
    size = 3,
    vjust = 1.1,
    hjust = 2
  )

ggsave(file.path(score_dir,"Score_Violin areagroup.pdf"),V, width=7, height=6)
ggsave(file.path(score_dir,"Score_Violin areagroup.png"),V, width=7, height=6)

##评分空间frature-----
Imagefeature_dir <- file.path(score_dir,"Imagefeature")
dir.create(Imagefeature_dir)

sample_label <- unique(xenium@meta.data$sample)
for (sample_id in sample_label) {
  # 智能匹配样本ID
  j <- names(xenium@images)[grepl(sample_id, names(xenium@images))]

      DefaultBoundary(xenium@images[[j]]) <- "segmentation"
      
      p <- ImageFeaturePlot(
        object = xenium,
        fov = j,
        nmols = 20000,       # 每个图像中显示的最大分子数（随机子抽样）
        border.color = "black",  # 细胞边界的颜色（深蓝色）
        border.size = 0.1,    # 细胞边界的线宽
        features = "portein_Score",
        max.cutoff = "q90",
        axes = FALSE,
        cols = c('#003366', "#E7475E"),
        size = 0.8,
        dark.background = F
      )
      
      ggsave(
        filename = file.path(Imagefeature_dir, paste0(sample_id, "_" , "spatial.png")),
        plot = p,
        width = 12,
        height = 12,
        dpi = 300
      )
      
      ggsave(
        filename = file.path(Imagefeature_dir, paste0(sample_id, "_" , "spatial.pdf")),
        plot = p,
        width = 12,
        height = 12
      )
      
      message("成功保存: ", sample_id)
}

#10、通路打分-----
library(clusterProfiler)
library(org.Mm.eg.db)
library("org.Hs.eg.db")
library(Seurat)

go_target_human <- list(
  # CC 细胞组分
  CC_extracellular_region = "GO:0005576",
  CC_extracellular_space = "GO:0005615",
  CC_plasma_membrane = "GO:0005886",
  CC_cell_periphery = "GO:0071944",
  # MF 分子功能
  MF_metal_ion_binding = "GO:0046872",
  MF_ion_binding = "GO:0043167",
  MF_transition_metal_binding = "GO:0046914",
  MF_cation_binding = "GO:0043169",
  MF_phosphodiester_hydrolase = "GO:0008081",
  MF_zinc_ion_binding = "GO:0008270",
  MF_hydrolase_activity = "GO:0016787",
  MF_phosphoester_hydrolase = "GO:0042578",
  # BP 生物学过程【客户重点：炎症负调控】
  BP_NTP_catabolic = "GO:0009143",
  BP_neg_reg_inflammation = "GO:0050728"
)


go_gene_sets <- list()
for (path_name in names(go_target_human)) {
  go_id <- go_target_human[[path_name]]
  # 提取GO对应的人类ENTREZID，转换为Symbol
  map_df <- select(
    org.Hs.eg.db,
    keys = go_id,
    keytype = "GOALL",
    columns = c("SYMBOL")
  )
  gene_vec <- unique(na.omit(map_df$SYMBOL))
  go_gene_sets[[path_name]] <- gene_vec
  cat(path_name, "通路基因总数：", length(gene_vec), "\n")
}

# 过滤仅保留Xenium面板存在的基因
go_filtered <- lapply(go_gene_sets, function(gene_list) {
  intersect(gene_list, rownames(xenium))
})
# 删除基因过少无法计算分数的通路
go_filtered <- go_filtered[sapply(go_filtered, length) >= 2]

xenium <- AddModuleScore(
  object = xenium,
  features = go_filtered,
  name = "PathScore_",
  ctrl = 100 # 随机选取100个背景基因校正表达偏差
)

score_cols <- grep("^PathScore_", colnames(xenium@meta.data), value = TRUE)

modules <- names(go_filtered) |> unique()
modules

# 定义新旧列名映射关系
new_names <- modules
old_names <- score_cols #paste0("PathScore_", 1:14)

# 直接修改metadata列名
colnames(xenium@meta.data)[match(old_names, colnames(xenium@meta.data))] <- new_names

pathway_dir <- file.path(result_dir,"pathway_result")
dir.create(pathway_dir)

##小提琴-----
for (i in modules) {
  V <- VlnPlot(
    xenium,
    features = i,
    group.by = "group",
    pt.size = 0,
    fill.by = "group"
  ) +
    geom_boxplot(width = 0.1, fill = "white") +
    theme_bw() +
    ggtitle("") +
    stat_compare_means(
      method = "wilcox.test",
      label = "p.signif",
      ref.group = "WT", # 以Con为对照组，自动和所有其他组比较
      hide.ns = F,
      size = 3,
      vjust = 1.1,
      hjust = 2
    ) #+
  # facet_grid(~modules,ncol = 3)

  ggsave(file.path(pathway_dir, paste0(i, "_Score_Violin.pdf")), V, width = 7, height = 6)
  ggsave(file.path(pathway_dir, paste0(i, "_Score_Violin.png")), V, width = 7, height = 6)
}


xenium$area_group %>% unique()
xenium$area_group <- factor(xenium$area_group,levels = c("WT_Tumor","MUT_Tumor"))



## 小提琴-----
for (i in modules) {
  V <- VlnPlot(
    xenium,
    features = i,
    group.by = "area_group",
    pt.size = 0,
    fill.by = "area_group"
  ) +
    geom_boxplot(width = 0.1, fill = "white") +
    theme_bw() +
    ggtitle("") +
    stat_compare_means(
      method = "wilcox.test",
      label = "p.signif",
      ref.group = "WT_Tumor", # 以Con为对照组，自动和所有其他组比较
      hide.ns = F,
      size = 3,
      vjust = 1.1,
      hjust = 2
    )

  ggsave(file.path(pathway_dir, paste0(i, "Score_Violin areagroup.pdf")), V, width = 7, height = 6)
  ggsave(file.path(pathway_dir, paste0(i, "Score_Violin areagroup.png")), V, width = 7, height = 6)

}



##通路空间featureplt----
Imagefeature_dir <- file.path(pathway_dir,"Imagefeature")
dir.create(Imagefeature_dir)

module_names <- new_names

for (sample_id in unique(xenium@meta.data$sample)) {
  # 智能匹配样本ID
  j <- names(xenium@images)[grepl(sample_id, names(xenium@images))]
  
  if (length(j) == 0) {
    warning("跳过未匹配样本: ", sample_id)
    next
  }
  
  for (module in module_names) {
    tryCatch({
      # 空值检查
      if (!is.null(xenium@images[[j]])) {
        DefaultBoundary(xenium@images[[j]]) <- "segmentation"
      }
      
      p <- ImageFeaturePlot(
        object = xenium,
        fov = j,
        nmols = 20000,       # 每个图像中显示的最大分子数（随机子抽样）
        border.color = "black",  # 细胞边界的颜色（深蓝色）
        border.size = 0.1,    # 细胞边界的线宽
        features = module,
        max.cutoff = "q90",
        axes = FALSE,
        cols = c('#003366', "#E7475E"),
        size = 0.8,
        dark.background = F
      )+
        ggtitle(paste(module)) +
        theme(
          plot.title = element_text(size = 14, face = "bold"),
          legend.position = "right",
          legend.title = element_text(size = 10)
        )
      
      ggsave(
        filename = file.path(Imagefeature_dir, paste0(sample_id, "_", module,"_spatial.png")),
        plot = p,
        width = 12,
        height = 12,
        dpi = 300
      )
      ggsave(
        filename = file.path(Imagefeature_dir, paste0(sample_id, "_", module,"_spatial.pdf")),
        plot = p,
        width = 12,
        height = 12
      )
      
      message("成功保存: ", sample_id, "_", module)
      
    }, error = function(e) {
      warning("失败[", sample_id, ":", "Module_",module,"]: ", e$message)
    })
  }
}
