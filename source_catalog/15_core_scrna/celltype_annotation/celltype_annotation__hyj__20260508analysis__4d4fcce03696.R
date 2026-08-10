# 1. TAM亚群细分与功能分型-----
#subtype.R中

# 2. TAM枢纽化的空间验证------
all_data$annotation %>% unique()

Macrophage <- readRDS("result/Macrophage_subtype/Macrophage.rds")


##映射——add亚类注释到all_data----
library(tibble)
library(data.table)
subtype <- rbind(#All_T[['subtype']] %>% rownames_to_column('cells'), 
                 Macrophage[['subtype']]%>% rownames_to_column('cells'))

meta.data.t <-
  subtype%>%
  right_join(
    all_data[['annotation']] %>% rownames_to_column('cells'),
    by = "cells"
  ) %>% 
  mutate(
    subtype = ifelse((as.character(annotation) %in% c(#'T cells',
                                                      "TAM")),
                     as.character(subtype), 
                     as.character(annotation)))%>%
  dplyr::select(cells, subtype) %>%
  column_to_rownames('cells')


all_data<-AddMetaData(all_data,meta.data.t['subtype'])
all_data$annotation %>% table()
all_data$subtype %>% table()

table(all_data$subtype, all_data$annotation, useNA = "ifany")

saveRDS(all_data, 'result/all_data_subtype.rds')

##空间meta构建--------
cells_boundary_ls <- readRDS(file.path("~/R_OUTPUT/data_for_analysis", "cells_boundary_ls.rds"))
cells_position_ls <- readRDS(file.path("~/R_OUTPUT/data_for_analysis", "cells_position_ls.rds"))

library(Seurat)
library(tidyverse)

source("spatial analysis.R")

all_data@meta.data %>% head()
all_data@meta.data$cell_id %>% head()
rownames(all_data@meta.data) %>% head()
#all_data@meta.data$cell_id <- paste0("TMA_",all_data@meta.data$cell_id)

table(all_data$group,all_data$sample)

# saveRDS(all_data,"RCTD/result_noEPI/all_data.rds")

###sample名称和slide名称需一致，slide名称是cells_position_ls各FOV的名字

cells_position_ls %>% head() 

# all_data_cn_meta
all_data_meta_pos <-
  cells_position_ls %>%
  imap(function(cells_position, slide){ # browser()
    
    all_data_meta <- all_data[[]]  
    
    proj_meta_sub <-
      all_data_meta %>%  
      dplyr::select(subtype, annotation, sample, "group") %>% 
      rownames_to_column("cell_id") %>%
      dplyr::filter(sample == slide)
    
    cells_position[, 1:3] %>%
      inner_join(proj_meta_sub, by = "cell_id") %>%  
      mutate(subtype = factor(subtype)) %>%
      dplyr::rename(subtype = subtype,
                    X = x_centroid,
                    Y = y_centroid,
                    Slide_ID = sample,
                    Group = group
      )
  }) %>%
  enframe() %>%
  unnest(cols = value) %>%
  dplyr::select(-name) %>%
  column_to_rownames("cell_id")

all_data_meta_pos %>% head()
all_data_meta_pos$Slide_ID %>% table()
all_data_meta_pos$Group %>% table()
all_data_meta_pos$subtype %>% table()

all_data_meta_pos$cell_id <- rownames(all_data_meta_pos)
all_data_meta_pos <- all_data_meta_pos %>%
  as.data.frame() %>%       
  remove_rownames()    

write.csv(all_data_meta_pos,file.path("result/","all_data_meta_pos.csv"))
table(all_data_meta_pos$Slide_ID,all_data_meta_pos$Group)


##最近邻分析 -------
nearest_neighbor_out_dir <- "result/2、nearest_neighbor_out_dir"
dir.create(nearest_neighbor_out_dir)

# 由于最近邻只需要最近的那个细胞的信息
# 因此不同于proximity analysis，有以下变化：
# classification调为NULL
# searchtype 邻近度使用radius
# k不需要太大，10即可，最近邻分析只需要第二个细胞的信息
# radius适当调大100um，防止部分细胞比较离散，距离其他细胞较远
# returntype与cell_enough_ratio_cutoff也没有意义了

# 参数介绍
# x,y 指定细胞的坐标
# group.by 一般是切片或者样本编号，代表分别计算最近邻关系
# classification NULL即可，返回每一个query细胞的细胞编号、距离、query和total细胞的id
# k 计算的最近邻的细胞数
# radius是只查看当前细胞的多少um内的细胞
# mpp 坐标与um的转换比例，xenium是1，celldive是0.325
res_nn_ls <-
  ProximutyAnalysis_NN_multi(
    all_data_meta_pos,
    query = all_data_meta_pos,
    x = "X",
    y = "Y",
    group.by = "Slide_ID",
    classification = NULL, # 若为NULL: 提前返回window的细胞近邻的list
    k = 10, # 'auto',
    searchtype = "radius",
    radius = 100, 
    mpp = 1, # !important, 1 for xenium
    cell_enough_ratio_cutoff = 20
  )

# res_nn_ls <- readRDS("data/res_nn_ls.rds")

res_nn_df <-
  res_nn_ls %>%
  map(function(x){ # browser()
    cell_idx  <- x$nn.idx[, 2]
    cell_idx[cell_idx == 0] <- NA # some cells have no nn cells
    cell_id   <- x$total_cell_id[cell_idx]
    cell_dists <- x$nn.dists[, 2]
    
    all_data_used <- all_data_meta_pos %>% column_to_rownames("cell_id")
    cell_anno <- all_data_used[cell_id, 'annotation']
    
    res <-
      data.frame(
        Others = cell_anno,
        distance = cell_dists
      ) %>%
      cbind(all_data_used[x$query_cell_id, c('annotation', 'Group', 'Slide_ID')]) %>%
      rownames_to_column("cell_id") %>%
      dplyr::rename(Target = annotation) 
    
    res
  }) %>%
  bind_rows()

if(any(is.na(res_nn_df$Others) )){
  message(str_glue("{sum(is.na(res_nn_df$Others))} cells have no nearest neighbor cell"))
  res_nn_df <- drop_na(res_nn_df)
}

###量化6种非肿瘤细胞向TAM邻近的数据（MUT vs WT）-----
res_nn_df_filter <- res_nn_df %>% 
  filter(Target == "TAM"& !Others %in% c("Tumor","TAM"))

# plot
p_nn_ls <-
  res_nn_df_filter %>%
  group_by(Target) %>%
  nest() %>%
  deframe() %>%
  imap(function(x, target){ # browser()
    x %>%
      ggplot(aes(x = distance, color = Group, fill = Group)) + 
      geom_density(alpha = 0.5) +
      facet_wrap(~ Others) +
      ggpubr::theme_pubclean() +
      labs(
        title = str_glue("Target cell: {target}"),
        x = "distance(um)"
      )
    
  })

# export plot
pdf(
  file.path(nearest_neighbor_out_dir, "nearest neighbor analysis.pdf"),
  height = 12,
  width  = 18
)
p_nn_ls %>% print()
dev.off()

res_nn_df$Target %>% table()
table(res_nn_df_filter$Target,res_nn_df_filter$Others)

### plot 箱线图-----
p_nn_ls <-
  res_nn_df_filter %>%
  group_by(Target) %>%
  nest() %>%
  deframe() %>%
  imap(function(x, target){  
    x %>%
      group_by(Slide_ID, Group, Others) %>%
      summarise(count = n(), .groups = "drop") %>%  
      mutate(freq = count / sum(count)) %>%
      
      ggplot(aes(x = Group, y = freq, fill = Group)) + 
      geom_boxplot(width = 0.6, outlier.size = 1) +  # 箱线图
      # stat_compare_means(
      #   comparisons = list(c("WT","MUT")),  
      #   method = "t.test",  
      #   label = "p.signif",  # p.signif, "p.value"
      #   p.adjust.method = "bonferroni"  
      # ) +
      
      facet_wrap(~ Others, 
                 scales = "free_y"
      ) +
      ggpubr::theme_pubclean() +
      labs(
        title = str_glue("Target cell: {target}"),
        y = "freq"
      )
  })

# export plot
pdf(
  file.path(nearest_neighbor_out_dir, "nearest_neighbor_analysis_boxplot_free_y.pdf"),
  height = 12,
  width  = 18
)
print(p_nn_ls)
dev.off()


##邻近度分析--------
proximity_out_dir <- "result/2、proximity_out_dir"
dir.create(proximity_out_dir)

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
#### 2.1 proximity analysis --------
res_pa <-
  ProximutyAnalysis_NN_multi(
    all_data_meta_pos,
    query = all_data_meta_pos,
    x = "X",
    y = "Y",
    group.by = "Slide_ID",
    classification = "annotation", # 若为NULL: 提前返回window的细胞近邻的list
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


# summarized by Group and celltype
# pa_groups <-
#   with(
#     all_data_meta_pos,
#     str_c(Classification1, Group, Slide_ID, sep = '..')
#   )


pa_groups <-
  all_data_meta_pos %>%
  group_by(Slide_ID) %>% # compute pa_groups per Slide_ID
  nest() %>%
  deframe() %>%
  imap(function(dat, Slide_ID){
    with(
      dat,
      str_c(annotation, Group, Slide_ID, sep = '..')
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

# p_res_pa <-
#   res_pa_df %>%
#   # dplyr::filter(Target == "Epithelial cells") %>%
#   ggplot(aes(x = Others, y = cell_props, fill = Group)) +
#   geom_boxplot(width = 0.35, outlier.colour = NA) +
#   geom_jitter(width = 0.25) +
#   facet_wrap(~ Target, scales = "free_y",ncol = 4)+
#   ggpubr::theme_pubclean() +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
# 
# # export plot
# pdf(
#   file.path(proximity_out_dir, "proximity analysis.pdf"),
#   height = 10,
#   width  = 20
# )
# p_res_pa %>% print()
# dev.off()

table(res_pa_df$Target,res_pa_df$Group)


all_colors <- c(
  "#A6CEE3", "#5B9EC9", "#2D82AF", "#7EBA98", "#98D277",
  "#52AF43", "#6F9E4C", "#DD9A88", "#FDBB69", "#FE982C",
  "#F78620", "#D9A295", "#feee7d", "#02c9c9", "#600473",
  "#B294C7", "#7D54A5", "#9E8099", "#F0EB99", "#DBB466",
  "#B15928"
)


p_res_pa <-
  res_pa_df %>%
  filter(Target == "TAM"& !Others %in% c("Tumor","TAM")) %>% 
  
  ggplot(aes(x = Others, y = cell_props, fill = Group
  )) +
  geom_boxplot(width = 0.35, outlier.colour = NA) +
  geom_jitter(width = 0.25) +
  facet_wrap(~ Target, scales = "free_y",ncol = 3)+
  ggpubr::theme_pubclean() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

# export plot
pdf(
  file.path(proximity_out_dir, paste0("TAM","_proximity analysis.pdf")),
  height = 8,
  width  = 8
)
p_res_pa %>% print()
dev.off()


# 3. TAM与基质细胞的互作网络-----
# 重点分析：Fibroblast、Pericyte、Endothelial
# 分析这三种基质细胞向TAM邻近的通讯特征
all_data$annotation %>% table()

cellchat_data <- all_data
cellchat_data <- JoinLayers(cellchat_data)

# library(future)
# plan("multisession", workers = 10, gc = TRUE)
# options(future.globals.maxSize = 50 * 1024^3) 
# cellchat_data <- NormalizeData(cellchat_data)

cellchat_4types_data <- cellchat_data[ ,cellchat_data$annotation %in% c("Endothelial","Fibroblast","Pericyte","TAM")]

cellchat_4types_data$Celltype <- cellchat_4types_data$annotation
cellchat_4types_data$Celltype <- as.character(cellchat_4types_data$Celltype)
cellchat_4types_data$Celltype %>% unique()
# cellchat_4types_data$Celltype <- gsub(" ","_",cellchat_4types_data$Celltype)
# cellchat_4types_data$Celltype %>% unique()

cellchat_4types_data$group%>% unique()
cellchat_4types_data$group <- as.character(cellchat_4types_data$group)

saveRDS(cellchat_4types_data, file = "/mnt/NexData_file/home/hyj/shell/FWSP20255039_fandandan/Cellchat_4types/cellchat_4types_data.rds")

## 投bash------
#/mnt/NexData_file/home/hyj/shell/FWSP20255039_fandandan/Cellchat_4types




# 4. TAM空间定位与血管生成关联-----


# 5. CellChat TAM枢纽通讯推断-----


# 6. TAM基因表达特征与通路富集-----

