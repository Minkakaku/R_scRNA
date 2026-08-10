cds <- readRDS(file.path(monocle_out_dir, "cds.rds"))
out_dir <- "Figure_2"
dir.create(out_dir)
states_color <- c(
  # ========== 左侧：未分化/早期（冷灰-蓝-紫，低饱和）==========
  "1"  = "#5B7C99",   # 灰蓝
  "2"  = "#7A8B69",   # 灰绿（与cell_colors的翠绿区分）
  "3"  = "#8B7D9E",   # 灰紫
  "4"  = "#6B8E9F",   # 青灰
  
  # ========== 中心：起点（中性浅色）==========
  "5"  = "#D4C5B5",   # 暖灰/米白，分化原点
  
  # ========== 右侧：分化进程（暖色-大地色，与cell_colors区分）==========
  "6"  = "#A67B5B",   # 浅褐
  "7"  = "#8B6F47",   # 驼色
  "8"  = "#B85C4E",   # 暗红（与cell_colors的深橙红区分）
  "9"  = "#C9A227",   # 暗金（与cell_colors的亮黄/金黄区分）
  "10" = "#6B5B73",   # 灰紫褐（与cell_colors的紫红/紫粉区分）
  "11" = "#4A6741"    # 深灰绿（与cell_colors的翠绿/青绿区分）
)
# some ploting
# cds <- orderCells(cds, root_state = 4)
plot_cell_trajectory(cds,color_by="Pseudotime", size=1,show_backbone=TRUE) 
plot_cell_trajectory(cds,color_by="annotation", size=1,show_backbone=TRUE) + scale_color_manual(values = col_plan_1)
plot_cell_trajectory(cds, color_by = "State",size=1,show_backbone=TRUE) + scale_color_manual(values = col_plan_2)
plot_cell_trajectory(cds, color_by = "State") + facet_wrap("~State", nrow = 1)
plot_complex_cell_trajectory(cds, x = 1, y = 2, color_by = "annotation")

# p1 <- plot_genes_jitter(cds[s.genes,], grouping = "State", color_by = "State")
# p2 <- plot_genes_violin(cds[s.genes,], grouping = "State", color_by = "State")
# p3 <- plot_genes_in_pseudotime(cds[s.genes,], color_by = "State")
# 
# p1=plot_cell_trajectory(cds, color_by = "CCL5")  + scale_color_gsea()

pdf(file.path(out_dir, "cell trajectory.pdf"), width = 8, height = 4)

plot_cell_trajectory(cds,color_by="Pseudotime", size=1, show_backbone=TRUE) + theme(legend.position = "right")

plot_cell_trajectory(cds,color_by="group", size=1,show_backbone=TRUE) +
  scale_color_manual(values = c("#D55E00","#0072B2")) +
  theme(legend.position = "right")

plot_cell_trajectory(cds,color_by="annotation", size=1,show_backbone=TRUE) + 
  scale_color_manual(values = cell_colors) +
  theme(legend.position = "right")

plot_cell_trajectory(cds, color_by = "State",cell_size =1,show_backbone=TRUE) +
  scale_color_manual(values = states_color) +
  theme(legend.position = "right")

dev.off()


pdf(file.path(out_dir, "cell trajectory - annotation - split by group.pdf"), width = 12, height = 4)
local({
  pData(cds)$group <- pData(cds)$group %>% factor(levels = c("C", "T"))
  plot_cell_trajectory(cds, color_by = "annotation",size=1,show_backbone=TRUE) +
    scale_color_manual(values = cell_colors) +
    facet_wrap(~ group) +
    theme(legend.position = "right")
})
dev.off()



pdf(file.path(out_dir, "cell trajectory - tree.pdf"), width = 28, height = 10)
plot_complex_cell_trajectory(cds, x = 1, y = 2, color_by = "annotation") +
  scale_color_manual(values = cell_colors) +
  facet_wrap(~ annotation,nrow = 1)
dev.off()

# cell barplot 
p_col_cell_state <-
  cds %>% 
  pData() %>%
  dplyr::select(annotation, State) %>%
  group_by(annotation, State) %>%
  mutate(annotation = as.character(annotation)) %>%
  mutate(annotation = factor(annotation, levels = c(
    "Syncytiotrophoblast Precursor", 
    "Villous Cytotrophoblast", 
    "Syncytiotrophoblast",   
    "Column Cytotrophoblast", 
    "Extravillous Trophoblast"
  ))) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = annotation, y = freq)) +
  geom_col(aes(fill = State)) +
  scale_fill_manual(values = states_color) +
  scale_y_continuous(labels = scales::percent) +
  coord_flip() +
  labs(x = NULL, y = NULL) +
  ggpubr::theme_pubr() +
  theme(
    legend.position = "right"
  )
ggsave(p_col_cell_state, filename = file.path(out_dir, "cell stat barplot.pdf"), height = 4, width = 8)

# state barplot
p_col_state_cell <-
  cds %>% 
  pData() %>%
  dplyr::select(annotation, State) %>%
  # 核心修改1：按 State 分组统计
  group_by(State, annotation) %>%
  mutate(annotation = as.character(annotation)) %>%
  mutate(annotation = factor(annotation, levels = c(
    "Syncytiotrophoblast Precursor", 
    "Villous Cytotrophoblast", 
    "Syncytiotrophoblast",   
    "Column Cytotrophoblast", 
    "Extravillous Trophoblast"
  ))) %>%
  summarise(count = n(), .groups = "drop") %>%
  # 核心修改2：按 State 计算比例（每个State内部的细胞类型占比）
  group_by(State) %>%
  mutate(freq = count/sum(count)) %>%
  
  # 核心修改3：x轴 = State，fill = 细胞类型annotation
  ggplot(aes(x = State, y = freq)) +
  geom_col(aes(fill = annotation)) +
  # 颜色映射改成细胞类型（如果你有对应配色直接替换）
  scale_fill_manual(values = cell_colors) + # 若无配色可删除此行
  scale_y_continuous(labels = scales::percent) +
  coord_flip() +
  labs(x = NULL, y = NULL, fill = "Cell Type") +
  ggpubr::theme_pubr() +
  theme(
    legend.position = "right"
  )
print(p_col_state_cell)

ggsave(p_col_state_cell, filename = file.path(out_dir, "stat cell barplot.pdf"), height = 4, width = 8)

#
p_col_state_group <-
  cds %>% 
  pData() %>%
  # 这里改成选择 State + 你的group列名
  dplyr::select(State, group) %>%  # 把 group 换成你真实的列名！
  
  # 按 State + group 分组计数
  group_by(State, group) %>%
  summarise(count = n(), .groups = "drop") %>%
  
  # 每个State内部算比例（核心）
  group_by(State) %>%
  mutate(freq = count / sum(count)) %>%
  
  # 绘图：x=State，填充=group
  ggplot(aes(x = State, y = freq)) +
  geom_col(aes(fill = group)) +
  scale_fill_manual(values = c("#D55E00","#0072B2")) + # 两个颜色
  scale_y_continuous(labels = scales::percent) +
  coord_flip() +
  labs(x = NULL, y = NULL, fill = "Group") +
  ggpubr::theme_pubr() +
  theme(legend.position = "right")
print(p_col_state_group)

ggsave(p_col_state_group, filename = file.path(out_dir, "stat group barplot.pdf"), height = 4, width = 6)


#
# de analysis
#这里是把排序基因（ordergene）提取出来做回归分析，来找它们是否跟拟时间有显著的关系
#如果不设置，就会用所有基因来做它们与拟时间的相关性


#####
BEAM_res <- read.csv(file.path(monocle_out_dir, "BEAM_res.csv"))
#选前100个基因可视化
BEAM_genes <- BEAM_res %>%
  top_n(n = 100, wt = -qval) %>%   # -qval 表示取最小的
  pull(gene_short_name) %>%
  as.character()
p <- plot_genes_branched_heatmap(
  cds[BEAM_genes,],  
  branch_labels = c("Extravillous Trophoblast", "Syncytiotrophoblast"),
  branch_point = 5, 
  num_clusters = 4, 
  show_rownames = T, 
  return_heatmap = T
)
ggsave(file.path(out_dir, "BEAM_heatmap_1.pdf"), p$ph_res, width = 6.5, height = 10)


######

#####alltype循环
all_cell_types <- names(cell_colors)

# 遍历每种细胞类型
for(target_cell_type in all_cell_types){
  
  message("正在绘制：", target_cell_type)
  
  # 强制锁定当前细胞类型和颜色（解决循环+imap冲突！）
  current_cell <- target_cell_type
  current_color <- cell_colors[[current_cell]]
  
  p_sp_ct_ls <- cells_position_ls %>%
    imap(function(coords, sam){
      
      # 强制内部使用锁定的变量，不受循环影响
      ct <- current_cell
      #cc <- current_color
      
      dat2 <- coords %>% left_join(plot_dat, by = "cell_id")
      
      # 高亮分组
      dat2$highlight <- ifelse(dat2[[annotation_col]] == ct, ct, "Other")
      
      dat2$highlight <- factor(
        dat2$highlight,
        levels = c(ct, "Other")
      )
      
      # 绘图
      p <- ggplot(dat2, aes(x = SPATIAL_1, y = SPATIAL_2)) +
        geom_tile(aes(fill = highlight,alpha = highlight), height = 10, width = 10) +
        theme_void() +
        coord_equal() +
        scale_y_reverse() +
        ##透明度
        scale_alpha_manual(
          values = setNames(c(1, 0.3), c(ct, "Other"))  # 目标不透明，其他15%透明
        ) +
        # 核心：强制命名颜色（解决ggplot不认变量问题）
        scale_fill_manual(
          values = setNames(c(current_color, "gray60"), c(ct, "Other")),
          guide = guide_legend(override.aes = list(alpha = c(1, 1)))
        ) +
        
        theme(
          panel.background = element_rect(fill = NA, color = NA),
          plot.title = element_text(hjust = 0.5)
        ) +
        ggtitle(sam)
      
      return(p)
    })
  
  # 输出PDF
  pdf(
    file.path(out_dir, paste0("spatial plot of ", target_cell_type, ".pdf")),
    height = 8, width = 8
  )
  walk(p_sp_ct_ls, print)
  dev.off()
}

message("全部完成！已输出每个细胞类型独立PDF！")

