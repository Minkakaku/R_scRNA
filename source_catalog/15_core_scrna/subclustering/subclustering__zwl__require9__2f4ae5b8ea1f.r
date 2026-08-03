setwd('/mnt/NexData_file/home/zwl/Project/SC2511076348/require9')
library(tidyr)
library(dplyr)
library(rstatix)
library(ggplot2)
library(ggpubr)
rm(list = ls())
gc()

################################################################################
######DM组和N组BMSC中KEGG
################################################################################
pathways <- c(
  "Wnt signaling pathway",
  "ECM-receptor interaction",
  "Focal adhesion",
  "PI3K-Akt signaling pathway",
  "TGF-beta signaling pathway",
  "MAPK signaling pathway",
  "Regulation of actin cytoskeleton",
  "TNF signaling pathway",
  "NF-kappa B signaling pathway",
  "Calcium signaling pathway",
  "HIF-1 signaling pathway",
  "Glycolysis / Gluconeogenesis",
  "Glutathione metabolism",
  "AGE-RAGE signaling pathway in diabetic complications",
  "Toll-like receptor signaling pathway"
)

up_data = read.table(file = '/mnt/NexData_file/home/zwl/Project/SC2511076348/require7/Fibro_py/BMSC/DM_vs_N/2.enrich.KEGG_up.csv',
                     sep = ',',header = TRUE)%>% dplyr::mutate(change = 'up')
down_data = read.table(file = '/mnt/NexData_file/home/zwl/Project/SC2511076348/require7/Fibro_py/BMSC/DM_vs_N/2.enrich.KEGG_down.csv',
                       sep = ',',header = TRUE)%>% dplyr::mutate(change = 'down')
merge_data = rbind(up_data,down_data)
filter_data1 = merge_data %>% 
  dplyr::filter(Description %in% pathways) %>%
  mutate(LogP = case_when(
    change == "up" ~ -log10(p.adjust),   # up 取负
    change == "down" ~ log10(p.adjust)   # down 不取负
  ))

filter_data <- filter_data1 %>%
  group_by(Description) %>%
  slice_min(p.adjust, n = 1) %>%  # 取最小p.adj
  ungroup()

# 3. 排序：LogP 从小到大 → 从上到下
filter_data <- filter_data %>% arrange(LogP)
filter_data$Description <- factor(filter_data$Description, levels = unique(filter_data$Description))
write.csv(file = '4.BMCS_care_pathway.csv',filter_data)

# 1. 绘图主体
p <- ggplot(filter_data, aes(x = LogP, y = Description, fill = change)) +
  geom_col() +
  scale_fill_manual(values = c("up"="#F08080", "down"="#87CEEB")) +
  theme_bw() +
  labs(x = "-log10(Adjusted P-value)")

# 2. 自定义主题
mytheme <- theme(
  legend.position = 'none',
  axis.text.y = element_blank(),
  axis.ticks.y = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.border = element_blank(),
  axis.line.x = element_line(color = 'grey60', size = 1.1),
  axis.text = element_text(size = 12)
)

p1 <- p + mytheme

# 3. 按正负拆分数据（用 filter_data！不是dt！）
up_data   <- filter_data[filter_data$LogP > 0, ]
down_data <- filter_data[filter_data$LogP < 0, ]

p2 <- p1 +
  geom_text(data = up_data,
            aes(x = -0.2, y = Description, label = Description),
            size = 3.5, hjust = 1) +
  geom_text(data = down_data,
            aes(x = 0.2, y = Description, label = Description),
            size = 3.5, hjust = 0)#右对齐

p2
ggsave(filename = paste0('4.BMCS_care_pathway.png'),w=10,h=6,p2)
ggsave(filename = paste0('4.BMCS_care_pathway.pdf'),w=10,h=6,p2)

################################################################################
######DM组和N组成骨细胞中KEGG
################################################################################
up_data = read.table(file = '/mnt/NexData_file/home/zwl/Project/SC2511076348/require7/Fibro_py/Osteoblast/DM_vs_N/2.enrich.KEGG_up.csv',
                     sep = ',',header = TRUE)%>% dplyr::mutate(change = 'up')
down_data = read.table(file = '/mnt/NexData_file/home/zwl/Project/SC2511076348/require7/Fibro_py/Osteoblast/DM_vs_N/2.enrich.KEGG_down.csv',
                       sep = ',',header = TRUE)%>% dplyr::mutate(change = 'down')
merge_data = rbind(up_data,down_data)
filter_data1 = merge_data %>% 
  dplyr::filter(Description %in% pathways) %>%
  mutate(LogP = case_when(
    change == "up" ~ -log10(p.adjust),   # up 取负
    change == "down" ~ log10(p.adjust)   # down 不取负
  ))

filter_data <- filter_data1 %>%
  group_by(Description) %>%
  slice_min(p.adjust, n = 1) %>%  # 取最小p.adj
  ungroup()

# 3. 排序：LogP 从小到大 → 从上到下
filter_data <- filter_data %>% arrange(LogP)
filter_data$Description <- factor(filter_data$Description, levels = unique(filter_data$Description))
write.csv(file = '4.Osteoblast_care_pathway.csv',filter_data)

# 1. 绘图主体
p <- ggplot(filter_data, aes(x = LogP, y = Description, fill = change)) +
  geom_col() +
  scale_fill_manual(values = c("up"="#F08080", "down"="#87CEEB")) +
  theme_bw() +
  labs(x = "-log10(Adjusted P-value)")

# 2. 自定义主题
mytheme <- theme(
  # legend.position = 'none',
  # axis.text.y = element_blank(),
  # axis.ticks.y = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.border = element_blank(),
  axis.line.x = element_line(color = 'grey60', size = 1.1),
  axis.text = element_text(size = 12)
)

p1 <- p + mytheme
p1

ggsave(filename = paste0('4.Osteoblast_care_pathway.png'),w=10,h=6,p1)
ggsave(filename = paste0('4.Osteoblast_care_pathway.pdf'),w=10,h=6,p1)

################################################################################
######DM组和N组BMSC中KEGG
################################################################################
pathways <- c(
  "Osteoclast differentiation",
  "NF-kappa B signaling pathway",
  "TNF signaling pathway",
  "MAPK signaling pathway",
  "Calcium signaling pathway",
  "Cytokine-cytokine receptor interaction",
  "Chemokine signaling pathway",
  "PI3K-Akt signaling pathway",
  "HIF-1 signaling pathway",
  "AGE-RAGE signaling pathway in diabetic complications",
  "Lysosome",
  
  # 脂质代谢（通用标准）
  "Fatty acid metabolism",
  "Glycerolipid metabolism",
  "Glycerophospholipid metabolism",
  "Sphingolipid metabolism",
  "Fatty acid biosynthesis",
  "Biosynthesis of unsaturated fatty acids",
  "Cholesterol metabolism",
  
  # 通用炎症通路（你要的 IL-17 + TLR）
  "IL-17 signaling pathway",
  "Toll-like receptor signaling pathway",
  "NOD-like receptor signaling pathway"
)

# A <- dplyr::bind_rows(up_data, up_data1)
up_data = read.table(file = '/mnt/NexData_file/home/zwl/Project/SC2511076348/require7/Mye_py/Osteoclast/DM_vs_N/2.enrich.KEGG_up.csv',
                     sep = ',',header = TRUE)%>% dplyr::mutate(change = 'up')
down_data = read.table(file = '/mnt/NexData_file/home/zwl/Project/SC2511076348/require7/Mye_py/Osteoclast/DM_vs_N/2.enrich.KEGG_down.csv',
                       sep = ',',header = TRUE)%>% dplyr::mutate(change = 'down')
merge_data = rbind(up_data,down_data)
filter_data1 = merge_data %>% 
  dplyr::filter(Description %in% pathways) %>%
  mutate(LogP = case_when(
    change == "up" ~ -log10(p.adjust),   # up 取负
    change == "down" ~ log10(p.adjust)   # down 不取负
  ))

filter_data <- filter_data1 %>%
  group_by(Description) %>%
  slice_min(p.adjust, n = 1) %>%  # 取最小p.adj
  ungroup()

# 3. 排序：LogP 从小到大 → 从上到下
filter_data <- filter_data %>% arrange(LogP)
filter_data$Description <- factor(filter_data$Description, levels = unique(filter_data$Description))
write.csv(file = '4.Osteoclast_care_pathway.csv',filter_data)

# 1. 绘图主体
p <- ggplot(filter_data, aes(x = LogP, y = Description, fill = change)) +
  geom_col() +
  scale_fill_manual(values = c("up"="#F08080", "down"="#87CEEB")) +
  theme_bw() +
  labs(x = "-log10(Adjusted P-value)")

# 2. 自定义主题
mytheme <- theme(
  legend.position = 'none',
  axis.text.y = element_blank(),
  axis.ticks.y = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.border = element_blank(),
  axis.line.x = element_line(color = 'grey60', size = 1.1),
  axis.text = element_text(size = 12)
)

p1 <- p + mytheme
p1
# 3. 按正负拆分数据（用 filter_data！不是dt！）
up_data   <- filter_data[filter_data$LogP > 0, ]
down_data <- filter_data[filter_data$LogP < 0, ]

p2 <- p1 +
  geom_text(data = up_data,
            aes(x = -0.2, y = Description, label = Description),
            size = 3.5, hjust = 1) +
  geom_text(data = down_data,
            aes(x = 0.2, y = Description, label = Description),
            size = 3.5, hjust = 0)#右对齐

p2
ggsave(filename = paste0('4.Osteoclast_care_pathway.png'),w=10,h=6,p2)
ggsave(filename = paste0('4.Osteoclast_care_pathway.pdf'),w=10,h=6,p2)

################################################################################
######DM组和N组两种巨噬细胞KEGG
################################################################################
pathways <- c(
  # 炎症核心通路
  "TNF signaling pathway",
  "Cytokine-cytokine receptor interaction",
  "Chemokine signaling pathway",
  "Toll-like receptor signaling pathway",
  "NOD-like receptor signaling pathway",
  
  # 脂代谢相关
  "Adipocytokine signaling pathway",
  "Fatty acid metabolism",
  
  # 糖代谢与应激
  "Glycolysis / Gluconeogenesis",
  
  # 缺氧与骨代谢
  "HIF-1 signaling pathway",
  "Osteoclast differentiation",
  
  # 信号与降解
  "MAPK signaling pathway",
  "Lysosome"
)

up_data = read.table(file = '/mnt/NexData_file/home/zwl/Project/SC2511076348/require7/Mye_py/Macro1/DM_vs_N/2.enrich.KEGG_up.csv',
                     sep = ',',header = TRUE)%>% dplyr::mutate(change = 'up')
down_data = read.table(file = '/mnt/NexData_file/home/zwl/Project/SC2511076348/require7/Mye_py/Macro1/DM_vs_N/2.enrich.KEGG_down.csv',
                       sep = ',',header = TRUE)%>% dplyr::mutate(change = 'down')
merge_data = rbind(up_data,down_data)
filter_data1 = merge_data %>% 
  dplyr::filter(Description %in% pathways) %>%
  mutate(LogP = case_when(
    change == "up" ~ -log10(p.adjust),   # up 取负
    change == "down" ~ log10(p.adjust)   # down 不取负
  ))

filter_data <- filter_data1 %>%
  group_by(Description) %>%
  slice_min(p.adjust, n = 1) %>%  # 取最小p.adj
  ungroup()

# 3. 排序：LogP 从小到大 → 从上到下
filter_data <- filter_data %>% arrange(LogP)
filter_data$Description <- factor(filter_data$Description, levels = unique(filter_data$Description))
write.csv(file = '5.Macro1_care_pathway.csv',filter_data)

# 1. 绘图主体
p <- ggplot(filter_data, aes(x = LogP, y = Description, fill = change)) +
  geom_col() +
  scale_fill_manual(values = c("up"="#F08080", "down"="#87CEEB")) +
  theme_bw() +
  labs(x = "-log10(Adjusted P-value)")

# 2. 自定义主题
mytheme <- theme(
  # legend.position = 'none',
  # axis.text.y = element_blank(),
  # axis.ticks.y = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.border = element_blank(),
  axis.line.x = element_line(color = 'grey60', size = 1.1),
  axis.text = element_text(size = 12)
)

p1 <- p + mytheme
p1
ggsave(filename = paste0('5.Macro1_care_pathway.png'),w=10,h=6,p1)
ggsave(filename = paste0('5.Macro1_care_pathway.pdf'),w=10,h=6,p1)


################################################################################
######DM组和N组两种巨噬细胞KEGG##################
################################################################################

up_data = read.table(file = '/mnt/NexData_file/home/zwl/Project/SC2511076348/require7/Mye_py/Macro2/DM_vs_N/2.enrich.KEGG_up.csv',
                     sep = ',',header = TRUE)%>% dplyr::mutate(change = 'up')
down_data = read.table(file = '/mnt/NexData_file/home/zwl/Project/SC2511076348/require7/Mye_py/Macro2/DM_vs_N/2.enrich.KEGG_down.csv',
                       sep = ',',header = TRUE)%>% dplyr::mutate(change = 'down')
merge_data = rbind(up_data,down_data)
filter_data1 = merge_data %>% 
  dplyr::filter(Description %in% pathways) %>%
  mutate(LogP = case_when(
    change == "up" ~ -log10(p.adjust),   # up 取负
    change == "down" ~ log10(p.adjust)   # down 不取负
  ))

filter_data <- filter_data1 %>%
  group_by(Description) %>%
  slice_min(p.adjust, n = 1) %>%  # 取最小p.adj
  ungroup()

# 3. 排序：LogP 从小到大 → 从上到下
filter_data <- filter_data %>% arrange(LogP)
filter_data$Description <- factor(filter_data$Description, levels = unique(filter_data$Description))
write.csv(file = '5.Macro2_care_pathway.csv',filter_data)

# 1. 绘图主体
p <- ggplot(filter_data, aes(x = LogP, y = Description, fill = change)) +
  geom_col() +
  scale_fill_manual(values = c("up"="#F08080", "down"="#87CEEB")) +
  theme_bw() +
  labs(x = "-log10(Adjusted P-value)")

# 2. 自定义主题
mytheme <- theme(
  # legend.position = 'none',
  # axis.text.y = element_blank(),
  # axis.ticks.y = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.border = element_blank(),
  axis.line.x = element_line(color = 'grey60', size = 1.1),
  axis.text = element_text(size = 12)
)

p1 <- p + mytheme
p1
ggsave(filename = paste0('5.Macro2_care_pathway.png'),w=10,h=6,p1)
ggsave(filename = paste0('5.Macro2_care_pathway.pdf'),w=10,h=6,p1)

################################################################################
######小提琴##################
################################################################################

metadata = read.csv('/mnt/NexData_file/home/zwl/Project/SC2511076348/require9/7.Fibro_score.csv',row.names = 1)
metadata_f = metadata %>%
  dplyr::filter(cell_type %in% c("BMSC","Osteoblast"))%>%
  dplyr::select(-c("GOBP_REGULATION_OF_RESPONSE_TO_OXIDATIVE_STRESS","GOBP_REGULATION_OF_CELLULAR_RESPONSE_TO_OXIDATIVE_STRESS"))
metadata_f$group = factor(metadata_f$group,levels = c('DM','N'))

for (i in colnames(metadata_f)[-(1:3)]) {
  print(i)
  for (j in unique(metadata_f$cell_type)) {
    metadata_sub = metadata_f[,c("sample","group","cell_type",i)] %>%
      dplyr::filter(cell_type == j)
    colnames(metadata_sub) = c("sample","group","cell_type",'exp')
    
    p_value = metadata_sub%>%
      group_by(cell_type) %>% 
      rstatix::wilcox_test(formula = exp~group)%>%
      add_significance() %>% add_xy_position(x = 'cell_type')
    
    p = ggplot(metadata_sub, aes(x = cell_type, y = exp, fill = group)) +
      geom_violin(position = position_dodge(0.9), width=0.8,alpha = 0.7) +  # 必加
      geom_boxplot(position = position_dodge(0.9), width=0.2, 
                   color="black", outlier.shape=NA) + # 必加
      geom_jitter(position = position_jitterdodge(jitter.width=0.1, dodge.width=0.9), 
                  color='black', size=0.1, alpha=0.5) +
      stat_pvalue_manual(p_value, label = "p.signif", hide.ns = TRUE) +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size=10))+
      labs(y=i)
    p
    ggsave(filename = paste0('7.1.Fibro_',j,'_',i,".png"),w=6,h=6,p)
    ggsave(filename = paste0('7.1.Fibro_',j,'_',i,".pdf"),w=6,h=6,p)
  }

}
################################################################################

metadata = read.csv('/mnt/NexData_file/home/zwl/Project/SC2511076348/require9/7.Myeloid_score.csv',row.names = 1)
metadata_f = metadata %>%
  dplyr::filter(cell_type %in% c("Macro1","Osteoclast","Macro2"))%>%
  dplyr::select(-c("GOBP_REGULATION_OF_RESPONSE_TO_OXIDATIVE_STRESS","GOBP_REGULATION_OF_CELLULAR_RESPONSE_TO_OXIDATIVE_STRESS"))
metadata_f$group = factor(metadata_f$group,levels = c('DM','N'))

for (i in colnames(metadata_f)[-(1:3)]) {
  print(i)
  for (j in unique(metadata_f$cell_type)) {
    metadata_sub = metadata_f[,c("sample","group","cell_type",i)] %>%
      dplyr::filter(cell_type == j)
    colnames(metadata_sub) = c("sample","group","cell_type",'exp')
    
    p_value = metadata_sub%>%
      group_by(cell_type) %>% 
      rstatix::wilcox_test(formula = exp~group)%>%
      add_significance() %>% add_xy_position(x = 'cell_type')
    
    p = ggplot(metadata_sub, aes(x = cell_type, y = exp, fill = group)) +
      geom_violin(position = position_dodge(0.9), width=0.8,alpha = 0.7) +  # 必加
      geom_boxplot(position = position_dodge(0.9), width=0.2, 
                   color="black", outlier.shape=NA) + # 必加
      geom_jitter(position = position_jitterdodge(jitter.width=0.1, dodge.width=0.9), 
                  color='black', size=0.1, alpha=0.5) +
      stat_pvalue_manual(p_value, label = "p.signif", hide.ns = TRUE) +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size=10))+
      labs(y=i)
    p
    ggsave(filename = paste0('7.3.Myeloid_',j,'_',i,".png"),w=6,h=6,p)
    ggsave(filename = paste0('7.3.Myeloid_',j,'_',i,".pdf"),w=6,h=6,p)
  }
  
}
################################################################################
######小提琴-改##################
################################################################################

metadata = read.csv('/mnt/NexData_file/home/zwl/Project/SC2511076348/require8/Fib_Mye_pathway_score.csv',row.names = 1)
metadata_f = metadata %>%
  dplyr::filter(cell_type %in% c("Macro1","Osteoclast","Macro2"))
metadata_f$group = factor(metadata_f$group,levels = c('DM','N'))

for (i in colnames(metadata_f)[-(1:4)]) {
  print(i)
  for (j in unique(metadata_f$cell_type)) {
    metadata_sub = metadata_f[,c("sample","group","celltype","cell_type",i)] %>%
      dplyr::filter(cell_type == j)
    colnames(metadata_sub) = c("sample","group","celltype","cell_type",'exp')
    
    p_value = metadata_sub%>%
      group_by(cell_type) %>% 
      rstatix::wilcox_test(formula = exp~group)%>%
      add_significance() %>% add_xy_position(x = 'cell_type')
    
    p = ggplot(metadata_sub, aes(x = cell_type, y = exp, fill = group)) +
      geom_violin(position = position_dodge(0.9), width=0.8,alpha = 0.7) +  # 必加
      geom_boxplot(position = position_dodge(0.9), width=0.2, 
                   color="black", outlier.shape=NA) + # 必加
      geom_jitter(position = position_jitterdodge(jitter.width=0.1, dodge.width=0.9), 
                  color='black', size=0.1, alpha=0.5) +
      stat_pvalue_manual(p_value, label = "p.signif", hide.ns = TRUE) +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size=10))+
      labs(y=i)
    p
    ggsave(filename = paste0('8.2.Myeloid_',j,'_',i,".png"),w=6,h=6,p)
    ggsave(filename = paste0('8.2.Myeloid_',j,'_',i,".pdf"),w=6,h=6,p)
  }
  
}
################################################################################

metadata = read.csv('/mnt/NexData_file/home/zwl/Project/SC2511076348/require9/7.Myeloid_score.csv',row.names = 1)
metadata_f = metadata %>%
  dplyr::filter(cell_type %in% c("BMSC","Osteoblast"))%>%
  dplyr::select(-c("GOBP_OSTEOCLAST_DIFFERENTIATION","GOBP_REGULATION_OF_OSTEOCLAST_DIFFERENTIATION"))
metadata_f$group = factor(metadata_f$group,levels = c('DM','N'))

for (i in colnames(metadata_f)[-(1:3)]) {
  print(i)
  for (j in unique(metadata_f$cell_type)) {
    metadata_sub = metadata_f[,c("sample","group","cell_type",i)] %>%
      dplyr::filter(cell_type == j)
    colnames(metadata_sub) = c("sample","group","cell_type",'exp')
    
    p_value = metadata_sub%>%
      group_by(cell_type) %>% 
      rstatix::wilcox_test(formula = exp~group)%>%
      add_significance() %>% add_xy_position(x = 'cell_type')
    
    p = ggplot(metadata_sub, aes(x = cell_type, y = exp, fill = group)) +
      geom_violin(position = position_dodge(0.9), width=0.8,alpha = 0.7) +  # 必加
      geom_boxplot(position = position_dodge(0.9), width=0.2, 
                   color="black", outlier.shape=NA) + # 必加
      geom_jitter(position = position_jitterdodge(jitter.width=0.1, dodge.width=0.9), 
                  color='black', size=0.1, alpha=0.5) +
      stat_pvalue_manual(p_value, label = "p.signif", hide.ns = TRUE) +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size=10))+
      labs(y=i)
    p
    ggsave(filename = paste0('8.1.Myeloid_',j,'_',i,".png"),w=6,h=6,p)
    ggsave(filename = paste0('8.1.Myeloid_',j,'_',i,".pdf"),w=6,h=6,p)
  }
  
}
################################################################################
##############相关性图############################
################################################################################
library(readxl)
obj = readRDS('/mnt/NexData_file/home/zwl/Project/SC2511076348/require7/Myeloids/annoraw.rds')
subdata = subset(obj,cell_type %in% 'Macro1')
target_genes = c('FABP4',"SPP1", "TREM2", "CTSB", "CTSL", "ACP5", "MMP7", "MMP9", "LPL", "APOE")
expr_matrix <- Seurat::GetAssayData(object = subdata, assay = "RNA", slot = "data")[target_genes, ]
expr_df <- as.data.frame(t(expr_matrix))
final_df <- cbind(expr_df, subdata@meta.data[, c("cell_type", "group")])

Blank_data = final_df %>% dplyr::filter(group =='DM')
setwd('/mnt/NexData_file/home/zwl/Project/SC2511076348/require9')
for(gene in target_genes[-1]){
  
  # 1. 计算相关性
  cor_res <- cor.test(Blank_data$FABP4, Blank_data[[gene]])
  
  # 2. 绘图
  gg <- ggplot(Blank_data, aes(x = FABP4, y = .data[[gene]])) +
    geom_jitter(size = 0.5, stroke = 0.5, shape = 16, alpha = 0.6) +
    geom_smooth(method = "lm", se = TRUE, color = "red3", formula = y ~ x) +
    annotate("text", x = Inf, y = Inf,
             label = sprintf("r = %.2f\np = %.2g", 
                             cor_res$estimate, 
                             cor_res$p.value),
             hjust = 1.1, vjust = 1.1, size = 5) +
    labs(x = "FABP4", y = gene) +
    theme_classic(base_size = 12)
  
  # 3. 保存图片（自动命名）
  ggsave(filename = paste0('11.Macro1_DM_FABP4_vs_', gene, '_cor_scatterplot.png'),
         dpi = 200, w=5, h=5, gg)
  ggsave(filename = paste0('11.Macro1_DM_FABP4_vs_', gene, '_cor_scatterplot.pdf'),
         dpi = 200, w=5, h=5, gg)
}

subdata = subset(obj,cell_type %in% 'Macro2')
target_genes = c('FABP4',"SPP1", "TREM2", "CTSB", "CTSL", "ACP5", "MMP7", "MMP9", "LPL", "APOE")
expr_matrix <- Seurat::GetAssayData(object = subdata, assay = "RNA", slot = "data")[target_genes, ]
expr_df <- as.data.frame(t(expr_matrix))
final_df <- cbind(expr_df, subdata@meta.data[, c("cell_type", "group")])

Blank_data = final_df %>% dplyr::filter(group =='DM')
setwd('/mnt/NexData_file/home/zwl/Project/SC2511076348/require9')
for(gene in target_genes[-1]){
  
  # 1. 计算相关性
  cor_res <- cor.test(Blank_data$FABP4, Blank_data[[gene]])
  
  # 2. 绘图
  gg <- ggplot(Blank_data, aes(x = FABP4, y = .data[[gene]])) +
    geom_jitter(size = 0.5, stroke = 0.5, shape = 16, alpha = 0.6) +
    geom_smooth(method = "lm", se = TRUE, color = "red3", formula = y ~ x) +
    annotate("text", x = Inf, y = Inf,
             label = sprintf("r = %.2f\np = %.2g", 
                             cor_res$estimate, 
                             cor_res$p.value),
             hjust = 1.1, vjust = 1.1, size = 5) +
    labs(x = "FABP4", y = gene) +
    theme_classic(base_size = 12)
  
  # 3. 保存图片（自动命名）
  ggsave(filename = paste0('11.Macro2_DM_FABP4_vs_', gene, '_cor_scatterplot.png'),
         dpi = 200, w=5, h=5, gg)
  ggsave(filename = paste0('11.Macro2_DM_FABP4_vs_', gene, '_cor_scatterplot.pdf'),
         dpi = 200, w=5, h=5, gg)
}

################################################################################
##############相关性图2############################
################################################################################
library(readxl)
obj = readRDS('/mnt/NexData_file/home/zwl/Project/SC2511076348/require7/Myeloids/annoraw.rds')
score = read.csv('/mnt/NexData_file/home/zwl/Project/SC2511076348/require8/Fib_Mye_pathway_score.csv',row.names = 1)
target_score = c('HALLMARK_INFLAMMATORY_RESPONSE',"GOBP_OSTEOCLAST_DIFFERENTIATION","GOBP_REGULATION_OF_OSTEOCLAST_DIFFERENTIATION")

subdata = subset(obj,cell_type %in% 'Macro1')
target_genes = c('FABP4',"IL1B","TNF","NLRP3","NFKBIA")
expr_matrix <- Seurat::GetAssayData(object = subdata, assay = "RNA", slot = "data")[target_genes, ]
expr_df <- as.data.frame(t(expr_matrix))
final_df <- cbind(expr_df, subdata@meta.data[, c("cell_type", "group")],score[rownames(expr_df),target_score])

Blank_data = final_df %>% dplyr::filter(group =='DM')
setwd('/mnt/NexData_file/home/zwl/Project/SC2511076348/require9')

for(gene in c(target_genes[-1],target_score)){
  
  # 1. 计算相关性
  cor_res <- cor.test(Blank_data$FABP4, Blank_data[[gene]])
  
  # 2. 绘图
  gg <- ggplot(Blank_data, aes(x = FABP4, y = .data[[gene]])) +
    geom_jitter(size = 0.5, stroke = 0.5, shape = 16, alpha = 0.6) +
    geom_smooth(method = "lm", se = TRUE, color = "red3", formula = y ~ x) +
    annotate("text", x = Inf, y = Inf,
             label = sprintf("r = %.2f\np = %.2g", 
                             cor_res$estimate, 
                             cor_res$p.value),
             hjust = 1.1, vjust = 1.1, size = 5) +
    labs(x = "FABP4", y = gene) +
    theme_classic(base_size = 12)
  
  # 3. 保存图片（自动命名）
  ggsave(filename = paste0('16.Macro1_DM_FABP4_vs_', gene, '_cor_scatterplot.png'),
         dpi = 200, w=5, h=5, gg)
  ggsave(filename = paste0('16.Macro1_DM_FABP4_vs_', gene, '_cor_scatterplot.pdf'),
         dpi = 200, w=5, h=5, gg)
}

subdata = subset(obj,cell_type %in% 'Macro2')
target_genes = c('FABP4',"IL1B","TNF","NLRP3","NFKBIA")
expr_matrix <- Seurat::GetAssayData(object = subdata, assay = "RNA", slot = "data")[target_genes, ]
expr_df <- as.data.frame(t(expr_matrix))
final_df <- cbind(expr_df, subdata@meta.data[, c("cell_type", "group")],score[rownames(expr_df),target_score])

Blank_data = final_df %>% dplyr::filter(group =='DM')
setwd('/mnt/NexData_file/home/zwl/Project/SC2511076348/require9')

for(gene in c(target_genes[-1],target_score)){
  
  # 1. 计算相关性
  cor_res <- cor.test(Blank_data$FABP4, Blank_data[[gene]])
  
  # 2. 绘图
  gg <- ggplot(Blank_data, aes(x = FABP4, y = .data[[gene]])) +
    geom_jitter(size = 0.5, stroke = 0.5, shape = 16, alpha = 0.6) +
    geom_smooth(method = "lm", se = TRUE, color = "red3", formula = y ~ x) +
    annotate("text", x = Inf, y = Inf,
             label = sprintf("r = %.2f\np = %.2g", 
                             cor_res$estimate, 
                             cor_res$p.value),
             hjust = 1.1, vjust = 1.1, size = 5) +
    labs(x = "FABP4", y = gene) +
    theme_classic(base_size = 12)
  
  # 3. 保存图片（自动命名）
  ggsave(filename = paste0('16.Macro2_DM_FABP4_vs_', gene, '_cor_scatterplot.png'),
         dpi = 200, w=5, h=5, gg)
  ggsave(filename = paste0('16.Macro2_DM_FABP4_vs_', gene, '_cor_scatterplot.pdf'),
         dpi = 200, w=5, h=5, gg)
}
