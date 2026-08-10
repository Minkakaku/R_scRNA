library(Seurat)
library(ggplot2)
library(ggpubr)

setwd('/mnt/NexData_file/home/wmx/Project/SC2604092832-01/shouhou4/')
# L2-3 IT-MT 细胞类型
down_gene = read.table(file = "/mnt/NexData_file/home/wmx/Project/SC2604092832-01/shouhou4/scenic/L2_3_IT_MT/res/3.commen_regulon_gene/0.commen_down_regulon_csv.txt")
up_gene = read.table(file = "/mnt/NexData_file/home/wmx/Project/SC2604092832-01/shouhou4/scenic/L2_3_IT_MT/res/3.commen_regulon_gene/0.commen_up_regulon_csv.txt")
scRNA=readRDS("/mnt/NexData_file/home/wmx/Project/SC2604092832-01/shouhou4/L2-3_IT-MT_anno_raw.rds")
unique(scRNA$cell_type)
scRNA$g_celltype = paste0(scRNA$cell_type,'_',scRNA$group)
outdir = '/mnt/NexData_file/home/wmx/Project/SC2604092832-01/shouhou4/scenic/L2_3_IT_MT/res/4.commen_regulon_score/'
dir.create(outdir, recursive = T)

All_gene = rbind(down_gene,up_gene)
for (gi in 1:length(All_gene$V1)) {
  genelist = unlist(strsplit(All_gene[gi,'V2'], split = ","))
  scRNA <- AddModuleScore(  object = scRNA,  features = list(genelist),  
                                  name = paste0(All_gene[gi,'V1'],' target gene score'))
  plot = scRNA@meta.data # %>% dplyr::filter(cell_type %in% c('L5/6 NP2','L2-3 IT-MT'))

my_comparisons <- list(
  c("CASP3_P80", "CT_P80"),
  c("CASP3_P42", "CT_P42")
)

y_col <- paste0(All_gene[gi, "V1"], " target gene score1")

y_max <- max(plot[[y_col]], na.rm = TRUE)
y_min <- min(plot[[y_col]], na.rm = TRUE)
y_range <- y_max - y_min
plot$group <- factor(plot$group, levels = c("CT_P42", "CASP3_P42", "CT_P80", "CASP3_P80"))
p <- ggplot(
  plot,
  aes(
    x = group,
    y = .data[[y_col]]
  )
) +
  geom_violin(
    aes(fill = group),
    alpha = 1,
    width = 0.8
  ) +
  geom_jitter(
    color = "black",
    size = 0.1,
    alpha = 0.5,
    width = 0.15
  ) +
  geom_boxplot(
    width = 0.2,
    color = "black",
    outlier.shape = NA
  ) +
  ggpubr::stat_compare_means(
    comparisons = my_comparisons,
    method = "t.test",
    label = "p.format",
    size = 4.5,
    label.y = c(
      y_max + 0.08 * y_range,
      y_max + 0.25 * y_range
    ),
    tip.length = 0.02,
    bracket.size = 0.6
  ) +
  coord_cartesian(
    ylim = c(
      y_min,
      y_max + 0.40 * y_range
    )
  )
  p
  ggsave(paste0(outdir,All_gene[gi,'V1'],' target gene score_volcano.pdf'), plot = p, width = 7, height = 6)
  ggsave(paste0(outdir,All_gene[gi,'V1'],' target gene score_volcano.png'), plot = p, width = 7, height = 6)
  
  plot$sample <- factor(plot$sample, levels = c("CT_1", "CT_2","CT_3","CASP3_1","CASP3_2","CASP3_3" ,"CT_P80_1","CT_P80_2","CT_P80_3", "CASP3_P80_1","CASP3_P80_2","CASP3_P80_3"))
  p2 = ggplot(plot, aes(x = sample,
                       y = .data[[paste0(All_gene[gi,'V1'],' target gene score1')]],
                       group = sample)) +
    geom_violin(aes(fill = group)) +
    geom_jitter(
      aes(group = sample),
      color = 'black',
      size = 0.1,
      alpha = 0.8
    )+
    geom_boxplot(width = 0.2,  color = "black", outlier.shape = NA) +
    labs(
      x = "sample",
      y = paste0(All_gene[gi,'V1'],' target gene score'),
      color = "Group"
    ) +
    theme_classic() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
  p2
  ggsave(paste0(outdir,All_gene[gi,'V1'],' target gene score_sample_volcano.pdf'), plot = p2, width = 10, height = 6)
  ggsave(paste0(outdir,All_gene[gi,'V1'],' target gene score_sample_volcano.png'), plot = p2, width = 10, height = 6)
  
}

# L2-3_IT_Cmss1+ 细胞类型
down_gene = read.table(file = "/mnt/NexData_file/home/wmx/Project/SC2604092832-01/shouhou4/scenic/L2_3_IT_Cmss1/res/3.commen_regulon_gene/0.commen_down_regulon_csv.txt")
up_gene = read.table(file = "/mnt/NexData_file/home/wmx/Project/SC2604092832-01/shouhou4/scenic/L2_3_IT_Cmss1/res/3.commen_regulon_gene/0.commen_up_regulon_csv.txt")
scRNA=readRDS("/mnt/NexData_file/home/wmx/Project/SC2604092832-01/shouhou4/L2-3_IT_Cmss1_anno_raw.rds")
unique(scRNA$cell_type)
scRNA$g_celltype = paste0(scRNA$cell_type,'_',scRNA$group)
outdir = '/mnt/NexData_file/home/wmx/Project/SC2604092832-01/shouhou4/scenic/L2_3_IT_Cmss1/res/4.commen_regulon_score/'
dir.create(outdir, recursive = T)

All_gene = rbind(down_gene,up_gene)
for (gi in 1:length(All_gene$V1)) {
  genelist = unlist(strsplit(All_gene[gi,'V2'], split = ","))
  scRNA <- AddModuleScore(  object = scRNA,  features = list(genelist),  
                                  name = paste0(All_gene[gi,'V1'],' target gene score'))
  plot = scRNA@meta.data # %>% dplyr::filter(cell_type %in% c('L5/6 NP2','L2-3 IT-MT'))

my_comparisons <- list(
  c("CASP3_P80", "CT_P80"),
  c("CASP3_P42", "CT_P42")
)

y_col <- paste0(All_gene[gi, "V1"], " target gene score1")

y_max <- max(plot[[y_col]], na.rm = TRUE)
y_min <- min(plot[[y_col]], na.rm = TRUE)
y_range <- y_max - y_min
plot$group <- factor(plot$group, levels = c("CT_P42", "CASP3_P42", "CT_P80", "CASP3_P80"))
p <- ggplot(
  plot,
  aes(
    x = group,
    y = .data[[y_col]]
  )
) +
  geom_violin(
    aes(fill = group),
    alpha = 1,
    width = 0.8
  ) +
  geom_jitter(
    color = "black",
    size = 0.1,
    alpha = 0.5,
    width = 0.15
  ) +
  geom_boxplot(
    width = 0.2,
    color = "black",
    outlier.shape = NA
  ) +
  ggpubr::stat_compare_means(
    comparisons = my_comparisons,
    method = "t.test",
    label = "p.format",
    size = 4.5,
    label.y = c(
      y_max + 0.08 * y_range,
      y_max + 0.25 * y_range
    ),
    tip.length = 0.02,
    bracket.size = 0.6
  ) +
  coord_cartesian(
    ylim = c(
      y_min,
      y_max + 0.40 * y_range
    )
  )
  p
  ggsave(paste0(outdir,All_gene[gi,'V1'],' target gene score_volcano.pdf'), plot = p, width = 7, height = 6)
  ggsave(paste0(outdir,All_gene[gi,'V1'],' target gene score_volcano.png'), plot = p, width = 7, height = 6)
  
  plot$group <- factor(plot$sample, levels = c("CT_1", "CT_2","CT_3","CASP3_1","CASP3_2","CASP3_3" ,"CT_P80_1","CT_P80_2","CT_P80_3", "CASP3_P80_1","CASP3_P80_2","CASP3_P80_3"))
  p2 = ggplot(plot, aes(x = sample,
                       y = .data[[paste0(All_gene[gi,'V1'],' target gene score1')]],
                       group = sample)) +
    geom_violin(aes(fill = group)) +
    geom_jitter(
      aes(group = sample),
      color = 'black',
      size = 0.1,
      alpha = 0.8
    )+
    geom_boxplot(width = 0.2,  color = "black", outlier.shape = NA) +
    labs(
      x = "sample",
      y = paste0(All_gene[gi,'V1'],' target gene score'),
      color = "Group"
    ) +
    theme_classic() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
  p2
  ggsave(paste0(outdir,All_gene[gi,'V1'],' target gene score_sample_volcano.pdf'), plot = p2, width = 10, height = 6)
  ggsave(paste0(outdir,All_gene[gi,'V1'],' target gene score_sample_volcano.png'), plot = p2, width = 10, height = 6)
  
}
