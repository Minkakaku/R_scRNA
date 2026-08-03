setwd('/mnt/NexData_file/home/zwl/Project/FW2024-285/')
library(scMetabolism)
library(ggplot2)
library(rsvd)
library(Seurat)
rm(list=ls());gc()
#exchange h5ad to rds############################################
library(sceasy)
library(Seurat)
library(reticulate)
# sceasy::convertFormat('/mnt/NexData_file/home/zwl/Project/FW2024-285/Mono_adata_raw.h5ad', from="anndata", to="seurat",
#                       outFile='Mono_adata_raw.rds')
obj <- readRDS(file = './Mono_adata_raw.rds')

# 删除不认识的基因
valid_genes <- !grepl("^ENSRNOG", rownames(obj))

# 保留符合条件的基因
obj_filtered <- obj[valid_genes, ]

# 查看过滤前后的基因数量
cat("过滤前基因数量:", nrow(obj), "\n")
cat("过滤后基因数量:", nrow(obj_filtered), "\n")

Seurat::DimPlot(obj_filtered,reduction = 'umap',group.by = c('cell_type','sample'))
####################################同源基因转换
library(homologene)
new_gene = mouse2human(rownames(obj_filtered))
counts <- GetAssayData(obj_filtered, slot = "counts")[new_gene$mouseGene,]
rownames(counts) <- new_gene$humanGene  # 更新基因名
rownames(counts) <- make.unique(rownames(counts), sep = "_")


countexp.Seurat<-sc.metabolism(countexp = counts, 
                               method = "AUCell",
                               imputation = F, 
                               ncores = 2,
                               metabolism.type = "KEGG")
saveRDS(file = './scMetabolism_KEGG.rds',compress = F,countexp.Seurat)

#只有代谢通路，不好用

####################################gene set####################################
forplot1 = read.csv('./Growth_ISG_geneset_score.csv',row.names = 1)
forplot1$group = factor(forplot1$group,levels = c('Blank','Control1','TTT'))
custom_colors <- brewer.pal(length(unique(forplot1$group)), "Set3")
names(custom_colors) <- unique(forplot1$group)
library(ggpubr)
library(rstatix)
pdf(file = './Growth_ISG_geneset_score.pdf',width = 4,height = 3)
for (i in colnames(forplot1)[-1]){
  forplot = forplot1[,c('group',i)]
  
  p_value = forplot%>%
    rstatix::t_test(formula = as.formula(paste0(i, " ~ group")))%>%
    add_significance() %>% add_xy_position(step.increase =0.4)
  
  gg = ggplot(forplot, aes_string('group', i)) +
    geom_violin(aes(fill = group),position = position_dodge(width = 0.8), scale = "width", alpha = 0.6, width = 0.6, size = 0.8) +
    geom_boxplot(aes(fill = group),position = position_dodge(width = 0.8), outlier.size = 0, width = 0.2, show.legend = FALSE) +
    #geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8), size=0.5,stroke=0.5, shape=16) +
    scale_color_manual(values = custom_colors)+
    scale_fill_manual(values = custom_colors)+
    stat_pvalue_manual(p_value, label = "p.adj.signif",hide.ns = T)+
    theme_classic()+labs(y=i,x=NULL)#+
  #theme(text = element_text(family = "TimesNewRoman"))
  print(gg)
  #ggsave(filename = './Pdcd1_Ifng_Cd69_CD8T_vlnplot.pdf',w=7,h=4,dpi = 200,gg)
}
dev.off()


######################计算相关性########################################
library(readxl)
forplot_O = read.csv('/mnt/NexData_file/home/zwl/Project/FW2024-285/mono_sc.csv',sep = '\t',row.names = 1)# 
all_data = cbind(forplot_O[rownames(forplot1),],forplot1[,-1]) 
Blank_data = all_data %>% dplyr::filter(group =='Blank')

cor_results <- list(
  ISG = cor.test(Blank_data$Oxidative.phosphorylation, Blank_data$ISG_score),
  growth = cor.test(Blank_data$Oxidative.phosphorylation, Blank_data$growth_score)
)

gg <- ggplot(Blank_data, aes(x = Oxidative.phosphorylation, y = ISG_score)) +
  geom_jitter(size = 0.5, stroke = 0.5, shape = 16, alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE, color = "red3", formula = y ~ x) +
  annotate("text", x = Inf, y = Inf,
           label = sprintf("r = %.2f\np = %.2g", 
                           cor_results$ISG$estimate, 
                           cor_results$ISG$p.value),
           hjust = 1.1, vjust = 1.1, size = 5) +
  theme_classic(base_size = 12)
ggsave(filename = 'Mono_Blank_Oxidative.phosphorylation vs ISG cor scatterplot.png',dpi = 200,w=5,h=5,gg)
ggsave(filename = 'Mono_Blank_Oxidative.phosphorylation vs ISG cor scatterplot.pdf',dpi = 200,w=5,h=5,gg)
gg_1 = ggplot(all_data, aes(x = Oxidative.phosphorylation, y = ISG_score,colour = group)) +
  geom_jitter(size = 0.8, stroke = 0.5, shape = 16, alpha = 1)+
  theme_classic(base_size = 12)
ggsave(filename = 'ALL_Oxidative.phosphorylation vs ISG cor scatterplot.png',dpi = 200,w=5,h=5,gg_1)
ggsave(filename = 'ALL_Oxidative.phosphorylation vs ISG cor scatterplot.pdf',dpi = 200,w=5,h=5,gg_1)

gg1 <- ggplot(Blank_data, aes(x = Oxidative.phosphorylation, y = growth_score)) +
  geom_jitter(size = 0.5, stroke = 0.5, shape = 16, alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE, color = "red", formula = y ~ x) +
  annotate("text", x = Inf, y = Inf,
           label = sprintf("r = %.2f\np = %.2g", 
                           cor_results$ISG$estimate, 
                           cor_results$ISG$p.value),
           hjust = 1.1, vjust = 1.1, size = 5) + 
  theme_classic(base_size = 12)
gg1
ggsave(filename = 'Mono_Blank_Oxidative.phosphorylation vs Growth cor scatterplot.png',dpi = 200,w=5,h=5,gg1)
ggsave(filename = 'Mono_Blank_Oxidative.phosphorylation vs Growth cor scatterplot.pdf',dpi = 200,w=5,h=5,gg1)
gg1_1 = ggplot(all_data, aes(x = Oxidative.phosphorylation, y = growth_score,colour = group)) +
  geom_jitter(size = 0.8, stroke = 0.5, shape = 16, alpha = 1)+
  theme_classic(base_size = 12)
ggsave(filename = 'ALL_Oxidative.phosphorylation vs Growth cor scatterplot.png',dpi = 200,w=5,h=5,gg1_1)
ggsave(filename = 'ALL_Oxidative.phosphorylation vs Growth cor scatterplot.pdf',dpi = 200,w=5,h=5,gg1_1)

