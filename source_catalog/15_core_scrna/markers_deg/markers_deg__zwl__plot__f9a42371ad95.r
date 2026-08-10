setwd('/mnt/NexData_file/home/zwl/Project/SC2511076348/require6/')
rm(list=ls())
gc()
metadata= read.csv('/mnt/NexData_file/home/zwl/Project/SC2511076348/require6/FABP4_exp.csv',head = T,row.names = 1)
cellname = c("Mono2","BMSC","Endothelial cells","Macro1","Mono1","Osteoblast","Osteoclast","Macro2")
for (i in cellname) {
  metadata_long = metadata %>% dplyr::filter(cell_type == i) %>%
    pivot_longer( cols =FABP4,
                  names_to = "gene",
                  values_to = "exp")
  # metadata_long$group = do::Replace(metadata_long$group,pattern = c('G1:Control','G2:Phage','G3:hvKp','G4:hvKp+Phage'))
  metadata_long$group = factor(metadata_long$group,levels = c('DM','N'))
  
  p_value = metadata_long%>%
    group_by(gene) %>% 
    rstatix::wilcox_test(formula = exp~group)%>%
    add_significance() %>% add_xy_position()
  
  
  p = ggplot(metadata_long%>%filter(gene=='FABP4'), aes(x = group, y = exp, colour = group)) +
    geom_violin(position = position_dodge(width = 0.8), scale = "width", alpha = 0.5, width = 0.6, size = 0.8) + 
    geom_jitter(position = position_jitterdodge(jitter.width = 0.6, dodge.width = 0.8), alpha = 0.5,size=0.5,stroke=0.5, shape=16) +
    geom_boxplot(position = position_dodge(width = 2),fill='grey', outlier.size = 0, width = 0.04, show.legend = FALSE) +
    stat_pvalue_manual(p_value%>%filter(gene=='FABP4'), label = "p.signif",hide.ns = T)+
    scale_color_manual(values = c("#F18F01","#2E86AB")  , guide = "none") + 
    theme(		axis.text.x  = element_text(angle = 0, hjust = 0.5,size = 5)) +
    theme_classic()+
    labs(y= paste0('FABP4 express in ',i))
  p
  ggsave(filename = paste0("./2.FABP4_gene_exp_",i,".png"),w=5,h=4,p)
  ggsave(filename = paste0("./2.FABP4_gene_exp_",i,".pdf"),w=5,h=4,p)
}

################monocle###########################
HSMM = readRDS(paste0('/mnt/NexData_file/home/zwl/Project/SC2511076348/require6/Fib_mono/monocle.rds'))
deg <- read.csv(paste0('/mnt/NexData_file/home/zwl/Project/SC2511076348/require6/Fib_mono/monocle.DEG.xls'),sep="\t")
rownames(deg) = deg$Genes
ordergene <- rownames(deg)[1:min(1000,nrow(deg))]
# 重新排序
plot_cell_trajectory(HSMM, color_by = "State")+
  facet_wrap(~cell_type,ncol = 4)
HSMM_myo <- orderCells(HSMM, root_state = 3)

# 重新画图
setwd('/mnt/NexData_file/home/zwl/Project/SC2511076348/require6/Fib_mono')
getwd()
clustercolname = 'cell_type'

p = plot_cell_trajectory(HSMM_myo,color_by="Pseudotime",size=1,show_backbone=TRUE)
ggsave("monocle.trajectory.pseudotime.pdf",p,width = 7,height = 7.5)
ggsave("monocle.trajectory.pseudotime.png",p,width = 7,height = 7.5)

p = plot_complex_cell_trajectory(HSMM_myo,x=1,y=2,color_by = clustercolname)
ggsave("monocle.trajectory.tree.pdf",p,width = 7,height = 7.5)
ggsave("monocle.trajectory.tree.png",p,width = 7,height = 7.5)
ncols = ceiling(sqrt(length(unique(HSMM_myo@phenoData@data[,clustercolname]))))
nrows = ceiling(length(unique(HSMM_myo@phenoData@data[,clustercolname]))/ncols)
p = plot_complex_cell_trajectory(HSMM_myo,x=1,y=2,color_by = clustercolname)+facet_wrap(vars(!!sym(clustercolname)),nrow=nrows)
ggsave("monocle.trajectory.tree.faceted.pdf",p,width = 4*ncols,height = 4.2*nrows)
ggsave("monocle.trajectory.tree.faceted.png",p,width = 4*ncols,height = 4.2*nrows)

pseudotime_df = HSMM_myo@phenoData@data[,c("sample",clustercolname,"Size_Factor","num_genes_expressed","Pseudotime","State")]
df = data.frame(Cells=rownames(pseudotime_df),pseudotime_df)
write.table(df,"3.monocle.metadata_reorder.xls",sep = "\t",quote = F,row.names = F)

saveRDS(file = '3.obj_reorder.rds',HSMM_myo,compress = F)

gs = c('CD36','PLIN2','LPL','PPARG','APOE','APOC1','CD9',
       'LGALS3','FABP4','FABP5','ABCA1','ABCG1')
p = plot_cell_trajectory(HSMM_myo,markers=gs,use_color_gradient=T)
  # scale_color_gradientn(colours = c('grey','#ff7f0e','firebrick'),
                        # values = c(0,0.5,1,1.5))
ggsave("3.gene_monocle.trajectory.tree.faceted.pdf",p,width = 9,height = 4)
ggsave("3.gene_monocle.trajectory.tree.faceted.png",p,width = 9,height = 4)
gs1 = gs[gs %in%rownames(HSMM_myo)]
plot_genes_in_pseudotime(HSMM_myo[gs1,],color_by = "cell_type",nrow= 4, ncol = NULL )

