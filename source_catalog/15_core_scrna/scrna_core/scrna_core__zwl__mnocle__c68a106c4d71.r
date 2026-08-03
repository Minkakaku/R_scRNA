setwd('/mnt/NexData_file/home/zwl/Project/other/SC2507294017_xzx/')
library(monocle)
rm(list = ls())
gc()
load("/mnt/NexData_file/home/zwl/Project/other/mono.Rdata")
p = plot_cell_trajectory(cds, 
                         color_by = "Pseudotime", 
                         size = 0.01, 
                         show_backbone = TRUE,
                         reduction_method = "DDRTree")
p

cnv <- read.csv("./gene_meta.csv",row.names = 1)
colnames(cnv)
cds@phenoData@data$SHMT1 <- cnv[rownames(cds@phenoData@data),'SHMT1']
cds@phenoData@data[['SHMT2']] <- cnv[rownames(cds@phenoData@data),'SHMT2']

for (i in colnames(cnv)) {
  cds@phenoData@data[[i]] <- cnv[rownames(cds@phenoData@data),i]
}

for (m in colnames(cnv)) {
  p = plot_cell_trajectory(cds,x=1,y=2,color_by = m,cell_size = .5)+
    scale_color_gradientn(colours = c('grey','#ff7f0e','firebrick'),
                          values = c(0,0.5,1,1.5))
  ggsave(filename = paste0("monocle_pseudotime_",m,".png"),p,h=7,w=10,dpi = 300)
}
p = plot_cell_trajectory(cds,x=1,y=2,color_by = 'SHMT2',cell_size = .5)+
  scale_color_gradientn(colours = c('grey','#ff7f0e','firebrick'),
                                               values = c(0,0.5,1,1.5))
p
ggsave(filename = "monocle_pseudotime_SHMT1.png",p,h=7,w=10,dpi = 300)