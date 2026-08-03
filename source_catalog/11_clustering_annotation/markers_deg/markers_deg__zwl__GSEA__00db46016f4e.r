# conda activate diffgene
#加载包
library(clusterProfiler)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(ggplot2)
library(ggrepel)
# library(ReactomePA)
library(enrichplot)
library(DOSE)
library(ComplexUpset)
library(UpSetR)
library(tidyr)
library(dplyr)
library(grDevices)
library(grid) 
rm(list=ls())
gc()



input = '/mnt/NexData_file/home/zwl/Project/SC2512157505/require3/Allcelltype_DEG/'
setwd(input)
output = input
files <- list.files(input, pattern = "*_All_expressed_genes.xls", recursive = TRUE, full.names = TRUE)


gmt_df <- clusterProfiler::read.gmt("/mnt/NexData_file/home/zwl/Project/FWAE20242338/mmu.pathway.gmt")
for(n in files){
  target_dirs <- file.path( basename(dirname(dirname(n))),basename(dirname(n)))
  print(target_dirs)
  setwd(output)
  dir.create(target_dirs,recursive = T)
  outdir = paste0(output,target_dirs,'/')
  
  data = read.table(file = n,sep = '\t',header = T)
  
  data_all_sort <- data %>% arrange(desc(logfoldchanges))
  geneList <- data_all_sort$logfoldchanges
  names(geneList) <- data_all_sort$names
  
  res <- GSEA(geneList,TERM2GENE = gmt_df,  
              minGSSize = 5,
              pvalueCutoff = 0.1,pAdjustMethod = 'BH',
              maxGSSize = 500)
  write.csv(file =paste0(outdir,"3.GSEA_res.csv"),res@result)
  num = length(res@result$ID)
  print(num)
  if (num >= 1) {
    # 最多画前10个
    n_plot <- min(num, 10)
    
    for (i in 1:n_plot) {
      p <- enrichplot::gseaplot2(res, geneSetID = i,
                                 color = 'green2',
                                 title = res$Description[i])
      
      cairo_pdf(filename = paste0(outdir,'3.',res@result[["ID"]][i],'_GSEA_plot.pdf'),
                width=8,height=6)
      print(p)
      grid.text(label = paste0("NES = ", round(res$NES[i],2), "\n", 
                               "p value = ", signif(res$pvalue[i], 2),"\n", 
                               "FDR = ", signif(res$p.adjust[i], 2)),
                x = unit(0.95,"npc"), y = unit(0.94,"npc"),
                just = c("right","top"), gp = gpar(fontsize=13,fontface="bold"))
      dev.off()
    }
  }else{
    print('数量小于1')
  }

}


