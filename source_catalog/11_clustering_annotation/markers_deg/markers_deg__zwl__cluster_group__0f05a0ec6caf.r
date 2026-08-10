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
rm(list=ls())
gc()

clus_enrichment <- function(gene_set,orgdb="org.Mm.eg.db",trend='ALL',keggdb = "mmu",outdir){
  
  
  if(length(gene_set) > 0){
    GO_ALL <- enrichGO(gene = gene_set,  #基因列表(转换的ID)
                       keyType = "SYMBOL",  #指定的基因ID类型，默认为ENTREZID
                       OrgDb=orgdb,  #物种对应的org包
                       ont = "ALL",   #CC细胞组件，MF分子功能，BP生物学过程
                       pvalueCutoff = 0.05,  #p值阈值
                       pAdjustMethod = "BH",  #多重假设检验校正方式
                       minGSSize = 3,   #注释的最小基因集，默认为10
                       maxGSSize = 500,  #注释的最大基因集，默认为500
                       qvalueCutoff = 1,  #q值阈值
                       readable = TRUE)  #基因ID转换为基因名
    write.table(GO_ALL, paste0(outdir, "1.GO_",trend,"_enrich.go.csv"), sep = ',', row.names = FALSE, quote = TRUE)
    
    if(!is.null(GO_ALL) && nrow(GO_ALL) > 0){
      dotplot(GO_ALL,  #GO富集分析结果
              x = "GeneRatio",  #横坐标,默认GeneRation,也可以为Count
              color = "p.adjust",  #右纵坐标,默认p.adjust,也可以为pvalue和qvalue
              showCategory = 10,  #展示前20个点，默认为10个
              size = NULL,  #点的大小
              title = "GO_dotplot",  #设置图片的标题
              split='ONTOLOGY',
              label_format = 80)+
        facet_grid(ONTOLOGY~.,scale="free")
      
      ggsave(paste0(outdir,"1.GO_",trend,"_dotplot.pdf"),width = 12,height = 9)
      ggsave(paste0(outdir,"1.GO_",trend,"_dotplot.png"),width = 12,height = 9)
      
      barplot(GO_ALL,  #GO富集分析结果
              x = "Count",  #横坐标,默认Count,也可以为GeneRation
              color = "p.adjust",  #右纵坐标,默认p.adjust,也可以为pvalue和qvalue
              showCategory = 10,  #展示前20个，默认为10个
              size = NULL,  
              title = "GO_barplot",  #设置图片的标题
              split='ONTOLOGY',
              label_format = 80)+
        facet_grid(ONTOLOGY~.,scale="free")
      
      ggsave(paste0(outdir,"1.GO_",trend,"_barplot.pdf"),width = 12,height = 9)
      ggsave(paste0(outdir,"1.GO_",trend,"_barplot.png"),width = 12,height = 9)
      
      # upset_data <- GO_ALL@result[1:min(10, nrow(GO_ALL@result)), ]
      # gene_go_long <- upset_data %>% select(Description, geneID) %>% separate_rows(geneID, sep = "/") %>% mutate(present = 1) %>% distinct()  # 去重
      # gene_go_wide <- gene_go_long %>% pivot_wider( id_cols = geneID, names_from = Description, values_from = present, values_fill = 0 )
      # upset_matrix <- as.data.frame(gene_go_wide[, -1])  # 去掉geneID列
      # rownames(upset_matrix) <- gene_go_wide$geneID
      # 
    # # 5. 绘制UpSet图
    #   pdf(paste0(outdir,"1.GO_",trend,"_upsetplot.pdf"), width = 20, height = 10)
    #   print(upset(upset_matrix, nsets = 10, 
    #               number.angles = 30, 
    #               point.size = 3.5, 
    #               line.size = 2, 
    #               mainbar.y.label = "Intersections Size", 
    #               sets.x.label = "Number of DEGs", 
    #               order.by = "freq",
    #               sets.bar.color = "#1f77b4",
    #               matrix.color = "#ff7f0e",
    #               text.scale = c(1.8, 1.6, 1.3, 1.0, 1.8, 1.3),
    #               mb.ratio = c(0.6, 0.4)))
    #   dev.off()
    #   png(paste0(outdir,"1.GO_",trend,"_upsetplot.png"), width = 20, height = 10, units = "in", res = 300)
    #   print(upset(upset_matrix, nsets = 10, 
    #               number.angles = 30, 
    #               point.size = 3.5, 
    #               line.size = 2, 
    #               mainbar.y.label = "Intersections Size", 
    #               sets.x.label = "Number of DEGs", 
    #               order.by = "freq",
    #               sets.bar.color = "#1f77b4",
    #               matrix.color = "#ff7f0e",
    #               text.scale = c(1.8, 1.6, 1.3, 1.0, 1.8, 1.3),
    #               mb.ratio = c(0.6, 0.4)))
    #   dev.off()
    #   
    }
    
    gene_symbol <- tryCatch({
      bitr(geneID = gene_set, fromType="SYMBOL", toType="ENTREZID", OrgDb=orgdb)
    }, error = function(e) {
      message("no gene")
      return(NULL)
    })
    
    # 如果转换失败，直接退出函数，不往下走
    if(is.null(gene_symbol)) return()
    gene_KEGG <- gene_symbol[,2]
    
    KEGG<- enrichKEGG(gene = gene_KEGG,   #基因列表(同GO) 
                      organism = keggdb,  #物种
                      keyType = "kegg",  #指定的基因ID类型，默认为kegg
                      minGSSize = 3, 
                      maxGSSize = 500,
                      pvalueCutoff = 0.05,  
                      pAdjustMethod = "BH",
                      qvalueCutoff = 1,
                      use_internal_data =T)
    
    id_map <- setNames(gene_symbol$SYMBOL, gene_symbol$ENTREZID)
    if(!is.null(KEGG) && nrow(KEGG) > 0){
      KEGG@result$geneID <- sapply(KEGG@result$geneID,
                                 function(x) paste(id_map[strsplit(x, "/")[[1]]], collapse = "/"))
      KEGG@result$Description <- sub(" - Mus musculus \\(house mouse\\)$", "", KEGG@result$Description)
    }
    write.table(KEGG, paste0(outdir,"2.enrich.KEGG_",trend,".csv"), sep = ',', row.names = FALSE, quote = TRUE)
    
    
    if(!is.null(KEGG) && nrow(KEGG) > 0){
      p <- dotplot(KEGG,  #GO富集分析结果
                   x = "GeneRatio",  #横坐标,默认GeneRation,也可以为Count
                   color = "p.adjust",  #右纵坐标,默认p.adjust,也可以为pvalue和qvalue
                   showCategory = 10,  #展示前20个点，默认为10个
                   size = NULL,  #点的大小
                   title = "KEGG_dotplot",  #设置图片的标题
                   label_format = 100)
      p + guides(colour = "Count")
      ggsave(paste0(outdir,"2.KEGG_",trend,"_dotplot.pdf"),width = 12,height = 5)
      ggsave(paste0(outdir,"2.KEGG_",trend,"_dotplot.png"),width = 12,height = 5)
      
      barplot(KEGG,  #GO富集分析结果
              x = "Count",  #横坐标,默认Count,也可以为GeneRation
              color = "p.adjust",  #右纵坐标,默认p.adjust,也可以为pvalue和qvalue
              showCategory = 10,  #展示前20个，默认为10个
              size = NULL,  
              title = "KEGG_barplot",  #设置图片的标题
              label_format = 80)
      
      ggsave(paste0(outdir,"2.KEGG_",trend,"_barplot.pdf"),width = 12,height = 5)
      ggsave(paste0(outdir,"2.KEGG_",trend,"_barplot.png"),width = 12,height = 5)
      
      # upset_data <- KEGG@result[1:min(10, nrow(KEGG@result)), ]
      # gene_go_long <- upset_data %>% select(Description, geneID) %>% separate_rows(geneID, sep = "/") %>% mutate(present = 1) %>% distinct()  # 去重
      # gene_go_wide <- gene_go_long %>% pivot_wider( id_cols = geneID, names_from = Description, values_from = present, values_fill = 0 )
      # upset_matrix <- as.data.frame(gene_go_wide[, -1])  # 去掉geneID列
      # rownames(upset_matrix) <- gene_go_wide$geneID
      # 
      # # 5. 绘制UpSet图
      # pdf(paste0(outdir,"2.KEGG_",trend,"_upsetplot.pdf"), width = 20, height = 10)
      # print(upset(upset_matrix, nsets = 10, 
      #             number.angles = 30, 
      #             point.size = 3.5, 
      #             line.size = 2, 
      #             mainbar.y.label = "Intersections Size", 
      #             sets.x.label = "Number of DEGs", 
      #             order.by = "freq",
      #             sets.bar.color = "#1f77b4",
      #             matrix.color = "#ff7f0e",
      #             text.scale = c(1.8, 1.6, 1.3, 1.0, 1.8, 1.3),
      #             mb.ratio = c(0.6, 0.4)))
      # dev.off()
      # png(paste0(outdir,"2.KEGG_",trend,"_upsetplot.png"), width = 20, height = 10, units = "in", res = 300)
      # print(upset(upset_matrix, nsets = 10, 
      #             number.angles = 30, 
      #             point.size = 3.5, 
      #             line.size = 2, 
      #             mainbar.y.label = "Intersections Size", 
      #             sets.x.label = "Number of DEGs", 
      #             order.by = "freq",
      #             sets.bar.color = "#1f77b4",
      #             matrix.color = "#ff7f0e",
      #             text.scale = c(1.8, 1.6, 1.3, 1.0, 1.8, 1.3),
      #             mb.ratio = c(0.6, 0.4)))
      # dev.off()
      
    }
    
    # Reactome <-enrichPathway(gene=gene_set,pvalueCutoff=0.05, readable=T)
    # Reactome<- enrichPathway(gene = gene_KEGG, #表示前景基因，即待富集的基因列表;[,1]表示对entrezid_all数据集的第1列进行处理
    #                      organism='mouse', #物种选择 "human", "rat", "mouse", "celegans", "yeast", "zebrafish", "fly".
    #                      minGSSize = 3, #默认 10
    #                      maxGSSize = 500, #默认 500
    #                      pAdjustMethod = "BH", # 指定多重假设检验矫正的方法，选项包含 "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
    #                      pvalueCutoff = 0.05,  #指定 p 值阈值（指定 1 以输出全部），默认为0.05
    #                      qvalueCutoff = 1, #指定 q 值阈值（指定 1 以输出全部），默认0.2
    #                      readable = T) # 是否将gene ID转换到 gene symbol
    # write.table(Reactome, paste0(outdir,cluster_diff_sub,"_",diff1,"-",diff2, "_enrich.Reactome.csv"), sep = ',', row.names = FALSE, quote = TRUE)
    
    # Reactome@result$Description <- sub(" - Mus musculus \\(house mouse\\)$", "", Reactome@result$Description)
    # if(!is.null(Reactome) && nrow(Reactome) > 0){
    # p <- dotplot(Reactome,  #Reactome富集分析结果
    #         x = "GeneRatio",  #横坐标,默认GeneRation,也可以为Count
    #         color = "p.adjust",  #右纵坐标,默认p.adjust,也可以为pvalue和qvalue
    #         showCategory = 10,  #展示前20个点，默认为10个
    #         size = NULL,  #点的大小
    #         title = "Reactome_dotplot",  #设置图片的标题
    # 				label_format = 100)
    # # p + guides(colour = "Count")
    # ggsave(paste0(outdir,cluster_diff_sub,"_",diff1,"-",diff2, ".Reactome_dotplot.pdf"),width = 12,height = 5)
    # ggsave(paste0(outdir,cluster_diff_sub,"_",diff1,"-",diff2, ".Reactome_dotplot.png"),width = 12,height = 5)
    
    # barplot(Reactome,  #Reactome富集分析结果
    #         x = "Count",  #横坐标,默认Count,也可以为GeneRation
    #         color = "p.adjust",  #右纵坐标,默认p.adjust,也可以为pvalue和qvalue
    #         showCategory = 10,  #展示前20个，默认为10个
    #         size = NULL,  
    #         title = "Reactome_barplot",  #设置图片的标题
    # 				label_format = 80)
    
    # ggsave(paste0(outdir,cluster_diff_sub,"_",diff1,"-",diff2, ".Reactome_barplot.pdf"),width = 12,height = 5)
    # ggsave(paste0(outdir,cluster_diff_sub,"_",diff1,"-",diff2, ".Reactome_barplot.png"),width = 12,height = 5)
    
    # upset_data <- Reactome@result[1:min(10, nrow(Reactome@result)), ]
    # gene_go_long <- upset_data %>% select(Description, geneID) %>% separate_rows(geneID, sep = "/") %>% mutate(present = 1) %>% distinct()  # 去重
    # gene_go_wide <- gene_go_long %>% pivot_wider( id_cols = geneID, names_from = Description, values_from = present, values_fill = 0 )
    # upset_matrix <- as.data.frame(gene_go_wide[, -1])  # 去掉geneID列
    # rownames(upset_matrix) <- gene_go_wide$geneID
    
    # pdf(paste0(outdir,cluster_diff_sub,"_",diff1,"-",diff2,".Reactome_upsetplot.pdf"), width = 20, height = 10)
    # print(upset(upset_matrix, nsets = 10, 
    #       number.angles = 30, 
    #       point.size = 3.5, 
    #       line.size = 2, 
    #       mainbar.y.label = "Intersections Size", 
    #       sets.x.label = "Number of DEGs", 
    # 			order.by = "freq",
    #       sets.bar.color = "#1f77b4",
    #       matrix.color = "#ff7f0e",
    #       text.scale = c(1.8, 1.6, 1.3, 1.0, 1.8, 1.3),
    #       mb.ratio = c(0.6, 0.4)))
    # dev.off()
    # png(paste0(outdir,cluster_diff_sub,"_",diff1,"-",diff2,".Reactome_upsetplot.png"), width = 20, height = 10, units = "in", res = 300)
    # print(upset(upset_matrix, nsets = 10, 
    #       number.angles = 30, 
    #       point.size = 3.5, 
    #       line.size = 2, 
    #       mainbar.y.label = "Intersections Size", 
    #       sets.x.label = "Number of DEGs", 
    # 			order.by = "freq",
    #       sets.bar.color = "#1f77b4",
    #       matrix.color = "#ff7f0e",
    #       text.scale = c(1.8, 1.6, 1.3, 1.0, 1.8, 1.3),
    #       mb.ratio = c(0.6, 0.4)))
    # dev.off()
    
    # 		}
  }
}


# setwd('/mnt/NexData_file/home/zwl/Project/FWSC20253010/require3/')
# input = '/mnt/NexData_file/home/zwl/Project/FWSC20253010/require3/cluster_groupDEG'
# output = '/mnt/NexData_file/home/zwl/Project/FWSC20253010/require3/cluster_groupenrich/'
setwd('/mnt/NexData_file/home/zwl/Project/SC2601060090/require1/T cells/celltype2/groupDEG')
input = '/mnt/NexData_file/home/zwl/Project/SC2601060090/require1/T cells/celltype2/groupDEG/'
output = '/mnt/NexData_file/home/zwl/Project/SC2601060090/require1/T cells/celltype2/groupDEG/'
files <- list.files(input, pattern = "*_differentially_expressed_genes.xls", recursive = TRUE, full.names = TRUE)


for(n in files[15:33]){
  target_dirs <- file.path( basename(dirname(dirname(n))),basename(dirname(n)))
  print(target_dirs)
  setwd(output)
  dir.create(target_dirs,recursive = T)
  outdir = paste0(output,target_dirs,'/')
  
  data = read.table(file = n,sep = '\t',header = T)
  # colnames(data) <- c("gene_name", "PValue", "log2FC", "pct_nz_group" , "pct_nz_reference", "FDR")
  up_gene = data %>% dplyr::filter(pvals_adj <0.05,logfoldchanges>0.25) %>% pull(names)
  if (length(up_gene)>0) {
    clus_enrichment(up_gene,orgdb="org.Mm.eg.db",trend='up',keggdb = "mmu",outdir)
  }
  
  
  down_gene = data %>% dplyr::filter(pvals_adj <0.05,logfoldchanges< -0.25) %>% pull(names)
  if (length(down_gene)>0) {
    clus_enrichment(down_gene,orgdb="org.Mm.eg.db",trend='down',keggdb = "mmu",outdir)
  }
}





