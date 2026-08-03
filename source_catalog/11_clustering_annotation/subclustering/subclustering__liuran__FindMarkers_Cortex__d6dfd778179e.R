setwd("/mnt/NexData_file/home/liuran/projects/3_Xenium/20260115_XE2511036190-01_mouse_brain")
# 读取数据
all_data <- readRDS(file.path("data_for_analysis", "all_data_ROI_annotation.rds"))


# FindMarkers 海马差异分析 --------
library(tidyverse)
library(Seurat)
library(patchwork)
library(stringr)
library('clusterProfiler')
library("org.Mm.eg.db")
library(RColorBrewer)
library(EnhancedVolcano)
library(dplyr) 

all_data_Cortex <- all_data[, all_data@meta.data$ROI == "Cortex"]

unique(all_data_Cortex$annotation)
unique(all_data_Cortex$group)

table(all_data_Cortex$annotation, all_data_Cortex$group)

# Date <- Sys.Date()
# deg_outdir <- str_glue("FindMarkers-{Date}")
# dir.create(deg_outdir)


All_type <- as.character(unique(all_data_Cortex$annotation))
print(All_type)
length(All_type)


All_type %>% 
  setNames(.,rep('all_data_Cortex', 9)) %>% 
  imap(function(celltype,name){ # browser()
    tryCatch({
      object <- get(name)
      Idents(object)<-'annotation'
      # print(celltype)
      # print(name)
      # return()
      # 验证亚型隔离
      cat("Processing:", celltype, "| Cells:", sum(Idents(object) == celltype), "\n")
      
      # 使用子集对象提高效率
      subset_obj <- subset(JoinLayers(object), idents = celltype)
      
      diffmarkers <-
        FindMarkers(
          subset_obj,
          ident.1 = 'CTE',
          ident.2 = 'NC',
          test.use = "wilcox",
          logfc.threshold = 0.25,
          min.pct = 0.1,
          group.by = 'group'
        )
      # ==== 新增：检查空结果 ====
      if (nrow(diffmarkers) == 0) {
        cat("No DEGs found for", celltype, ". Skipping enrichment analysis.\n")
        write.csv(data.frame(), str_glue('ROI_analysis/FindMarkers_Cortex/{name}/{celltype}/NO_DEGs.csv'))
        return(NULL)
      }
      
      diffmarkers_up <- diffmarkers %>% filter(avg_log2FC > 0) %>% arrange(desc(avg_log2FC))
      diffmarkers_down <- diffmarkers %>% filter(avg_log2FC < 0) %>% arrange(desc(avg_log2FC))
      
      dir.create(file.path(str_glue('ROI_analysis/FindMarkers_Cortex/{name}/{celltype}/')), recursive = TRUE)
      write.csv(diffmarkers,file.path(str_glue('ROI_analysis/FindMarkers_Cortex/{name}/{celltype}/{celltype} group diffmarker.csv')))
      write.csv(diffmarkers_up,file.path(str_glue('ROI_analysis/FindMarkers_Cortex/{name}/{celltype}/{celltype} group diffmarker up.csv')))
      write.csv(diffmarkers_down,file.path(str_glue('ROI_analysis/FindMarkers_Cortex/{name}/{celltype}/{celltype} group diffmarker down.csv')))
      
      # GO up
      enrichGO_up_ls <- enrichGO(
        rownames(diffmarkers_up),
        OrgDb = org.Mm.eg.db, 
        ont='ALL',
        pAdjustMethod = 'BH',
        pvalueCutoff = 0.05, 
        qvalueCutoff = 0.05,
        keyType = 'SYMBOL')
      write.csv(enrichGO_up_ls@result,file = file.path(str_glue('ROI_analysis/FindMarkers_Cortex/{name}/{celltype}/{celltype} Go up enriment table.csv')),row.names = FALSE)
      
      enrichGO_up_dot <- enrichplot::dotplot(
        enrichGO_up_ls,
        showCategory = 10,
        label_format = 50)
      ggsave(
        enrichGO_up_dot,
        filename = file.path(str_glue('ROI_analysis/FindMarkers_Cortex/{name}/{celltype}/{celltype} Go up dotplot.pdf')),
        height = 6,
        width = 8)
      enrichGO_up_bar <- barplot(
        enrichGO_up_ls,
        showCategory = 10,
        label_format = 50)
      ggsave(
        enrichGO_up_bar,
        filename = file.path(str_glue('ROI_analysis/FindMarkers_Cortex/{name}/{celltype}/{celltype} Go up barplot.pdf')),
        height = 6,
        width = 8)
      
      
      # GO down
      enrichGO_down_ls <- enrichGO(
        rownames(diffmarkers_down),
        OrgDb = org.Mm.eg.db, 
        ont='ALL',
        pAdjustMethod = 'BH',
        pvalueCutoff = 0.05, 
        qvalueCutoff = 0.05,
        keyType = 'SYMBOL')
      write.csv(enrichGO_down_ls@result,file = file.path(str_glue('ROI_analysis/FindMarkers_Cortex/{name}/{celltype}/{celltype} Go down enriment table.csv')),row.names = FALSE)
      
      enrichGO_down_dot <- enrichplot::dotplot(
        enrichGO_down_ls,
        showCategory = 10,
        label_format = 50)
      ggsave(
        enrichGO_down_dot,
        filename = file.path(str_glue('ROI_analysis/FindMarkers_Cortex/{name}/{celltype}/{celltype} Go down dotplot.pdf')),
        height = 6,
        width = 8)
      enrichGO_down_bar <- barplot(
        enrichGO_down_ls,
        showCategory = 10,
        label_format = 50)
      ggsave(
        enrichGO_down_bar,
        filename = file.path(str_glue('ROI_analysis/FindMarkers_Cortex/{name}/{celltype}/{celltype} Go down barplot.pdf')),
        height = 6,
        width = 8)
      
      # KEGG up
      de_up_genes_id_for_enrichmenet <- 
        clusterProfiler::bitr(
          rownames(diffmarkers_up),
          OrgDb = "org.Mm.eg.db",
          fromType = "SYMBOL",
          toType = "ENTREZID"
        ) %>%
        .$ENTREZID
      enrichKEGG_up_ls <- enrichKEGG(
        de_up_genes_id_for_enrichmenet,
        organism = "mmu", 
        keyType = "ncbi-geneid",
        # keyType = "kegg", 
        pvalueCutoff = 0.05, 
        pAdjustMethod = "BH", 
        minGSSize = 10,
        maxGSSize = 500, 
        qvalueCutoff = 0.2,
        use_internal_data = FALSE
        # use_internal_data = TRUE
      )
      enrichKEGG_up_ls <- setReadable(enrichKEGG_up_ls, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
      write.csv(enrichKEGG_up_ls@result,file = file.path(str_glue('ROI_analysis/FindMarkers_Cortex/{name}/{celltype}/{celltype} KEGG up enriment table.csv')),row.names = FALSE)
      
      enrichKEGG_up_dot <- enrichplot::dotplot(
        enrichKEGG_up_ls,
        showCategory = 10,
        label_format = 50)
      enrichKEGG_up_bar <- barplot(
        enrichKEGG_up_ls,
        showCategory = 10,
        label_format = 50)
      ggsave(
        enrichKEGG_up_dot,
        filename = file.path(str_glue('ROI_analysis/FindMarkers_Cortex/{name}/{celltype}/{celltype} KEGG up dotplot.pdf')),
        height = 6,
        width = 8)
      ggsave(
        enrichKEGG_up_bar,
        filename = file.path(str_glue('ROI_analysis/FindMarkers_Cortex/{name}/{celltype}/{celltype} KEGG up barplot.pdf')),
        height = 6,
        width = 8)
      
      
      # KEGG down
      de_down_genes_id_for_enrichmenet <- 
        clusterProfiler::bitr(
          rownames(diffmarkers_down),
          OrgDb = "org.Mm.eg.db",
          fromType = "SYMBOL",
          toType = "ENTREZID"
        ) %>%
        .$ENTREZID
      enrichKEGG_down_ls <- enrichKEGG(
        de_down_genes_id_for_enrichmenet,
        organism = "mmu", 
        keyType = "ncbi-geneid",
        # keyType = "kegg", 
        pvalueCutoff = 0.05, 
        pAdjustMethod = "BH", 
        minGSSize = 10,
        maxGSSize = 500, 
        qvalueCutoff = 0.2,
        use_internal_data = FALSE
        # use_internal_data = TRUE
      )
      enrichKEGG_down_ls <- setReadable(enrichKEGG_down_ls, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
      write.csv(enrichKEGG_down_ls@result,file = file.path(str_glue('ROI_analysis/FindMarkers_Cortex/{name}/{celltype}/{celltype} KEGG down enriment table.csv')),row.names = FALSE)
      
      enrichKEGG_down_dot <- enrichplot::dotplot(
        enrichKEGG_down_ls,
        showCategory = 10,
        label_format = 50)
      enrichKEGG_down_bar <- barplot(
        enrichKEGG_down_ls,
        showCategory = 10,
        label_format = 50)
      ggsave(
        enrichKEGG_down_dot,
        filename = file.path(str_glue('ROI_analysis/FindMarkers_Cortex/{name}/{celltype}/{celltype} KEGG down dotplot.pdf')),
        height = 6,
        width = 8)
      ggsave(
        enrichKEGG_down_bar,
        filename = file.path(str_glue('ROI_analysis/FindMarkers_Cortex/{name}/{celltype}/{celltype} KEGG down barplot.pdf')),
        height = 6,
        width = 8)
      
      
      top <- rbind(diffmarkers %>% top_n(n = 25, wt = avg_log2FC),diffmarkers %>% top_n(n = -25, wt = avg_log2FC))
      p<-DoHeatmap(object, label = F,
                   features = rownames(top),
                   group.by = "group",
                   disp.max= 1,
                   disp.min= -1,
                   assay = 'RNA')+ #设置组别颜色
        scale_fill_gradientn(colors = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100))
      pdf(file.path(str_glue("ROI_analysis/FindMarkers_Cortex/{name}/{celltype}/diffgene_heatmap.pdf")), height = 8, width =10)
      print(p)
      dev.off()
      
      p_volcano <- EnhancedVolcano(diffmarkers,
                                   lab = rownames(diffmarkers),
                                   x = 'avg_log2FC',
                                   y = 'p_val_adj',
                                   title = str_glue('{celltype} volcano plot'),
                                   pCutoff = 10e-2,
                                   FCcutoff = 0,
                                   pointSize = 3.0,
                                   labSize = 6.0)
      ggsave(p_volcano,
             filename = file.path(str_glue("ROI_analysis/FindMarkers_Cortex/{name}/{celltype}/diffgene valcano.pdf")),
             height = 6,
             width = 8)
    }, error = function(e) {
      # ===== 错误处理代码 ===== #
      message("\nERROR in ", celltype, ": ", e$message)
      # 创建专门存放错误对象的目录
      error_dir <- "ROI_analysis/FindMarkers_Cortex/all_data_error_objects/"
      dir.create(error_dir, showWarnings = FALSE, recursive = TRUE)
      # 保存引发错误的数据对象
      saveRDS(
        list(celltype = celltype, 
             error_msg = e$message, 
             subset_obj = if(exists("subset_obj")) subset_obj else NULL),
        file.path(error_dir, paste0(celltype, "_error.rds"))
      )
      return(NULL)
    })
  })


