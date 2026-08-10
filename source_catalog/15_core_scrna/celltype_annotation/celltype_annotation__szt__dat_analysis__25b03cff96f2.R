library(tidyverse)
library(Seurat)

all_data <- readRDS("data/adata_CN15.rds")

all_data$project <- all_data$sample
all_data$annotation <- all_data$celltype

all_data <- all_data %>% NormalizeData()
all_data <- all_data %>% RenameCells(add.cell.id = 'cell')
all_data$cell_id <- colnames(all_data)


out_dir <- "de and enrichment analysis"
dir.create(out_dir)


# 1 major group de analysis ---------
Idents(all_data_cn15) <- "group"
all_marker_grp <- 
  FindAllMarkers(
    all_data_cn15,
    return.thresh = 1.1
  )

#  export
all_marker_grp %>%
  openxlsx::write.xlsx(
    file = file.path(out_dir, "1.2 de analysis table.xlsx"),
    asTable = TRUE
  )


# 2 enrichment analysis -----------

all_geneid_marker_grp_ls <-
    all_marker_grp %>%
    dplyr::filter(avg_log2FC > 0.25, p_val_adj < 0.05, pct.1 > 0.1) %>%
    group_by(cluster) %>%
    nest %>%
    deframe() %>%
    map(~{
      clusterProfiler::bitr(
        .x$gene,
        fromType = "SYMBOL",
        toType   = "ENTREZID",
        OrgDb    = "org.Hs.eg.db"
      ) %>%
        .$ENTREZID
    })


all_enrich_grp_ls <- 
  all_geneid_marker_grp_ls %>%
  map(function(x){
    clusterProfiler::enrichGO(
      gene = x,
      ont = "BP",
      OrgDb = "org.Hs.eg.db",
      pvalueCutoff  = 1,  
      pAdjustMethod = "BH",
      qvalueCutoff  = 1  
    )
  })

# enrichment 
all_enrich_grp_ls <-
  all_enrich_grp_ls %>%
  map(function(x){ # browser()
    xx <- x@result$geneID %>% str_split("/")
    
    all_genes <- unlist(xx)
    gene_db <- 
      clusterProfiler::bitr(
        all_genes,
        fromType = "ENTREZID",
        toType   = "SYMBOL",
        OrgDb    = "org.Hs.eg.db"
      ) 
    
    xx2 <-
      xx %>%
      map_chr(function(g){
        str_c(gene_db$SYMBOL[match(g, gene_db$ENTREZID)], collapse = ",")
      })
    x@result$geneSymbol <- xx2
    x
  })


# enrichment plot
p_all_enrich_grp_ls <-
  all_enrich_grp_ls %>%
  imap(function(enrich, name){ # browser()
    enrichplot::dotplot(
      enrich,
      x = "GeneRatio",
      color = "p.adjust",
      showCategory = 10,
      label_format = 40
    ) + 
      # scale_fill_gradient2(low = "gray", high = "red") +
      ggtitle(name)
  })


#  export
all_enrich_grp_ls %>%
  map('result') %>%
  openxlsx::write.xlsx(
    file = file.path(out_dir, "1.2 enrichment analysis table.xlsx"),
    asTable = TRUE
  )
  


pdf(
  file.path(out_dir, "1.2 de enrichment analysis dotplot.pdf"), 
  height = 6.5, 
  width  = 8
)
p_all_enrich_grp_ls %>%  walk(print)
dev.off()
