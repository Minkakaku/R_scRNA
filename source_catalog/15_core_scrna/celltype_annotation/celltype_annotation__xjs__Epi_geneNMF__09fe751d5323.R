library(Seurat)
library(tidyverse)
library(GeneNMF) # https://carmonalab.r-universe.dev/GeneNMF/doc/manual.html
library(msigdbr)
library(fgsea)
library(AUCell)

out_dir <- 'Epi_geneNMF'
dir.create(out_dir)

Epithelial <- readRDS("~/CN_analysis_epi/all_data_CN_epis.rds")
Epithelial$annotation %>% table()
Epithelial$CN2 %>% table()


# runNMF -----------
ndim <- 15
Epithelial <- FindVariableFeatures(Epithelial) # nfeatures = length(rownames(Epithelial))
Epithelial <- runNMF(JoinLayers(Epithelial), k = ndim)
Epithelial@reductions$NMF


Epithelial <- RunUMAP(Epithelial, 
                      reduction = "NMF", 
                      dims=1:ndim, 
                      reduction.name = "NMF_UMAP", 
                      reduction.key = "nmfUMAP_")
DimPlot(Epithelial, reduction = "NMF_UMAP", 
        group.by = "sample", label=T) + 
  theme(aspect.ratio = 1,axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank()) +
  ggtitle("NMF UMAP") + NoLegend()
ggsave(file.path(out_dir, 'Dimplot_NMF_UMAP_sample.pdf'), width = 6, height = 6)




# find MPs ------------
seu.list <- SplitObject(Epithelial, split.by = "sample")
geneNMF.programs <- multiNMF(seu.list, 
                             assay="RNA", 
                             slot="data", 
                             k=4:6, # 4:9
                             nfeatures = 2000) # 官方是1000 length(rownames(Epithelial))

# 可以将多个样本和数量k中识别出的基因程序组合成元程序（MPs），即在不同 NMF 运行中稳健识别的共识程序。在此，我们将定义10个 MPs
geneNMF.metaprograms <- getMetaPrograms(geneNMF.programs,
                                        nMP=10, # 10
                                        weight.explained = 0.7,
                                        # min.confidence = 0.7,
                                        max.genes=100)


plotMetaPrograms(geneNMF.metaprograms)
geneNMF.metaprograms$metaprograms.metrics
# sampleCoverage(MP 在多少样本中被检测到）；silhouette(一个MP中的个体程序相对于其他MP中的程序有多相似——数值越高越好）
# sampleCoverage  silhouette meanSimilarity numberGenes numberPrograms
# MP1       1.0000000 -0.07623899          0.042           7            683
# MP2       0.7894737  0.20352527          0.319         100            171
# MP3       0.5789474  0.22618637          0.324          75            119
# MP4       0.5263158  0.76662423          0.795          16             89
# MP5       0.4473684  0.83709234          0.854         100             93
# MP6       0.3947368  0.28203423          0.366         100             83
# MP7       0.3947368  0.23379610          0.286          37             99
# MP8       0.3684211  0.18437231          0.304          83             68
# MP9       0.2368421  0.56715948          0.617          24             44
# MP10      0.2105263  0.74063832          0.776          41             33

# sampleCoverage  silhouette meanSimilarity numberGenes numberPrograms
# MP1      1.0000000 -0.03962427          0.061           2           1003
# MP2      0.7894737  0.77766671          0.802           2            135
# MP3      0.7368421  0.32923025          0.437          33            155
# MP4      0.4736842  0.79222722          0.816           2            106
# MP5      0.3947368  0.43987872          0.542          58             83

# sampleCoverage silhouette meanSimilarity numberGenes numberPrograms
# MP1      0.7894737  0.7430198          0.802           2            135
# MP2      0.7368421  0.3274550          0.437          33            155
# MP3      0.6842105  0.1678602          0.285          24            168
# MP4      0.6578947  0.1470367          0.270          98            116
# MP5      0.4736842  0.7919742          0.816           2            106
# MP6      0.4736842  0.3850460          0.463          12             74
# MP7      0.3947368  0.7131272          0.769           3             39
# MP8      0.3947368  0.4209957          0.542          58             83
# MP9      0.2368421  0.7858019          0.812          10             35

geneNMF_metaprograms_filtered <- 
  dropMetaPrograms(
    geneNMF.metaprograms, 
    dropMP = c("MP1")
  )

ph <- plotMetaPrograms(geneNMF_metaprograms_filtered)
ggsave(file.path(out_dir, 'Heatmap_MPs_10.pdf'), ph, width = 7.5, height = 6)




# 查看MP中的驱动基因 ------------
lapply(geneNMF_metaprograms_filtered$metaprograms.genes, head)
saveRDS(geneNMF_metaprograms_filtered$metaprograms.genes, file.path(out_dir, 'MP_genes.rds'))




# 对基因programs进行GSEA分析 -----------
runGSEA <- function (genes, universe = NULL, category = "H", subcategory = NULL, 
                     species = "Homo sapiens", custom.db = NULL, pval.thr = 0.05) 
{
  if (!requireNamespace("fgsea", quietly = TRUE) | !requireNamespace("msigdbr", 
                                                                     quietly = TRUE)) {
    stop("Function 'runGSEA' requires the 'fgsea' and 'msigdbr' packages.\n            Please install them.", 
         call. = FALSE)
  }
  if (any(duplicated(genes))) {
    genes <- genes[!duplicated(genes)]
  }
  if (!is.null(custom.db)) {
    if (!is.list(custom.db) || is.null(names(custom.db))) 
      stop("custom.db should be a named list", call. = FALSE)
    DB_list <- custom.db
  }
  else {
    msig_df <- msigdbr::msigdbr(species = species, category = category, 
                                subcategory = subcategory)
    DB_list <- split(x = msig_df$gene_symbol, f = msig_df$gs_name)
  }
  fgRes <- fgsea::fora(pathways = DB_list, genes = genes, universe = universe)
  fgRes <- fgRes[fgRes$pval <= pval.thr, ]
  return(fgRes)
}


MP_HALLMARK <- lapply(geneNMF_metaprograms_filtered$metaprograms.genes, function(program) {
  runGSEA(program, universe=rownames(Epithelial), species = "Homo sapiens", category = "H", pval.thr = 1.1 )
})
saveRDS(MP_HALLMARK, file.path(out_dir, 'MP_HALLMARK.rds'))

MP_KEGG <- lapply(geneNMF_metaprograms_filtered$metaprograms.genes, function(program) {
  runGSEA(program, universe=rownames(Epithelial), species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG", pval.thr = 1.1 )
})
saveRDS(MP_KEGG, file.path(out_dir, 'MP_KEGG.rds'))


MP_GOBP <- lapply(geneNMF_metaprograms_filtered$metaprograms.genes, function(program) {
  runGSEA(program, universe=rownames(Epithelial), species = "Homo sapiens", category = "C5", subcategory = "GO:BP", pval.thr = 1.1 )
})
saveRDS(MP_GOBP, file.path(out_dir, 'MP_GOBP.rds'))

MP_GOCC <- lapply(geneNMF_metaprograms_filtered$metaprograms.genes, function(program) {
  runGSEA(program, universe=rownames(Epithelial), species = "Homo sapiens", category = "C5", subcategory = "GO:CC", pval.thr = 1.1 )
})
saveRDS(MP_GOCC, file.path(out_dir, 'MP_GOCC.rds'))

MP_GOMF <- lapply(geneNMF_metaprograms_filtered$metaprograms.genes, function(program) {
  runGSEA(program, universe=rownames(Epithelial), species = "Homo sapiens", category = "C5", subcategory = "GO:MF", pval.thr = 1.1 )
})
saveRDS(MP_GOMF, file.path(out_dir, 'MP_GOMF.rds'))

MP_GO <- 
  MP_GOBP %>% 
  imap(function(enrich, MP){
    bind_rows(enrich, MP_GOCC[[MP]], MP_GOMF[[MP]])
  })

MP_GO  
  
  


## plot ----------
MP_enrich_ls <- list(MP_HALLMARK, MP_KEGG, MP_GO)
names(MP_enrich_ls) <- c('enrichment_HALLMARK', 'enrichment_KEGG', 'enrichment_GO')


p_enrich_dotplot_ls <- 
  MP_enrich_ls %>% 
  imap(function(MP_enrich, enrichment_db){
    MP_enrich %>% 
      imap(function(enrich, name) {
        # 计算 GeneRatio 和 -log10(padj)
        enrich <- enrich %>% 
          filter(padj < 0.05) %>% 
          mutate(
            GeneRatio = overlap / size,
            log10p = -log10(padj)
          ) %>%
          # 按 p.adjust 排序取前10条（最小 padj 的 Top10）
          arrange(padj) %>%
          slice_head(n = 10) %>%
          # 让 pathway 按 GeneRatio 升序排列（显示时从上到下递增，可通过 scale_y_discrete 调整）
          mutate(pathway = factor(pathway, levels = pathway[order(GeneRatio)]))
        
        ggplot(enrich, aes(x = GeneRatio, y = pathway)) +
          geom_point(aes(size = GeneRatio, color = log10p)) +
          scale_color_gradient(low = "blue", high = "red", name = "-log10(p.adjust)") +
          scale_size_continuous(range = c(2, 8), name = "GeneRatio") +
          labs(x = "GeneRatio", y = NULL, title = name) +
          theme_bw() +
          theme(
            axis.text.y = element_text(size = 8),
            plot.title = element_text(hjust = 0.5)
          )
      })
  })

p_enrich_dotplot_ls %>% 
  iwalk(function(p_enrich_dotplot, enrichtype){
    pdf(file.path(out_dir, str_glue('MP enrichment analysis dotplot - {enrichtype}.pdf')), width = 10, height = 6)
    p_enrich_dotplot %>% walk(print)
    dev.off()
  })







# score ----------------
mp.genes <- geneNMF_metaprograms_filtered$metaprograms.genes
Epithelial <- AddModuleScore(Epithelial, 
                            features = mp.genes, 
                            assay="RNA", 
                            ncores=4, 
                            name = 'tmp')
Epithelial@meta.data %>% head
Epithelial$MP2 <- NULL
Epithelial$MP3 <- NULL
Epithelial$MP4 <- NULL
Epithelial$MP5 <- NULL
Epithelial[[]] <- Epithelial[[]] %>%
  dplyr::rename_with(.cols = starts_with("tmp"), 
                     .fn = ~ set_names(names(mp.genes)))
Epithelial@meta.data %>% head

VlnPlot(Epithelial, features=names(mp.genes), 
        group.by = "CN2",
        pt.size = 0, ncol=5) & 
  theme(axis.text.x = element_text(angle = 90))
ggsave(file.path(out_dir, 'Epithelial_CN_MP_10.pdf'), width = 15, height = 10)

saveRDS(Epithelial, file.path(out_dir, 'Epithelial.rds'))
