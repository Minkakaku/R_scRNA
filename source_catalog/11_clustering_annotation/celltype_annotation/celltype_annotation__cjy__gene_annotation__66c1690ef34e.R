# R script to perform gene annotation using local Ensembl database based on R package ensembldb

if (length(snakemake@log) > 0) {
    log <- file(snakemake@log[[1]], open="wt")
    sink(log)
    sink(log, type="message")
}

library(rtracklayer)
library(dplyr)

# longestgff <- rtracklayer::import('longest.gff3')
# longestgff <- as.data.frame(longestgff)
# longestgff.use <- longestgff %>% filter(type == "CDS") %>%
# dplyr::select(c("protein_id"))
# longest.proteinid <- unique(longestgff.use$protein_id)


gff <- rtracklayer::import(snakemake@params[["gff3"]])
gff <- as.data.frame(gff)
gff.use <- gff %>% filter(type == "gene") %>%
dplyr::select(c("gene_id","seqnames","start","end","strand","Name","biotype","description")) %>% 
arrange(gene_id) %>%
dplyr::rename(Gene = gene_id) %>%
dplyr::rename(chromesome = seqnames) %>%
dplyr::rename(symbol_gtf = Name) %>%
dplyr::rename(description_gtf = description)
# write.table(gff.use,"gff3_annotation.xls",sep = "\t", quote = FALSE, row.names = FALSE)

# library(AnnotationHub)
# library(ensembldb)

# # local Ensembl database used for converting gene ID to gene name and gene annotation
# ensdbfile <- snakemake@params[["ensdb"]]

# counts table
df <- read.table(snakemake@input[["counts"]], sep='\t', header=1, check.names=FALSE)

df_annotation1 = merge(df, gff.use, by = "Gene", all.x = TRUE)

eggnog_anno <- read.table(snakemake@input[["eggnog_mapper_annotation"]], sep='\t', header=1, check.names=FALSE, quote = "")
eggnog_anno.use <- eggnog_anno %>%
dplyr::rename(Gene = gene_id) %>%
dplyr::rename(seed_ortholog_eggnog = seed_ortholog) %>%
dplyr::rename(symbol_eggnog = Preferred_name) %>%
dplyr::rename(description_eggnog = Description) %>%
dplyr::rename(COG_category_eggnog = COG_category) %>%
dplyr::rename(GOs_eggnog = GOs) %>%
dplyr::rename(EC_eggnog = EC) %>%
dplyr::rename(KEGG_ko_eggnog = KEGG_ko) %>%
dplyr::rename(KEGG_Pathway_eggnog = KEGG_Pathway) %>%
dplyr::rename(PFAMs_eggnog = PFAMs)
# edb <- EnsDb(ensdbfile)
# g2g <- AnnotationDbi::select(
#   edb,
#   keys = df$Gene,
#   keytype = "GENEID",
#   columns = c(
#     "GENEID",
#     "GENENAME",
#     "SEQNAME",
#     "GENESEQSTART",
#     "GENESEQEND",
#     "SEQSTRAND",
#     "GENEBIOTYPE",
#     "DESCRIPTION"
#   ),
#   multiVals = "first"
# )

df_annotation2 = merge(df_annotation1, eggnog_anno.use, by = "Gene", all.x = TRUE)

# # Since 'multiVals = "first"' parameter doesn't work in AnnotationDbi::select (ONLY works for AnnotationDbi::mapIds)
# # Manully make sure 1:1 annotation mapping for each gene
# g2g <- g2g[!duplicated(g2g$GENEID),]

# colnames(g2g) <- c("ensembl_gene_id", "external_gene_name",
#                    "gene_chr", "gene_start", "gene_end", "gene_strand",
#                    "gene_biotype", "gene_description")
# g2g$gene_strand <- ifelse(g2g$gene_strand == 1, "+", ifelse(g2g$gene_strand == -1, "-", ""))

# annotated <- merge(df, g2g, by.x="Gene", by.y="ensembl_gene_id")
# annotated$external_gene_name <- ifelse(annotated$external_gene_name == '', annotated$Gene, annotated$external_gene_name)
# annotated <- annotated %>%
#   dplyr::rename(gene_name = external_gene_name)

write.table(df_annotation2, snakemake@output[["annotation"]], sep='\t', row.names = FALSE, quote = FALSE)

