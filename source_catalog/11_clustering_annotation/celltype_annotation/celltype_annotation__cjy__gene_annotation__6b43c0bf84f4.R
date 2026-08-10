# R script to perform gene annotation using local Ensembl database based on R package ensembldb

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(AnnotationHub)
library(ensembldb)

# local Ensembl database used for converting gene ID to gene name and gene annotation
ensdbfile <- snakemake@params[["ensdb"]]

# counts table
df <- read.table(snakemake@input[["counts"]], sep='\t', header=1, check.names=FALSE)

edb <- EnsDb(ensdbfile)
g2g <- AnnotationDbi::select(
  edb,
  keys = df$Gene,
  keytype = "GENEID",
  columns = c(
    "GENEID",
    "GENENAME",
    "SEQNAME",
    "GENESEQSTART",
    "GENESEQEND",
    "SEQSTRAND",
    "GENEBIOTYPE",
    "DESCRIPTION"
  ),
  multiVals = "first"
)

# Since 'multiVals = "first"' parameter doesn't work in AnnotationDbi::select (ONLY works for AnnotationDbi::mapIds)
# Manully make sure 1:1 annotation mapping for each gene
g2g <- g2g[!duplicated(g2g$GENEID),]

colnames(g2g) <- c("ensembl_gene_id", "external_gene_name",
                   "gene_chr", "gene_start", "gene_end", "gene_strand",
                   "gene_biotype", "gene_description")
g2g$gene_strand <- ifelse(g2g$gene_strand == 1, "+", ifelse(g2g$gene_strand == -1, "-", ""))

annotated <- merge(df, g2g, by.x="Gene", by.y="ensembl_gene_id")
annotated$external_gene_name <- ifelse(annotated$external_gene_name == '', annotated$Gene, annotated$external_gene_name)
annotated <- annotated %>%
  dplyr::rename(gene_name = external_gene_name)

write.table(annotated, snakemake@output[["annotation"]], sep='\t', row.names=F)


