library(data.table)
library(dplyr)

# 读取注释文件
message("读取蛋白注释文件...")
data.go=fread(paste(snakemake@params[["db_dir"]],'GO_annotation.tsv',sep="/"))
data.kegg=fread(paste(snakemake@params[["db_dir"]],'KEGG_annotation.tsv',sep="/"))
data.pfam=fread(paste(snakemake@params[["db_dir"]],'Pfam_annotation.tsv',sep="/"))
data.nog=fread(paste(snakemake@params[["db_dir"]],'NOG_annotation.tsv',sep="/"))

# 读取差异表达蛋白数据
message("读取差异表达蛋白数据...")
diff.data <- fread(snakemake@input[["dep_results"]])
diff.list=subset(diff.data, Regulation!='unchanged')
bp.diff=merge(subset(data.go,Category=='Biological process'),diff.list,by="Protein accession",allow.cartesian=TRUE)
mf.diff=merge(subset(data.go,Category=='Molecular function'),diff.list,by="Protein accession",allow.cartesian=TRUE)
cc.diff=merge(subset(data.go,Category=='Cellular component'),diff.list,by="Protein accession",allow.cartesian=TRUE)
kegg.diff=merge(data.kegg,diff.list,by="Protein accession",allow.cartesian=TRUE)
pfam.diff=merge(data.pfam,diff.list,by="Protein accession",allow.cartesian=TRUE)
kog.diff=merge(data.nog,diff.list,by="Protein accession",allow.cartesian=TRUE)

# 写出差异蛋白注释结果文件
message("保存注释结果...")
fwrite(bp.diff,file=(snakemake@output[["bp_annotation"]]),sep="\t")
fwrite(cc.diff,file=snakemake@output[["cc_annotation"]],sep="\t")
fwrite(mf.diff,file=snakemake@output[["mf_annotation"]],sep="\t")
fwrite(kegg.diff,file=snakemake@output[["kegg_annotation"]],sep="\t")
fwrite(pfam.diff,file=snakemake@output[["pfam_annotation"]],sep="\t")
fwrite(kog.diff,file=snakemake@output[["kog_annotation"]],sep="\t")
