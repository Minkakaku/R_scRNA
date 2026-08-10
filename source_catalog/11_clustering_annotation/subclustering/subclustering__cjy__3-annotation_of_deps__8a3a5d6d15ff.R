# log_con <- file(snakemake@log[[1]], "w")
# sink(log_con, type = "message")

# 差异表达蛋白注释脚本
library(data.table)
library(dplyr)

sub.go.list=c("GO:0009507","GO:0005737","GO:0005856","GO:0005783","GO:0005576","GO:0005794","GO:0005739","GO:0016604","GO:0005777","GO:0005886","GO:0005774")

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
bp.diff=merge(subset(data.go,Category=='Biological process'),diff.list,by="Protein accession")
mf.diff=merge(subset(data.go,Category=='Molecular function'),diff.list,by="Protein accession")
cc.diff=merge(subset(data.go,Category=='Cellular component'),diff.list,by="Protein accession")
kegg.diff=merge(data.kegg,diff.list,by="Protein accession")
pfam.diff=merge(data.pfam,diff.list,by="Protein accession")
kog.diff=merge(data.nog,diff.list,by="Protein accession")
subcell.diff=merge(subset(data.go,ID %in% sub.go.list,select=c("Terms","Protein accession")),diff.list,by="Protein accession") %>%
  dplyr::rename(Subcellular_Location="Terms")

# 写出差异蛋白注释结果文件
message("保存注释结果...")
fwrite(bp.diff,file=(snakemake@output[["bp_annotation"]]),sep="\t")
fwrite(cc.diff,file=snakemake@output[["cc_annotation"]],sep="\t")
fwrite(mf.diff,file=snakemake@output[["mf_annotation"]],sep="\t")
fwrite(kegg.diff,file=snakemake@output[["kegg_annotation"]],sep="\t")
fwrite(pfam.diff,file=snakemake@output[["pfam_annotation"]],sep="\t")
fwrite(kog.diff,file=snakemake@output[["kog_annotation"]],sep="\t")
fwrite(subcell.diff,file=snakemake@output[["subcell_annotation"]],sep="\t")

# message("差异表达蛋白注释完成！")
# sink(type = "message")
# close(log_con)