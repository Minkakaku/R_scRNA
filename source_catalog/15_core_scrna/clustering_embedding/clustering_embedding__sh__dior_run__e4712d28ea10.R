.libPaths("/home/wzb/R/x86_64-pc-linux-gnu-library/4.2/")
library(dior)
library(Seurat)
library(getopt)
# library(loupeR)
command = matrix(c('hf','r',1,"character",
                   'outpath','o',1,"character",
                   'filename','f',1,"character",
                   'select','s',2,"character"
                   )
                 ,byrow=TRUE, ncol=4)
args = getopt(command)

# args[1]h5文件输入
# args[2]输出文件夹
# args[3]rds文件输出名
print(args$select)
print(args$hf)
print(c(args$hf,args$outpath,args$filename))
a = dior::read_h5(args$hf)
if(is.null(args$select)){
  select = null}else{
  select = args$select
  }
if(is.null(select)){
  if("cell_type" %in% colnames(a@meta.data)){
    Idents(a) = a$cell_type
  }else if ("leiden" %in% colnames(a@meta.data)) {
    Idents(a) = a$leiden
  }else if ("louvain" %in% colnames(a@meta.data)) {
    Idents(a) = a$louvain
  }
}else{  
  if("cell_type" == select){
    Idents(a) = a$cell_type
  }else if ("leiden" == select) {
    Idents(a) = a$leiden
  }else if ("louvain" == select) {
    Idents(a) = a$louvain
  }}

# a = NormalizeData(a) 
a = FindVariableFeatures(a,assay = 'RNA')
setwd(args$outpath)
saveRDS(a, args$filename)
# create_loupe_from_seurat(a,output_dir =args$outpath)