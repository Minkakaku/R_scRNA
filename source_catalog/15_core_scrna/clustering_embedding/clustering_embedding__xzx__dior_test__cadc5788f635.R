library(reticulate)
use_condaenv("/mnt/NexData_file/home/xzx/envs/scDIOR/",required = T)

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
#a = FindVariableFeatures(a,assay = 'RNA')
#a@assays[["counts"]] <- NULL
setwd(args$outpath)
saveRDS(a, args$filename)
# create_loupe_from_seurat(a,output_dir =args$outpath)


################################################################################
setwd("/mnt/NexData_file/home/xzx/project/SC2508154388/report/")
sc <- dior::read_h5("./anno_fin.h5")
min(sc@assays[["counts"]]@counts@x)
max(sc@assays[["counts"]]@counts@x)


# 1. 提取目标向量（非零元素）
counts_x <- sc@assays[["RNA"]]@counts@x

# 2. 查看向量的存储类型（判断是否原生整数型）
cat("向量存储类型：", class(counts_x), "\n")  # 输出 integer 或 numeric

# 3. 核心验证：所有元素是否为整数（数值上无小数部分）
# 方法1：判断元素与自身取整后是否相等（最直接）
all_integers <- all(counts_x == as.integer(counts_x))
# 方法2：检查小数部分是否全为0（避免浮点数精度问题，推荐）
all_integers2 <- all(floor(counts_x) == counts_x)  # floor() 向下取整



sc <- readRDS("./h5_trans.rds")
GetAssayData(sc, assay = 'RNA',slot = 'counts')
