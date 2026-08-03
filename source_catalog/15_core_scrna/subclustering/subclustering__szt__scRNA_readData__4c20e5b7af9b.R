
##==准备seurat对象列表==##
library(Seurat)
packageVersion("Seurat")
library(ggplot2)
script_path <- normalizePath(dirname(rstudioapi::getSourceEditorContext()$path))
setwd(script_path)

#############################################################
#修改barcode,读取样本
dir <- dir("../GMU_scData/")
sample_name <- dir
dir <- paste0("../GMU_scData/", dir)

# 原始数据通过read10X读取
for (i in 1:length(sample_name)) {
  
  #i <-1 
  assign(sample_name[i], Read10X(data.dir = dir[i]))
  tmp <- get(sample_name[i])
  print(sample_name[i])
  print(nchar(tmp@Dimnames[[2]][1])) #length of barcodes
  colnames(tmp) <- paste0(sample_name[i], "_", colnames(tmp))
  assign(sample_name[i], tmp)

  print(dim(get(sample_name[i])))#dim of orignal data
}

rm(tmp)
#############################################################
#创建SeuratObject
for (i in 1:length(sample_name)) {
  #i <- 1
  assign(x = paste0("sdata.", sample_name[i]), 
         value = CreateSeuratObject(get(sample_name[i]), project=sample_name[i], min.cells = 5, min.features = 500))

  #计算线粒体相关基因含量
  assign(x = paste0("sdata.", sample_name[i]),
         value = PercentageFeatureSet(get(paste0("sdata.", sample_name[i])), "^MT-", col.name = "percent_mito"))

  #计算核糖体相关基因含量
  assign(x = paste0("sdata.", sample_name[i]), 
         value = PercentageFeatureSet(get(paste0("sdata.", sample_name[i])), "^RP[SL]", col.name = "percent_ribo"))
  #----------------------------------------------------------
  #数据过滤筛选
  data.filt <- get(paste0("sdata.", sample_name[i]))
  print(dim(data.filt))
  
  selected_mito <- WhichCells(data.filt, expression = percent_mito < 15)
  data.filt <- subset(data.filt, cells = selected_mito)
  
  selected_ribo <- WhichCells(data.filt, expression = percent_ribo > 5)
  data.filt <- subset(data.filt, cells = selected_ribo)
  
  selected_RNA <- WhichCells(data.filt, expression = nCount_RNA < 20000)
  data.filt <- subset(data.filt, cells = selected_RNA)
  
  print(dim(data.filt))
  gc()
  #--------------------------------------------------------
  #过滤线粒体基因
  data.filt <- data.filt[!grepl("^MT-", rownames(data.filt)), ]
  
  #过滤核糖体基因
  data.filt <- data.filt[!grepl('^RP[SL]', rownames(data.filt)), ]
  
  #过滤血红蛋白基因
  data.filt <- data.filt[!grepl("^HB[^(P)]", rownames(data.filt)), ]
  
  print(dim(data.filt))
  saveRDS(data.filt, file = paste0("sdata.", sample_name[i], ".rds"))
}
#############################################################
#删除不用的数据对象
for (i in 1:length(sample_name)) {
  rm(list = sample_name[i])
  rm(list = paste0("sdata.", sample_name[i]))
}
#run garbage collect to free up memory
gc()

