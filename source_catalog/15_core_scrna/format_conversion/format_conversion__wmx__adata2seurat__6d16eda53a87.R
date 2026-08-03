#!/usr/bin/env Rscript
# options(encoding = "UTF-8")
# conda activate /mnt/NexData_file/share/scRNA-tools/envs/pip_SCENICR/
# debug:zhangwenlu
library(argparser)
library(reticulate)
library(sceasy)
use_condaenv('/mnt/NexData_file/share/scRNA-tools/envs/pip_SCENICR/bin/python3.8')

argv <- arg_parser('')
argv <- add_argument(argv,"--h5ad", help="input h5ad")
argv <- add_argument(argv,"--rds", help="output rds file")
argv <- parse_args(argv)

h5ad <- argv$h5ad
rds <- argv$rds
sceasy::convertFormat(h5ad,from='anndata',to='seurat',main_layer = "counts",outFile=rds)
