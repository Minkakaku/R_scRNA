library(dplyr)
library(getopt)
library(preprocessCore)#装，安装方法写在注释了，bioconductor安装很慢很慢
library(parallel)
library(e1071)#要装
library(tidyr)#装
library(tidyverse)#要装
CoreAlg <- function(X, y){
  #try different values of nu
  svn_itor <- 3
  
  res <- function(i){
    if(i==1){nus <- 0.25}
    if(i==2){nus <- 0.5}
    if(i==3){nus <- 0.75}
    model<-e1071::svm(X,y,type="nu-regression",kernel="linear",nu=nus,scale=F)
    model
  }
  
  if(Sys.info()['sysname'] == 'Windows') out <- parallel::mclapply(1:svn_itor, res, mc.cores=1) else out <- parallel::mclapply(1:svn_itor, res, mc.cores=3)
  
  nusvm <- rep(0,svn_itor)
  corrv <- rep(0,svn_itor)
  
  #do cibersort
  t <- 1
  while(t <= svn_itor) {
    weights = t(out[[t]]$coefs) %*% out[[t]]$SV#coefs系数乘以训练标签,SV是生成的支持向量，%*%表示矩阵相乘，转置coefs矩阵
    weights[which(weights<0)]<-0
    w<-weights/sum(weights)
    u <- sweep(X,MARGIN=2,w,'*')
    k <- apply(u, 1, sum)
    nusvm[t] <- sqrt((mean((k - y)^2)))
    corrv[t] <- cor(k, y)
    # 计算相关性和协方差和方差 应该是算相关性 
    t <- t + 1
  }
  #总体上是 下面doperm函数随机选取yr，随机挑选对应细胞和表达量，然后跑支持向量机，得到随机的nusvm，选一组最合适的model，一共运行1000次，再排序筛选找到最好的结
  
  #pick best model
  rmses <- nusvm
  mn <- which.min(rmses)
  print(mn)
  if(length(mn)==0){
    newList = NULL
    print("yr")
  }else{
    
    model <- out[[mn]]#算不出mn就返回空值，doperm把这次运行跳过，随机选取有可能会出错
  
  #get and normalize coefficients
  q <- t(model$coefs) %*% model$SV
  q[which(q<0)]<-0
  w <- (q/sum(q))
  
  mix_rmse <- rmses[mn]
  mix_r <- corrv[mn]
  
  newList <- list("w" = w, "mix_rmse" = mix_rmse, "mix_r" = mix_r)
  }
}

#' do permutations
#' @param perm Number of permutations
#' @param X cell-specific gene expression
#' @param y mixed expression per sample
#' @export
doPerm <- function(perm, X, Y){
  itor <- 1
  Ylist <- as.list(data.matrix(Y))
  dist <- matrix()
  
  while(itor <= perm){
    
    
    #random mixture
    yr <- as.numeric(Ylist[sample(length(Ylist),dim(X)[1])])
    
    #standardize mixture
    yr <- (yr - mean(yr)) / sd(yr)
    
    #run CIBERSORT core algorithm
    result <- CoreAlg(X, yr)
    
    if(is.null(result)){
      print(paste("error:",itor,sep=""))#打印第几次运行报错
      next()#跳过
    }
    
    mix_r <- result$mix_r
    
    #store correlation
    if(itor == 1) {dist <- mix_r}
    else {dist <- rbind(dist, mix_r)}
    
    itor <- itor + 1
    print(paste("now:",itor,sep=""))
  }
  newList <- list("dist" = dist)
}

#' Main functions
#' @param sig_matrix file path to gene expression from isolated cells
#' @param mixture_file heterogenous mixed expression
#' @param perm Number of permutations
#' @param QN Perform quantile normalization or not (TRUE/FALSE)
#' @export
CIBERSORT <- function(sig_matrix, mixture_file,mark_matrix=NULL,perm=0, QN=TRUE){
  #read in data
  X <- read.table(sig_matrix,header=T,sep="\t",row.names=1,check.names=F)
  if(!is.null(mark_matrix)){
    MARK <- read.table(mark_matrix,header=T,sep="\t",check.names=F)
    MARK = unique(MARK[,1])
    gene1=intersect(rownames(X),MARK)#跟去过重复的markergenes取交集
    X<-X[gene1,]#把行名在向量gene1里面的行保存下来
  }
  Y <- read.table(mixture_file, header=T, sep="\t", check.names=F)
  if(max(X) < 50) {
    X <- exp(1)^X 
    X<-X-1
  }
  Y <- Y[!duplicated(Y[,1]),]#把重复的基因都删了(把第一列重复的基因所在的行删了)
  rownames(Y)<-Y[,1]#数据要处理一下
  Y<-Y[,-1]
  if(max(Y) < 50) {
    Y <- exp(1)^Y 
         Y<-Y-1
         }
  #判断是否取过对数
  #数据的标化永远是+1取log的
  #order
  
  X <- data.matrix(X)
  Y <- data.matrix(Y)
  
  X <- X[order(rownames(X)),]
  Y <- Y[order(rownames(Y)),]
  
  P <- perm #number of permutations
  #quantile normalization of mixture file
  if(QN == TRUE){
    tmpc <- colnames(Y)
    tmpr <- rownames(Y)
    Y <- preprocessCore::normalize.quantiles(Y)
    colnames(Y) <- tmpc
    rownames(Y) <- tmpr
  }
  #intersect genes
  Xgns <- row.names(X)
  Ygns <- row.names(Y)
  YintX <- Ygns %in% Xgns
  Y <- Y[YintX,]
  XintY <- Xgns %in% row.names(Y)
  X <- X[XintY,]
  
  #standardize sig matrix
  X <- (X - mean(X)) / sd(as.vector(X))
  
  #empirical null distribution of correlation coefficients
  if(P > 0) {nulldist <- sort(doPerm(P, X, Y)$dist)}
   
  #print(nulldist)
  
  header <- c('Mixture',colnames(X),"P-value","Correlation","RMSE")
  #print(header)
  
  output <- matrix()
  itor <- 1
  mixtures <- dim(Y)[2]
  pval <- 9999
  
  #iterate through mixtures
  while(itor <= mixtures){
    
    y <- Y[,itor]
    
    #standardize mixture
    y <- (y - mean(y)) / sd(y)
    
    #run SVR core algorithm
    result <- CoreAlg(X, y)
    
    if(is.null(result)){
      print(paste("error:",itor,sep=""))#打印第几次运行报错
      next()#跳过
    }
    
    #get results
    w <- result$w
    mix_r <- result$mix_r
    mix_rmse <- result$mix_rmse
    
    #calculate p-value
    if(P > 0) {pval <- 1 - (which.min(abs(nulldist - mix_r)) / length(nulldist))}
    
    #print output
    out <- c(colnames(Y)[itor],w,pval,mix_r,mix_rmse)
    if(itor == 1) {output <- out}
    else {output <- rbind(output, out)}
    
    itor <- itor + 1
    
  }
  
  #save results
  write.table(rbind(header,output), file="CIBERSORT-Results.txt", sep="\t", row.names=F, col.names=F, quote=F)
  
  #return matrix object containing all results
  obj <- rbind(header,output)
  obj <- obj[,-1]
  obj <- obj[-1,]
  obj <- matrix(as.numeric(unlist(obj)),nrow=nrow(obj))
  rownames(obj) <- colnames(Y)
  colnames(obj) <- c(colnames(X),"P-value","Correlation","RMSE")
  save(obj,file = 'output_obj.Rdata')
}   #记得修改完函数重新调用一下

options(stringsAsFactors = FALSE)
options(future.globals.maxSize= 30*1024*1024*1024)
command = matrix(c('sig_matrix','s',1,"character",#每个参数传入都有变量名，0不需要写入，1必须写入东西，2都可以（类似于：传matrix就必须是1，因为要matrix=）
                'mixture_file','m',1,"character",
                'mark_matrix','k',2,"character",
                'QN','q',1,"logical",
                'perm','r',1,"integer",
                'path','o',1,"character")
                ,byrow=TRUE, ncol=4)
args=getopt(command)#调用传参矩阵
sig_matrix = args$sig_matrix
mixture_file = args$mixture_file
mark_matrix = args$mark_matrix
QN = args$QN
perm = args$perm
path = args$path#声明变量
#到这位置都是传参，
#写脚本的时候，函数可以不写进来，函数传在其他位置
sig_matrix="/lvdata/wzb/scRNA/FW2024-089/2024.05.30/Cibersort/gene_cluster.txt"
mixture_file="/lvdata/wzb/scRNA/FW2024-089/2024.05.30/Cibersort/GSE134347_mix.txt"
perm=100
QN=F
mark_matrix = NULL
path = '/lvdata/wzb/scRNA/FW2024-089/2024.05.30/Cibersort'

setwd(path)
CIBERSORT(sig_matrix,mixture_file,mark_matrix,perm, QN)
#CIBERSORT(sig_matrix,mixture_file="mixtures.txt",perm, QN)

cibersort_raw <- read.table("CIBERSORT-Results.txt",header = T,sep = '\t',check.names=F) %>% rename("Patients" = "Mixture") %>% select(-c("P-value","Correlation","RMSE"))
cibersort_tidy <- cibersort_raw %>% remove_rownames() %>% column_to_rownames("Patients")
flag <- apply(cibersort_tidy,2,function(x) sum(x == 0) <  dim(cibersort_tidy)[1]/2)
cibersort_tidy <- cibersort_tidy[,which(flag)] %>% as.matrix()
spotsname<-rownames(cibersort_tidy)
cibersort_tidy<-cbind(spotsname,cibersort_tidy)
write.table(cibersort_tidy,"result_matrix.txt",sep="\t",row.names = F,quote = F)
#' CIBERSORT R script v1.03 (last updated 07-10-2015)
#' Note: Signature matrix construction is not currently available; use java version for full functionality.
#' Author: Aaron M. Newman, Stanford University (amnewman@stanford.edu)
#' Requirements:
#'       R v3.0 or later. (dependencies below might not work properly with earlier versions)
#'       install.packages('e1071')
#'       install.pacakges('parallel')
#'       install.packages('preprocessCore')
#'       if preprocessCore is not available in the repositories you have selected, run the following:
#'            BiocManager::install(version="3.10")
#'            BiocManager::install('preprocessCore')
#' Windows users using the R GUI may need to Run as Administrator to install or update packages.
#' This script uses 3 parallel processes.  Since Windows does not support forking, this script will run
#' single-threaded in Windows.
#'
#' Usage:
#'       Navigate to directory containing R script
#'
#'   In R:
#'       source('CIBERSORT.R')
#'       results <- CIBERSORT('sig_matrix_file.txt','mixture_file.txt', perm, QN)
#'
#'       Options:
#'       i)  perm = No. permutations; set to >=100 to calculate p-values (default = 0)
#'       ii) QN = Quantile normalization of input mixture (default = TRUE)
#'
#' Input: signature matrix and mixture file, formatted as specified at http://cibersort.stanford.edu/tutorial.php
#' Output: matrix object containing all results and tabular data written to disk 'CIBERSORT-Results.txt'
#' License: http://cibersort.stanford.edu/CIBERSORT_License.txt
#' Core algorithm
#' @param X cell-specific gene expression
#' @param y mixed expression per sample
#' @export
