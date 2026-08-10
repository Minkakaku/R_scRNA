library(data.table)
library(dplyr)
library(tidyverse)
library(seqinr)
library(gggenomes)
library(Biostrings)

setwd("C:/Users/css/Desktop/蛋白组学工作/FW2024-351-2例细菌外泌体LFQ/database/")

# aa.percentage --去重复
fasta <- readDNAStringSet("Helicobacter_combined_20240507.fasta") # 见脚本 get_AA.py
t = names(fasta) 
tail(t)
head(t)
t1 = t[grep("\\|", t)]
t2 = sapply(strsplit(t1, split="\\|"), function(x) x[[2]])
which(duplicated(t2)) 

# NOG annotation
nog_annotation <- fread("Helicobacter_combined_20240507/NOG_annotation.tsv")
head(nog_annotation)
COG <- fread("4_Annotation/4_Annotation/COG/gene_to_COG.xls") %>% select(Gene_ID, Code, Name, COG_ID, COG_name)
COG$Name = gsub(" ;","",COG$Name)
colnames(COG) <- colnames(nog_annotation)
nog_annotation1 <- rbind(nog_annotation, COG)

write.table(nog_annotation1, file = "NOG_annotation.tsv",sep = "\t",quote = F,row.names = F)

# Pfam annotation
pfam_annotation <- fread("Helicobacter_combined_20240507/Pfam_annotation.tsv")
pfam <- fread("4_Annotation/4_Annotation/blast_best_hit/PFAM_blast.out.infor.best.xls")

t <- sapply(str_split(pfam$Sbjct, "\\."), function(x) x[[1]])
t1 <- sapply(str_split(t, ","), function(x) x[[1]])
pfam$Terms <- sapply(str_split(t, ", "), function(x) x[[3]])
pfam$pfam_ids = sapply(str_split(t1, " "), function(x) x[[2]])
pfam$pfam_ids = gsub("pfam","PF", pfam$pfam_ids)
pfam = subset(pfam, select=c(Query_ID, pfam_ids, Terms))
colnames(pfam) = colnames(pfam_annotation)

pfam_annotation = rbind(pfam_annotation, pfam) 
write.table(pfam_annotation, file = "Pfam_annotation.tsv",sep = "\t",quote = F,row.names = F)

# product 
product <- fread("Helicobacter_combined_20240507/product.tsv")
swiss <- fread("4_Annotation/4_Annotation/blast_best_hit/Swissprot_blast.out.infor.best.xls")
swiss$Source = sapply(str_split(swiss$Sbjct, "\\|"), function(x) x[[1]])

string = rbind(fread("string_mapping_1.tsv"), fread("string_mapping_2.tsv")) %>% select(queryItem, stringId, preferredName)
colnames(string) <- c("Protein accession", "Cross-reference (STRING)","Gene name")
  
tmp1 = left_join(product[grep("PROKKA",product$Source),] %>% select(Source,`Protein description`,  Length) %>% dplyr::rename("Protein accession" =Source), subset(string, select = c("Protein accession", "Gene name"))) %>% 
  left_join(subset(swiss,select =c("Query_ID","Source")) %>% dplyr::rename("Protein accession"=Query_ID)) %>% select(`Protein accession`, `Gene name`,`Protein description`,Source, Length)
product = rbind(product[-grep("PROKKA",product$Source),], tmp1)
product$`Gene name`[which(is.na(product$`Gene name`))] = "---"

write.table(product, file = "product.tsv",sep = "\t",quote = F,row.names = F)

# KEGG annotation 
library(RCurl)
library(stringr)
library(XML)
library(remotes)
library(downloader)

kegg.anno = fread("Helicobacter_combined_20240507/KEGG_annotation.tsv")
kegg = fread("4_Annotation/4_Annotation/KEGG/kegg_annot1.xls")

kegg$`KEGG ID`=paste0("hhe:", kegg$KO)
kegg = kegg %>% separate_rows(pathway, sep = ";")
kegg$`Pathway ID`=sapply(str_split(kegg$pathway, ","), function(x) x[[1]])
kegg$Terms=sapply(str_split(kegg$pathway, ","), function(x) x[[2]])

API <- "https://rest.kegg.jp"
# Download KO list 
download(paste0(API, "/list/ko"),"ko.txt")
ko <- read.table("ko.txt", sep = "\t", stringsAsFactors = F, quote="")
colnames(ko) <- c("KEGG category","KEGG description")

kegg = merge(kegg, ko, by.x= "KO",by.y = "KEGG category")
kegg = kegg %>% select(Gene_ID, `KEGG ID`, KO, `KEGG description`, `Pathway ID`, Terms)
colnames(kegg) = c("Protein accession", "KEGG ID", "KEGG category", "KEGG description", "Pathway ID","Terms")
# kegg$`Pathway ID` = gsub("ko","", kegg$`Pathway ID`)

# Download pathway list
download(paste0(API, "/get/br:br08901/json"), "pathway.json")
path_list <- rjson::fromJSON(file = "pathway.json")
path_list <- path_list$children
pathway_dat <- data.frame()
for (i in 1:length(path_list)){
  level1 <- path_list[[i]]$name
  tmp <- path_list[[i]]$children
  level2 <- unlist(lapply(tmp, '[',1))
  
  data2 <- data.frame()
  for (j in 1:length(tmp)){
    tmp2 <- tmp[[j]]$children 
    level3 <- unlist(lapply(tmp2, '[',1))
    data <- data.frame(level1 = level1, level2 = level2[j], level3 = level3, row.names = NULL)
    data2 <- rbind(data2, data)
  }
  pathway_dat <- rbind(pathway_dat, data2)
}
pathway_dat$ID <- str_extract(pathway_dat$level3,"\\d+")
pathway_dat$level3 <- gsub("[[:digit:]]", "",pathway_dat$level3)
pathway_dat <- pathway_dat %>% dplyr::select(ID, level3, level1, level2)
colnames(pathway_dat) <- c("Pathway ID","Terms","Pathway level 1", "Pathway level 2")
pathway_dat$Terms = trimws(pathway_dat$Terms,which = "both") 
pathway_dat$`Pathway ID` = paste0("ko", pathway_dat$`Pathway ID`)

kegg = left_join(kegg, subset(pathway_dat, select=c("Pathway ID","Pathway level 1", "Pathway level 2"))) 
kegg.anno = rbind(kegg.anno, kegg)
write.table(kegg.anno , file = "KEGG_annotation.tsv",sep = "\t",quote = F,row.names = F)

# GO annotation 
go.anno = fread("Helicobacter_combined_20240507/GO_annotation.tsv")

go1 = fread("4_Annotation/4_Annotation/GO/full_GO_annot.xls") %>% select(Query_ID, Ontology)
colnames(go1) = c("Protein accession","ID")
go2 = fread("4_Annotation/4_Annotation/GO/full_GO_annot_all_level.xls") %>% select(Level, GO, Term, Ontology)
colnames(go2) = c("GO levels","ID","Terms","Category")
go = merge(go1, go2, by="ID") %>% select(`Protein accession`, Category, ID, `GO levels`, Terms)
go$Category <- gsub("_"," ", go$Category)
library(Hmisc)
go$Category <- capitalize(go$Category)
go <- rbind(go.anno, go)
write.table(go , file = "GO_annotation.tsv",sep = "\t",quote = F,row.names = F)

# annotation 
annotation <- fread("Helicobacter_combined_20240507/annotation.tsv")

product.tmp <- product %>% select(`Protein accession`,`Protein description`,`Gene name`) 
go.tmp <- go %>% group_by(`Protein accession`) %>% mutate(`Gene ontology IDs` = paste(ID, collapse ="; ")) %>% distinct(`Protein accession`,`Gene ontology IDs`)
kegg.tmp <- kegg.anno %>% group_by(`Protein accession`) %>% mutate(`Cross-reference (KEGG)` = paste(`KEGG ID`, collapse = ";")) %>% distinct(`Protein accession`,`Cross-reference (KEGG)`)
nog.tmp <- nog_annotation1 %>% group_by(`Protein accession`) %>% mutate(`Cross-reference (eggNOG)` = paste(`eggNOG ID`, collapse = ";")) %>% distinct(`Protein accession`,`Cross-reference (eggNOG)`)
pfam.tmp <- pfam_annotation %>% group_by(`Protein accession`) %>% mutate(`Cross-reference (Pfam)` = paste(`Pfam ID`, collapse = ";")) %>% distinct(`Protein accession`,`Cross-reference (Pfam)`)
string.tmp <- string %>% group_by(`Protein accession`) %>% mutate(`Cross-reference (STRING)`= paste(`Cross-reference (STRING)`, collapse = ";")) %>% distinct(`Protein accession`,`Cross-reference (STRING)`)
string.tmp$`Cross-reference (STRING)` = paste0(string.tmp$`Cross-reference (STRING)`, ";")
string.tmp = rbind(subset(annotation, select = c(Entry, `Cross-reference (STRING)`)) %>% dplyr::rename("Protein accession"=Entry), string.tmp)

temp = left_join(product.tmp, go.tmp)
temp = left_join(temp, kegg.tmp)
temp = left_join(temp, nog.tmp)
temp = left_join(temp, pfam.tmp)
temp = left_join(temp, string.tmp)
colnames(temp) = colnames(annotation)
write.table(annotation, "annotation.tsv", sep = "\t",quote = F,row.names = F)

# KEGG image download 
##----------- KEGG pathway image & html & xml download -------------##
library(downloader)

pathways <- unique(kegg$`Pathway ID`)
mydir <- "C:/Users/css/Desktop/蛋白组学工作/FW2024-351-2例细菌外泌体LFQ/database/"
for (i in 1:length(pathways)){
  API <- "http://rest.kegg.jp/"
  ID <- paste0("hhe",pathways[i])
  filename <- paste0("kegg.dir/",ID)
  
  ## 下载KEGG pathway图片和xml版本
  download(paste0(API,"get/",ID,"/image"), paste0(filename,".png"), mode = "wb")
  download(paste0(API,"get/",ID,"/kgml"), paste0(filename,".xml"), mode = "wb")
  
  # 下载KEGG pathway的html网页版
  download(paste0("http://www.kegg.jp/pathway/", ID), paste0(filename, ".html"), mode = "wb") 
  hsa <- readLines(paste0(filename, ".html"))
  hsa <- gsub(" href=\"/", " href=\"https://www.kegg.jp/", hsa)
  
  hsa <- gsub(paste0("/kegg/pathway/hhe/",ID, ".png"), paste0(mydir,"/",filename, ".png"), hsa)
  # hsa <- gsub("/kegg/pathway/hsa/hsa04810.png", "https://rest.kegg.jp/get/hsa04810/image", hsa)
  write.table(hsa, paste0(filename, ".html"), quote = F , col.names = F, row.names = F)
}

mydir = "C:/Users/css/Desktop/蛋白组学工作/FW2024-351-2例细菌外泌体LFQ/database/kegg.dir"
get_kegg_image <- function(annotation,API = "https://rest.kegg.jp"){
  pathway <- unique(annotation$`Pathway ID`)
  # if(dir.exists("DB/kegg.dir")){
  #   unlink("DB/kegg.dir",recursive = T)
  # }
  # dir.create("DB/kegg.dir")
  for(i in 1:length(pathway)){
    if(pathway[i] == "ko00473"|pathway[i] == "ko00471"|pathway[i]=="ko00312"){
      next
    }
    download(paste(API,"get",pathway[i],"image",sep = "/"),paste("kegg.dir/",pathway[i],".png",sep = ""),mode = "wb")
    download(paste(API,"get",pathway[i],"kgml",sep = "/"),paste("kegg.dir/",pathway[i],".xml",sep = ""),mode = "w")
    download(paste("http://kegg.jp/pathway",pathway[i],sep = "/"),paste("kegg.dir/",pathway[i],".html",sep = ""),mode = "w")
    html <- readLines(paste("kegg.dir/",pathway[i],".html",sep = ""))
    html <- gsub(" href=\"/"," href=\"https://www.kegg.jp/",html)
    html <- gsub("/kegg/pathway/hhe","image",html)
    write.table(html,paste("kegg.dir/",pathway[i],".html",sep = ""),quote = F,row.names = F,col.names = F)
  }
}

get_kegg_image(kegg,API = "https://rest.kegg.jp")

evidence <- fread("C:/Users/css/Desktop/蛋白组学工作/FW2024-351-2例细菌外泌体LFQ/Q2_FW2024_351_Combine_Proteins.txt")
