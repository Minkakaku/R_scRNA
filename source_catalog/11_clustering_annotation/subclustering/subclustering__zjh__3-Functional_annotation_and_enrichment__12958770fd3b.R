library(data.table)
library(dplyr)
library(tidyverse)
library(ggsci)
library(ggh4x)
library(RColorBrewer)
library(corto)
library(ggrepel)
library(ggplot2)
library(readxl)
library(igraph)
library(ggraph)
library(ggpubr) 
library(reshape2)
############################## Parameters ######################################
# workflow = "Metabolomics" # Must be one of "Metabolomics", "Lipidomics" 
# organism = "hsa" # Must be organism IDs from KEGG 
# sample_info = "input/samples.csv"
# comparison = "input/comparison.txt"
# db_dir = "DB/HMDB_20250114"
# results_dir = "Results"

diff_file = snakemake@input[["dep_result"]]
workflow = snakemake@params[["workflow"]]
organism = snakemake@params[["organism"]]
sample_info = snakemake@input[["sample_info"]]
db_dir1 = snakemake@params[["db_dir"]]
results_dir = snakemake@params[["results_dir"]]

############################ Functions defined #################################
diff.classy.plot1=function(data, font = 14,color.pal="Set3",  angle = 75, margin =15, margin2 = 10, type ="Metabolomics", label.font=4){
  if(type=="Metabolomics"){
    n = brewer.pal.info[color.pal,"maxcolors"]
    if (length(unique(data$HMDB_SuperClass))>n){
      colors = colorRampPalette(ggsci.colors[[color.pal]][1:n])(length(unique(data$HMDB_SuperClass)))
    }else{
      colors = ggsci.colors[[color.pal]][1:length(unique(data$HMDB_SuperClass))]
    }
    names(colors)=unique(data$HMDB_SuperClass)
    data=data %>% arrange(-Freq)
    data$HMDB_SuperClass = factor(data$HMDB_SuperClass, levels=unique(data$HMDB_SuperClass))
    p1 = ggplot(data,aes(reorder(HMDB_Class,-Freq),Freq,fill=HMDB_SuperClass))+geom_col()+facet_grid(.~HMDB_SuperClass,scales = 'free',space = 'free')+
      scale_fill_manual(values = colors)+
      theme_classic(base_size = font)+scale_y_continuous(expand = expansion(mult = c(0,0.1)))+
      theme(strip.text = element_blank(),axis.text.x=element_text(angle = angle,hjust=1),axis.text = element_text(color='black'), 
            plot.margin=unit(c(margin2,5,5,margin),'mm'),
            legend.background=element_blank(), legend.box.background=element_blank())+
      geom_text(aes(label=Freq),size=label.font, vjust=0)+labs(x = "HMDB Category",y="Number of compounds")+
      theme(legend.position ="top",legend.title.position="top")+
      guides(fill=guide_legend(ncol=3))
  }else{
    colors = ggsci.colors[[color.pal]][1:length(unique(data$Categories))]
    names(colors)=unique(data$Categories)
    p1 = ggplot(data, aes(reorder(Name,Freq),Freq,fill=Categories))+geom_col()+
      facet_grid(.~Categories,scales = 'free',space = 'free')+scale_fill_manual(values = colors)+theme_classic(base_size = font)+
      scale_y_continuous(expand = expansion(mult = c(0,0.1)))+
      theme(strip.text = element_blank(),axis.text.x=element_text(angle = angle,hjust=1),axis.text = element_text(color='black'), 
            plot.margin=unit(c(margin2,5,5,margin),'mm'))+
      geom_text(aes(label=Freq),size=label.font, vjust=0)+theme(legend.background=element_blank(),
                                                                legend.box.background=element_blank())+labs(x = "Lipid category",y="Number of compounds")+
      theme(legend.position ="top",legend.title.position="top")+
      guides(fill=guide_legend(ncol = 2))
  }
  return(p1)
}
diff.classy.plot2=function(data, font=14, color.pal="Set2", angle = 60, margin = 15, margin2 = 10,label.font=4){
  colors = ggsci.colors[[color.pal]][1:length(unique(data$Pathway.level.1))]
  names(colors)=unique(data$Pathway.level.1)
  p2 = ggplot(data,
              aes(reorder(Pathway.level.2,Freq),Freq,fill=Pathway.level.1))+geom_col()+
    facet_grid(.~Pathway.level.1,scales = 'free',space = 'free')+scale_fill_manual(values = colors)+theme_classic(base_size = 14)+
    scale_y_continuous(expand = expansion(mult = c(0,0.1)))+
    theme(strip.text = element_blank(),axis.text.x=element_text(angle = angle,hjust=1),axis.text = element_text(color='black'), 
          plot.margin=unit(c(margin2,5,5,margin),'mm'))+
    geom_text(aes(label=Freq),vjust=0,size=label.font)+theme(legend.background=element_blank(), 
                                                             legend.box.background=element_blank())+labs(x = "KEGG pathway",y="Number of compounds", fill ="Pathway Level 1")+
    theme(legend.position ="top", legend.title.position="top")+
    guides(fill=guide_legend(ncol = 2))
  return(p2)
}
ggsci.colors <- list(
  "Set1"=brewer.pal(9,"Set1"),
  "Set2"=brewer.pal(8,"Set2"),
  "Set3"=brewer.pal(12,"Set3"),
  "Pastel1"=brewer.pal(9,"Pastel1"),
  "Pastel2"=brewer.pal(8,"Pastel2"),
  "Paired"=brewer.pal(12,"Paired"),
  "Dark2"=brewer.pal(8,"Dark2"),
  "Accent"=brewer.pal(8,"Accent"),
  "npg"=pal_npg()(10),
  "aaas"=pal_aaas()(10),
  "nejm"=pal_nejm()(8),
  "lancet"=pal_lancet()(9),
  "jama"=pal_jama()(7),
  "jco"=pal_jco()(10),
  "ucscgb"=pal_ucscgb()(26),
  "d3_category10"=pal_d3("category10")(10),
  "d3_category20"=pal_d3("category20")(20),
  "d3_category20c"=pal_d3("category20c")(20),
  "igv"=pal_igv()(51),
  "uchicago"=pal_uchicago()(9),
  "rickandmorty"=pal_rickandmorty()(12),
  "simpsons"=pal_simpsons()(16)
)
enrich.test=function(x){a=x[1]; b=x[2]-a; c=x[3]-a; d=x[4]-a-b-c; p.value=fisher.test(matrix(c(a,b,c,d),nrow = 2),alternative = "greater")$p.value}
enrich.demo=function(go.data,diff.list){
  go.data.diff=subset(go.data,`Metabolite name` %in% diff.list)
  mapping=data.frame(table(go.data.diff$Terms))
  background=data.frame(table(go.data$Terms))
  go.data.state=merge(mapping,background,by="Var1")
  colnames(go.data.state)=c("Terms","mapping","background")
  go.data.state$`All mapping`=length(unique(go.data.diff$`Metabolite name`))
  go.data.state$`All background`=length(unique(go.data$`Metabolite name`))
  go.data.state$`Fold enrich`=round(go.data.state$mapping/go.data.state$background/go.data.state$`All mapping`*go.data.state$`All background`,2)
  go.data.state=subset(go.data.state,mapping>1)
  go.data.state$p.value=signif(apply(go.data.state[2:5],1,enrich.test),3)
  go.data.state$fdr=signif(p.adjust(go.data.state$p.value),3)
  go.data.state=subset(go.data.state,select=c("Terms","mapping","Fold enrich","p.value","fdr")) %>% merge(subset(go.data.diff, select=c("Terms","Metabolite name")), by="Terms") %>% 
    group_by(Terms, mapping, `Fold enrich`,p.value, fdr) %>% summarize(Metabolites = paste(`Metabolite name`, collapse = " ; "))
  return(go.data.state)
}
gseaplot=function(t){
  temp.data=data.frame(`Runing Enrichment Score`=t$running_score,reflist=t$reflist,inSet=t$inSet,Rank=rank(-t$reflist))
  if(t$nes<0){line.color='darkgreen'} else{line.color='red'}
  ggarrange(
    ggplot(temp.data)+geom_col(aes(Rank,reflist),color="steelblue")+theme_classic()+
      theme(axis.text.x=element_blank(),
            axis.title.x=element_blank(),
            axis.ticks.x=element_blank(),
            axis.line.x = element_blank())+
      annotate(geom="text",color=line.color,label=sprintf("P value=%.0e; ES=%.3f; NES=%.3f",t$p.value,t$es,t$nes),x = nrow(temp.data)/2,
               y=min(temp.data$reflist)/2)+
      labs(title = t$name,y="Log2fc")+theme(plot.title = element_text(hjust = 0.5,size =10)),
    ggplot(temp.data)+geom_line(aes(Rank,t$running_score),color=line.color)+geom_point(data=temp.data %>% filter(inSet==1),aes(x=Rank,y=0),shape=3,size=5)+
      theme_classic()+geom_hline(yintercept = 0)+labs(y="Running Enrichment Score"),
    nrow=2,heights = c(1,2),align = 'v'
  )
}
getpath=function(obo,term){
  p=obo$parents[[term]]
  if(length(p)>0){
    out.data=data.frame(source=p,target=term)
    out=rbindlist(lapply(obo$parents[[term]],function(x) getpath(obo,x)))
    return(rbind(out.data,out))
  }else{
    return(NULL)
  }
}
obo.plot=function(enrich.out, obo, col.1 = "blue", col.2="red"){
  enrich.out %>% slice_min(p.value,n = 5)
  network=rbindlist(lapply(enrich.out$ID[enrich.out$p.value<0.05],function(x) getpath(obo,x))) %>% distinct_all()
  node=data.frame(ID=unique(c(network$source,network$target))) %>% merge(enrich.out,by="ID",all.x=T)
  ig=graph_from_data_frame(network,directed = T,vertices = node)
  ggraph(ig)+geom_edge_link(color='black',arrow = arrow(length = unit(0.02, "npc")))+
    scale_fill_gradientn(colours =c(col.1,col.2),limits=c(0,5),oob=scales::squish)+theme_void()+
    geom_node_point(aes(size=`Fold enrich`,fill=-log10(p.value)),shape=21)+
    geom_node_text(aes(label=stringr::str_wrap(Terms, width=20)),size=4,repel = T)+
    labs(title = enrich.out$compare[1])+theme(plot.title = element_text(face = 'bold'))
}

######################## Functional enrichment  ################################
set.seed(123)
compare.name=sub(sub(diff_file,pattern=".*of ",replacement=""),pattern=".xls$",replacement="")
if(workflow=="Metabolomics"){
  print(db_dir1)
  if(db_dir1 != "-"){db_dir = db_dir1}else{db_dir = "/mnt/NexData_file/cosmicdarahub/Bioinformatics/Database/Metabolism/DB/HMDB_20250114"}
  db_dir = gsub("/home/ztt/data/","/mnt/NexData_file/cosmicdarahub/Bioinformatics/Database/",db_dir)
  db = sub(".*\\/","",db_dir)
  
  data.diff = fread(diff_file)
  
  kegg.data=fread(paste0(db_dir,"/KEGG_compound_pathway_annotation.tsv")) %>%
    filter(`Pathway ID` %in% sub(fread(paste0("https://rest.kegg.jp//list/pathway/",organism),header = F)$V1,pattern=organism,replacement="map"))
  diff.data = separate_rows(data.diff,cpd_id,sep = "; ")
  
  if(db == "HMDB_20250114"){
    data.kegg=merge(diff.data,kegg.data,by=c("cpd_id","HMDB ID")) %>% distinct() %>% dplyr::select(`Pathway ID`,Terms,`Pathway level 1`,`Pathway level 2`,everything())
  }else{
    data.kegg=merge(diff.data,kegg.data,by=c("cpd_id")) %>% distinct() %>% dplyr::select(`Pathway ID`,Terms,`Pathway level 1`,`Pathway level 2`,everything())
  }
  
  if(length(data.kegg$`Metabolite name`[data.kegg$regulation!='unchanged'])>0){
    kegg.enrich=enrich.demo(data.kegg %>% filter(`Pathway level 1`=='Metabolism') %>% 
                              distinct(Terms,`Metabolite name`),data.diff$`Metabolite name`[data.diff$regulation != "unchanged"])
  }else{
    kegg.enrich = NULL
  }
  data.kegg.sub=data.kegg %>% group_by(Terms) %>% filter(length(Terms)>5,`Pathway level 1`=='Metabolism') %>% 
    ungroup() %>% select(!cpd_id) %>% distinct()
  data.kegg.anno=data.kegg %>% group_by(Terms) %>% summarise(Metabolites=paste(unique(`Metabolite name`),collapse = "; "))
  
  if(!is.null(kegg.enrich)){
    if(nrow(kegg.enrich)>30){
      kegg.tmp=kegg.enrich %>% as.data.frame() %>% slice_min(p.value, n=30, with_ties = FALSE)
    }else{kegg.tmp =kegg.enrich}
    xmax1=as.integer(max(kegg.tmp$`Fold enrich`))+1
    xmin1=as.integer(min(kegg.tmp$`Fold enrich`))-0.3
    enrich.plot1=ggplot(kegg.tmp,aes(`Fold enrich`,reorder(Terms,`Fold enrich`)))+
      geom_point(aes(size=mapping,fill=p.value),shape=21)+
      theme_test(base_size=14)+scale_size(range=c(3,6))+
      scale_fill_material("red",reverse = T)+ 
      xlim(xmin1-0.05*xmin1, xmax1) +
      geom_text(aes(label=case_when(p.value<0.01~'**',p.value<0.05~'*')),hjust=2,vjust=0.8,size=5) + 
      labs(y = "KEGG pathway",x="Fold enrich")
  }else{
    enrich.plot1=ggplot()+
      annotate("text",x=0.5,y=0.5,label="No terms enriched.",hjust=0.5,vjust=0.5,size=6)+
      xlim(0,1)+ylim(0,1)+
      theme_void()
  }
  
  #---------MSEA enrichment -----------#
  set.seed(123)
  gsea.res=lapply(split(data.kegg.sub,data.kegg.sub$`Pathway ID`),function(x){t=gsea(setNames(data.diff$log2fc,data.diff$`Metabolite name`),x$`Metabolite name`,method = 'permutation',np = 500);t$name=x$Terms[1];t$ID=x$`Pathway ID`[1];t$N=sum(t$inSet);return(t)})
  gsea.enrich = rbindlist(lapply(gsea.res, function(x) {
    data.frame(
      Terms = x$name,
      NES = signif(x$nes, 3),
      p.value = x$p.value,
      Count = x$N
    )
  }), idcol = "Pathway ID") %>%
    merge(data.kegg.anno, by = "Terms")
  tmp.data <- gsea.enrich %>% filter(p.value < 0.05)
  if(!is.null(gsea.enrich)){
    xmax2=as.integer(max(-log10(gsea.enrich$p.value)))+1
    xmin2=0
    enrich.plot2 = ggplot(gsea.enrich %>% filter(p.value<0.05),aes(-log10(p.value),reorder(Terms,-log10(p.value)),fill=NES))+
      geom_point(aes(size=Count),shape=21)+
      scale_fill_gradientn(colors = c("#00A087FF","white","#E64B35FF"),limits=c(-2,2),oob=scales::squish)+
      theme_test(base_size=14)+guides(color='none')+ xlim(xmin2, xmax2) +
      geom_text(aes(label=case_when(p.value<0.01~'**',p.value<0.05~'*')),hjust=2,vjust=0.8,size=5) + 
      labs(y = "Terms",x="-log10(p.value)")
  }else{
    enrich.plot2=ggplot()+
      annotate("text",x=0.5,y=0.5,label="No terms enriched.",hjust=0.5,vjust=0.5,size=6)+
      xlim(0,1)+ylim(0,1)+
      theme_void()
  }
  if(nrow(tmp.data)>0){
    res.plot=lapply(gsea.res[sapply(gsea.res,function(x) x$p.value<0.05)],gseaplot)
    ncol=4
    nrow=ceiling(length(res.plot)/ncol)
    enrich.plot3=ggarrange(plotlist = res.plot,ncol = ncol,nrow=nrow) 
  }else{
    res.plot = NULL
    nrow = 1
    enrich.plot3=ggplot()+
      annotate("text",x=0.5,y=0.5,label="No terms enriched.",hjust=0.5,vjust=0.5,size=6)+
      xlim(0,1)+ylim(0,1)+
      theme_void()
  }
  
  # Diff compound classification
  classify.data1= data.diff %>% filter(regulation != "unchanged") %>% 
    group_by(HMDB_SuperClass, HMDB_Class) %>% 
    summarize(Freq= length(unique(`Metabolite name`)), Metabolites = paste(`Metabolite name`, collapse = " ; ")) %>% ungroup()
  if(nrow(classify.data1)>0){
    comclass.plot1=diff.classy.plot1(classify.data1 %>% slice_max(Freq,n=40),font=14,color.pal="Set3",type=workflow, angle= 75, margin =16, margin2=10, label.font=4)
  }else{
    comclass.plot1=ggplot()+theme(panel.background = element_blank())
  }
  
  classify.data2=data.kegg %>% filter(regulation!="unchanged") %>% 
    group_by(`Pathway level 1`,`Pathway level 2`) %>% 
    summarize(Freq=length(unique(`Metabolite name`)), HMDBs= paste(`HMDB ID`, collapse = ";"), Metabolites=paste(`Metabolite name`, collapse = " ; ")) %>% ungroup()
  if(nrow(classify.data2)>0){
    colnames(classify.data2)=gsub(" ",".", colnames(classify.data2))
    comclass.plot2=diff.classy.plot2(classify.data2 %>% slice_max(Freq,n=40),font=14,color.pal="Set2", angle= 60, margin =18,margin2=14,label.font=4)
  }else{
    comclass.plot2=ggplot()+
      annotate("text",x=0.5,y=0.5,label="No pathway mapping.",hjust=0.5,vjust=0.5,size=6)+
      xlim(0,1)+ylim(0,1)+
      theme_void()
  }
}

if(workflow=="Lipidomics"){
  if(db_dir1 != "-"){db_dir = db_dir1}else{db_dir = "/mnt/NexData_file/cosmicdarahub/Bioinformatics/Database/Metabolism/DB/LION_20250116"}
  db_dir = gsub("/home/ztt/data/","/mnt/NexData_file/cosmicdarahub/Bioinformatics/Database/",db_dir)
  sample = fread(sample_info)
  if(all(sample$Sample == "-") == T){
    sample = sample %>% mutate(Sample = ID)
  } else {
    sample = sample 
  }
  meta.data=read_excel(paste(results_dir,"1-MS_identified_summary/Table 1.2 MS identified lipids detail_normalized.xlsx",sep="/"))
  data.diff = fread(diff_file)
  
  obo=readRDS(paste0(db_dir, "/obo.RDS"))
  values = list()
  #----------LIMAPS enrichment-----------#
  limap.enrich=enrich.demo(data.diff %>% dplyr::rename(Terms="Name"), data.diff$`Metabolite name`[data.diff$regulation != "unchanged"])
  
  #----------LION enrichment-------------#
  LION.annot=rbindlist(lapply(split(meta.data, meta.data$`Metabolite name`), function(x){
    id=unique(getpath(obo,x$`LION ID`[1])$source)
    data.frame(ID=id, Name=obo$name[id])
  }),idcol="Metabolite name", fill = T)
  LION.enrich=enrich.demo(LION.annot %>% dplyr::rename(Terms ="Name"), 
                          data.diff$`Metabolite name`[data.diff$regulation != "unchanged"]) %>% 
    merge(LION.annot %>% distinct(Name, ID), by.x = "Terms", by.y="Name")

  s = unlist(str_split(compare.name,"-vs-"))

  if(nrow(limap.enrich)>0){
    xmax1=as.integer(max(limap.enrich$`Fold enrich`)+1)
    xmin1=as.integer(min(limap.enrich$`Fold enrich`))
    enrich.plot1=ggplot(limap.enrich,aes(`Fold enrich`,reorder(Terms,`Fold enrich`)))+
      geom_point(aes(size=mapping,fill=p.value),shape=21)+ 
      theme_test(base_size = 14)+scale_size(range=c(3,6))+
      scale_fill_material("red",reverse = T)+ 
      geom_text(aes(label=case_when(p.value<0.01~'**',p.value<0.05~'*')),hjust=2,vjust=0.7,size=5)+
      labs(y = "LIMAPS Name",x="Fold enrich")+ 
      xlim(xmin1-0.05*xmin1, xmax1)
  }else{
      enrich.plot1 = ggplot()+theme(panel.background = element_blank())
  }
  tmp = LION.enrich %>% filter(p.value<0.05)
  
  if(nrow(tmp)>0){
    xmax2=as.integer(max(LION.enrich$`Fold enrich`)+1)
    xmin2=as.integer(min(LION.enrich$`Fold enrich`))
    enrich.plot2=ggplot(LION.enrich %>% filter(p.value<0.05) %>% slice_min(p.value,n=30),aes(`Fold enrich`,reorder(Terms,`Fold enrich`)))+
      geom_point(aes(size=mapping,fill=p.value),shape=21)+
      theme_test(base_size=14)+scale_size(range=c(3,6))+
      scale_fill_material("red",reverse = T)+
      geom_text(aes(label=case_when(p.value<0.0001~'***',p.value<0.01~'**',p.value<0.05~'*')),hjust=1.5,vjust=0.7,size=5)+
      labs(y = "Lipid ontology",x="Fold enrich")+
      xlim(xmin2-0.05*xmin2, xmax2)
    enrich.plot3=obo.plot(LION.enrich, obo=obo)
  }else{
    enrich.plot2 = ggplot()+theme(panel.background = element_blank())
    enrich.plot3 = ggplot()+theme(panel.background = element_blank())
  }
  # Diff compound classification
  classify.data1=data.diff %>% filter(regulation!='unchanged') %>% ungroup() %>% group_by(Categories, Name) %>% summarize(Freq=length(unique(`Metabolite name`)), Metabolites = paste(`Metabolite name`, collapse = " ; "))
  comclass.plot1=diff.classy.plot1(classify.data1 %>% slice_max(Freq,n=40),font=14,color.pal="Set3",type=workflow, angle= 60, margin =10,margin2=10)
  classify.data2=data.diff %>% melt(measure.vars = sample$Sample[sample$Group %in% s],variable.name = "Sample") %>% 
    merge(sample, by="Sample") %>%  filter(Group %in% s) %>% 
    filter(regulation !="unchanged") %>% group_by(Name, Group) %>%
    summarize(mean=mean(value), ymin=mean_se(value)[[2]], ymax=mean_se(value)[[3]])
  comclass.plot2 = classify.data2 %>% ggplot(aes(x=Name, y = mean, fill = Group))+
    geom_bar(stat="identity",position="dodge", color = "black") +
    geom_errorbar(aes(x=Name,ymin=ymin, ymax=ymax), stat="identity", position=position_dodge2(width = 0.4, padding = 0.5)) +
    theme_bw(base_size = 14) +
    scale_fill_manual(values = c("#71D0F5FF","#197EC0FF")) + 
    theme(axis.text.x =element_text(angle = 60,hjust = 1),plot.margin=unit(c(5,5,5,10),'mm')) +
    labs(x= "Category", y = "Lipid Abundance")
}

########################## Export Results #####################################
temp.dir = file.path(results_dir,"4-Functional Annotation and Enrichment")
if(!dir.exists(temp.dir)) dir.create(path = temp.dir)


if(workflow=="Metabolomics"){
  dir.create(paste0(temp.dir, "/msea_ind"))
  
  ggsave(comclass.plot1,filename = paste0(temp.dir,"/Figure 4a. Diff Metabolite Classification by Category of ",compare.name,".pdf"), width= 12, height=7)
  ggsave(comclass.plot1,filename = paste0(temp.dir,"/Figure 4a. Diff Metabolite Classification by Category of ",compare.name,".png"), width= 12, height=7)
  ggsave(comclass.plot2,filename = paste0(temp.dir,"/Figure 4b. Diff Metabolite Classification by Pathway of ",compare.name,".pdf"), width= 10, height=7)
  ggsave(comclass.plot2,filename = paste0(temp.dir,"/Figure 4b. Diff Metabolite Classification by Pathway of ",compare.name,".png"), width= 10, height=7)

  fwrite(data.kegg,file=paste0(temp.dir, "/Table 4a.KEGG Annotation of ", compare.name,".xls"),sep ="\t")
  if(is.null(kegg.enrich)){
    fwrite(
      data.table(
        Terms=character(),
        mapping=integer(),
        `Fold enrich`=numeric(),
        p.value=numeric(),
        fdr=numeric(),
        Metabolites=character()
      ),
      file=paste0(temp.dir, "/Table 4b.KEGG Enrichment of ",compare.name,".xls"),
      sep="\t"
    )
  }else{
    fwrite(kegg.enrich,file=paste0(temp.dir, "/Table 4b.KEGG Enrichment of ",compare.name,".xls"),sep ="\t")
  }
  ggsave(enrich.plot1,filename = paste0(temp.dir,"/Figure 4c. KEGG enrichment of ",compare.name,".pdf"), width= 8, height=6)
  ggsave(enrich.plot1,filename = paste0(temp.dir,"/Figure 4c. KEGG enrichment of ",compare.name,".png"), width= 8, height=6)
  fwrite(gsea.enrich,file=paste0(temp.dir, "/Table 4c.MSEA Enrichment of ",compare.name,".xls"),sep ="\t")
  ggsave(enrich.plot2,filename = paste0(temp.dir,"/Figure 4d. MSEA enrichment1 of ",compare.name,".pdf"), width= 8, height=6)
  ggsave(enrich.plot2,filename = paste0(temp.dir,"/Figure 4d. MSEA enrichment1 of ",compare.name,".png"), width= 8, height=6)
  ggsave(enrich.plot3,filename = paste0(temp.dir,"/Figure 4e. MSEA enrichment2 of ",compare.name,".pdf"), width= 12, height=3*nrow)
  ggsave(enrich.plot3,filename = paste0(temp.dir,"/Figure 4e. MSEA enrichment2 of ",compare.name,".png"), width= 12, height=3*nrow,bg = "white")
  
  dir.create(paste0(temp.dir, "/msea_ind/",compare.name))
  saveRDS(gsea.res, file = paste0(temp.dir, "/MSEA res of ",compare.name, ".Rds"))
  if(!is.null(res.plot)){
    lapply(names(res.plot),function(y){
      p = plot(res.plot[[y]])
      ggsave(paste0(temp.dir, "/msea_ind/",compare.name,"/MSEA Enrichment plot for ",y,".pdf"),plot = p, width =8, height = 8)
      ggsave(paste0(temp.dir, "/msea_ind/",compare.name,"/MSEA Enrichment plot for ",y,".png"),plot = p, width =8, height = 8)
    })
  }else{
    p = ggplot()+theme(panel.background = element_blank())
    ggsave(paste0(temp.dir, "/msea_ind/",compare.name,"/MSEA Enrichment plot.pdf"),plot = p, width =8, height = 8)
    ggsave(paste0(temp.dir, "/msea_ind/",compare.name,"/MSEA Enrichment plot.png"),plot = p, width =8, height = 8)
  }
}
if(workflow=="Lipidomics"){
  fwrite(limap.enrich,file=paste0(temp.dir, "/Table 4a.LIPID MAPS Enrichment of ",compare.name,".xls"),sep ="\t")
  fwrite(LION.enrich,file=paste0(temp.dir, "/Table 4b.LION Enrichment of ",compare.name,".xls"),sep ="\t")
  
  ggsave(comclass.plot1,filename = paste0(temp.dir,"/Figure 4a. Diff Lipid Classification by Category of ",compare.name,".pdf"), width= 10, height=6)
  ggsave(comclass.plot1,filename = paste0(temp.dir,"/Figure 4a. Diff Lipid Classification by Category of ",compare.name,".png"), width= 10, height=6)
  ggsave(comclass.plot2,filename = paste0(temp.dir,"/Figure 4b. Diff Lipid Abundance by Category of ",compare.name,".pdf"), width= 10, height=6)
  ggsave(comclass.plot2,filename = paste0(temp.dir,"/Figure 4b. Diff Lipid Abundance by Category of ",compare.name,".png"), width= 10, height=6)
  ggsave(enrich.plot1,filename = paste0(temp.dir,"/Figure 4c. LIPID MAPS enrichment of ",compare.name,".pdf"), width= 8, height=6)
  ggsave(enrich.plot1,filename = paste0(temp.dir,"/Figure 4c. LIPID MAPS enrichment of ",compare.name,".png"), width= 8, height=6)
  ggsave(enrich.plot2,filename = paste0(temp.dir,"/Figure 4d. LION enrichment1 of ",compare.name,".pdf"), width= 8, height=6)
  ggsave(enrich.plot2,filename = paste0(temp.dir,"/Figure 4d. LION enrichment1 of ",compare.name,".png"), width= 8, height=6)
  ggsave(enrich.plot3,filename = paste0(temp.dir,"/Figure 4e. LION enrichment2 of ",compare.name,".pdf"), width= 11, height=11)
  ggsave(enrich.plot3,filename = paste0(temp.dir,"/Figure 4e. LION enrichment2 of ",compare.name,".png"), width= 11, height=11, bg="white")
  
}
