.libPaths("/home/wzb/R/x86_64-pc-linux-gnu-library/4.2")
library(dplyr)
library(ggplot2)
library(ggsci)

args=commandArgs(T)
# args[1] 传入数据metadata
# args[2] 定义列
# args[3] 分组列
# args[4] 输出地址
# args[5] 定義分組的顔色列表
args[1] = as.character(args[1])
args[2] = as.character(args[2])
args[3] = as.character(args[3])
args[4] = as.character(args[4])
args[5] = as.character(args[5])
args[6] = as.character(args[6])
color_man = NULL
# args[1] = '/lvdata/wzb/scRNA/FW2024-225_04/Recluster/NK/celltype/group/meta.xls'
# args[2] = 'cell_type'
# args[3] = 'group'
# args[4] = '/lvdata/wzb/scRNA/FW2024-225_04/Recluster/NK/celltype/group'
# args[5] = '/lvdata/wzb/scRNA/FW2024-225_04/order.xls'
# args[6] = '/lvdata/wzb/scRNA/FW2024-225_04/Recluster/NK/celltype/group/cell_typecolor_dict.xls'

data = read.table(args[1], sep = '\t', header = T, row.names = 1, check.names = F)
if (args[5] != "None"){
  group_order = read.table(args[5],sep="\t",header= T)[,1,drop=FALSE]
  }else{group_order = NULL}
cellnum <- table(data[,args[2]], data[,args[3]])
cell.prop <- as.data.frame(prop.table(cellnum))
colnames(cell.prop)<-c("Celltype","Group", "Proportion")
write.table(cell.prop,paste0(args[4],"/cell_proportion.tsv"),sep = "\t",quote = F,row.names = F)
if (!is.null(group_order)){
  cell.prop$Group = factor(cell.prop$Group,levels= group_order[,1])
  }

n = length(unique(cell.prop$Celltype))
if(n<10){
  mycolor = pal_npg()(10)
} else {
  mycolor = colorRampPalette(colors = pal_npg()(10))(n)
}
if (args[6] != 'None'){
  cell_type_order = read.csv(args[6],sep="\t",header= T)
  if(ncol(cell_type_order)==1){
  cell.prop$Celltype = factor(cell.prop$Celltype,levels= cell_type_order[,1])
  }
  if(ncol(cell_type_order)==2){
  cell.prop$Celltype = factor(cell.prop$Celltype,levels= cell_type_order[,1])
  color_man = cell_type_order[,2]
  }
  }

if(!is.null(color_man)){
  print(color_man)
  p = ggplot(cell.prop,aes(Group, Proportion, fill=Celltype))+
    geom_bar(stat="identity",position="fill")+
    ggtitle("cell proportion")+
    theme_classic()+
    scale_fill_manual(values=color_man)+
    theme(
      # panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),
          axis.text.x=element_text(color='black', angle=45, hjust=1))
    guides(fill=guide_legend(title=NULL))


}else{
  print('None color')
  p = ggplot(cell.prop,aes(Group, Proportion, fill=Celltype))+
  geom_bar(stat="identity",position="fill")+
  ggtitle("cell proportion")+
  theme_classic()+
  scale_fill_manual(values=rev(mycolor))+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),
        axis.text.x=element_text(color='black', angle=45, hjust=1))
  guides(fill=guide_legend(title=NULL))


}
  p1 = p+
  # facet_wrap(~Group,ncol=4)+
    coord_polar(theta="y",direction=-1,clip="off")+   theme(axis.text.x=element_blank(),axis.line = element_blank()) 
  # p2 = ggplot(cell.prop,aes(Proportion, Celltype, fill=Group))+
  # geom_bar(stat="identity",position="stack")+
  # ggtitle("cell proportion")+
  # theme_classic()+
  # scale_fill_manual(values=rev(mycolor))+
  # theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),
  #       axis.text.x=element_text(color='black', angle=45, hjust=1))
  # guides(fill=guide_legend(title=NULL))
  cell.prop_normalized <- cell.prop %>%
  group_by(Group) %>%
  mutate(Proportion_normalized = Proportion / sum(Proportion))
  cell.prop_normalized$Celltype = factor(cell.prop_normalized$Celltype,levels= rev(levels(cell.prop_normalized$Celltype)))
  # 绘制归一化的堆积条形图
  p2 = ggplot(cell.prop_normalized, aes(Proportion_normalized,Celltype , fill=Group)) +
    geom_bar(stat="identity", position="fill") +
    ggtitle("cell proportion") +
    theme_classic() +
    scale_fill_manual(values=rev(mycolor)) +
    theme(panel.border = element_rect(fill=NA, color="black", size=0.5, linetype="solid"),
          axis.text.x=element_text(color='black', angle=90, hjust=1)) +
    guides(fill=guide_legend(title=NULL))+geom_vline(xintercept  = c(0.25, 0.5, 0.75), linetype = "dashed", color = "blue", size = 1)
    # guides(fill=guide_legend(title=NULL))+geom_vline(xintercept  = c(0.33, 0.66), linetype = "dashed", color = "blue", size = 1)
    # guides(fill=guide_legend(title=NULL))+geom_vline(xintercept  = c(0.5), linetype = "dashed", color = "blue", size = 1)
ggsave(paste(args[4],"cell_proportion_cluster.png",sep = "/"), plot=p, width = 8, height = 6)
ggsave(paste(args[4],"cell_proportion_cluster.pdf",sep = "/"), plot=p, width = 8, height = 6)
ggsave(paste(args[4],"cell_proportion_cluster2.png",sep = "/"), plot=p2, width = 8, height = 6)
ggsave(paste(args[4],"cell_proportion_cluster2.pdf",sep = "/"), plot=p2, width = 8, height = 6)
ggsave(paste(args[4],"cell_proportion_pie_cluster.png",sep = "/"), plot=p1, width = 8, height = 6)
ggsave(paste(args[4],"cell_proportion_pie_cluster.pdf",sep = "/"), plot=p1, width = 8, height = 6)


