#!/usr/bin/env Rscript
#args = commandArgs(trailingOnly=TRUE) #loregic_out_file celltype.state (e.g: Ex.AD)
rm(list=ls())

library(tidyr)
library(Loregic)
library(dplyr)

celltypes=c("Mic.AD","Oli.AD","Ex.AD","In.AD","Mic.Ctrl","Oli.Ctrl","Ex.Ctrl","In.Ctrl")

for (k in 1:length(celltypes))
{

cellname=paste("Loregic",celltypes[k],sep="/")
cellname=paste(cellname,"loregic.out",sep=".")

loregicOut=read.table(cellname, header=T, sep="\t")

tmp=loregicOut %>% unite(x, c(RF1,RF2,target), sep="_")
rownames(tmp)=tmp$x
tmp=tmp[,-1]
mat=tmp[apply(tmp[,], 1, function(x) !all(x==0)),] # remove gates with no scores
colnames(mat)=c("T=0","AND","RF1*~RF2","RF1","~RF1*RF2","RF2","XOR","OR","NOR","XNOR","~RF2","RF1+~RF2","~RF1","~RF1+RF2","NAND","T=1")

l=apply(mat,1,function(x) which(x==max(x))) #store max of each row in a list

DF=data.frame(matrix(NA, nrow = 1, ncol = 16))
colnames(DF)=colnames(mat)

#for each list, check of the number of items in more than 2; if yes then gate-inconsistent
for(i in 1:length(l))
{
    nGates=length(l[[i]])
    if(nGates<2)
    {
      DF=rbind(DF,mat[i,]) #store all gate consistent in a DF df
    }
}

DF=DF[-1,] #stores all gate consistent triplets
tmp=DF
tmp$triplet=rownames(tmp)
tmp=tmp %>% separate(triplet, c("RF1","RF2","target"),sep="_")
rownames(tmp)=c()
gate_consistent_trips=tmp

#count occurences per gate
#df=data.frame(matrix(NA, nrow = 1, ncol = 1))
#colnames(df)=c("Count")
df=data.frame(Count=NULL)
for (j in 1:ncol(DF))
{
  count=as.data.frame(dim(DF[which(DF[,j]>0),])[1])
  rownames(count)=colnames(DF)[j]
  colnames(count)="Count"
  df=rbind(df,count)
}

df$cell=celltypes[k]
df$gate=rownames(df)
name=paste(celltypes[k],"logics.txt",sep=".")
assign(name,df)

}

 df=rbind(Ex.AD.logics.txt,In.AD.logics.txt)
 df=rbind(df,Mic.AD.logics.txt)
 df=rbind(df,Oli.AD.logics.txt)

 df2=rbind(Ex.Ctrl.logics.txt,In.Ctrl.logics.txt)
 df2=rbind(df2,Mic.Ctrl.logics.txt)
 df2=rbind(df2,Oli.Ctrl.logics.txt)

 df=rbind(df,df2)
df$state=ifelse(df$cell %like% "AD","AD","Ctrl")
df$cell=gsub("\\.AD","",df$cell)
df$cell=gsub("\\.Ctrl","",df$cell)
 p=ggplot(df,aes(x=gate,y=Count,fill=state))+geom_bar(stat="identity",position="dodge")+facet_wrap(~cell,ncol=1)+
 theme(axis.text.x=element_text(angle=90))

 ggsave(p,filename="Figures/p.logics.pdf", device="pdf",width=3,height=3,units="in")
