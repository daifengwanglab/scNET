#!/usr/bin/env Rscript
rm(list=ls())
source('~/work/scNET-devel/scripts/load_libraries.R')
source('~/work/scNET-devel/scripts/functions_for_network_analysis.R')

datadir=c("logics/motifs/1000_rands_with_100_sampling/FFL_members_full_list")

celltypes=c("Mic.AD","Oli.AD","Ex.AD","In.AD","Mic.Ctrl","Oli.Ctrl","Ex.Ctrl","In.Ctrl")

for (k in 1:length(celltypes))
{

cellname=paste(celltypes[k],"gate_consistent_FLL_trips.txt",sep=".")

name=paste(datadir,cellname,sep="/")

loregicOut=read.table(name, header=T, sep="\t")
colnames(loregicOut)=c("T=0","AND","RF1*~RF2","~RF1*RF2","XOR","OR","NOR","XNOR","~RF2","RF1+~RF2","~RF1","~RF1+RF2","NAND","T=1","TF1","TF2","TF.target","pvalue")

#colnames(mat)=c("T=0","AND","RF1*~RF2","RF1","~RF1*RF2","RF2","XOR","OR","NOR","XNOR","~RF2","RF1+~RF2","~RF1","~RF1+RF2","NAND","T=1")


tmp=loregicOut %>% unite(x, c(TF1,TF2,TF.target), sep="_")
rownames(tmp)=tmp$x
tmp$x=NULL
#tmp=tmp[tmp$pvalue <= 0.01,]
tmp$pvalue=NULL
tmp$ct=celltypes[k]
name=paste(celltypes[k],"ffl_logics.pval_lt_pt1",sep=".")
assign(name, tmp)
}


 df=rbind(Ex.AD.ffl_logics.pval_lt_pt1,In.AD.ffl_logics.pval_lt_pt1)
 df=rbind(df,Mic.AD.ffl_logics.pval_lt_pt1)
 df=rbind(df,Oli.AD.ffl_logics.pval_lt_pt1)

 df2=rbind(Ex.Ctrl.ffl_logics.pval_lt_pt1,In.Ctrl.ffl_logics.pval_lt_pt1)
 df2=rbind(df2,Mic.Ctrl.ffl_logics.pval_lt_pt1)
 df2=rbind(df2,Oli.Ctrl.ffl_logics.pval_lt_pt1)

 df=rbind(df,df2)
 df=melt(df)
df$condition=ifelse(df$ct %like% "AD","AD","Ctrl")
df$ct=gsub("\\.AD","",df$ct)
df$ct=gsub("\\.Ctrl","",df$ct)
 p=ggplot(df,aes(x=as.character(variable),fill=condition))+
 geom_bar(stat="identity")+facet_wrap(~ct,ncol=1)+
 scale_fill_manual(values=c("AD"=npgcolors[1],"Ctrl"=npgcolors[2]))+
 theme_bw(base_size=12) + labs(y="Count",x="logic gates")+
 theme(axis.text.x=element_text(angle=90))+theme(legend.position = "top")

 ggsave(p,filename="Figures/p.logics.pdf", device="pdf",width=3,height=3,units="in")
