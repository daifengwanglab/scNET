#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE) #network_file celltype (e.g: Ex1)

library(Loregic)
library(dplyr)
library(ArrayBin)
#library(tidyverse)

load("~/work/scNET_manuscript/data/Ting_Jin/AD_MIT/MIT_imputed_gexpr.rdata")

celltypes=c("Mic","Oli","Ex","In")

allmat=gexpr_AD
table(rownames(allmat)==allmat$TAG)
rownames(allmat)=allmat$new.broad.cell.type


for (i in 1:length(celltypes))
{
  tag=celltypes[i]
  gexpr=allmat[rownames(allmat) %like% tag, ]
  #gexpr = t(gexpr[,log10(colSums(gexpr)+1)> 1])
  gexpr=t(gexpr)
  variances = apply(X=gexpr, MARGIN=2, FUN=var)
  sorted = sort(variances, decreasing=TRUE, index.return=TRUE)
  sorted.indx=sorted$ix[1:100]
  gexpr.highvariance =   gexpr[, sorted.indx]
  name=paste(tag,"gexpr_highvar",sep=".")
  assign(name,gexpr.highvariance)
}

mat1=merge(Mic.gexpr_highvar,Oli.gexpr_highvar,by="row.names")
mat2=merge(Ex.gexpr_highvar,In.gexpr_highvar,by="row.names")
mat=merge(mat1,mat2,by="Row.names")
rownames(mat)=mat$Row.names
mat$Row.names=NULL

#################
net=read.table(args[1], header=T, sep="\t", stringsAsFactors=FALSE)
#data=read.table("data/scGRN_openchrom_en/Mic_openchrom_en_GRN.csv", header=T, sep=",", stringsAsFactors=FALSE)
net=net[net$mse<0.1 & net$abs_coef > 0.01,]

net=net %>% arrange(-abs_coef)
net=net[,c("TF","TG","abs_coef")]
net= net %>% distinct()
net=net[,1:2]
nodes=append(unique(net$TF),unique(net$TG))


indx=match(nodes,rownames(mat))
indx=indx[!is.na(indx)]
nodes.ct.gexpr=mat[indx,]

nodes.ct.gexpr.bin=binarize.array(nodes.ct.gexpr)
ct.loregic.out=loregic(net,nodes.ct.gexpr.bin)

filename=paste(args[2],"loregic.out",sep=".")
write.table(ct.loregic.out,file=filename, col.names=TRUE, row.names=FALSE, sep="\t", quote=F)

filename=paste(args[2],"gexpr.bin.mat",sep=".")
write.table(data.frame("Gene"=rownames(nodes.ct.gexpr.bin),nodes.ct.gexpr.bin),file=paste("results",filename,sep="/"), row.names=FALSE,quote=F, sep="\t")
