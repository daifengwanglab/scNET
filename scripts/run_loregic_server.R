#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE) #network_file celltype (e.g: Ex1)

library(Loregic)
library(dplyr)
library(ArrayBin)
#library(tidyverse)

load("data/gematrix/gexpr_All.rdata")

gexpr=t(gexpr)

#data=read.table(args[1], header=T, sep=",", stringsAsFactors=FALSE)
data=read.table("data/scGRN_openchrom_en/Mic_openchrom_en_GRN.csv", header=T, sep=",", stringsAsFactors=FALSE)
data=data %>% arrange(-abs_coef)
data=data[,c("TF","TG","abs_coef")]
tmp= data %>% distinct()
tmp=tmp[1:50000,]
data=tmp[,1:2]
nodes=append(unique(data$TF),unique(data$TG))

ct.gexpr=dplyr::select(as.data.frame(gexpr),contains("Mic"))
#ct.gexpr=dplyr::select(as.data.frame(gexpr),contains(args[2]))
indx=match(nodes,rownames(ct.gexpr))
indx=indx[!is.na(indx)]
nodes.ct.gexpr=ct.gexpr[indx,]

#calculate variances in each sample
variances = apply(X=nodes.ct.gexpr, MARGIN=2, FUN=var)
sorted = sort(variances, decreasing=TRUE, index.return=TRUE)
sorted.indx=sorted$ix[1:10]
gexpr.highvariance = nodes.ct.gexpr[, sorted.indx]
print(dim(gexpr.highvariance))

nodes.ct.gexpr.bin=binarize.array(gexpr.highvariance)
ct.loregic.out=loregic(data,nodes.ct.gexpr.bin)

filename=paste(args[2],"loregic.out",sep=".")
write.table(ct.loregic.out,file=paste("results",filename,sep="/"), col.names=TRUE, row.names=FALSE, sep="\t", quote=F)
filename=paste(args[2],"gexpr.bin.mat",sep=".")
write.table(data.frame("Gene"=rownames(nodes.ct.gexpr.bin),nodes.ct.gexpr),file=paste("results",filename,sep="/"), row.names=FALSE,quote=F, sep="\t")
