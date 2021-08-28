
library(dplyr)


args = commandArgs(trailingOnly=TRUE) #network_file celltype(e.g: Ex1)

data=read.table(args[1], header=T, sep="\t",stringsAsFactors=FALSE)
#data=read.table("~/scGRN/MIT_AD_analysis/MIT.AD.Mic.grn.txt",header=T, , sep="\t",stringsAsFactors=FALSE)
net=data[data$mse<0.1 & data$abs_coef > 0.01,]
TF.net=net[,c("TF","TG","mse")]
TF.net=distinct(TF.net)
TF.net=TF.net[TF.net$TG %in% TF.net$TF,]

TF.net=TF.net[,1:2]

entrez=read.table('~/scGRN/data/entrez_to_hgnc_mapping.txt',header=F)
colnames(entrez)=c("symbol","id")

new = TF.net  # create a copy of df
# using lapply, loop over columns and match values to the look up table. store in "new".
new[] = lapply(TF.net, function(x) entrez$id[match(x, entrez$symbol)])
new=distinct(new)
new$w=1

filename=paste(args[2],"TF-TF.entrez.txt",sep=".")
write.table(new,file=filename,sep="\t",col.names=F,row.names=F,quote=F)
