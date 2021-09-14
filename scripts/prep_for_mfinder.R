library(dplyr)


args = commandArgs(trailingOnly=TRUE) #network_file celltype(e.g: Ex1)

data=read.table(args[1], header=T, sep="\t",stringsAsFactors=FALSE)
#data=read.table("~/scGRN/MIT_AD_analysis/MIT.AD.Mic.grn.txt",header=T, , sep="\t",stringsAsFactors=FALSE)
net=data[data$mse<0.1 & data$abs_coef > 0.01,]
TF.net=net[,c("TF","TG","mse")]
TF.net=distinct(TF.net)
TF.net=TF.net[TF.net$TG %in% TF.net$TF,]
TF.net=TF.net[TF.net$mse >= mean((TF.net$mse)),]

TF.net=TF.net[,1:2]

tmp=as.data.frame(unique(append(TF.net$TF,TF.net$TG)))
tmp$ID=seq(1:nrow(tmp))
colnames(tmp)=c("symbol","id")
filename=paste(args[2],"node_symbol-integar.key.txt",sep=".")
write.table(tmp,file=filename,sep="\t",col.names=T,row.names=F,quote=F)

new = TF.net  # create a copy of df
# using lapply, loop over columns and match values to the look up table. store in "new".
new[] = lapply(TF.net, function(x) tmp$id[match(x, tmp$symbol)])
new=distinct(new)
new$w=1

filename=paste(args[2],"TF-TF.entrez.txt",sep=".")
write.table(new,file=filename,sep="\t",col.names=F,row.names=F,quote=F)
