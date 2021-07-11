args = commandArgs(trailingOnly=TRUE) #network_file celltype(e.g: Ex1) state

#args=c("~/scGRN/MIT_AD_analysis/Ex/MIT.AD.Ex.grn.txt","Ex","AD")

library(Loregic)
library(dplyr)
library(ArrayBin)
library(data.table)

load("/ua/cgupta8/scGRN/data/gematrix/MIT_imputed_gexpr.rdata")

if (args[3] == "AD")
{
AD_cells.tmp=AD_cells[,c("TAG","broad.cell.type")]
AD_cells.tmp$new.broad.cell.type = sub('[.]', '.', make.names(AD_cells.tmp$broad.cell.type, unique=TRUE))
gexpr_AD.tmp=gexpr_AD
gexpr_AD.tmp$newrows=rownames(gexpr_AD.tmp)
table=setDT(gexpr_AD.tmp)
lookup=setDT(AD_cells.tmp)
setkey(gexpr_AD.tmp, "newrows")
setkey(AD_cells.tmp, "TAG")
joined=gexpr_AD.tmp[AD_cells.tmp]
joined.df=as.data.frame(joined)
rownames(joined.df)=joined.df$new.broad.cell.type
joined.df=joined.df[,1:ncol(gexpr_AD)]
}

if (args[3] == "Ctrl")
{
CTL_cells.tmp=CTL_cells[,c("TAG","broad.cell.type")]
CTL_cells.tmp$new.broad.cell.type = sub('[.]', '.', make.names(CTL_cells.tmp$broad.cell.type, unique=TRUE))
gexpr_CTL.tmp=gexpr_CTL
gexpr_CTL.tmp$newrows=rownames(gexpr_CTL.tmp)
table=setDT(gexpr_CTL.tmp)
lookup=setDT(CTL_cells.tmp)
setkey(gexpr_CTL.tmp, "newrows")
setkey(CTL_cells.tmp, "TAG")
joined=gexpr_CTL.tmp[CTL_cells.tmp]
joined.df=as.data.frame(joined)
rownames(joined.df)=joined.df$new.broad.cell.type
joined.df=joined.df[,1:ncol(gexpr_CTL)]

}

cell.joined.df=joined.df[rownames(joined.df) %like% args[2], ]
gexpr= t(cell.joined.df[,log10(colSums(cell.joined.df)+1)> 1])


data=read.table(args[1], header=T, sep="\t",stringsAsFactors=FALSE)
data=data[data$mse<0.1 & data$abs_coef > 0.01,]
TF.net1=data[data$TFbs %in% "promoter", ]
TF.net2=data[data$TFbs %in% "both", ]
TF.net=rbind(TF.net1,TF.net2)
TF.net=TF.net[,c("TF","TG","mse")]
TF.net=distinct(TF.net)
TF.net=TF.net[,1:2]

nodes=append(unique(TF.net$TF),unique(TF.net$TG))

rm(data,TF.net1,TF.net2)

indx=match(nodes,rownames(gexpr))
indx=indx[!is.na(indx)]
nodes.ct.gexpr=gexpr[indx,]


nodes.ct.gexpr.bin=binarize.array(nodes.ct.gexpr[,1:100])

statetag=paste(args[2],args[3],sep=".")
filename=paste(statetag,"gexpr.bin.mat",sep=".")

write.table(data.frame("Gene"=rownames(nodes.ct.gexpr.bin),nodes.ct.gexpr.bin),file=filename, row.names=FALSE,quote=F, sep="\t")

print("running loregic")
ct.loregic.out=loregic(TF.net,nodes.ct.gexpr.bin)
print("finished loregic")

statetag=paste(args[2],args[3],sep=".")

filename=paste(statetag,"loregic.out",sep=".")

write.table(ct.loregic.out,file=filename, col.names=TRUE, row.names=FALSE, sep="\t", quote=F)
