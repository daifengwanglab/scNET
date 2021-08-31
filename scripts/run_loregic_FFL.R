data=read.table('MIC.Ctrl.FFL.txt')
entrez=read.table('~/scGRN/data/entrez_to_hgnc_mapping.txt',header=F)
 colnames(entrez)=c("symbol","id")
 colnames(data)=c("RF1","RF2","target")
 new=data
 new[] = lapply(data, function(x) entrez$symbol[match(x, entrez$id)])
 library(dplyr)

 new=distinct(new)

 library(Loregic)
  library(dplyr)
 library(ArrayBin)
 library(data.table)

 load("/ua/cgupta8/scGRN/data/gematrix/MIT_imputed_gexpr.rdata")

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


 cell.joined.df=joined.df[rownames(joined.df) %like% "Mic", ]
 gexpr= t(cell.joined.df[,log10(colSums(cell.joined.df)+1)> 1])



 nodes=append(unique(as.character(new$RF1)),unique(as.character(new$RF2)))
nodes=append(unique(nodes),unique(as.character(new$target)))
nodes=unique(nodes)
indx=match(nodes,rownames(gexpr))
 indx=indx[!is.na(indx)]
 nodes.ct.gexpr=gexpr[indx,]
 nodes.ct.gexpr.bin=binarize.array(nodes.ct.gexpr[,1:100])
ct.loregic.out=loregic(new,nodes.ct.gexpr.bin)
write.table(data.frame("Gene"=rownames(nodes.ct.gexpr.bin),nodes.ct.gexpr.bin),file="Mic.ffl.ctrl.gexpr.bin.mat", row.names=FALSE,quote=F, sep= "\t")
write.table(ct.loregic.out,file="Mic.ffl.ctrl.loregic.out", col.names=TRUE, row.names=FALSE, sep="\t", quote=F)
