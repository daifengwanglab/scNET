#source('')


#make PR NetworksMatrix
pr.list=list(ex.pr,in.pr,oligo.pr,micro.pr)
df=Reduce(function(x, y) merge(x, y, all=TRUE, by="gene"), pr.list, accumulate=FALSE)
df[is.na(df)]=0

library(pheatmap)
library(viridis)
dat=read.table("pr.mat-Human_GO_bp_with_GO_iea_symbol.fdr-based.zscore.mat",header=T,row.names=1)
dat=dat[,3:6]
mat.new=dat

paletteLength=1000
myColor <- colorRampPalette(c("yellow", "white", "blue"))(paletteLength)
colorv=viridis(direction=-1,n=1000)
myBreaks <- c(seq(min(mat.new), 0, length.out=ceiling(paletteLength/2) + 1),
+               seq(max(mat.new)/paletteLength, max(mat.new), length.out=floor(paletteLength/2)))

col=viridis(n=100)
pheatmap(mat.new,show_rownames=TRUE,filename=NA,col=colorv,cluster_cols=FALSE,border_color = "black")


pheatmap(dat,cellwidth=10,cellheight=10,show_rownames=TRUE,filename="hubscore.pdf")
