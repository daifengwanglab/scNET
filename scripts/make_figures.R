#source('')


#make PR NetworksMatrix
pr.list=list(ex.pr,in.pr,oligo.pr,micro.pr)
df=Reduce(function(x, y) merge(x, y, all=TRUE, by="gene"), pr.list, accumulate=FALSE)
df[is.na(df)]=0


dat=read.table("pr.scaled.mat-Human_GO_bp_with_GO_iea_symbol.mt3lt500.fdr-based.zscore.mat",header=T,row.names=1)
dat=dat[,3:6]
pheatmap(dat,cellwidth=10,cellheight=10,show_rownames=TRUE,filename="hubscore.pdf")
