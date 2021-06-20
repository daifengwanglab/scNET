
#########################################

#############################
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



library(ggplot2)
library(plyr) # might be not needed here anyway it is a must-have package I think in R
library(reshape2) # to "melt" your dataset
library (scales) # it has a "rescale" function which is needed in heatmaps
library(RColorBrewer) # for convenience of heatmap colors, it reflects your mood sometimes
nba <- read.csv("http://datasets.flowingdata.com/ppg2008.csv")
nba <- as.data.frame(cor(nba[2:ncol(nba)])) # convert the matrix correlations to a dataframe
nba.m <- data.frame(row=rownames(nba),nba) # create a column called "row"
rownames(nba) <- NULL #get rid of row names
nba <- melt(nba)
nba.m$value<-cut(nba.m$value,breaks=c(-1,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1),include.lowest=TRUE,label=c("(-0.75,-1)","(-0.5,-0.75)","(-0.25,-0.5)","(0,-0.25)","(0,0.25)","(0.25,0.5)","(0.5,0.75)","(0.75,1)")) # this can be customized to put the correlations in categories using the "cut" function with appropriate labels to show them in the legend, this column now would be discrete and not continuous
nba.m$row <- factor(nba.m$row, levels=rev(unique(as.character(nba.m$variable)))) # reorder the "row" column which would be used as the x axis in the plot after converting it to a factor and ordered now
#now plotting
ggplot(nba.m, aes(row, variable)) +
geom_tile(aes(fill=value),colour="black") +
scale_fill_brewer(palette = "RdYlGn",name="Correlation")  # here comes the RColorBrewer package, now if you ask me why did you choose this palette colour I would say look at your battery charge indicator of your mobile for example your shaver, won't be red when gets low? and back to green when charged? This was the inspiration to choose this colour set.


color.scheme <- rev(brewer.pal(8,"Greys")) # generate the color scheme to use
pdf(file="clustet_toy3.pdf")
heatmap.2(kmedoids.cor, Rowv = NULL, Colv = NULL,
          dendrogram = "none",
          trace = "none", density.info = "none",
          col = color.scheme, key = FALSE,
          labRow = FALSE, labCol = FALSE)


dev.off()

dat=read.table('loregic.out')
colnames(dat)=c("T=0","AND","RF1*~RF2","RF1","~RF1*RF2","RF2","XOR","OR","NOR","XNOR","~RF2","RF1+~RF2","~RF1","~RF1+RF2","NAND","T=1")
melt=melt(dat, rownames(dat))

library(ochRe)
ggplot(melt, aes(y=ct,fill = variable,x=value))+geom_bar(stat = "identity") + scale_fill_manual(values=values$emu_woman_paired) + theme_minimal()
