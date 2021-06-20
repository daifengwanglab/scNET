

rm(list=ls())
source('~/work/scNET-devel/scripts/load_libraries.R')
source('~/work/scNET-devel/scripts/functions_for_network_analysis.R')


cent.mat=read.table("results/centralities/centrality_matrix.En-networks.all-cellTypes.txt", header=T)
markers=read.table("~/work/scNET_manuscript/data/marker_genes.txt", header=T)
markers$Cluster=gsub("Oligo","Oli",markers$Cluster)
markers$Cluster=gsub("Microglia","Mic",markers$Cluster)

upset.tbl= data.frame(matrix(ncol = 2, nrow = 0))
colnames(upset.tbl)=c("Gene","Cluster")

#function
get_top_cent.genes = function(cent.df,n) #eg n=0.3 for 30%
{
  name=colnames(cent.df)[2]
  tag=strsplit(name,"[.]")
  colnames(cent.df)=c("gene","score")
  cent.df=cent.df[order(-cent.df$score), , drop = FALSE]
  topn=round(dim(cent.df)[1]*n)
  cent.df=cent.df[1:topn,]
  cent.df=cent.df[cent.df$score > 0 ,]
  cent.df$score=name
  colnames(cent.df)=c("Gene","Cluster")
  #ct.markers=markers[markers$Cluster==tag[[1]][1],]
  #tbl=rbind(cent.df,ct.markers)
  #tbl
  cent.df
}


for (i in 2:ncol(cent.mat))
{
  ct=paste(colnames(cent.mat)[i])
  ct=paste(ct,"tbl",sep=".")
  df=cent.mat[,c(1,i)]
  assign(ct, get_top_cent.genes(df,0.01))
  upset.tbl=rbind(upset.tbl,get(ct))
  tbl=get(ct)
  ct=paste(colnames(cent.mat)[i])
  ct=paste(ct,"list",sep=".")
  assign(ct,c(tbl$Gene))
}


upset.tbl$score=1
upset.tbl=distinct(upset.tbl)


write.table(upset.tbl[,1:2], file="results/gsea/En.centrality.genesets",col.names=FALSE,row.names=FALSE, sep="\t",quote=F)

#Run jaccrard in perl
rm(list=ls())

get_lower_tri<-function(cormat){
   cormat[upper.tri(cormat)] <- NA
   return(cormat)
 }
 # Get upper triangle of the correlation matrix
 get_upper_tri <- function(cormat){
   cormat[lower.tri(cormat)]<- NA
   return(cormat)
 }
 reorder_cormat <- function(cormat){

  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
  }

#source http://www.sthda.com/english/wiki/ggplot2-quick-correlation-matrix-heatmap-r-software-and-data-visualization

dat=read.table("/Users/chiraggupta/work/scNET_manuscript/results/gsea/Results/degree_in.En.for_heatmap", header=T)
dat$Factor=gsub(".betweenness","",dat$Factor)
dat$GeneSet=gsub(".betweenness","",dat$GeneSet)
tmp=acast(dat, Factor~GeneSet, value.var="Jaccard" )
tmp[is.na(tmp)] <- 0
mat=tmp
#mat <- reorder_cormat(mat)
upper_tri <- get_upper_tri(mat)
melted_cormat <- melt(upper_tri, na.rm = TRUE)
p=ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
 geom_tile(color = "white")+
 scale_fill_gradient2(low = "white", high = "black", mid = "white",
   midpoint = median(dat$Jaccard), limit = c(min(dat$Jaccard),max(dat$Jaccard)), space = "Lab",
   name="Jaccard\nIndex") +
  theme_minimal()+
 theme(axis.text.x = element_text(angle = 45, vjust = 1,
    size = 12, hjust = 1))+
 coord_fixed()+labs(x="",y="") +

theme(
  panel.grid.major = element_blank(),
  panel.border = element_blank(),
  panel.background = element_blank(),
  axis.ticks = element_blank(),
  legend.justification = c(1, 0),
  legend.position = c(0.4, 0.7),
  legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                title.position = "top", title.hjust = 0.5))

pdf(file="En.central_overlap_betweenness.pdf")
p
dev.off()
