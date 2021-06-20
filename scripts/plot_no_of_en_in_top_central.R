
#Dont run from command line; last part not executed.
rm(list=ls())
source('~/work/scNET-devel/scripts/load_libraries.R')
source('~/work/scNET-devel/scripts/functions_for_network_analysis.R')


cent.mat=read.table("results/centralities/centrality_matrix.En-networks.all-cellTypes.txt", header=T)

En.count.tbl= data.frame(matrix(ncol = 2, nrow = 0))
colnames(En.count.tbl)=c("Cluster","Count")


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
  tbl=get(ct)
  all.cent=dim(tbl)[1]
  tbl=tbl %>% filter(Gene %like% "chr")
  en.cent=dim(tbl)[1]
  frac=(en.cent/all.cent)*100
  ct=paste(colnames(cent.mat)[i])
  tmp=data.frame(Cluster=ct,Count=frac)
  En.count.tbl=rbind(En.count.tbl,tmp)
}
