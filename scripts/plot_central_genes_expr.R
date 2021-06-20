

rm(list=ls())
source('~/work/scNET-devel/scripts/load_libraries.R')
source('~/work/scNET-devel/scripts/functions_for_network_analysis.R')

load("data/gematrix/gexpr_All.rdata")

cent.mat=read.table("centrality_matrix.all-cellTypes.txt", header=T)


#function
get_top_cent.genes = function(cent.df,n) #eg n=0.1 for top 10%
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
  ct.markers=markers[markers$Cluster==tag[[1]][1],]
  tbl=rbind(cent.df,ct.markers)
  tbl
}


for (i in 2:ncol(cent.mat))
{
  ct=paste(colnames(cent.mat)[i])
  ct=paste(ct,"tbl",sep=".")
  df=cent.mat[,c(1,i)]
  assign(ct, get_top_cent.genes(df,0.1))
  upset.tbl=rbind(upset.tbl,get(ct))
}

upset.tbl$score=1
upset.tbl=distinct(upset.tbl)

Mic=upset.tbl[upset.tbl$Cluster %like% "Mic",]
tmp=acast(Mic, Gene~Cluster, value.var="score" )
tmp[is.na(tmp)] <- 0
m1=make_comb_mat(tmp)
m1=m1[comb_size(m1) >= 4]
mic.upset=UpSet(m1)

Oli=upset.tbl[upset.tbl$Cluster %like% "Oli",]
tmp=acast(Oli, Gene~Cluster, value.var="score" )
tmp[is.na(tmp)] <- 0
m1=make_comb_mat(tmp)
m1=m1[comb_size(m1) >= 4]
oli.upset=UpSet(m1)

Ex1=upset.tbl[upset.tbl$Cluster %like% "Ex1",]
tmp=acast(Ex1, Gene~Cluster, value.var="score" )
tmp[is.na(tmp)] <- 0
m1=make_comb_mat(tmp)
m1=m1[comb_size(m1) >= 4]
Ex1.upset=UpSet(m1)
