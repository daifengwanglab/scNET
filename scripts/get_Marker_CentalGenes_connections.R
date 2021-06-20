# Regulatory network format;: 3 columns (TF target score)
rm(list=ls())
source('~/work/scNET-devel/scripts/load_libraries.R')
source('~/work/scNET-devel/scripts/read_data.R')
source('~/work/scNET-devel/scripts/functions_for_network_analysis.R')

cent.mat=read.table("centrality_matrix.all-cellTypes.txt", header=T)
betweenness.mat=dplyr::select(cent.mat,contains(c("gene","betweenness")))
cent.mat=betweenness.mat


markers=read.table("~/work/scNET_manuscript/data/marker_genes.txt", header=T)
markers$Cluster=gsub("Oligo","Oli",markers$Cluster)
markers$Cluster=gsub("Microglia","Mic",markers$Cluster)


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
  ct.markers=markers[markers$Cluster==tag[[1]][1],]
  tbl=rbind(cent.df,ct.markers)
  tbl
}

meanShortPathDf=data.frame(mean="",ct="")
for (i in 2:ncol(cent.mat))
{
  ct=paste(colnames(cent.mat)[i])
  print(ct)
  ct=paste(ct,"tbl",sep=".")
  df=cent.mat[,c(1,i)]
  assign(ct, get_top_cent.genes(df,0.01))
  df=get(ct)
  ct=levels(factor(df$Cluster))[1]
  cent=levels(factor(df$Cluster))[2]
  gs.markers=df[df$Cluster==ct,]$Gene
  gs.cent=df[df$Cluster==cent,]$Gene
  pattern=paste(ct,"_openchrom_en_GRN.network",sep="")
  net=get(pattern)
  net1=net[net$TFbs %in% "promoter", ]
  net2=net[net$TFbs %in% "both", ]
  net=rbind(net1,net2)
  net=net[,c("TF","TG","abs_coef")]
  net=distinct(net)
  net.igraph=graph_from_data_frame(net, directed = TRUE, vertices = NULL)
  distMatrix = shortest.paths(net.igraph, v=V(net.igraph),to=V(net.igraph))
  indx1=match(gs.markers,colnames(distMatrix))
  indx1=indx1[!is.na(indx1)]
  indx2=match(gs.cent,rownames(distMatrix))
  indx2=indx2[!is.na(indx2)]
  centGenes.markers.shortestPaths=distMatrix[indx2,indx1]
  mean=mean(centGenes.markers.shortestPaths)
  meanDF=data.frame(mean,ct)
  no.of.CentGenes=dim(centGenes.markers.shortestPaths)[1]
  start_time <- Sys.time()
  for (j in 1:100)
  {
    random.cent.genes=sample(unique(net$TG),no.of.CentGenes)
    indx2=match(random.cent.genes,rownames(distMatrix))
    indx2=indx2[!is.na(indx2)]
    RandomCentGenes.markers.shortestPaths=distMatrix[indx2,indx1]
    tmp=data.frame(mean(RandomCentGenes.markers.shortestPaths),paste("random",ct,sep="."))
    colnames(tmp)=c("mean","ct")
    meanDF=rbind(meanDF,tmp)
    end_time <- Sys.time()
    print(end_time - start_time)
  }
  meanShortPathDf=rbind(meanShortPathDf,meanDF)
}

tmp=meanShortPathDf[meanShortPathDf$ct %like% "Ex1",]
p=ggplot(tmp[2:nrow(tmp),], aes(x=as.numeric(mean)))+ geom_density()+geom_vline(xintercept=as.numeric(tmp[1,]$mean), size=1.5, color="red")
