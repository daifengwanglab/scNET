
#https://igraph.org/r/doc/triad_census.html


#guide to claculate for every node, number of motifs it participates in
#https://stackoverflow.com/questions/12374534/how-to-mine-for-motifs-in-r-with-igraph
# Regulatory network format;: 3 columns (TF target score)
rm(list=ls())
source('~/work/scNET-devel/scripts/load_libraries.R')
source('~/work/scNET-devel/scripts/read_data.R')
source('~/work/scNET-devel/scripts/functions_for_network_analysis.R')

#source('../scNET-devel/scripts/convert_ct_nets_to_igraph.R')

nets=ls(pattern="*\\.network")
for(i in 1:length(nets))
{
  name=nets[i]
  tag=gsub(".network","",name)
  name=paste(tag,"motifs",sep=".")
  net=as.data.frame(lapply(nets[i],get))
  TF.net=net[,c("TF","TG","mse")]
  TF.net=distinct(TF.net)
  net.igraph=graph_from_data_frame(TF.net, directed = TRUE, vertices = NULL)
  assign(name, igraph::motifs(net.igraph))
}

motif_list=ls(pattern="*.motifs")
for(i in 1:length(motif_list))
{
  name=motif_list[i]
  tag=gsub("_counts.motifs","",name)
  name=paste(tag,"motif_counts",sep=".")
  count=as.data.frame(lapply(motif_list[i],get))
  colnames(count)=c("Count")
  #count=count[,c("MotifId","Count")]
  count$celltype=tag
  count$motifID=rownames(count)
  assign(name, count)
}

count_list=ls(pattern="*.motif_counts")
df=as.data.frame(lapply(count_list[1],get))
for(i in 2:length(count_list))
{
  df=rbind(df,as.data.frame(lapply(count_list[i],get)))
}
tmp=xtabs(Count~motifID+celltype, data=df)
tmp=tmp[rowSums(tmp[])>0,]
df=melt(tmp)

#remove motifID 7
#df=df[df$motifID !=7,]
#df=df[df$motifID !=3,]

p=ggplot(df,aes(x=value,y=celltype, fill=as.character(motifID))) + geom_bar(position="fill", stat="identity")
