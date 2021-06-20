rm(list=ls())
source('~/work/scNET-devel/scripts/load_libraries.R')
source('~/work/scNET-devel/scripts/read_data.R')
source('~/work/scNET-devel/scripts/functions_for_network_analysis.R')

nets=ls(pattern="*_openchrom_en_GRN.network")
for(i in 1:length(nets))
{
  name=nets[i]
  tag=gsub("_openchrom_en_GRN.network","",name)
  name=paste(tag,"coregnet.igraph",sep=".")
  net=as.data.frame(lapply(nets[i],get))
  TF.net1=net[net$TFbs %in% "promoter", ]
  TF.net2=net[net$TFbs %in% "both", ]
  TF.net=rbind(TF.net1,TF.net2)
  TF.net=net[,c("TF","TG","mse")]
  TF.net=distinct(TF.net)
  assign(name, get_coregnet_graph(TF.net,0.1,tag))
}
