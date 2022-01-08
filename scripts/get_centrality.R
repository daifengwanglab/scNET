

source('../scripts/load_libraries.R')
source('../scripts/read_data.R')
source('../scripts/functions_for_network_analysis.R')



celltypes=c("Mic","Oli","Ex","In")

###############################################
#calculate centralities
centralities=c("betweenness","degree_out","degree_in")

for (k in 1:length(centralities))
{
  c=centralities[k]
  nets=ls(pattern="*\\.network")
  for(i in 1:length(nets))
  {
    name=nets[i]
    tag=gsub(".network","",name)
    name=paste(tag,c,sep=".")
    net=as.data.frame(lapply(nets[i],get))
    TF.net=net[,c("TF","TG","mse")]
    TF.net=distinct(TF.net)
    assign(name, get_centrality(TF.net,c,tag))
  }
}

rm(list=nets)

for(i in 1:length(celltypes))
{
  tag=celltypes[i]
  pattern=paste(tag,"\\.",sep="")
  list=ls(pattern=pattern)
  data = Reduce(function(x, y) merge(x, y, all=T), lapply(list,get), accumulate=F)
  name=paste(tag,"centrality_matrix",sep=".")
  assign(name,data)
}

pattern="*.centrality_matrix"
list=ls(pattern=pattern)

old=get(list[1])
for (i in 2:length(list))
{
  old=merge(old, get(list[i]), by="gene", all=TRUE)
}

write.table(old, file="results/centrality_matrix.all-cellTypes.txt",col.names=TRUE,row.names=FALSE, sep="\t",quote=F)
