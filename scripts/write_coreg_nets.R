rm(list=ls())
source('~/work/scNET-devel/scripts/load_libraries.R')
source('~/work/scNET-devel/scripts/read_data.R')
source('~/work/scNET-devel/scripts/functions_for_network_analysis.R')


nets=ls(pattern="*\\.network")

th=0.2

for(i in 1:length(nets))
{
      name=nets[i]
      celltype=gsub(".network","",name)
      net=as.data.frame(lapply(nets[i],get))
      net=net[,c("TF","TG","abs_coef")]
      net=distinct(net)
      name=paste(name,"coreg.txt",sep=".")
      assign(name, get_coregnet_graph(net,th))
      write.table(get(name),file=name,col.names=T, row.names=F, sep="\t",quote=F)
}
