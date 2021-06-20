library(igraph)

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE) #realNet #number_of_random_nets #outdir
args=("real_networks/Ex1_openchrom_en_GRN.modified.for_snap.txt","Ex1",2,"random_networks/")
realnet=read.table(args[1], header=F, sep="\t")

net.igraph=graph_from_data_frame(realnet, directed=TRUE, vertices=NULL)

for (i in 1:args[2])
{
  randnet1=net.igraph %>% rewire(keeping_degseq(niter = 100))
  tmp=as.data.frame(get.edgelist(randnet1))
  head(tmp)
  outfilename=paste("random",i,sep=".")
  outfilename=paste(args[3],outfilename,sep="/")
  write.table(tmp, file=outfilename,col.names=FALSE,row.names=FALSE, sep="\t",quote=F)
}
