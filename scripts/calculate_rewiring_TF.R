rm(list=ls())

source('~/work/scNET-devel/scripts/load_libraries.R')
source('~/work/scNET-devel/scripts/read_data.R')
source('~/work/scNET-devel/scripts/functions_for_network_analysis.R')


celltypes=c("Mic","Oli","Ex","In")

###############################################


  for(i in 1:length(celltypes))
  {
    TF.rewiring.table=data.frame(TF=NULL,score=NULL,cell=NULL)
    pattern=paste(celltypes[i],"network",sep=".")
    nets=ls(pattern=pattern)
    net=as.data.frame(lapply(nets[1],get))
    TF.net=net[,c("TF","TG","mse")]
    TF.net1=distinct(TF.net)
    net=as.data.frame(lapply(nets[2],get))
    TF.net=net[,c("TF","TG","mse")]
    TF.net2=distinct(TF.net)
    allTFs=union(TF.net1$TF,TF.net2$TF)
    for(j in 1:length(allTFs))
    {
        qTF.net1.targets=TF.net1[TF.net1$TF %in% allTFs[j],]
        qTF.net2.targets=TF.net2[TF.net2$TF %in% allTFs[j],]
        intersection=length(intersect(qTF.net1.targets$TG,qTF.net2.targets$TG))
        union=length(union(qTF.net1.targets$TG,qTF.net2.targets$TG))
        rs=1-(intersection/union)
        tbl=data.frame(TF=allTFs[j],score=rs,cell=celltypes[i])
        TF.rewiring.table=rbind(TF.rewiring.table,tbl)
    }
    name=paste(celltypes[i],"rewiring.tbl",sep=".")
    assign(name,TF.rewiring.table)
  }
  source('~/work/scNET-devel/scripts/calculate_h_metric.R')
