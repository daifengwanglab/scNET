
#read disease data

  EdgeDensity.tbl=data.frame(celltype=NULL, Score=NULL, Metric=NULL, state=NULL)
  Transitivity.tbl=data.frame(celltype=NULL, Score=NULL, Metric=NULL,state=NULL)
  AvgDeg.tbl=data.frame(celltype=NULL, Score=NULL, Metric=NULL, state=NULL)

nets=ls(pattern="*_GENIE3_GRN.network")

for(i in 1:length(nets))
  {
    name=nets[i]
    tag=gsub("_GENIE3_GRN.network","",name)
    TF.net=as.data.frame(lapply(nets[i],get))
    TF.net=TF.net[order(-TF.net$weight),]
    TF.net=TF.net[1:50000,]
    TF.net=distinct(TF.net)
    net.igraph=graph_from_data_frame(TF.net, directed = TRUE, vertices = NULL)

    ed=edge_density(net.igraph, loops = FALSE)
    df=data.frame(celltype=tag, Score=ed, Metric="Edge Density",state="AD")
    EdgeDensity.tbl=rbind(EdgeDensity.tbl,df)

    trans = transitivity(net.igraph,weights = E(net.igraph)$weight)
    df=data.frame(celltype=tag, Score=trans, Metric="Transitivity",state="AD")
    Transitivity.tbl=rbind(Transitivity.tbl,df)

    deg=igraph::degree(net.igraph)
    avgdeg=sum(deg)/length(deg)
    df=data.frame(celltype=tag, Score=avgdeg, Metric="Avg.Degree", state="AD")
    AvgDeg.tbl=rbind(AvgDeg.tbl,df)
  }

 topology.AD.tbl=rbind(EdgeDensity.tbl,Transitivity.tbl,AvgDeg.tbl)


topology.tbl.subset=topology.tbl[topology.tbl$celltype %in% c("Ex1","In1a","Mic","Oli"),]
