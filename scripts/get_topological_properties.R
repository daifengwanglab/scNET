# Regulatory network format;: 3 columns (TF target score)

source('~/work/scNET-devel/scripts/functions_for_network_analysis.R')

  EdgeDensity.tbl=data.frame(celltype=NULL, Score=NULL, Metric=NULL)
  Transitivity.tbl=data.frame(celltype=NULL, Score=NULL, Metric=NULL)
  AvgDeg.tbl=data.frame(celltype=NULL, Score=NULL, Metric=NULL)

nets=ls(pattern="*\\.network")
for(i in 1:length(nets))
  {
    name=nets[i]
    tag=gsub("_openchrom_en_GRN.network","",name)
    net=as.data.frame(lapply(nets[i],get))
    TF.net1=net[net$TFbs %in% "promoter", ]
    TF.net2=net[net$TFbs %in% "both", ]
    TF.net=rbind(TF.net1,TF.net2)
    TF.net=net[,c("TF","TG","mse")]
    TF.net=distinct(TF.net)
    net.igraph=graph_from_data_frame(TF.net, directed = TRUE, vertices = NULL)

    ed=edge_density(net.igraph, loops = FALSE)
    df=data.frame(celltype=tag, Score=ed, Metric="Edge Density")
    EdgeDensity.tbl=rbind(EdgeDensity.tbl,df)

    trans = transitivity(net.igraph,weights = E(net.igraph)$weight)
    df=data.frame(celltype=tag, Score=trans, Metric="Transitivity")
    Transitivity.tbl=rbind(Transitivity.tbl,df)

    deg=igraph::degree(net.igraph)
    avgdeg=sum(deg)/length(deg)
    df=data.frame(celltype=tag, Score=avgdeg, Metric="Avg.Degree")
    AvgDeg.tbl=rbind(AvgDeg.tbl,df)

    #randomize using a simple method; optimize later
    #random_ed = c()
    #random_tr = c()
    #random_deg =c()
    #for(j in 1:5)
    #{
      #random=TF.net
      #random[] <- lapply(random, sample)
    #  vl_random_net = sample_degseq(deg, method = "vl")
    #  random_ed = c(random_ed, edge_density(vl_random_net))
    #}
  }

topology.tbl=rbind(EdgeDensity.tbl,Transitivity.tbl,AvgDeg.tbl)

topology.tbl$state = ifelse(topology.tbl$celltype %like% "AD","AD", "Healthy")

topology.tbl$ct = ifelse(topology.tbl$celltype %like% "Ex","Ex",
                  ifelse(topology.tbl$celltype %like% "In","In",
                  ifelse(topology.tbl$celltype %like% "Mic","Mic", "Oli")))

topology.tbl$celltype=gsub(".network","",topology.tbl$celltype)
multifacet= ggplot(topology.tbl, aes(x=Score, y=factor(celltype)))+
  geom_bar(stat="identity",width=.5,position = position_dodge(width = 0.1))+
  facet_grid(ct ~ Metric, scales="free")+labs(x="Score",y="cells")
pdf(file="topology.pdf")
multifacet
dev.off()

write.table(topology.tbl, file="globalTopology.txt",col.names=TRUE,row.names=FALSE, sep="\t",quote=F)

#source('~/work/scNet-devel/scripts/disease_topology.R')
#topology.tbl.subset=topology.tbl[topology.tbl$celltype %in% c("Ex1","In1a","Mic","Oli"),]
#table=rbind(topology.tbl.subset,topology.AD.tbl)
#multifacet=ggplot(table, aes(x=Score, y=celltype))+geom_bar(stat="identity")+facet_grid(state ~ Metric , scales="free")
#healthy=ggplot(topology.tbl.subset, aes(x=Score, y=celltype))+geom_bar(stat="identity")+facet_grid(~ Metric, scales="free")
#disease=ggplot(topology.AD.tbl, aes(x=Score, y=celltype))+geom_bar(stat="identity")+facet_grid(~ Metric, scales="free")
