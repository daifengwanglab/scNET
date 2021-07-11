# Regulatory network format;: 3 columns (TF target score)

source('~/work/scNET-devel/scripts/functions_for_network_analysis.R')

  EdgeDensity.tbl=data.frame(celltype=NULL, Score=NULL, Metric=NULL)
  Transitivity.tbl=data.frame(celltype=NULL, Score=NULL, Metric=NULL)
  AvgDeg.tbl=data.frame(celltype=NULL, Score=NULL, Metric=NULL)

nets=ls(pattern="*\\.network")
for(i in 1:length(nets))
  {
    name=nets[i]
    #tag=gsub("_openchrom_en_GRN.network","",name)
    tag=gsub(".network","",name)
    net=as.data.frame(lapply(nets[i],get))
    #TF.net1=net[net$TFbs %in% "promoter", ]
    #TF.net2=net[net$TFbs %in% "both", ]
    #TF.net=rbind(TF.net1,TF.net2)
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
  }

topology.tbl=rbind(EdgeDensity.tbl,Transitivity.tbl,AvgDeg.tbl)

topology.tbl$state = ifelse(topology.tbl$celltype %like% "AD","AD", "Healthy")

topology.tbl$ct = ifelse(topology.tbl$celltype %like% "Ex","Ex",
                  ifelse(topology.tbl$celltype %like% "In","In",
                  ifelse(topology.tbl$celltype %like% "Mic","Mic", "Oli")))

topology.tbl$celltype=gsub(".network","",topology.tbl$celltype)


p.topology=ggplot(topology.tbl, aes(x=Score, y=factor(celltype), fill=state))+
  geom_bar(stat="identity",position = "dodge")+
  scale_fill_npg() +
  facet_grid(ct ~ Metric, scales="free")+labs(x="Score",y="state") +
  theme_classic()+theme(text = element_text(size = 8),axis.text.x = element_text(angle = 90))


#ggsave(p.topology,file="Figures/topology.pdf",width = 4, height = 3, dpi = 300, units = "in", device='tiff')


write.table(topology.tbl, file="globalTopology.txt",col.names=TRUE,row.names=FALSE, sep="\t",quote=F)
