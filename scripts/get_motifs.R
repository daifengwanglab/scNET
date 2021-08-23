
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
motif.enrichment.tbl=data.frame(Motif=NULL, zscore=NULL,cell=NULL,state=NULL)
for(i in 1:length(nets))
{
  name=nets[i]
  tag=gsub(".network","",name)
  name=paste(tag,"motifs",sep=".")
  net=as.data.frame(lapply(nets[i],get))
  TF.net=net[,c("TF","TG","mse")]
  TF.net=distinct(TF.net)
  TF.net=TF.net[TF.net$TG %in% TF.net$TF,]
  TF.net.igraph=graph_from_data_frame(TF.net, directed = TRUE, vertices = NULL)
  #net.igraph=graph_from_data_frame(TF.net, directed = TRUE, vertices = NULL)
  motifs=igraph::motifs(TF.net.igraph)
  motifs.count=as.data.frame(motifs)
  colnames(motifs.count)=c("Count")
  motifs.count$Motif=rownames(motifs.count)
  motifs.count=motifs.count[,2:1]
  for(j in 1:1000)
  {
    TF.net.r=TF.net
    TF.net.r$TF=sample(TF.net.r$TF)
    TF.net.rand=graph_from_data_frame(TF.net.r, directed = TRUE, vertices = NULL)
  #  TF.net.rand=erdos.renyi.game(length(V(TF.net.igraph)),length(E(TF.net.igraph)), type="gnm",directed = TRUE)
    r.motifs=igraph::motifs(TF.net.rand)
    r.motifs.count=as.data.frame(r.motifs)
    colnames(r.motifs.count)=paste("Count_random",j,sep=".")
    r.motifs.count$Motif=rownames(r.motifs.count)
    r.motifs.count=r.motifs.count[,2:1]
    motifs.count=inner_join(motifs.count,r.motifs.count)
  }
  rand=motifs.count[,3:ncol(motifs.count)]
  motifs.count=motifs.count[,1:2]
  z.df=data.frame(Motif=NULL, zscore=NULL)
  for (k in 1:nrow(rand))
  {
    sd=sd(rand[k,])
    mean=rowMeans(rand[k,],na.rm=T)
    z=(motifs.count[k,2] - mean)/sd
    m.z.df=data.frame(Motif=k,zscore=z)
    z.df=rbind(z.df,m.z.df)
  }
  z.df$cell=tag
  motif.enrichment.tbl=rbind(motif.enrichment.tbl,z.df)

}
 motif.enrichment.tbl$state=ifelse(motif.enrichment.tbl$cell %like% "AD","AD","Ctrl")
  motif.enrichment.tbl$cell=gsub("Ctrl.","",motif.enrichment.tbl$cell)
  motif.enrichment.tbl$cell=gsub("AD.","",motif.enrichment.tbl$cell)

p=ggplot(motif.enrichment.tbl,aes(x=as.character(Motif),y=zscore,fill=state))+geom_bar(stat="identity",position="dodge")+facet_wrap(~cell,ncol=1)

ggsave(p,filename="Figures/p.motifs.pdf", device="pdf",width=4,height=6,units="in")
