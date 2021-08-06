

source('~/work/scNET-devel/scripts/load_libraries.R')
source('~/work/scNET-devel/scripts/read_data.R')
source('~/work/scNET-devel/scripts/functions_for_network_analysis.R')


celltypes=c("Mic","Oli","Ex","In")

###############################################

  nets=ls(pattern="*\\.network")
  pcounter.tbl=data.frame(pcounter=NULL, cell=NULL)

  for(i in 1:length(nets))
  {
    name=nets[i]
    tag=gsub(".network","",name)
    name=paste(tag,"hie.ht",sep=".")
    distTable=data.frame(gene=NULL,h=NULL)
    distTblName=paste(name, "randomDist",sep=".")
    net=as.data.frame(lapply(nets[i],get))
    TF.net=net[,c("TF","TG","mse")]
    TF.net=distinct(TF.net)
    TF.net=TF.net[TF.net$TG %in% TF.net$TF,]
    TF.net.igraph=graph_from_data_frame(TF.net, directed = TRUE, vertices = NULL)
    out.deg=get_centrality(TF.net,"degree_out",tag)
    in.deg=get_centrality(TF.net,"degree_in",tag)
    df=merge(out.deg,in.deg,by="gene",all=TRUE)
    colnames(df)=c("gene","out.d","in.d")
    #df$o_minus_i=df$out.d-df$in.d
    #df$o_plus_i=df$out.d+df$in.d
    df$h=(df$out.d-df$in.d)/(df$out.d+df$in.d)
    df=df[,c("gene","h")]
    newhcolname=paste(tag,"h",sep=".")
    colnames(df)=c("gene",newhcolname)
    assign(name,df)
    pcounter=0
    for (j in 1:10)
    {
      TF.net.rand=erdos.renyi.game(length(V(TF.net.igraph)),length(E(TF.net.igraph)), type="gnm",directed = TRUE)
      TF.net.rand.df=as.data.frame(get.edgelist(TF.net.rand))
      colnames(TF.net.rand.df)=c("TF","TG")
      out.deg=get_centrality(TF.net.rand.df,"degree_out",tag)
      in.deg=get_centrality(TF.net.rand.df,"degree_in",tag)
      df.r=merge(out.deg,in.deg,by="gene",all=TRUE)
      colnames(df.r)=c("gene","out.d","in.d")
      df.r$h=(df.r$out.d-df.r$in.d)/(df.r$out.d+df.r$in.d)
      df.r=df.r[,c("gene","h")]
      distTable=rbind(distTable,df.r)
      assign(distTblName,distTable)
      pval=ks.test(df.r$h,df[,2],alternative="less")$p.value
      if(pval < 0.05) {pcounter=pcounter+1}
    }
    tmp=data.frame(pcounter=pcounter,cell=tag)
    pcounter.tbl=rbind(pcounter.tbl,tmp)
  }

  #plot rand dist of all ct AD networks
  rand.all.cts.dist=data.frame(gene=NULL,h=NULL,cell=NULL,network=NULL)
  for (i in 1:length(celltypes))
  {
    pattern=paste(celltypes[i],".hie.ht.randomDist",sep="")
    pattern=paste("AD",pattern,sep=".")
    df=get(pattern)
    df$cell=celltypes[i]
    df$network="random"
    rand.all.cts.dist=rbind(rand.all.cts.dist,df)
  }

  AD.all.cts.dist=data.frame(gene=NULL,h=NULL,cell=NULL,network=NULL)
  for (i in 1:length(celltypes))
  {
    pattern=paste(celltypes[i],".hie.ht",sep="")
    pattern=paste("AD",pattern,sep=".")
    df=get(pattern)
    colnames(df)=c("gene","h")
    df$cell=celltypes[i]
    df$network="AD"
    AD.all.cts.dist=rbind(AD.all.cts.dist,df)
  }

  Ctrl.all.cts.dist=data.frame(gene=NULL,h=NULL,cell=NULL,network=NULL)
  for (i in 1:length(celltypes))
  {
    pattern=paste(celltypes[i],".hie.ht",sep="")
    pattern=paste("Ctrl",pattern,sep=".")
    df=get(pattern)
    colnames(df)=c("gene","h")
    df$cell=celltypes[i]
    df$network="Ctrl"
    Ctrl.all.cts.dist=rbind(Ctrl.all.cts.dist,df)
  }


  p.rand.all.ct.h.dist=ggplot(rand.all.cts.dist, aes(x=h))+geom_histogram()+
  labs(x="hierarchy height",y="Count")+ facet_wrap(~cell) +ggtitle("Random Networks")+
  annotate("rect", xmin=-1,xmax=-0.33,ymin=-0,ymax=Inf, alpha=0.2, fill="red")+
  annotate("rect", xmin=-0.33,xmax=0.33,ymin=-0,ymax=Inf, alpha=0.2, fill="green")+
  annotate("rect",xmin=0.33,xmax=1,ymin=-0,ymax=Inf, alpha=0.2, fill="blue")+
  theme(text = element_text(size = 10)) +
  theme_bw(base_size=12)

  pdf(file="Figures/p.rand.all.ct.h.dist.pdf")
  p.rand.all.ct.h.dist
  dev.off()

  p.AD.all.ct.h.dist=ggplot(AD.all.cts.dist, aes(x=h))+geom_histogram()+
  labs(x="hierarchy height",y="Count")+ facet_wrap(~cell) +ggtitle("AD Networks")+
  annotate("rect", xmin=-1,xmax=-0.33,ymin=-0,ymax=Inf, alpha=0.2, fill="red")+
  annotate("rect", xmin=-0.33,xmax=0.33,ymin=-0,ymax=Inf, alpha=0.2, fill="green")+
  annotate("rect",xmin=0.33,xmax=1,ymin=-0,ymax=Inf, alpha=0.2, fill="blue")+
  theme(text = element_text(size = 10)) +
  theme_bw(base_size=12)
  pdf(file="Figures/p.AD.all.ct.h.dist.pdf")
  p.AD.all.ct.h.dist
  dev.off()

  p.Ctrl.all.ct.h.dist=ggplot(Ctrl.all.cts.dist, aes(x=h))+geom_histogram()+
  labs(x="hierarchy height",y="Count")+ facet_wrap(~cell) +ggtitle("Control Networks")+
  annotate("rect", xmin=-1,xmax=-0.33,ymin=-0,ymax=Inf, alpha=0.2, fill="red")+
  annotate("rect", xmin=-0.33,xmax=0.33,ymin=-0,ymax=Inf, alpha=0.2, fill="green")+
  annotate("rect",xmin=0.33,xmax=1,ymin=-0,ymax=Inf, alpha=0.2, fill="blue")+
  theme(text = element_text(size = 10)) +
  theme_bw(base_size=12)
  pdf(file="Figures/p.Ctrl.all.ct.h.dist.pdf")
  p.Ctrl.all.ct.h.dist
  dev.off()



#TF net sankey
meanFCTbl=data.frame(Cell=NULL,Top=NULL, Middle=NULL, Bottom=NULL)
no.of.enh.tbl=data.frame(Cell=NULL,Promoter=NULL, Enhancer=NULL, Level=NULL, totalTF=NULL)
meanRWscore.tbl=data.frame(Cell=NULL,Top=NULL, Middle=NULL, Bottom=NULL)
for (i in 1:length(celltypes))
{
    pattern=paste(celltypes[i],"hie.ht",sep=".")
    list=ls(pattern=pattern)
    list=Filter(function(x) !any(grepl(".randomDist", x)), list)
    df1=get(list[1])
    df2=get(list[2])
    df=merge(df1,df2,by="gene")
    all.deg=read_xlsx("~/work/scNET_manuscript/data/gematrix/Diff.Exp.Genes.DataS2.MIT.xlsx",sheet=celltypes[i],skip=1)[,1:9]
    ct.deg=as.data.frame(all.deg[,c(1,5)])
    colnames(ct.deg)=c("gene","FC")
    ct.deg=ct.deg[ct.deg$gene %in% df$gene,]
    df.deg.hei=merge(ct.deg,df,by="gene")
    df.deg.hei$AD = cut(df.deg.hei[,3],c(-Inf,-0.33,0.33,1),labels=c("Bottom","Middle","Top"))
    df.deg.hei$Control = cut(df.deg.hei[,4],c(-Inf,-0.33,0.33,1),labels=c("Bottom","Middle","Top"))
    Topmean=mean(df.deg.hei$FC[df.deg.hei$AD=="Top"])
    Middlemean=mean(df.deg.hei$FC[df.deg.hei$AD=="Middle"])
    Bottomean=mean(df.deg.hei$FC[df.deg.hei$AD=="Bottom"])
    fctbl=data.frame(Cell=celltypes[i],Top=Topmean,Middle=Middlemean,Bottom=Bottomean)
    meanFCTbl=rbind(meanFCTbl,fctbl)
    network.name=paste("AD",celltypes[i],sep=".")
    network.name=paste(network.name,"network",sep=".")
    net=get(network.name)

    ### AD networks
    genes=df.deg.hei[df.deg.hei$AD=="Top",]$gene
    genes.prom=net[net$TF %in% genes & net$TFbs=='promoter',]
    TopNprom=length(unique(genes.prom$promoter))/length(unique(genes.prom$TF))
    genes.enh=net[net$TF %in% genes & net$TFbs=='enhancer',]
    TopNenh=length(unique(genes.enh$enhancer))/length(unique(genes.enh$TF))
    tbl=data.frame(Cell=celltypes[i],Promoter=TopNprom,Enhancer=TopNenh, Level="Top",totalTF=length(unique(genes)))
    no.of.enh.tbl=rbind(no.of.enh.tbl,tbl)

    genes=df.deg.hei[df.deg.hei$AD=="Middle",]$gene
    genes.prom=net[net$TF %in% genes & net$TFbs=='promoter',]
    TopNprom=length(unique(genes.prom$promoter))/length(unique(genes.prom$TF))
    genes.enh=net[net$TF %in% genes & net$TFbs=='enhancer',]
    TopNenh=length(unique(genes.enh$enhancer))/length(unique(genes.enh$TF))
    tbl=data.frame(Cell=celltypes[i],Promoter=TopNprom,Enhancer=TopNenh, Level="Middle",totalTF=length(unique(genes)))
    no.of.enh.tbl=rbind(no.of.enh.tbl,tbl)

    genes=df.deg.hei[df.deg.hei$AD=="Bottom",]$gene
    genes.prom=net[net$TF %in% genes & net$TFbs=='promoter',]
    TopNprom=length(unique(genes.prom$promoter))/length(unique(genes.prom$TF))
    genes.enh=net[net$TF %in% genes & net$TFbs=='enhancer',]
    TopNenh=length(unique(genes.enh$enhancer))/length(unique(genes.enh$TF))
    tbl=data.frame(Cell=celltypes[i],Promoter=TopNprom,Enhancer=TopNenh, Level="Bottom",totalTF=length(unique(genes)))
    no.of.enh.tbl=rbind(no.of.enh.tbl,tbl)


    #rewiring
    TF.rewiring.table=data.frame(gene=NULL,RWscore=NULL,cell=NULL)
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
        tbl=data.frame(gene=allTFs[j],RWscore=rs,cell=celltypes[i])
        TF.rewiring.table=rbind(TF.rewiring.table,tbl)
    }
    name=paste(celltypes[i],"rewiring.tbl",sep=".")
    assign(name,TF.rewiring.table)
    df=merge(df.deg.hei,TF.rewiring.table,by="gene")

    TopRWmean=mean(df$RWscore[df$AD=="Top"])
    MiddleRWmean=mean(df$RWscore[df$AD=="Middle"])
    BottomRWmean=mean(df$RWscore[df$AD=="Bottom"])
    RWtbl=data.frame(Cell=celltypes[i],Top=TopRWmean,Middle=MiddleRWmean,Bottom=BottomRWmean)
    meanRWscore.tbl=rbind(meanRWscore.tbl,RWtbl)


  #   plot sankey
  #  p.df = df.deg.hei %>% make_long(colnames(df.deg.hei[,5:6]))
  #  p=ggplot(p.df, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = factor(node), label = node)) +
  #    geom_sankey(flow.alpha = .6,
  #                node.color = "gray30") +
  #    scale_fill_manual(values = c("red", "green", "blue"), name="Hierarchy Levels") +
  #    labs(x=NULL)+
  #    theme_sankey(base_size = 12) +
  #    ggtitle(celltypes[i])+theme(plot.title = element_text(hjust = 0.30,vjust=0.9))
    name=paste("p",pattern,sep=".")
    name=paste(name,".sankey.TFNet",sep="")
    name=paste(name,"pdf",sep=".")
    name=paste("Figures",name,sep="/")
  #  ggsave(p,filename=name, device="pdf" )
}


npgcolors=pal_npg("nrc", alpha = 1)(10)
plotData=melt(meanRWscore.tbl)
p.rewiring=ggplot(plotData,aes(color=Cell, y=value, x=reorder(variable))) +
  geom_point(size=2)+labs(y="Average rewiring score",x="Levels")+
  coord_flip()+
  scale_color_manual(values=c("In"=npgcolors[1],"Ex"=npgcolors[2],"Mic"=npgcolors[3],"Oli"=npgcolors[4]))+
 theme_bw(base_size=12)+theme(axis.text.x=element_text(angle=90))+theme(legend.position = "none")
ggsave(p.rewiring,filename="Figures/p.rewiring.pdf", device="pdf",width=3,height=3,units="in")


plotData=melt(meanFCTbl)
p.foldchange=ggplot(plotData,aes(color=Cell,size=2, y=value, x=reorder(variable))) +
geom_point(size=2)+ labs(y="Average fold change",x="Levels")+
coord_flip()+
  scale_color_manual(values=c("In"=npgcolors[1],"Ex"=npgcolors[2],"Mic"=npgcolors[3],"Oli"=npgcolors[4]))+
  theme_bw(base_size=12) +
  theme(axis.text.x=element_text(angle=90))+theme(legend.position = "none")

ggsave(p.foldchange,filename="Figures/p.foldchange.pdf", device="pdf",width=3,height=3,units="in")


p.promoter=ggplot(no.of.enh.tbl,aes(color=Cell, y=Promoter, x=Level)) +
  geom_point(size=2)+ labs(y="Average # of promoters targeted",x="Levels")+
  coord_flip()+
  scale_color_manual(values=c("In"=npgcolors[1],"Ex"=npgcolors[2],"Mic"=npgcolors[3],"Oli"=npgcolors[4])) +
  theme_bw(base_size=12) +
  theme(axis.text.x=element_text(angle=90))+theme(legend.position = "none")
ggsave(p.promoter,filename="Figures/p.promoter.pdf", device="pdf",width=3.2,height=3,units="in" )

p.enhancer=ggplot(no.of.enh.tbl,aes(color=Cell, y=Enhancer, x=Level)) +
  geom_point(size=2)+
  coord_flip()+
  scale_color_manual(values=c("In"=npgcolors[1],"Ex"=npgcolors[2],"Mic"=npgcolors[3],"Oli"=npgcolors[4]))+
  theme_bw(base_size=12) + labs(y="Average # of enhancers targeted",x="Levels")+
  theme(axis.text.x=element_text(angle=90))+theme(legend.position = "none")
ggsave(p.enhancer,filename="Figures/p.enhancer.pdf", device="pdf",width=3.2,height=3,units="in" )

#legend
p.ct.legend=ggplot(no.of.enh.tbl,aes(color=Cell, y=Enhancer, x=Level)) +
  geom_point(size=2)+
  coord_flip()+
  scale_color_manual(values=c("In"=npgcolors[1],"Ex"=npgcolors[2],"Mic"=npgcolors[3],"Oli"=npgcolors[4]))+
  theme_bw(base_size=12) + labs(y="Average # of enhancers targeted",x="Levels")+
  theme(axis.text.x=element_text(angle=90))
ggsave(p.ct.legend,filename="Figures/p.ct.legend.pdf", device="pdf",width=3.2,height=3,units="in" )
