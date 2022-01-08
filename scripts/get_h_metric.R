
source('../scripts/load_libraries.R')
source('../scripts/read_data.R')
source('../scripts/functions_for_network_analysis.R')


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

  hei.height.matrix=rbind(Ctrl.all.cts.dist,AD.all.cts.dist)
  hei.height.matrix$level=cut(hei.height.matrix$h,c(-Inf,-0.33,0.33,1),labels=c("Bottom","Middle","Top"))
  write.table(hei.height.matrix,file="results/hei.height.matrix",sep="\t",col.names=T, row.names=F, quote=F)
