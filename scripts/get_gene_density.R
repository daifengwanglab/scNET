rm(list=ls())
source('~/work/scNET-devel/scripts/load_libraries.R')
source('~/work/scNET-devel/scripts/read_data.R')
source('~/work/scNET-devel/scripts/functions_for_network_analysis.R')

nets=ls(pattern="*\\.network")
for(i in 1:length(nets))
{
  name=nets[i]
  celltype=gsub(".network","",name)
  filename=paste(celltype,"density.tbl",sep=".")
  net=as.data.frame(lapply(nets[i],get))
  TF.net1=net[net$TFbs %in% "promoter", ]
  TF.net2=net[net$TFbs %in% "both", ]
  TF.net=rbind(TF.net1,TF.net2)
  TF.net=net[,c("TF","TG","mse")]
  TF.net=distinct(TF.net)
  name=paste(celltype,"coreg.net",sep=".")
  df= get_coregnet_graph(TF.net,0)
  nodes=unique(append(unique(df$g1), unique(df$g2)))
  print(paste(celltype,length(nodes),sep=":"))
  node.density.df=data.frame("gene"=NULL, "density"=NULL)
  for(j in 1:length(nodes))
  {
    node.net = df[df$g1 %in% nodes[j] | df$g2 %in% nodes[j],]
    net.igraph=graph_from_data_frame(node.net, directed = FALSE, vertices = NULL)
    net.igraph = set_edge_attr(net.igraph, "weight", value= node.net$jaccard)
    density=edge_density(net.igraph)
    tmp=data.frame("gene"=nodes[j], "density"=density)
    node.density.df=rbind(node.density.df,tmp)
  }
  colnames(node.density.df)=c("gene",celltype)
  assign(filename,node.density.df)
}

list=ls(pattern=".density.tbl")
celltypes=c("Mic","Oli","Ex","In")

data = Reduce(function(x, y) merge(x, y, all=T), lapply(list,get), accumulate=F)
data[is.na(data)]=0
rownames(data)=data$gene
data$gene=NULL
min=min(data[data > 0])
data[data == 0]=(min*0.01)

for (i in 1: length(celltypes))
{
  df=data[,colnames(data) %like% celltypes[i],]
  df = transform(df, lfc = log2(df[,colnames(df)%like% "AD"]/ df[,colnames(df)%like% "Ctrl"]))
  colnames(df)=c(colnames(df)[1],colnames(df)[2],paste(celltypes[i],"lfc",sep="_"))
  df$gene=rownames(df)
  rownames(df)=NULL
  filename=paste(celltypes[i],"density.lfc.mat",sep=".")
  assign(filename,df)
}

list=ls(pattern="density.lfc.mat")
data = Reduce(function(x, y) merge(x, y, all=T), lapply(list,get), accumulate=F)
rownames(data)=data$gene
data$gene=NULL
data=data[,colnames(data) %like% "_lfc"]
