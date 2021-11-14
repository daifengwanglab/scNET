rm(list=ls())
source('~/work/scNET-devel/scripts/load_libraries.R')
source('~/work/scNET-devel/scripts/read_data.R')
source('~/work/scNET-devel/scripts/functions_for_network_analysis.R')


nets=ls(pattern="*\\.network")

th=0.2
msize=30

for(i in 1:length(nets))
{
      name=nets[i]
      celltype=gsub(".network","",name)
      net=as.data.frame(lapply(nets[i],get))
      net=net[,c("TF","TG","abs_coef")]
      net=distinct(net)
      name=paste(name,"JI.coreg.mat",sep=".")
      assign(name, find_target_pairs_matrix(net))
      mat=get(name)
      mat[mat < th] <-0 #remove edges with less than th% overlap
      name=gsub(".mat",".modules", name)
      tag=paste("JI",th,sep="_")
      name=paste(name,tag,sep=".")
      tag=paste("ModSize",msize,sep="_")
      name=paste(name,tag,sep=".")
      assign(name, detect_modules(mat,msize))
}


#plotting module 2 genes (mod 2 GO BP phospolipids)

mod2=AD.Mic.network.JI.coreg.modules.JI_0.2.ModSize_30[AD.Mic.network.JI.coreg.modules.JI_0.2.ModSize_30$moduleID %in% "2",]$gene
net=AD.Mic.network
net=net[,c("TF","TG","abs_coef")]
net=distinct(net)

node.attr=as.data.frame(unique(net$TF))
node.attr$type="TF"
colnames(node.attr)=c("Gene","Type")
write.table(node.attr,file="mod2.coregnet.attr",row.names=F,col.names=T,sep="\t",quote=FALSE)

mat=find_target_pairs_matrix(net)
mat[mat < 0.2] <-0
indx.c=match(mod2,colnames(mat))
indx.r=match(mod2,rownames(mat))
mod2.coregnet.mat=mat[indx.r,indx.c]
g=graph.adjacency(mod2.coregnet.mat,weighted=TRUE)
df <- get.data.frame(igraph::simplify(g,remove.multiple = TRUE, remove.loops = TRUE))

#select top edges to plot
df.filtered=df[df$weight > 0.4,]
write.table(df.filtered,file="mod2.coregnet.dat",row.names=F,col.names=T,sep="\t",quote=FALSE)
