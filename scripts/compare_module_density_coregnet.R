rm(list=ls())
source('~/work/scNET-devel/scripts/load_libraries.R')
source('~/work/scNET-devel/scripts/read_data.R')
source('~/work/scNET-devel/scripts/functions_for_network_analysis.R')

#fig 2a in https://www.nature.com/articles/s12276-020-00528-0#Sec2
#fig 1c in https://www.nature.com/articles/s41467-019-10591-5/figures/1

#Read GO data
data=GSA.read.gmt('/Users/chiraggupta/work/scNET_manuscript/genome/genesets/BaderLab/Human_GO_bp_with_GO_iea_symbol.gmt')
genesets=data$genesets
names(genesets)=data$geneset.descriptions

#number of random nets needed
nrandnets=1
set.seed(123)

nets=ls(pattern="*\\.network")
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
  mat[mat < 0.5] <-0 #remove edges with less than 10% overlap
  name=gsub(".mat",".modules", name)
  assign(name, detect_modules(mat))
}

density.no.Modules.perturbed=data.frame(cell=NULL,pos=NULL,neg=NULL)
Modules.df=data.frame(Modulename=NULL,cell=NULL)

list=ls(pattern="*\\.modules")
list=Filter(function(x) !any(grepl("AD", x)), list)
for(i in 1:length(list))
{
  Module.density.tbl=data.frame("Modulename"=NULL,"Ctrl"=NULL,"AD"=NULL,lfc=NULL)
  cell=gsub(".network.JI.coreg.modules","",list[i])
  cell=gsub("Ctrl.","",cell)
  df=get(list[i])
  allmodules=unique(df$moduleID)
  for (j in 1:length(allmodules))
  {
    module=df[df$moduleID %in% allmodules[j],]$gene

    matname=paste("Ctrl",cell,sep=".")
    matname=paste(matname,".network.JI.coreg.mat",sep="")
    mat=get(matname)
    indx.c=match(module,colnames(mat))
    indx.r=match(module,rownames(mat))
    indx.c=indx.c[!is.na(indx.c)]
    indx.r=indx.r[!is.na(indx.r)]
    module.mat=mat[indx.r,indx.c]
    g=graph_from_adjacency_matrix(module.mat,weighted=TRUE, diag=FALSE, mode='undirected')
    df.2= get.data.frame(igraph::simplify(g,remove.multiple = TRUE, remove.loops = TRUE))
    e.density.Ctrl=sum(df.2$weight)/length(unique(module))

    matname=paste("AD",cell,sep=".")
    matname=paste(matname,".network.JI.coreg.mat",sep="")
    mat=get(matname)
    indx.c=match(module,colnames(mat))
    indx.r=match(module,rownames(mat))
    indx.c=indx.c[!is.na(indx.c)]
    indx.r=indx.r[!is.na(indx.r)]
    module.mat=mat[indx.r,indx.c]
    g=graph_from_adjacency_matrix(module.mat,weighted=TRUE, diag=FALSE, mode='undirected')
    df.2= get.data.frame(igraph::simplify(g,remove.multiple = TRUE, remove.loops = TRUE))
    e.density.AD=sum(df.2$weight)/length(unique(module))

    df.3=data.frame("Modulename"=allmodules[j],"Ctrl"=e.density.Ctrl,"AD"=e.density.AD)
    df.3$lfc = log2(df.3[,colnames(df.3)%like% "AD"]/ df.3[,colnames(df.3)%like% "Ctrl"])
    Module.density.tbl=rbind(Module.density.tbl,df.3)
  }
  colnames(Module.density.tbl)=c("Modulename","Ctrl","AD","lfc")
  Module.density.tbl$Cell=cell
  name=paste(cell,"Module.density.tbl",sep=".")

  d=Module.density.tbl
  d$abslfc=abs(d$lfc)
  d=d[order(-d$lfc),]
  assign(name,d)
  pos=nrow(d[d$lfc > 0,])
  neg=nrow(d[d$lfc < 0,])
  if(pos > 0 | neg > 0)
  {
    tmp.df=data.frame(Modulename=rownames(d[d$abslfc > 0,]), cell= cell)
    Modules.df=rbind(Modules.df,tmp.df)
  }
  tmp.df=data.frame(cell=cell,pos= pos,neg=neg)
  density.no.Modules.perturbed=rbind(density.no.Modules.perturbed,tmp.df)
  name=paste(cell,"module.density.lfc.tbl",sep=".")
  d$cell=cell
  d=d[,c("lfc","cell")]
  assign(name,d)

}

df_to_plot=melt(density.no.Modules.perturbed)


npgcolors=pal_npg("nrc", alpha = 1)(10)
p.no_density_Modules=ggplot(df_to_plot,aes(x=cell,y=value)) +
  geom_col(aes(fill=variable),width=0.5)+
  labs(y="# of Modules",x="Cell types")+
  scale_fill_manual(values=c("pos"=npgcolors[5],"neg"=npgcolors[6]),name=expression(Delta ~ "coregulation" ),labels=c("gain","loss"))+
 theme_bw(base_size=12)+theme(legend.position="top")
#ggsave(p.no_density_Modules,filename="Figures/p.no_density_Modules.pdf", device="pdf",width=3,height=3,units="in")



#box plot
for_boxplot=rbind(Mic.module.density.lfc.tbl,Oli.module.density.lfc.tbl)
for_boxplot=rbind(for_boxplot,Ex.module.density.lfc.tbl)
for_boxplot=rbind(for_boxplot,In.module.density.lfc.tbl)

npgcolors=pal_npg("nrc", alpha = 1)(10)
p.lfc_density_module_boxplot=ggplot(for_boxplot,aes(x=cell,y=lfc)) +
  geom_boxplot()+
  labs(y="change in module edge density",x="Cell types")+
 theme_bw(base_size=12)+theme(legend.position="top")
#ggsave(p.lfc_density_module_boxplot,filename="Figures/p.lfc_density_module_boxplot.pdf", device="pdf",width=3,height=3,units="in")



#plot synapse assembly subnets

#extract synapse assembly genes from GO
synapse=as.data.frame(genesets["synapse assembly"])$synapse.assembly

indx.c=match(synapse,colnames(Ctrl.Mic.network.JI.coreg.mat))
indx.c=indx.c[!is.na(indx.c)]
indx.r=match(synapse,rownames(Ctrl.Mic.network.JI.coreg.mat))
indx.r=indx.r[!is.na(indx.r)]

synapse.coregnet.mat=Ctrl.Mic.network.JI.coreg.mat[indx.r,indx.c]
g=graph.adjacency(synapse.coregnet.mat,weighted=TRUE)
df <- get.data.frame(igraph::simplify(g,remove.multiple = TRUE, remove.loops = TRUE))
#df=df[df$weight >= 0.3,]
synapse.coregnet.df=df
#write.table(synapse.coregnet.df,file="synapse_assembly.mic.ctrl.coregnet.dat",row.names=F,
#col.names=T,sep="\t",quote=FALSE)

#select synapse assembly genes in AD mic net
indx.c=match(synapse,colnames(AD.Mic.network.JI.coreg.mat))
indx.c=indx.c[!is.na(indx.c)]

indx.r=match(synapse,rownames(AD.Mic.network.JI.coreg.mat))
indx.r=indx.r[!is.na(indx.r)]

synapse.coregnet.mat=AD.Mic.network.JI.coreg.mat[indx.r,indx.c]
g=graph.adjacency(synapse.coregnet.mat,weighted=TRUE)
df <- get.data.frame(igraph::simplify(g,remove.multiple = TRUE, remove.loops = TRUE))
#df=df[df$weight >= 0.3,]
synapse.coregnet.df=df
#write.table(synapse.coregnet.df,file="synapse_assembly.mic.AD.coregnet.dat",row.names=F,
#col.names=T,sep="\t",quote=FALSE)
