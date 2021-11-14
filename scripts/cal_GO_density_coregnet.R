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


celltypes=c("Mic","Oli","Ex","In")

nets=ls(pattern="*\\.network")
for(i in 1:length(nets))
{
  name=nets[i]
  celltype=gsub(".network","",name)
  GO.density.tbl=data.frame("GOname"=NULL,"edgedensity"=NULL)
  net=as.data.frame(lapply(nets[i],get))
  net=net[,c("TF","TG","abs_coef")]
  net=distinct(net)
  mat=find_target_pairs_matrix(net)
  #mat[mat < 0.1] <-0 #remove edges with less than 10% overlap
  for (j in 2:length(genesets))
  {
    gsname=names(genesets)[j]
    if(gsname %in% "GO ID not found in OBO file ontology definitions")
    {
      next
    }
    if(length(genesets[[j]]) < 10 | length(genesets[[j]]) > 300 )
    {
      next
    }
    print(paste(gsname,length(genesets[[j]])),sep=":")
    gs=genesets[[j]]
    indx1=match(gs, rownames(mat))
    indx2=match(gs,colnames(mat))
    indx1=indx1[!is.na(indx1)]
    indx2=indx2[!is.na(indx2)]
    if(length(indx1) > 10 & length(indx1) < 300)
    {
      gs.mat=mat[indx1,indx2]
      g=graph_from_adjacency_matrix(gs.mat, weighted=TRUE, diag=FALSE, mode='undirected')
      df= get.data.frame(igraph::simplify(g,remove.multiple = TRUE, remove.loops = TRUE))
      tot=length(unique(gs))
      e.density=sum(df$weight)/(tot*(tot-1)/2)
      #e.density=edge_density(g, loops=FALSE)
      df=data.frame("GOname"=gsname,"edgedensity"=e.density)
      GO.density.tbl=rbind(GO.density.tbl,df)
    }
    else
    {
      df=data.frame("GOname"=gsname,"edgedensity"=0)
      GO.density.tbl=rbind(GO.density.tbl,df)
    }
  }
  colnames(GO.density.tbl)=c("GOname",celltype)
  name=paste(celltype,"GO.density.tbl",sep=".")
  assign(name,GO.density.tbl)
}

#for stacked barplot and boxplot
GOBP_terms.df=data.frame(GOname=NULL,cell=NULL)
density.no.BP.perturbed=data.frame(cell=NULL,pos=NULL,neg=NULL)
for (i in 1:length(celltypes))
{
  c.net=paste("Ctrl",celltypes[i],sep=".")
  c.net=paste(c.net, "GO.density.tbl",sep=".")
  a.net=paste("AD",celltypes[i],sep=".")
  a.net=paste(a.net, "GO.density.tbl",sep=".")
  df= merge(get(c.net),get(a.net),by="GOname")
  rownames(df)=df$GOname
  df$GOname=NULL
  df=df[rowSums(df) >0,]
  min=min(df[df > 0])
  df[df == 0]=(min*0.01)
  d = transform(df, lfc = log2(df[,colnames(df)%like% "AD"]/ df[,colnames(df)%like% "Ctrl"]))
  d$lfc=d$lfc
  d$abslfc=abs(d$lfc)
  d=d[order(-d$lfc),]
  name=paste(celltypes[i],"genesets.density.tbl",sep=".")
  assign(name,d)
  pos=nrow(d[d$lfc > 1,])
  neg=nrow(d[d$lfc < -1,])
  if(pos > 0 | neg > 0)
  {
    tmp.df=data.frame(GOname=rownames(d[d$abslfc > 1,]), cell=celltypes[i])
    GOBP_terms.df=rbind(GOBP_terms.df,tmp.df)
  }
  tmp.df=data.frame(cell=celltypes[i],pos= pos,neg=neg)
  density.no.BP.perturbed=rbind(density.no.BP.perturbed,tmp.df)
  name=paste(celltypes[i],"genesets.density.lfc.tbl",sep=".")
  d$cell=celltypes[i]
  d=d[,c("lfc","cell")]
  assign(name,d)
}

df_to_plot=melt(density.no.BP.perturbed)

npgcolors=pal_npg("nrc", alpha = 1)(10)
p.no_density_GOBP=ggplot(df_to_plot,aes(x=cell,y=value)) +
  geom_col(aes(fill=variable),width=0.5)+
  labs(y="# of GO BP",x="Cell types")+
  scale_fill_manual(values=c("pos"=npgcolors[4],"neg"=npgcolors[8]),name=expression(Delta ~ "coherence" ),labels=c("gain","loss"))+
 theme_bw(base_size=12)+theme(legend.position="top")
ggsave(p.no_density_GOBP,filename="Figures/p.no_density_GOBP.pdf", device="pdf",width=2.5,height=2,units="in")


#box plot
for_boxplot=rbind(Mic.genesets.density.lfc.tbl,Oli.genesets.density.lfc.tbl)
for_boxplot=rbind(for_boxplot,Ex.genesets.density.lfc.tbl)
for_boxplot=rbind(for_boxplot,In.genesets.density.lfc.tbl)

npgcolors=pal_npg("nrc", alpha = 1)(10)
p.lfc_density_GOBP_boxplot=ggplot(for_boxplot,aes(x=cell,y=lfc)) +
  geom_boxplot()+
  labs(y="change in GO edge density",x="Cell types")+
 theme_bw(base_size=12)+theme(legend.position="top")
#ggsave(p.lfc_density_GOBP_boxplot,filename="Figures/p.lfc_density_GOBP_boxplot.pdf", device="pdf",width=3,height=3,units="in")

#if using migsdb
#ggsave(p.lfc_density_GOBP_boxplot,filename="Figures/p.lfc_density_migsdb_boxplot.pdf", device="pdf",width=3,height=3,units="in")


#tmp1=merge(Ex.genesets.density.tbl,In.genesets.density.tbl,by="GOname")
#tmp2=merge(Mic.genesets.density.tbl,Oli.genesets.density.tbl,by="GOname")
#df=merge(tmp1,tmp2,by="GOname")
#rownames(df)=df$GOname
#df$GOname=NULL
#df=df[rowSums(df) >0,]
