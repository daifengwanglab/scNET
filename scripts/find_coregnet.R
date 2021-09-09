rm(list=ls())
source('~/work/scNET-devel/scripts/load_libraries.R')
source('~/work/scNET-devel/scripts/read_data.R')
source('~/work/scNET-devel/scripts/functions_for_network_analysis.R')

#fig 2a in https://www.nature.com/articles/s12276-020-00528-0#Sec2
#fig 1c in https://www.nature.com/articles/s41467-019-10591-5/figures/1

#Read GO data
data=GSA.read.gmt('~/work/scNET_manuscript/genome/genesets/GO_annotations-9606-inferred-allev.gmt')
genesets=data$genesets
names(genesets)=data$geneset.descriptions

#number of random nets needed
nrandnets=1
set.seed(123)

nets=ls(pattern="*\\.network")

th=0.2
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
  mat[mat < th] <-0 #remove edges with less than 10% overlap
  name=gsub(".mat",".modules", name)
  assign(name, detect_modules(mat))
  #for (j in 1:nrandnets)
  #{
  #  net.rand.df=net
#    net.rand.df$TF=sample(net.rand.df$TG)
#    name=paste(nets[i],"randomCoregNet",sep=".")
#    name=gsub(".network.JI.coreg.mat","",name)
#    name=paste(name, j,sep=".")
#    assign(name, find_target_pairs_matrix(net.rand.df))
#    mat=get(name)
#    mat[mat < 0.1] <-0 #remove edges with less than 10% overlap
#    name=paste(name,".modules", sep="")
#    assign(name, detect_modules(mat))
#  }
}


##################
##enrichment analysis
###########################
#GO enrichment
diff.cent.enrich.tbl=data.frame("label"=NULL,"pval"=NULL,"fdr"=NULL,"signature"=NULL,"geneset"=NULL,
"overlap"=NULL,"background"=NULL,"hits"=NULL,"cell"=NULL,"module"=NULL)

module.enrich.tbl=data.frame("cell"=NULL,"total"=NULL,"annotated"=NULL)
#genesets <- msigdb_gsets("Homo sapiens", "C2", "CP:KEGG", clean=TRUE)


list=ls(pattern="*\\.modules")
for (i in 1:length(list))
{
  total=0
  annotated=0
  totalBP=0
  name=list[i]
  name=gsub(".network.JI.coreg.modules","",name)
  df=get(list[i])
  modnames=unique(factor(df$moduleID))
#  universe=df$gene
  total=length(modnames)
  print(paste(name,length(modnames),sep=":"))
  for (j in 1:length(modnames))
  {
    mod.genes=df[df$moduleID %in% modnames[j],]$gene
    hyp_obj = hypeR(mod.genes, genesets,fdr=0.01)
    hyp_df =  hyp_obj$data
    if(nrow(hyp_df) > 0)
    {
      annotated=annotated+1
      print(paste(modnames[j],nrow(hyp_df),sep=":"))
      hyp_df$cell=name
      hyp_df$module=modnames[j]
    }
    diff.cent.enrich.tbl=rbind(diff.cent.enrich.tbl,  hyp_df)
    newtbl=data.frame("cell"=name,"total"=total,"annotated"=annotated)
  }
  module.enrich.tbl=rbind(module.enrich.tbl,newtbl)
}

#split random and real into two dfs
#rand.module.enrich.tbl=module.enrich.tbl[module.enrich.tbl$cell %like% "randomCoregNet",]
#module.enrich.tbl=module.enrich.tbl[!module.enrich.tbl$cell %like% "randomCoregNet",]

#rand.diff.cent.enrich.tbl=diff.cent.enrich.tbl[diff.cent.enrich.tbl$cell %like%  "randomCoregNet",]
#diff.cent.enrich.tbl=diff.cent.enrich.tbl[!diff.cent.enrich.tbl$cell %like%  "randomCoregNet",]

#count total BP per net
ct.totalBP=data.frame("cell"=NULL,"TotalBP"=NULL)
list=unique(diff.cent.enrich.tbl$cell)
for (i in 1:length(list))
{
  df=diff.cent.enrich.tbl[diff.cent.enrich.tbl$cell %in%  list[i],]
  tmp.df=data.frame("cell"=list[i],"TotalBP"=length(unique(df$label)))
  ct.totalBP=rbind(ct.totalBP,tmp.df)
}

#count total BP per random net
#rand.ct.totalBP=data.frame("cell"=NULL,"Total"=NULL)
#for (i in 1:nrow(ct.totalBP))
#{
#  count=0
#  tag=ct.totalBP[i,]$cell
#  for (j in 1:nrandnets)
#  {
#    name=paste(tag, "network.randomCoregNet",sep=".")
#    name=paste(name,j,sep=".")
#    df=rand.diff.cent.enrich.tbl[rand.diff.cent.enrich.tbl$cell %like%  name,]
#    if (length(unique(df$label)) >= ct.totalBP[i,]$TotalBP)
#    {
#      count=count+1
#    }
#  }
#  tmp.df=data.frame("cell"=ct.totalBP[i,]$cell,Total=count
#  rand.ct.totalBP=rbind(rand.ct.totalBP,tmp.df)
#}
