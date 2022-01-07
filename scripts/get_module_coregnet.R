rm(list=ls())
source('../scripts/read_data.R')
source('../scripts/functions_for_network_analysis.R')
source('../scripts/load_libraries.R')


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

all.modules=data.frame("gene"=NULL, "moduleID"=NULL, "celltype"=NULL, "condition"=NULL)
list=ls(pattern=pattern)

for (i in 1:length(list))
{
  name=list[i]
  name=gsub(".network.JI.coreg.modules.JI_0.2.ModSize_30","",name)
  tag=strsplit(name,"\\.")[[1]]
  condition=tag[1]
  cell=tag[2]
  df=get(list[i])
  df$cell=cell
  df$condition=condition
  all.modules=rbind(all.modules,df)
}

all.modules$module=paste(all.modules$cell,all.modules$condition,sep="_")
all.modules$module=paste(all.modules$moduleID,all.modules$module,sep="_")
all.modules=all.modules[,c("gene","module")]

background=as.data.frame(unique(all.modules$gene))
colnames(background)=""

genesets=as.data.frame(unique(all.modules$module))
colnames(genesets)="genesets"
genesets$name=genesets$genesets
genesets=unique(genesets)


#for paper
write.table(all.modules,file="~/work/scNET_manuscript/AD_MIT/supp_data/all.modules.txt",row.names=F,col.names=T,sep="\t",quote=FALSE)
