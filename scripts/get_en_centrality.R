
# Regulatory network format;: 3 columns (TF target score)
rm(list=ls())
source('~/work/scNET-devel/scripts/load_libraries.R')
source('~/work/scNET-devel/scripts/read_data.R')
source('~/work/scNET-devel/scripts/functions_for_network_analysis.R')


celltypes=c("Ex1","Ex2","Ex3e","Ex4","Ex5b","Ex6a","Ex6b","Ex8","Ex9",
"In1a","In1b","In1c","In3","In4a","In4b","In6a","In6b","In7","In8",
"Mic","Oli")


centralities=c("betweenness","degree_out","degree_in")

for (k in 1:length(centralities))
{
  c=centralities[k]
  nets=ls(pattern="*_openchrom_en_GRN.network")
  for(i in 1:length(nets))
  {
    name=nets[i]
    tag=gsub("_openchrom_en_GRN.network",".Enhancer",name)
    name=paste(tag,c,sep=".")
    net=as.data.frame(lapply(nets[i],get))

    #select only TF-TG interactions through promoters
    TF.net1=net[net$TFbs %in% "promoter", ]
    TF.net2=net[net$TFbs %in% "both", ]
    TF.net=rbind(TF.net1,TF.net2)
    TF.net=net[,c("TF","TG","mse")]
    TF.net=distinct(TF.net)

    #En acts as both target as well as RFs
    En.net=net[net$TFbs %in% "enhancer", ]
    TF.En.net=En.net[,c("TF","enhancer","mse")]
    En.TG.net=En.net[,c("enhancer","TG","mse")]
    colnames(En.TG.net)=colnames(TF.En.net)=colnames(TF.net)
    En.net=rbind(TF.En.net,En.TG.net)
    En.net=distinct(En.net)

    mixed.net=rbind(TF.net,En.net)
    assign(name, get_centrality(mixed.net,c,tag))
    }
}

for(i in 1:length(celltypes))
{
  tag=celltypes[i]
  pattern=paste(tag,"\\.",sep="")
  list=ls(pattern=pattern)
  tmp=merge(get(list[1]),get(list[2]))
  tmp=merge(tmp,get(list[3]))
  name=paste(tag,"centrality_matrix",sep=".")
  assign(name,tmp)
}

pattern="*.centrality_matrix"
list=ls(pattern=pattern)

old=get(list[1])
for (i in 2:length(list))
{
  old=merge(old, get(list[i]), by="gene", all=TRUE)
}

write.table(old, file="/Users/chiraggupta/work/scNET_manuscript/results/centralities/centrality_matrix.En-networks.all-cellTypes.txt",col.names=TRUE,row.names=FALSE, sep="\t",quote=F)
