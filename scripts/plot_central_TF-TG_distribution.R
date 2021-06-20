# Regulatory network format;: 3 columns (TF target score)
rm(list=ls())
source('~/work/scNET-devel/scripts/load_libraries.R')
source('~/work/scNET-devel/scripts/read_data.R')
source('~/work/scNET-devel/scripts/functions_for_network_analysis.R')

celltypes=c("Ex1","Ex2","Ex3e","Ex4","Ex5b","Ex6a","Ex6b","Ex8","Ex9",
"In1a","In1b","In1c","In3","In4a","In4b","In6a","In6b","In7","In8",
"Mic","Oli")


centralities=c("betweenness","degree_out","degree_in")

mean_indegree=function(ct.centrality_matrix)
{
  tmp=ct.centrality_matrix[,3]
  return(mean(tmp))
}

mean_outdegree=function(ct.centrality_matrix)
{
  tmp=ct.centrality_matrix %>% filter(ct.centrality_matrix[,4] > 0)
  tmp=tmp[,4]
  return(mean(tmp))
}

no.TF.in.TG.tbl= data.frame(matrix(ncol = 2, nrow = 0))
colnames(no.TF.in.TG.tbl)=c("Cluster","percentTFasTargets")

for (k in 1:length(centralities))
{
  c=centralities[k]
  nets=ls(pattern="*_openchrom_en_GRN.network")
  for(i in 1:length(nets))
  {
    name=nets[i]
    tag=gsub("_openchrom_en_GRN.network","",name)
    name=paste(tag,c,sep=".")
    net=as.data.frame(lapply(nets[i],get))
    net1=net[net$TFbs %in% "promoter", ]
    net2=net[net$TFbs %in% "both", ]
    net=rbind(net1,net2)
    net=net[,c("TF","TG","abs_coef")]
    net=distinct(net)
    TF=unique(net$TF)
    TG=unique(net$TG)
    nTFs=length(intersect(TG,TF))/length(TF)
    nTFs=nTFs*100
    no.TF.in.TG.tbl=rbind(no.TF.in.TG.tbl,data.frame(Cluster=tag,percentTFasTargets=nTFs))
    assign(name, get_centrality(net,c,tag))
    }
}

mean_indegree.tbl= data.frame(matrix(ncol = 2, nrow = 0))
colnames(mean_indegree.tbl)=c("Cluster","in_deg_mean")

mean_outdegree.tbl= data.frame(matrix(ncol = 2, nrow = 0))
colnames(mean_outdegree.tbl)=c("Cluster","out_deg_mean")

for(i in 1:length(celltypes))
{
  tag=celltypes[i]
  pattern=paste(tag,"\\.",sep="")
  list=ls(pattern=pattern)
  tmp=merge(get(list[1]),get(list[2]))
  tmp=merge(tmp,get(list[3]))
  name=paste(tag,"centrality_matrix",sep=".")
  assign(name,tmp)
  mean_in=mean_indegree(get(name))
  mean_indegree.tbl=rbind(mean_indegree.tbl,data.frame(Cluster=tag,in_deg_mean=mean_in))
  mean_out=mean_outdegree(get(name))
  mean_outdegree.tbl=rbind(mean_outdegree.tbl,data.frame(Cluster=tag,out_deg_mean=mean_out))
}

no.TF.in.TG.tbl=distinct(no.TF.in.TG.tbl[order(no.TF.in.TG.tbl$percentTFasTargets),])
mean_outdegree.tbl=mean_outdegree.tbl[order(mean_outdegree.tbl$out_deg_mean),]
mean_indegree.tbl=mean_indegree.tbl[order(mean_indegree.tbl$in_deg_mean),]

pdf(file="meanIndegree.pdf")
ggplot(mean_indegree.tbl, aes(y=in_deg_mean,x=reorder(Cluster,-in_deg_mean))) +
 geom_bar(stat="identity")+labs(x="Cluster",y="In-degree (mean)") +
 geom_hline(yintercept=mean(mean_indegree.tbl$in_deg_mean))
dev.off()

pdf(file="meanOutdegree.pdf")
ggplot(mean_outdegree.tbl, aes(y=out_deg_mean,x=reorder(Cluster,-out_deg_mean))) +
 geom_bar(stat="identity")+labs(x="Cluster",y="Out-degree (mean)") +
 geom_hline(yintercept=mean(mean_outdegree.tbl$out_deg_mean))
dev.off()

pdf(file="Frac_TFTG.pdf")
ggplot(no.TF.in.TG.tbl, aes(y=percentTFasTargets,x=reorder(Cluster,percentTFasTargets))) +
geom_bar(stat="identity")+labs(x="Cluster",y="TFs with incoming edges (%)") +
geom_hline(yintercept=mean(no.TF.in.TG.tbl$percentTFasTargets))
dev.off()
