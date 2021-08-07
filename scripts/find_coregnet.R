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
  name=gsub(".mat",".modules", name)
  assign(name, detect_modules(mat))
  net.igraph=graph_from_data_frame(net, directed = FALSE, vertices = NULL)
  for (j in 1:10)
  {
    #net.rand=erdos.renyi.game(length(V(net.igraph)),length(E(net.igraph)), type="gnm",directed = TRUE)
    net.rand.df=as.data.table(net)[, lapply(.SD, sample)]

  }
}


##################
##enrichment analysis

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
    hyp_obj = hypeR(mod.genes, genesets,fdr=0.05)
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

ct.totalBP=data.frame("cell"=NULL,"TotalBP"=NULL)
list=unique(diff.cent.enrich.tbl$cell)
for (i in 1:length(list))
{
  df=diff.cent.enrich.tbl[diff.cent.enrich.tbl$cell %in%  list[i],]
  tmp.df=data.frame("cell"=list[i],"TotalBP"=length(unique(df$label)))
  ct.totalBP=rbind(ct.totalBP,tmp.df)
}

##########################################
#cross compare modules
list=ls(pattern="*\\.modules")
for (i in 1:length(list))
{
  name=list[i]
  name=gsub(".network.JI.coreg.modules","",name)
  df=get(list[i])
  df$moduleID=paste(name,df$moduleID,sep="_")
  df=df[!(df$moduleID %like% "_0"), ] #remove singletons (module 0)
  assign(list[i],df)
}

data = Reduce(function(x, y) merge(x, y, all=T), lapply(list,get), accumulate=F)
data$score=1
colnames(data)=c("TF","target","score")
m=acast(data, TF~target, value.var="score")
m=t(m)
m[is.na(m)]=0 #set NA =0
#find cardinalities
# Find that paper and add reference
i12 = m %*% t(m)
s = diag(i12) %*% matrix(1, ncol = length(diag(i12)))
u12 = s + t(s) - i12
jacc= i12/u12
jacc[jacc > 0.5] = 1
diag(jacc)=0
pheatmap(jacc, color=rev(inferno(10)), border=FALSE)
################################################################

#plotting
npgcolors=pal_npg("nrc", alpha = 1)(10)

list=unique(diff.cent.enrich.tbl$cell)
for (i in 1:length(list))
{
  df=diff.cent.enrich.tbl[diff.cent.enrich.tbl$cell %in% list[i],]
  wordcloudtble=as.data.frame(table(df$label))
  wordcloudtble$Var1=gsub("_"," ",wordcloudtble$Var1)
  filename=paste(list[i],"wordcloud.pdf", sep=".")
  pdf(file=filename,width=10, height=10)
  wordcloud(words=wordcloudtble$Var1, freq=wordcloudtble$Freq,  min.freq = 1, random.order=FALSE, rot.per=0, fixed.asp=TRUE,colors=npgcolors)
  dev.off()
}


#plot no. of module vs annotated
module.enrich.tbl$state=ifelse(module.enrich.tbl$cell %like% "AD","AD","Healthy")
module.enrich.tbl$cell=gsub("AD.","",gsub("Ctrl.","",module.enrich.tbl$cell))
df2=module.enrich.tbl
df2$fraction=df2$annotated/df2$total
p.barplot=ggplot(df2, aes(x=cell,y=annotated,fill=cell)) +geom_bar(stat="identity", position=position_dodge()) +
facet_wrap(~state)
