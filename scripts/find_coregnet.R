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
  mat[mat < 0.1] <-0 #remove edges with less than 10% overlap
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

#DO
diff.cent.enrich.tbl=data.frame("label"=NULL,"pval"=NULL,"fdr"=NULL,"signature"=NULL,"geneset"=NULL,
"overlap"=NULL,"background"=NULL,"hits"=NULL,"cell"=NULL,"module"=NULL)

module.DO.enrich.tbl=data.frame("cell"=NULL,"total"=NULL,"annotated"=NULL)

source('~/work/scNET_manuscript/get_api.R')

#DO enrichment
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
    hyp_obj = disease_enrichment( entities =mod.genes, vocabulary = "HGNC", database = "ALL")
    hyp_df =  hyp_obj@qresult[hyp_obj@qresult$FDR<0.01, c("Description", "FDR", "Ratio",  "BgRatio","shared_symbol")]
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
  module.DO.enrich.tbl=rbind(module.DO.enrich.tbl,newtbl)
}
################################
#extract ALZ modules based on DO enrichment
alz=diff.cent.enrich.tbl[diff.cent.enrich.tbl$Description %like% "Alzheimer",]


#ctrl has module 2 with alz genes signif enriched; extract the hits
shared_symbol=alz[alz$cell %in% "Ctrl.Mic",]
shared_symbol=as.data.frame(shared_symbol[,c("shared_symbol")])
colnames(shared_symbol)=c("alz.hits.ctrl.module")
shared_symbol=data.frame(alz.hits.ctrl.module = unlist(strsplit(as.character(shared_symbol$alz.hits.ctrl.module), ";")))
alz.risk.ctrl.module.genes=unique(shared_symbol$alz.hits.ctrl.module)

#select module 2 alz risk genes in ctrl mic net
indx.c=match(alz.risk.ctrl.module.genes,colnames(Ctrl.Mic.network.JI.coreg.mat))
indx.r=match(alz.risk.ctrl.module.genes,rownames(Ctrl.Mic.network.JI.coreg.mat))
ctrl.mic.alz.genes.coregnet.mat=Ctrl.Mic.network.JI.coreg.mat[indx.r,indx.c]
g=graph.adjacency(ctrl.mic.alz.genes.coregnet.mat,weighted=TRUE)
df <- get.data.frame(igraph::simplify(g,remove.multiple = TRUE, remove.loops = TRUE))
df=df[df$weight >= 0.3,]
ctrl.mic.alz.genes.coregnet.df=df
write.table(ctrl.mic.alz.genes.coregnet.df,file="ctrl.mod2.mic.alz.genes.coregnet.dat",row.names=F,
col.names=T,sep="\t",quote=FALSE)

#select alz risk genes in AD mic modules
shared_symbol=alz[alz$cell %in% "AD.Mic",]
shared_symbol=as.data.frame(shared_symbol[,c("shared_symbol")])
colnames(shared_symbol)=c("alz.hits.AD.module")
shared_symbol=data.frame(alz.hits.AD.module = unlist(strsplit(as.character(shared_symbol$alz.hits.AD.module), ";")))
alz.risk.AD.module.genes=unique(shared_symbol$alz.hits.AD.module)

indx.c=match(alz.risk.AD.module.genes,colnames(AD.Mic.network.JI.coreg.mat))
indx.r=match(alz.risk.AD.module.genes,rownames(AD.Mic.network.JI.coreg.mat))
AD.mic.alz.genes.coregnet.mat=AD.Mic.network.JI.coreg.mat[indx.r,indx.c]
g=graph.adjacency(AD.mic.alz.genes.coregnet.mat,weighted=TRUE)
df <- get.data.frame(igraph::simplify(g,remove.multiple = TRUE, remove.loops = TRUE))
df=df[df$weight >= 0.3,]
AD.mic.alz.genes.coregnet.df=df
write.table(AD.mic.alz.genes.coregnet.df,file="AD.mic.alz.genes.coregnet.dat",row.names=F,
col.names=T,sep="\t",quote=FALSE)



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
rand.ct.totalBP=data.frame("cell"=NULL,"Total"=NULL)
for (i in 1:nrow(ct.totalBP))
{
  count=0
  tag=ct.totalBP[i,]$cell
  for (j in 1:nrandnets)
  {
    name=paste(tag, "network.randomCoregNet",sep=".")
    name=paste(name,j,sep=".")
    df=rand.diff.cent.enrich.tbl[rand.diff.cent.enrich.tbl$cell %like%  name,]
    if (length(unique(df$label)) >= ct.totalBP[i,]$TotalBP)
    {
      count=count+1
    }
  }
  tmp.df=data.frame("cell"=ct.totalBP[i,]$cell,Total=count)
  rand.ct.totalBP=rbind(rand.ct.totalBP,tmp.df)
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

npgcolors=pal_npg("nrc", alpha = 1)(10)
tmp=data.frame(cell=colnames(jacc))
tmp$celltype=ifelse(tmp$cell %like% "Ex","Ex",ifelse(tmp$cell %like% "In","In",ifelse(tmp$cell %like% "Mic","Mic","Oli")))
tmp$state=ifelse(tmp$cell %like% "AD","AD","Control")

rowannot=rowAnnotation(Celltype=tmp$celltype,State=tmp$state,
col = list(Celltype = c("In" = npgcolors[1], "Ex" = npgcolors[2], "Mic"=npgcolors[3],"Oli"=npgcolors[4]),
      State=c("AD"=npgcolors[5],"Control"=npgcolors[6]),width = unit(2, "mm")
  )
)
column_ha=HeatmapAnnotation(Celltype=tmp$celltype,State=tmp$state,
col = list(Celltype = c("In" = npgcolors[1], "Ex" = npgcolors[2], "Mic"=npgcolors[3],"Oli"=npgcolors[4]),
      State=c("AD"=npgcolors[5],"Control"=npgcolors[6])
  ))
heatmap=Heatmap(jacc,col=viridis(3),show_column_names = FALSE)+rowannot

pdf(file="Figures/module_cross_heatmpa.pdf",width=5,height=4)
heatmap
dev.off()

####cell type specific module compare
celltypes=c("Mic","Oli","Ex","In")
for (i in 1:length(celltypes))
{
  ct.df=data[data$target %like% celltypes[i],]
  m=acast(ct.df, TF~target, value.var="score")
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
  tmp=data.frame(cell=colnames(jacc))
  tmp$celltype=ifelse(tmp$cell %like% "Ex","Ex",ifelse(tmp$cell %like% "In","In",ifelse(tmp$cell %like% "Mic","Mic","Oli")))
  tmp$state=ifelse(tmp$cell %like% "AD","AD","Control")


  column_ha=HeatmapAnnotation(State=tmp$state,
    col = list(State=c("AD"=npgcolors[1],"Control"=npgcolors[2]),width = unit(2, "mm"))
  )

  rowannot=rowAnnotation(State=tmp$state,
    col = list(State=c("AD"=npgcolors[1],"Control"=npgcolors[7]),width = unit(2, "mm"))
  )
  name=paste(celltypes[i],"module_compare_htmap",sep=".")
  heatmap=Heatmap(jacc,col=viridis(3),show_column_names = FALSE,bottom_annotation=column_ha)+rowannot
  assign(name, heatmap)
}

pdf(file="Figures/EX.module_compare.heatmap.pdf",width=5,height=4)
Ex.module_compare_htmap
dev.off()

pdf(file="Figures/In.module_compare.heatmap.pdf",width=5,height=4)
In.module_compare_htmap
dev.off()

pdf(file="Figures/Mic.module_compare.heatmap.pdf",width=5,height=4)
Mic.module_compare_htmap
dev.off()

pdf(file="Figures/Oli.module_compare.heatmap.pdf",width=5,height=4)
Oli.module_compare_htmap
dev.off()

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
