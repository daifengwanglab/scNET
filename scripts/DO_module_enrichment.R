

#copy this to the end of 'find_coregnet.R'
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
#prep for fig 5 subnets
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


#select module 2 alz risk genes in AD mic net
indx.c=match(alz.risk.ctrl.module.genes,colnames(AD.Mic.network.JI.coreg.mat))
indx.r=match(alz.risk.ctrl.module.genes,rownames(AD.Mic.network.JI.coreg.mat))
AD.mic.alz.genes.coregnet.mat=AD.Mic.network.JI.coreg.mat[indx.r,indx.c]
g=graph.adjacency(AD.mic.alz.genes.coregnet.mat,weighted=TRUE)
df <- get.data.frame(igraph::simplify(g,remove.multiple = TRUE, remove.loops = TRUE))
df=df[df$weight >= 0.3,]
AD.mic.alz.genes.coregnet.df=df
write.table(AD.mic.alz.genes.coregnet.df,file="AD.ctrl.mod2.mic.alz.genes.coregnet.dat",row.names=F,
col.names=T,sep="\t",quote=FALSE)



#fold change values
all.deg=read_xlsx("~/work/scNET_manuscript/data/gematrix/Diff.Exp.Genes.DataS2.MIT.xlsx",sheet="Mic",skip=1)[,1:9]
mic.deg=as.data.frame(all.deg[,c(1,5)])
colnames(mic.deg)=c("gene","FC")
attr.tbl=left_join(attr.tbl,mic.deg)
write.table(attr.tbl,file="attr.txt",row.names=F,
col.names=T,sep="\t",quote=FALSE)



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



#Heatmapp

paletteLength=100
myColor <- colorRampPalette(c("yellow", "white", "blue"))(paletteLength)
myBreaks <- c(seq(-2, 0, length.out=ceiling(paletteLength/2) + 1),
+ seq(max(mat.new)/paletteLength, 2, length.out=floor(paletteLength/2)))


p=pheatmap(mat.new,cellwidth=12,cellheight=10,show_rownames=TRUE,
main="Module DE",
 breaks=myBreaks,
 cluster_cols=FALSE,
 cluster_rows=TRUE,
 gaps_col = 3,
 border_color = "grey60",
 angle_col=c("45")
 )

pdf(file="Module-FC_page.pdf")
p
dev.off()
