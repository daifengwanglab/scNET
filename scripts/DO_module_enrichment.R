

#copy this to the end of 'find_coregnet.R'


##################
##enrichment analysis
###########################
#GO enrichment
diff.cent.enrich.tbl=data.frame("label"=NULL,"pval"=NULL,"fdr"=NULL,"signature"=NULL,"geneset"=NULL,
"overlap"=NULL,"background"=NULL,"hits"=NULL,"cell"=NULL,"module"=NULL,"th"=NULL,"msize"=NULL)

module.enrich.tbl=data.frame("cell"=NULL,"total"=NULL,"annotated"=NULL, "th"=NULL,"msize"=NULL)

list=ls(pattern="*\\.modules")
list=Filter(function(x) !any(grepl("JI_0.5", x)), list)
list=Filter(function(x) !any(grepl("JI_0.6", x)), list)
list=Filter(function(x) !any(grepl("JI_0.7", x)), list)
list=Filter(function(x) !any(grepl("JI_0.8", x)), list)
list=Filter(function(x) !any(grepl("JI_0.9", x)), list)


for (i in 1:length(list))
{
  total=0
  annotated=0
  totalBP=0
  df=get(list[i])
  df=df[as.character(df$moduleID) != "0",]
  name=list[i]
  name=gsub(".network.JI.coreg.modules.","",name)
  name=gsub("JI_"," ",gsub(".ModSize_"," ",name))
  th=sapply(strsplit(name, " "), head, 3)[2]
  msize=sapply(strsplit(name, " "), head, 3)[3]
  name=sapply(strsplit(name, " "), head, 3)[1]
  ct=ifelse(name %like% "Ex","Ex",ifelse(name %like% "In","In",ifelse(name %like% "Mic","Mic","Oli")))
  cond=ifelse(name %like% "AD","AD","Ctrl")
  modnames=unique(factor(df$moduleID))
#  universe=df$gene
  total=length(modnames)
  message=paste(name,total,sep=":")
  message=paste(message,th,sep=":")
  message=paste(message,msize,sep=":")
  print (message)
  if (total > 0)
  {
    for (j in 1:length(modnames))
    {
      mod.genes=df[df$moduleID %in% modnames[j],]$gene
      hyp_obj = hypeR(mod.genes, genesets,fdr=0.05)
      hyp_df =  hyp_obj$data
      if(nrow(hyp_df) > 0)
      {
        annotated=annotated+1
        message=paste(name,modnames[j],sep="_ModName:")
        message=paste(message,nrow(hyp_df),sep="_totBP:")
        print (message)
        hyp_df$cell=name
        hyp_df$module=modnames[j]
        hyp_df$th=th
        hyp_df$msize=msize
      }
      diff.cent.enrich.tbl=rbind(diff.cent.enrich.tbl,  hyp_df)
      newtbl=data.frame("cell"=name,"total"=total,"annotated"=annotated,"th"=th,"msize"=msize)
    }
    module.enrich.tbl=rbind(module.enrich.tbl,newtbl)
  }
}
GO.diff.cent.enrich.tbl=diff.cent.enrich.tbl



module.enrich.tbl$PercAnnot=(module.enrich.tbl$annotated/module.enrich.tbl$total)*100
module.enrich.tbl$condition=ifelse(module.enrich.tbl$cell %like% "AD","AD","Ctrl")
module.enrich.tbl$cell=gsub("AD.","",module.enrich.tbl$cell)
module.enrich.tbl$cell=gsub("Ctrl.","",module.enrich.tbl$cell)


p3=ggplot(module.enrich.tbl,aes(x=cell,y=PercAnnot,fill=condition))+
geom_bar(stat="identity",position="dodge")+facet_grid(th~as.numeric(msize))+
scale_fill_manual(values=c("AD"=npgcolors[1],"Ctrl"=npgcolors[2]))+
theme_bw(base_size=12) + labs(y="Fraction of predicted modules\n with GO annotations",x="cell types")+
theme(legend.position = "top")+theme(axis.text.x=element_text(angle=90))

ggsave(p3,filename="../Figures/p.frac_mod_annot.pdf", device="pdf",width=6,height=6,units="in")



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

#fold change values
all.deg=read_xlsx("~/work/scNET_manuscript/data/gematrix/Diff.Exp.Genes.DataS2.MIT.xlsx",sheet="Mic",skip=1)[,1:9]
mic.deg=as.data.frame(all.deg[,c(1,5)])
colnames(mic.deg)=c("gene","FC")
attr.tbl=left_join(attr.tbl,mic.deg)
write.table(attr.tbl,file="attr.txt",row.names=F,
col.names=T,sep="\t",quote=FALSE)



##########################################
#cross compare modules

tm=c(0.2)
msize=c(30)

for(i in 1:length(nets))
{
  for (j in 1:length(th))
  {
    for (k in 1:length(msize))
    {
      name=nets[i]
      celltype=gsub(".network","",name)
      net=as.data.frame(lapply(nets[i],get))
      net=net[,c("TF","TG","abs_coef")]
      net=distinct(net)
      name=paste(name,"JI.coreg.mat",sep=".")
      assign(name, find_target_pairs_matrix(net))
      mat=get(name)
      mat[mat < th[j]] <-0 #remove edges with less than th% overlap
      name=gsub(".mat",".modules", name)
      tag=paste("JI",th[j],sep="_")
      name=paste(name,tag,sep=".")
      tag=paste("ModSize",msize[k],sep="_")
      name=paste(name,tag,sep=".")
      assign(name, detect_modules(mat,msize[k]))
    }
  }
}


pattern="*.network.JI.coreg.modules.JI_0.2.ModSize_30"
list=ls(pattern=pattern)

for (i in 1:length(list))
{
  name=list[i]
  name=gsub(".network.JI.coreg.modules","",name)
  df=get(list[i])
  df=df[!(df$moduleID %like% "_0"), ] #remove singletons (module 0)
  df$moduleID=paste(name,df$moduleID,sep="_")
  assign(list[i],df)
}

data = Reduce(function(x, y) merge(x, y, all=T), lapply(list,get), accumulate=F)
data$score=1
colnames(data)=c("TF","target","score")
data$target=gsub(".JI_0.2.ModSize_30","",data$target)
m=acast(data, TF~target, value.var="score")
m=t(m)
m[is.na(m)]=0 #set NA =0
#find cardinalities
# Find that paper and add reference
i12 = m %*% t(m)
s = diag(i12) %*% matrix(1, ncol = length(diag(i12)))
u12 = s + t(s) - i12
jacc= i12/u12
#jacc[jacc > 0.5] = 1
diag(jacc)=0

npgcolors=pal_npg("nrc", alpha = 0)(10)
tmp=data.frame(cell=colnames(jacc))
tmp$celltype=ifelse(tmp$cell %like% "Ex","Ex",ifelse(tmp$cell %like% "In","In",ifelse(tmp$cell %like% "Mic","Mic","Oli")))
tmp$state=ifelse(tmp$cell %like% "AD","AD","Ctrl")


rowannot=rowAnnotation(Celltype=tmp$celltype,State=tmp$state,
 col = list(Celltype = c("In" = npgcolors[1], "Ex" = npgcolors[2], "Mic"=npgcolors[3],"Oli"=npgcolors[4]),
       State=c("AD"=npgcolors[5],"Ctrl"=npgcolors[6]),width = unit(2, "mm")
   )
 )

column_ha=HeatmapAnnotation(Celltype=tmp$celltype,State=tmp$state,
col = list(Celltype = c("In"=npgcolors[1],"Ex"=npgcolors[2],"Mic"=npgcolors[3],"Oli"=npgcolors[4]),
      State=c("AD"=npgcolors[5],"Ctrl"=npgcolors[6])
  ))
heatmap=Heatmap(jacc,col=viridis(100),top_annotation = column_ha)+rowannot

pdf(file="../Figures/module_cross_heatmap.pdf",width=5,height=5)
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
