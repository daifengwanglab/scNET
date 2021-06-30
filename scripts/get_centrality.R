
source('~/work/scNET-devel/scripts/read_data.R')
source('~/work/scNET-devel/scripts/functions_for_network_analysis.R')
source('~/work/scNET-devel/scripts/load_libraries.R')

my_entrez_gene_info=read.table('~/work/scNET_manuscript/genome/hgncSymbols_to_entrez.human.txt', header=T,sep="\t")
colnames(my_entrez_gene_info)=c("entrezID","gene")


#celltypes=c("Ex1","Ex2","Ex3e","Ex4","Ex5b","Ex6a","Ex6b","Ex8","Ex9",
#"In1a","In1b","In1c","In3","In4a","In4b","In6a","In6b","In7","In8",
#"Mic","Oli")


celltypes=c("Mic","Oli","Ex","In")

centralities=c("betweenness","degree_out","degree_in")

for (k in 1:length(centralities))
{
  c=centralities[k]
  nets=ls(pattern="*\\.network")
  for(i in 1:length(nets))
  {
    name=nets[i]
    tag=gsub(".network","",name)
    name=paste(tag,c,sep=".")
    net=as.data.frame(lapply(nets[i],get))
    TF.net1=net[net$TFbs %in% "promoter", ]
    TF.net2=net[net$TFbs %in% "both", ]
    TF.net=rbind(TF.net1,TF.net2)
    TF.net=net[,c("TF","TG","mse")]
    TF.net=distinct(TF.net)
    assign(name, get_centrality(TF.net,c,tag))
  }
}

rm(list=nets)

for(i in 1:length(celltypes))
{
  tag=celltypes[i]
  pattern=paste(tag,"\\.",sep="")
  list=ls(pattern=pattern)
  data = Reduce(function(x, y) merge(x, y, all=T), lapply(list,get), accumulate=F)
  #tmp=merge(get(list[1]),get(list[2]))
  #tmp=merge(tmp,get(list[3]))
  name=paste(tag,"centrality_matrix",sep=".")
  assign(name,data)
}


pattern="*.centrality_matrix"
list=ls(pattern=pattern)

old=get(list[1])
for (i in 2:length(list))
{
  old=merge(old, get(list[i]), by="gene", all=TRUE)
}

#write.table(old, file="centrality_matrix.all-cellTypes.txt",col.names=TRUE,row.names=FALSE, sep="\t",quote=F)

###plotting

cent.mat=old
rownames(cent.mat)=cent.mat$gene
cent.mat$gene=NULL
cent.mat[is.na(cent.mat)]=0

#standardize
cent.mat.scaled=apply(cent.mat, MARGIN = 2, FUN = function(X) (X - min(X))/diff(range(X)))
colnames(cent.mat.scaled)=gsub(".betweenness","",gsub(".degree_out","",gsub(".degree_in","",colnames(cent.mat.scaled))))

TF.df=as.data.frame(rownames(cent.mat.scaled))
colnames(TF.df)=c("Gene")
#use real TF list
TF.df$TF=ifelse(TF.df$Gene %in% TF.net$TF, "Yes","No")
#row_ha = rowAnnotation(isTF = TF.df$TF,col = list(bar = c("Yes" = "red", "No" = "green", "c" = "blue")) )


markers=read.table("~/work/scNET_manuscript/data/marker_genes.txt", header=T)
markers$Cluster=gsub("Oligo","Oli",markers$Cluster)
markers$Cluster=gsub("Microglia","Mic",markers$Cluster)

markers.df=as.data.frame(rownames(cent.mat.scaled))
colnames(markers.df)="Gene"
markers.df$marker=ifelse(markers.df$Gene %in% markers$Gene, "Yes","No")


library(circlize)

col_fun = colorRamp2(c(0, 1), c("white", "red"))

#NOT WORKING
#https://github.com/jokergoo/ComplexHeatmap/issues/148
vip_colors = c("red","blue")
nr_TF=sort(unique(TF.df$TF))
names(vip_colors) = nr_TF
colz = list(TF = vip_colors)


#pdf(file="centrality.heatmap.pdf")

p=Heatmap(cent.mat.scaled, col=col_fun, column_dend_reorder = FALSE,
column_split = rep(c("betweenness","in degree","out degree"), 8),
  column_gap = unit(5, "mm"),border = TRUE,row_names_gp = gpar(fontsize = 5),name="centrality score") +
  rowAnnotation(isTF = TF.df$TF, isMarker=markers.df$marker,hei = anno_density(as.matrix(mic.ad), type = "heatmap", width = unit(6, "cm")))

#draw(p)
#dev.off()


#calculate fold change in centrality

cent.mat=old
rownames(cent.mat)=cent.mat$gene
cent.mat$gene=NULL
cent.mat[is.na(cent.mat)]=0
cent.mat.scaled=apply(cent.mat, MARGIN = 2, FUN = function(X) (X - min(X))/diff(range(X)))


#enrichment analysis

data=GSA.read.gmt('~/work/scNET_manuscript/genome/genesets/GO_annotations-9606-inferred-allev.gmt')
genesets=data$genesets
names(genesets)=data$geneset.descriptions


diff.cent.enrich.tbl=data.frame("label"=NULL,"pval"=NULL,"fdr"=NULL,"signature"=NULL,"geneset"=NULL,
"overlap"=NULL,"background"=NULL,"hits"=NULL,"cell"=NULL )

for(i in 1:length(celltypes))
{
  tag=celltypes[i]
  df=cent.mat.scaled[,colnames(cent.mat.scaled) %like% tag]
  df=df[,colnames(df) %like% "degree_in"]
  colnames(df)[1]=ifelse(colnames(df)[1] %like% "AD","AD","Ctrl")
  colnames(df)[2]=ifelse(colnames(df)[2] %like% "AD","AD","Ctrl")
  min=min(df[df > 0])
  df[df == 0]<-(min*0.01)
  d = transform(df, lfc = log2(df[,colnames(df)%like% "AD"]/ df[,colnames(df)%like% "Ctrl"]))
  d$lfc=abs(d$lfc)
  d=d[order(-d$lfc),]
#  colnames(d)=tag
  name=paste(tag,"degreeIn.df",sep=".")
  assign(name,d)
  signature=rownames(d[d$lfc>0,])
  bkgrnd=rownames(d)
  genesets= msigdb_gsets("Homo sapiens", "C2", "CP:KEGG")
  #genesets=msigdb_gsets("Homo sapiens", "H", "")
  hyp_obj = hypeR(signature, genesets, background=rownames(d),fdr=0.1)
  hyp_df =  hyp_obj$data
  name=paste(tag,"degreeIn.GO.hyperGeo",sep=".")
  assign(name, hyp_obj)
  if(nrow(hyp_df) > 0)
  {
    hyp_df$cell = tag
    diff.cent.enrich.tbl=rbind(diff.cent.enrich.tbl,  hyp_df)
  }
}

data=diff.cent.enrich.tbl
data$label=gsub("_"," ",data$label)
p=ggplot(data, aes(y=label, x=cell)) + geom_point(aes(size=-log10(fdr))) + theme(axis.text.x=element_text(angle=45,vjust=1,hjust=1))
#pdf(file="Diff.centrality.kegg.enrichment.pdf")
#p
#dev.off()


#scatter plots
pattern="*.degreeIn.df"
list=ls(pattern=pattern)


for(i in 1:length(list))
{

  tag=gsub(".degreeIn.df","",list[i])
  p=plot_scatter(get(list[i]),tag,3)
  name=paste(tag,"scatter.pdf",sep=".")
  ggsave(name, p, device="pdf")
}
### automate this
#tmp=cbind(Ex.degreeIn.df,In.degreeIn.df)
#tmp=cbind(tmp,Mic.degreeIn.df)
#tmp=cbind(tmp,Oli.degreeIn.df)

#write.table(data.frame("Gene"=rownames(tmp),tmp),"Differential_inDegree.mat", row.names=FALSE, sep="\t", quote=F)
