# Regulatory network format;: 3 columns (TF target score)
rm(list=ls())
source('~/work/scNET-devel/scripts/load_libraries.R')
source('~/work/scNET-devel/scripts/functions_for_network_analysis.R')
my_entrez_gene_info=read.table('~/work/scNET_manuscript/genome/hgncSymbols_to_entrez.human.txt', header=T,sep="\t")
colnames(my_entrez_gene_info)=c("entrezID","gene")

cent.mat=read.table("centrality_matrix.all-cellTypes.txt", header=T)

for (i in 2:ncol(cent.mat))
{

  name=paste(colnames(cent.mat)[i])
  name=paste(name,"GOBPgseaRes",sep=".")
  df=cent.mat[,c(1,i)]
  assign(name, gsea(df,tag))
}


centralities=c("betweenness","degree_out","degree_in")
for (i in 1:length(centralities))
{
  c=centralities[i]
  pattern=paste(c,"GOBPgseaRes",sep=".")
  listGseaDfs=ls(pattern=pattern)
  df=as.data.frame(lapply(listGseaDfs[1],get))
  for(i in 2:length(listGseaDfs))
  {
    df=rbind(df,as.data.frame(lapply(listGseaDfs[i],get)))
  }
  data=df[complete.cases(df), ]
  #mat=xtabs(-1*log(KS)~celltype+Term, data=data)
  data$KS=-1*log(data$KS)
  mat=acast(data, celltype~Term, value.var="KS")
  mat[is.na(mat)]=0 #set NA =0
  filename=paste(c,"GOBPgsea.pdf",sep=".")
  p=Heatmap(mat,  row_split = 4,  row_title = c,column_names_gp = gpar(fontsize = 8))
  pdf(file=filename)
  p
  dev.off()

}
