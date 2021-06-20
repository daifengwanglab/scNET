
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
  name=paste(name,"KEGGgseaRes",sep=".")
  df=cent.mat[,c(1,i)]
  assign(name, TopKEGGgsea(df,0.1,my_entrez_gene_info))
}


centralities=c("betweenness","degree_out","degree_in")
for (i in 1:length(centralities))
{
  c=centralities[i]
  pattern=paste(c,"KEGGgseaRes",sep=".")
  listGseaDfs=ls(pattern=pattern)
  df=as.data.frame(lapply(listGseaDfs[1],get))
  for(i in 2:length(listGseaDfs))
  {
    df=rbind(df,as.data.frame(lapply(listGseaDfs[i],get)))
  }
  data=df[complete.cases(df), ]
  #mat=xtabs(qvalue~ct+Description, data=data)
  mat=acast(data, ct~Description, value.var="qvalue")
  mat[is.na(mat)]=0 #set NA =0
  mat_breaks <- seq(min(mat), max(mat), length.out = 10)
  filename=paste(c,"KEGGgseaRes.pdf",sep=".")
  mainTitle=c
  pheatmap::pheatmap(
    mat               = t(mat),
    color             = inferno(length(mat_breaks) - 1),
    breaks            = mat_breaks,
    border_color      = NA,
    show_colnames     = TRUE,
    show_rownames     = TRUE,
    drop_levels       = TRUE,
    fontsize          = 10,
    main              = mainTitle,
    file=filename
  )
}
