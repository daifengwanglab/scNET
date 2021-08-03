
source('~/work/scNET-devel/scripts/functions_for_network_analysis.R')
source('~/work/scNET-devel/scripts/load_libraries.R')


#normalize cent mat
cent.mat=read.table('centralities/centrality_matrix.all-cellTypes.txt', header=T) #first col is gene name
rownames(cent.mat)=cent.mat$gene
cent.mat$gene=NULL
cent.mat[is.na(cent.mat)]=0
cent.mat.scaled=apply(cent.mat, MARGIN = 2, FUN = function(X) (X - min(X))/diff(range(X)))

#GO enrichment analysis
GOdata=GSA.read.gmt('~/work/scNET_manuscript/genome/genesets/GO_annotations-9606-inferred-allev.gmt')
genesets=GOdata$genesets
names(genesets)=GOdata$geneset.descriptions

celltypes=c("Mic","Oli","Ex","In")
diff.cent.enrich.tbl=data.frame("label"=NULL,"pval"=NULL,"fdr"=NULL,"signature"=NULL,"geneset"=NULL,
"overlap"=NULL,"background"=NULL,"hits"=NULL,"cell"=NULL,"centrality"=NULL)

for(i in 1:length(celltypes))
{
  tag=celltypes[i]
  df=cent.mat.scaled[,colnames(cent.mat.scaled) %like% tag]
  df=df[,colnames(df) %like% "degree_in"]
  colnames(df)[1]=ifelse(colnames(df)[1] %like% "AD","AD","Ctrl")
  colnames(df)[2]=ifelse(colnames(df)[2] %like% "AD","AD","Ctrl")
  min=min(df[df > 0])
  df[df == 0]=(min*0.01)
  d = transform(df, lfc = log2(df[,colnames(df)%like% "AD"]/ df[,colnames(df)%like% "Ctrl"]))
  d$lfc=abs(d$lfc)
  d=d[order(-d$lfc),]
  name=paste(tag,"degreeIn.df",sep=".")
  assign(name,d)
  signature=rownames(d[d$lfc>0,])
  bkgrnd=rownames(d)
  genesets=genesets
  hyp_obj = hypeR(signature, genesets, background=rownames(cent.mat.scaled),fdr=0.05)
  hyp_df =  hyp_obj$data
  if(nrow(hyp_df) > 0)
  {
    hyp_df$cell = tag
    hyp_df$centrality="degreeIn"
    diff.cent.enrich.tbl=rbind(diff.cent.enrich.tbl,  hyp_df)
  }
}

data=diff.cent.enrich.tbl
data$label=gsub("_"," ",data$label)
p.GO.degreeIn.diffCent=ggplot(data, aes(y=label, x=cell)) + geom_point(aes(size=-log10(fdr))) +
labs(x="cell type",y="Biological process")+ggtitle(expression(Delta ~ "In degree" ))+
theme(text = element_text(size = 10)) +
theme_bw(base_size=12)
pdf(file="Figures/p.GO.degreeIn.diffCent.pdf")
p.GO.degreeIn.diffCent
dev.off()




diff.cent.enrich.tbl=data.frame("label"=NULL,"pval"=NULL,"fdr"=NULL,"signature"=NULL,"geneset"=NULL,
"overlap"=NULL,"background"=NULL,"hits"=NULL,"cell"=NULL,"centrality"=NULL)
for(i in 1:length(celltypes))
{
  tag=celltypes[i]
  df=cent.mat.scaled[,colnames(cent.mat.scaled) %like% tag]
  df=df[,colnames(df) %like% "degree_out"]
  colnames(df)[1]=ifelse(colnames(df)[1] %like% "AD","AD","Ctrl")
  colnames(df)[2]=ifelse(colnames(df)[2] %like% "AD","AD","Ctrl")
  min=min(df[df > 0])
  df[df == 0]=(min*0.01)
  d = transform(df, lfc = log2(df[,colnames(df)%like% "AD"]/ df[,colnames(df)%like% "Ctrl"]))
  d$lfc=abs(d$lfc)
  d=d[order(-d$lfc),]
  name=paste(tag,"degreeOut.df",sep=".")
  assign(name,d)
  signature=rownames(d[d$lfc>2,])
  bkgrnd=rownames(d)
  genesets=genesets
  hyp_obj = hypeR(signature, genesets, background=rownames(cent.mat.scaled),fdr=0.01)
  hyp_df =  hyp_obj$data
  if(nrow(hyp_df) > 0)
  {
    hyp_df$cell = tag
    hyp_df$centrality="degreeOut"
    diff.cent.enrich.tbl=rbind(diff.cent.enrich.tbl,  hyp_df)
  }
}

data=diff.cent.enrich.tbl
data$label=gsub("_"," ",data$label)
p.GO.degreeOut.diffCent=ggplot(data, aes(y=label, x=cell)) + geom_point(aes(size=-log10(fdr))) +
labs(x="cell type",y="Biological process")+ggtitle(expression(Delta ~ "Out degree" ))+
theme(text = element_text(size = 10)) +
theme_bw(base_size=12)
pdf(file="Figures/p.GO.degreeOut.diffCent.pdf")
p.GO.degreeOut.diffCent
dev.off()



diff.cent.enrich.tbl=data.frame("label"=NULL,"pval"=NULL,"fdr"=NULL,"signature"=NULL,"geneset"=NULL,
"overlap"=NULL,"background"=NULL,"hits"=NULL,"cell"=NULL,"centrality"=NULL)
for(i in 1:length(celltypes))
{
  tag=celltypes[i]
  df=cent.mat.scaled[,colnames(cent.mat.scaled) %like% tag]
  df=df[,colnames(df) %like% "betweenness"]
  colnames(df)[1]=ifelse(colnames(df)[1] %like% "AD","AD","Ctrl")
  colnames(df)[2]=ifelse(colnames(df)[2] %like% "AD","AD","Ctrl")
  min=min(df[df > 0])
  df[df == 0]=(min*0.01)
  d = transform(df, lfc = log2(df[,colnames(df)%like% "AD"]/ df[,colnames(df)%like% "Ctrl"]))
  d$lfc=abs(d$lfc)
  d=d[order(-d$lfc),]
  name=paste(tag,"betweenness.df",sep=".")
  assign(name,d)
  signature=rownames(d[d$lfc>0,])
  bkgrnd=rownames(d)
  genesets=genesets
  hyp_obj = hypeR(signature, genesets, background=rownames(cent.mat.scaled),fdr=0.01)
  hyp_df =  hyp_obj$data
  if(nrow(hyp_df) > 0)
  {
    hyp_df$cell = tag
    hyp_df$centrality="betweenness"
    diff.cent.enrich.tbl=rbind(diff.cent.enrich.tbl,  hyp_df)
  }
}

data=diff.cent.enrich.tbl
data$label=gsub("_"," ",data$label)
p.GO.bet.diffCent=ggplot(data, aes(y=label, x=cell)) + geom_point(aes(size=-log10(fdr))) +
labs(x="cell type",y="Biological process")+ggtitle(expression(Delta ~ "Betweenness" ))+
theme(text = element_text(size = 10)) +
theme_bw(base_size=12)
pdf(file="Figures/p.GO.bet.diffCent.pdf")
p.GO.bet.diffCent
dev.off()



#scatter plots

df.for_scatter=data.frame(AD=NULL,Ctrl=NULL, lfc=NULL, cell=NULL, cent=NULL)
pattern="*.degreeIn.df"
list=ls(pattern=pattern)

for(i in 1:length(list))
{
  cell=gsub(".df","",list[i])
  cell=gsub(".degreeIn","",cell)
  d=get(list[i])
  d$Name=rownames(d)
  d$lfc=log2(d[,1]/d[,2])
  d$cell=cell
  d$cent="In degree"
  df.for_scatter=rbind(df.for_scatter,d)

  #p=plot_scatter(get(list[i]),tag,1)
  #name=paste(tag,"scatter.pdf",sep=".")
  #name=paste("Figures",name,sep="/")
  #ggsave(name, p, device="pdf")
}

pattern="*.degreeOut.df"
list=ls(pattern=pattern)

for(i in 1:length(list))
{
cell=gsub(".df","",list[i])
cell=gsub(".degreeOut","",cell)
d=get(list[i])
d$Name=rownames(d)
d$lfc=log2(d[,1]/d[,2])
d$cell=cell
d$cent="Out degree"
df.for_scatter=rbind(df.for_scatter,d)

  #tag=gsub(".df","",list[i])
  #p=plot_scatter(get(list[i]),tag,2)
  #name=paste(tag,"scatter.pdf",sep=".")
  #name=paste("Figures",name,sep="/")
  #ggsave(name, p, device="pdf")
}

pattern="*.betweenness.df"
list=ls(pattern=pattern)

for(i in 1:length(list))
{
cell=gsub(".df","",list[i])
cell=gsub(".betweenness","",cell)
d=get(list[i])
d$Name=rownames(d)
d$lfc=log2(d[,1]/d[,2])
d$cell=cell
d$cent="Betweenness"
df.for_scatter=rbind(df.for_scatter,d)
}


p.df.for_scatter=ggplot(df.for_scatter, aes(x= AD, y = Ctrl)) +
geom_point(color = "black",
             size = 1, alpha = 0.8) + ggtitle("")+ facet_grid(cell ~ cent)+
             labs(x="Normalized centrality in AD", y="Normalized centrality in control")+
						 theme(text = element_text(size = 10)) +
						 theme_bw(base_size=12)


pdf(file="Figures/p.cent.scatter.pdf")
p.df.for_scatter
dev.off()





#plot feature matrix

TF.df=as.data.frame(rownames(cent.mat.scaled))
colnames(TF.df)=c("Gene")
#use real TF list
TF.df$TF=ifelse(TF.df$Gene %in% TF.net$TF, "TF","Target")

#for heirarchy
data.hei=read.table("~/work/scNET_manuscript/AD_MIT/heirarchy/all_genes_heirarchy_levs.txt", header=T)
data.hei=as.data.frame(lapply(data.hei, function(y) gsub("Lev1", "Bottom", y)))
data.hei=as.data.frame(lapply(data.hei, function(y) gsub("Lev6", "Top", y)))
data.hei=as.data.frame(lapply(data.hei, function(y) gsub("Lev2", "Middle", y)))
data.hei=as.data.frame(lapply(data.hei, function(y) gsub("Lev3", "Middle", y)))
data.hei=as.data.frame(lapply(data.hei, function(y) gsub("Lev4", "Middle", y)))
data.hei=as.data.frame(lapply(data.hei, function(y) gsub("Lev5", "Middle", y)))


markers=read.table("~/work/scNET_manuscript/data/marker_genes.txt", header=T)
markers$Cluster=gsub("Oligo","Oli",markers$Cluster)
markers$Cluster=gsub("Microglia","Mic",markers$Cluster)

markers.df=as.data.frame(rownames(cent.mat.scaled))
colnames(markers.df)="Gene"
markers.df$marker=ifelse(markers.df$Gene %in% markers$Gene, "Yes","No")


library(circlize)

npgcolors=pal_npg("nrc", alpha = 1)(10)
troncolors=pal_tron("legacy", alpha = 1)(7)

col_fun = c(troncolors[5], npgcolors[8])


colnames(cent.mat.scaled)=gsub(".betweenness","",gsub(".degree_out","",gsub(".degree_in","",colnames(cent.mat.scaled))))


p.heatmap=Heatmap(cent.mat.scaled, col=col_fun, column_dend_reorder = FALSE,show_row_dend=FALSE,
column_split = rep(c("betweenness","in degree","out degree"), 8),
  column_gap = unit(5, "mm"),border = TRUE,row_names_gp = gpar(fontsize = 5),name="centrality score") +
  rowAnnotation(isTF = TF.df$TF,
  hei.Ex=data.hei$Ex.AD,
  hei.In=data.hei$In.AD,
  hei.Oli=data.hei$Oli.AD,
  hei.Mic=data.hei$Mic.AD,
  width = unit(6, "cm"),border = TRUE,
  col = list(isTF = c("TF" = npgcolors[9], "Target" = npgcolors[10]),
        hei.Ex=c("Bottom"=npgcolors[3],"Top"=npgcolors[4],"Middle"=npgcolors[5],"None"="white"),
        hei.In=c("Bottom"=npgcolors[3],"Top"=npgcolors[4],"Middle"=npgcolors[5],"None"="white"),
        hei.Oli=c("Bottom"=npgcolors[3],"Top"=npgcolors[4],"Middle"=npgcolors[5],"None"="white"),
        hei.Mic=c("Bottom"=npgcolors[3],"Top"=npgcolors[4],"Middle"=npgcolors[5],"None"="white")
    )
  )
#p.heatmap


#pdf(file="Figures/centrality_heatmap.pdf")
#draw(p.heatmap, heatmap_legend_side = "left", annotation_legend_side = "bottom")
#dev.off()

##########################################################
