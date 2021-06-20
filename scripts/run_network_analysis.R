
rm(list=ls())

source('~/work/scNET-devel/scripts/load_libraries.R')
source('~/work/scNET-devel/scripts/read_data.R')
source('~/work/scNET-devel/scripts/functions_for_network_analysis.R')
source('~/work/scNET-devel/scripts/main.R')


result_dir="~/work/scNET-devel/results"

celltypes=c("Ex1","Ex2","Ex3e","Ex4","Ex5b","Ex6a","Ex6b","Ex8","Ex9",
"In1a","In1b","In1c","In3","In4a","In4b","In6a","In6b","In7","In8",
"Mic","Oli")


#clean dir
rm(list=ls(pattern="*.scNET.out"))

#network analysis
for(i in 1:length(celltypes))
{
  tag=celltypes[i]
  outname=paste(celltypes[i],"scNET.out",sep=".")
  assign(outname, scNET(lapply(ls(pattern=celltypes[i]), get),celltypes[i],nedges))
}

#write data to disk
for(i in 1:length(celltypes))
{
  m=paste(celltypes[i],".scNET.out",sep="")
  tag=celltypes[i]
  modules=get(m)$modules
  modules$moduleID=paste(tag,modules$moduleID,sep="_")
  name=paste(celltypes[i],"modules",sep="_")
  assign(name,modules)
  filename=paste(name,"genesets",sep=".")
  write.table(get(name), file=paste(result_dir, filename, sep="/"), row.names=FALSE, col.names=FALSE, sep="\t", quote=F)
  filename=paste(name,"genesets.bg",sep=".")
  write.table(as.data.frame(unique(get(name)$gene)),file=paste(result_dir, filename, sep="/"),row.names=FALSE, col.names=FALSE, sep="\t", quote=F)
  filename=paste(name,"genesets.desc",sep=".")
  df=as.data.frame(unique(get(name)$moduleID))
  colnames(df)=c("id")
  df$name=df$id
  write.table(df,file=paste(result_dir, filename, sep="/"),row.names=FALSE, col.names=FALSE, sep="\t", quote=F)
}


rm(df)
df <- data.frame(GO.ID=as.Date(character()),
                 Term=character(),
                 KS=numeric(),
                celltype=character(),
                 stringsAsFactors=FALSE)

for(i in 1:length(celltypes))
{
  m=paste(celltypes[i],".scNET.out",sep="")
  tag=celltypes[i]
  pr=get(m)$degree
  gs.df=gsea(pr,2,celltypes[i])
  df=rbind(df,gs.df)
}

df.all=df
#Then run again with threshold.


GOBP.degree=df.all %>% filter (Term %in% df$Term)
tmp=acast(GOBP.degree, Term~celltype, value.var="KS" )
tmp[is.na(tmp)]=0
GOBP.degree=tmp


GOBP.pr=df.all %>% filter (Term %in% df$Term)
tmp=acast(GOBP.pr, Term~celltype, value.var="KS" )
tmp[is.na(tmp)]=0
GOBP.pr=tmp

ggplot(GOBP.pr,mapping=aes(x=celltype, y=Term, fill=-log(KS))) + geom_tile()



library(reshape2)
library(tidyr)
library(pheatmap)

dat=read.table("module_overlap.tbl")
colnames(dat)=c("x","y","z")

tmp=acast(dat, x~y, value.var="z")
tmp[is.na(tmp)]=0
annotation1=as.data.frame(colnames(tmp))
annotation2=as.data.frame(rownames(tmp))
 colnames(annotation1)=c("module")
 colnames(annotation2)=c("module")
annotation=rbind(annotation1,annotation2)
list=unique(annotation$module)
annotation=as.data.frame(list)
annotation$celltype=annotation$list
annotation=separate(annotation, celltype, into=c("a","b"))[,1:2]
colnames(annotation)=c("module","celltype")

pheatmap(mat.new,cellwidth=2,cellheight=2,show_rownames=FALSE,show_colnames=FALSE,filename=NA)



#loregic
#load ge matrix

for(i in 1:length(celltypes))
{
  tmp=as.data.frame(gexpr) %>% dplyr::select(starts_with(celltypes[i]))
  filename=paste(celltypes[i],"gexpr.avg",sep=".")
  data=as.data.frame(rowSums(tmp)/dim(tmp)[2])
  colnames(data)=celltypes[i]
  assign(filename, data)
}
gexpr.avg=bind_cols(Ex1.gexpr.avg,Ex2.gexpr.avg,Ex3e.gexpr.avg,Ex4.gexpr.avg,
Ex5b.gexpr.avg,Ex6a.gexpr.avg,Ex6b.gexpr.avg,Ex8.gexpr.avg,
Ex9.gexpr.avg,In1a.gexpr.avg,In1b.gexpr.avg,In1c.gexpr.avg,
In3.gexpr.avg,In4a.gexpr.avg,In4b.gexpr.avg,In6a.gexpr.avg,
In6b.gexpr.avg,In7.gexpr.avg,In8.gexpr.avg,Mic.gexpr.avg,
Oli.gexpr.avg)

gexpr.avg.bin=binarize.array(gexpr.avg)

gexpr.avg.scaled=scale(gexpr.avg)
gexpr.avg.scaled.bin=binarize.array(gexpr.avg.scaled)
source('scripts/run_loregic.R')



#x
ex.target.graph=get_target_graph(ex.consensus,0.5)
ex.igraph.cluster=cluster_fast_greedy(ex.target.graph)
ex.modularity=modularity(ex.target.graph, membership(ex.igraph.cluster))
ex.density=edge_density(ex.target.graph)



#create a list of networks for each cell type
#Celltype1:
ex.list=list(ex_genie,ex_grnbst,ex_pidc)
ex.consensus=ara(ex.list,"ex",nedges) #consensus network
#ex.loregic.out=edgelist2triplets(ex.consensus[,1:2])
ex.pr=get_centrality(ex.consensus,"pr","ex")
ex.betweenness=get_centrality(ex.consensus,"betweenness","ex")
ex.hub_score=get_centrality(ex.consensus,"hub_score","ex")
#ex.closeness=get_centrality(ex.consensus,"closeness","ex")
#ex.loregic.pr=calculate_triplet_hubScores(ex.loregic.out,ex.pr)
ex.target_pairs=find_target_pairs_matrix(ex.consensus,0.5) #0.5 is arbitrary; test a few values
ex.modules=detect_modules(ex.target_pairs)

ex.target.graph=get_target_graph(ex.consensus,0.5)
ex.igraph.cluster=cluster_fast_greedy(ex.target.graph)
ex.modularity=modularity(ex.target.graph, membership(ex.igraph.cluster))
ex.density=edge_density(ex.target.graph)

###############################################

#Celltype2:
in.list=list(in_genie,in_grnbst,in_pidc)
in.consensus=ara(in.list,"in",nedges) #consensus network
#in.loregic.out=edgelist2triplets(in.consensus[,1:2])
in.pr=get_centrality(in.consensus,"pr","in")
in.betweenness=get_centrality(in.consensus,"betweenness","in")
in.hub_score=get_centrality(in.consensus,"hub_score","in")
#in.closeness=get_centrality(in.consensus,"closeness","in")
#in.loregic.pr=calculate_triplet_hubScores(in.loregic.out,in.pr)
in.target_pairs=find_target_pairs_matrix(in.consensus,0.5)
in.modules=detect_modules(in.target_pairs)


in.target.graph=get_target_graph(in.consensus,0.5)
in.igraph.cluster=cluster_fast_greedy(in.target.graph)
in.modularity=modularity(in.target.graph, membership(in.igraph.cluster))
in.density=edge_density(in.target.graph)



############
oligo.list=list(oligo_genie,oligo_grnbst,oligo_pidc)
oligo.consensus=ara(oligo.list,"oligo",nedges) #consensus network
#oligo.loregic.out=edgelist2triplets(oligo.consensus[,1:2])
oligo.pr=get_centrality(oligo.consensus,"pr","oligo")
oligo.betweenness=get_centrality(oligo.consensus,"betweenness","oligo")
oligo.hub_score=get_centrality(oligo.consensus,"hub_score","oligo")
#oligo.closeness=get_centrality(oligo.consensus,"closeness","oligo")
#oligo.loregic.pr=calculate_triplet_hubScores(oligo.loregic.out,oligo.pr)
oligo.target_pairs=find_target_pairs_matrix(oligo.consensus,0.5)
oligo.modules=detect_modules(oligo.target_pairs)

oligo.target.graph=get_target_graph(oligo.consensus,0.5)
oligo.igraph.cluster=cluster_fast_greedy(oligo.target.graph)
oligo.modularity=modularity(oligo.target.graph, membership(oligo.igraph.cluster))

oligo.density=edge_density(oligo.target.graph)

####################
micro.list=list(micro_genie,micro_grnbst,micro_pidc)
micro.consensus=ara(micro.list,"micro",nedges) #consensus network
#micro.loregic.out=edgelist2triplets(micro.consensus[,1:2])
micro.pr=get_centrality(micro.consensus,"pr","micro")
micro.betweenness=get_centrality(micro.consensus,"betweenness","micro")
micro.hub_score=get_centrality(micro.consensus,"hub_score","micro")
#micro.closeness=get_centrality(micro.consensus,"closeness","micro")
#micro.loregic.pr=calculate_triplet_hubScores(micro.loregic.out,micro.pr)
micro.target_pairs=find_target_pairs_matrix(micro.consensus,0.5)
micro.modules=detect_modules(micro.target_pairs)


micro.target.graph=get_target_graph(micro.consensus,0.5)
micro.igraph.cluster=cluster_fast_greedy(micro.target.graph)
micro.modularity=modularity(micro.target.graph, membership(micro.igraph.cluster))

micro.density=edge_density(micro.target.graph)

mean_distance(micro.target.graph, directed=FALSE)
mean_distance(oligo.target.graph, directed=FALSE)
mean_distance(in.target.graph, directed=FALSE)
mean_distance(ex.target.graph, directed=FALSE)

#######combine centralities for functional enrichment
merge_mat=function(matrix)
{
matrix[is.na(matrix)]=0
rownames(matrix)=matrix$gene
matrix=matrix[,2:5]
matrix=scale(matrix)
rownames(matrix)=gsub(" ","",rownames(matrix))

return(matrix)
}


closeness.all=merge(in.closeness,ex.closeness, all=TRUE) %>% merge(oligo.closeness, all=TRUE) %>% merge(micro.closeness, all=TRUE)
closeness.all=merge_mat(closeness.all)
write.table(data.frame("Gene"=rownames(closeness.all),closeness.all),file="results/closeness.mat", sep="\t", quote=FALSE, row.names=FALSE)

pr.all=merge(in.pr,ex.pr, all=TRUE) %>% merge(oligo.pr, all=TRUE) %>% merge(micro.pr, all=TRUE)
pr.all=merge_mat(pr.all)
write.table(data.frame("Gene"=rownames(pr.all),pr.all),file="results/pr.mat", sep="\t", quote=FALSE, row.names=FALSE)

hub_score.all=merge(in.hub_score,ex.hub_score, all=TRUE) %>% merge(oligo.hub_score, all=TRUE) %>% merge(micro.hub_score, all=TRUE)


centrality.matrix=merge(closeness.all,pr.all, all=TRUE, by="gene") %>%  merge(hub_score.all,all=TRUE,by="gene")
centrality.matrix[is.na(centrality.matrix)]=0
rownames(centrality.matrix)=centrality.matrix$gene
centrality.matrix=centrality.matrix[,2:13]
centrality.matrix=scale(centrality.matrix)
rownames(centrality.matrix)=gsub(" ","",rownames(centrality.matrix))

write.table(centrality.matrix, file="results/centrality_matrix.mat", sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)



####
#Store results
#write.table(ex.loregic.out, file="results/ex.loregic.triplets.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
#write.table(in.loregic.out, file="results/in.loregic.triplets.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
#write.table(micro.loregic.out, file="results/micro.loregic.triplets.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
#write.table(oligo.loregic.out, file="results/oligo.loregic.triplets.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

write.table(oligo.modules, file="results/oligo.WGCNA.consensus.modules.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(micro.modules, file="results/micro.WGCNA.consensus.modules.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(ex.modules, file="results/ex.WGCNA.consensus.modules.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(in.modules, file="results/in.WGCNA.consensus.modules.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

write.table(oligo.consensus, file="results/oligo.consensus.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(micro.consensus, file="results/micro.consensus.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(ex.consensus, file="results/ex.consensus.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(in.consensus, file="results/in.consensus.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
####
