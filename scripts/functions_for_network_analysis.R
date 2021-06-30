

process_networks = function(list, celltypeTag, th)
{
	newList=list()
	for (i in 1:length(list))
	{
		network=rank_network_edges(as.data.frame(list[i]),th,	celltypeTag)
		colnames(network)=c("edge",paste("rank",i,sep="."))
		name=paste("Network", i, sep = ".")
		newList[[name]]=network
	}
	newList
}

rank_network_edges = function(net,th,tag) #network and threshold for # of edges to report
{
	#rank edges of each network from 1 to n
	colnames(net)=c("regulatoryGene","targetGene","weight")
	net= net %>% filter(regulatoryGene %in% TFs$V1)
	if(nrow(net) < th)
	{
		th=nrow(net)
	}
	net$edge=paste(net$regulatoryGene,"-",net$targetGene)
	net=net[,3:4]
	net=net[,c(2,1)]
	df=net %>% arrange(-weight) %>% mutate(rank = dense_rank(-weight))
	df=df[1:th,c(1,3)]
	colnames(df)=c("edge",paste(tag,"rank",sep="_"))
	df
}

combine_network_dfs = function(list)
{
	#merge all dfs columwise; NAs if the edge is absent
	df=Reduce(function(x, y) merge(x, y, all=TRUE), list, accumulate=FALSE)
	#replace NAs with a value that is equal to the highest rank +1 for each network
	l=colMax(df)
	for(i in 2:ncol(df))
	{
		n=l[i]
		vars=paste("rank",i-1,sep=".")

		#courtesy https://bit.ly/3pL5Fut
		df=df %>%
			mutate_at(.vars = vars, .funs = funs(ifelse(is.na(.), n, .)))
	}
	rownames(df)=df$edge
	df=df[,2:ncol(df)]
	df[] <- lapply(df, function(x) as.numeric(as.character(x)))
	df
}


colMax = function(data) sapply(data, max, na.rm = TRUE)

ara=function(list, celltypeTag, th) #no. of edges to return
{
	list=list
	tag=celltypeTag
	th=th
	NetList=process_networks(list,celltypeTag,nedges)
	NetworksMatrix=combine_network_dfs(NetList)
	matrix=NetworksMatrix
	ncol=ncol(matrix)
	matrix$mean = rowMeans(matrix)
	matrix=matrix %>% mutate (meanRank = dense_rank(mean)) %>% arrange(meanRank)
	matrix=matrix[1:nedges,]
	matrix$edge=rownames(matrix)
	matrix=matrix[,c("edge","meanRank")]
	rownames(matrix)=c()
	matrix=matrix %>% separate(edge, c("TF","target"),sep="-")
}

get_centrality=function(net, type, tag)
{
	net.igraph=graph_from_data_frame(net, directed = TRUE, vertices = NULL)
	if(type=="pr")
	{
		cent=page_rank(net.igraph)
		cent=as.data.frame(cent$vector)
		colnames(cent)=c("scores")
		cent=cent[order(-cent$scores), , drop = FALSE]
		colnames(cent)=paste(tag,type,sep=".")
		cent$gene=rownames(cent)
		rownames(cent)=c()
		return(cent)
	}
	else if (type=="closeness")
	{
		cent=closeness(net.igraph)
		cent=as.data.frame(cent)
		colnames(cent)=c("scores")
		cent=cent[order(-cent$scores), , drop = FALSE]
		colnames(cent)=paste(tag,type,sep=".")
		cent$gene=rownames(cent)
		rownames(cent)=c()
		return(cent)
	}
	else if (type=="betweenness")
	{
		cent=betweenness(net.igraph, directed=TRUE)
		cent=as.data.frame(cent)
		colnames(cent)=c("scores")
		cent=cent[order(-cent$scores), , drop = FALSE]
		colnames(cent)=paste(tag,type,sep=".")
		cent$gene=rownames(cent)
		rownames(cent)=c()
		return(cent)
	}
	else if (type=="hub_score")
	{
		cent=hub_score(net.igraph)
		cent=as.data.frame(cent$vector)
		colnames(cent)=c("scores")
		cent=cent[order(-cent$scores), , drop = FALSE]
		colnames(cent)=paste(tag,type,sep=".")
		cent$gene=rownames(cent)
		rownames(cent)=c()
		return(cent)
	}
	else if (type=="degree_in")
	{
		cent=	as.data.frame(igraph::degree(net.igraph,mode="in"))
		colnames(cent)=c("scores")
		cent=cent[order(-cent$scores), , drop = FALSE]
		colnames(cent)=paste(tag,type,sep=".")
		cent$gene=rownames(cent)
		rownames(cent)=c()
		return(cent)
	}
	else if (type=="degree_out")
	{
		cent=	as.data.frame(igraph::degree(net.igraph,mode="out"))
		colnames(cent)=c("scores")
		cent=cent[order(-cent$scores), , drop = FALSE]
		colnames(cent)=paste(tag,type,sep=".")
		cent$gene=rownames(cent)
		rownames(cent)=c()
		return(cent)
	}
	else if (type=="components")
	{
	 #no of connected components
		nCC=components(net.igraph)$no
		return(nCC)
	}
}

calculate_triplet_hubScores=function(loregicOut,hubtbl)
{
	colnames(hubtbl)=c("scores","gene")
	new=loregicOut
	new[] <- hubtbl$scores[match(unlist(loregicOut), hubtbl$gene)] #(reference: https://stackoverflow.com/questions/35636315/replace-values-in-a-dataframe-based-on-lookup-table)
	new$RF1=as.numeric(new$RF1)
	new$RF2=as.numeric(new$RF2)
	new$target=as.numeric(new$target)
	new$Sum=rowSums(new)
	new$triplet=paste(loregicOut$RF1,loregicOut$RF2,loregicOut$target, sep="-")
	new=new %>% arrange(-Sum)
	new
}

find_target_pairs_matrix=function(net,th) #network and JI threshold
{
	colnames(net)=c("TF","target","score")
	m=acast(net, TF~target, value.var="score")
	m=t(m)
	m[is.na(m)]=0 #set NA =0
	#find cardinalities
	# Find that paper and add reference
	i12 = m %*% t(m)
	s = diag(i12) %*% matrix(1, ncol = length(diag(i12)))
	u12 = s + t(s) - i12
	jacc= i12/u12
	#genes_jaccard_dist.dat=melt(as.matrix(jacc))
	#target_pairs=genes_jaccard_dist.dat[genes_jaccard_dist.dat[, ncol(genes_jaccard_dist.dat)]>th,]
	#colnames(target_pairs)=c("gene1","gene2","Jaccard")
	#target_pairs$Jaccard=1
	#target_pairs
	jacc[jacc < th] = 0
	jacc[jacc >= th] = 1
	diag(jacc)=1 #required for TOMsimilarity
	jacc
}

get_coregnet_graph=function(net,th,tag) #network, JI threshold, tag
{
	colnames(net)=c("TF","target","score")
	m=acast(net, TF~target, value.var="score")
	m=t(m)
	m[is.na(m)]=0 #set NA =0
	#find cardinalities
	# Find that paper and add reference
	i12 = m %*% t(m)
	s = diag(i12) %*% matrix(1, ncol = length(diag(i12)))
	u12 = s + t(s) - i12
	jacc= i12/u12
	genes_jaccard_dist.dat=melt(as.matrix(jacc))
	target_pairs=genes_jaccard_dist.dat[genes_jaccard_dist.dat[, ncol(genes_jaccard_dist.dat)]>th,]
	colnames(target_pairs)=c("gene1","gene2","Jaccard")
	target_pairs.igraph=graph_from_data_frame(target_pairs, directed=FALSE)
	target_pairs.igraph=simplify(target_pairs.igraph)
	target_pairs.igraph
}

calculate_geneset_density = function(net,geneset) #fullNetwork, listOfQuerygenes
{
	sub=induced.subgraph(net,vids=geneset)
	density=edge_density(sub)
	density
}


detect_modules = function(matrix)
{
	#ref: https://support.bioconductor.org/p/102857/
	#ref:
	TOM = TOMsimilarity(matrix,TOMType="unsigned");
	dissTOM = 1-TOM
	# Call the hierarchical clustering function
	geneTree = flashClust(as.dist(dissTOM),method="average");
	minModuleSize = 100;
	# Module identification using dynamic tree cut:
	dynamicMods = cutreeDynamic(dendro = geneTree,  method="tree", minClusterSize = minModuleSize)
	modules=cbind(as.data.frame(dynamicMods),rownames(matrix))
	modules=modules[,c(2,1)]
	colnames(modules)=c("gene","moduleID")
	modules
}

gsea = function(inputdf,tag)
{
tag=colnames(df)[2]
colnames(inputdf)= c("gene", "score")
inputdf=na.omit(inputdf)
genes=inputdf$score
names(genes)=inputdf$gene
selection <- function(allScore){ return(allScore > 0)}
allGO2genes <- annFUN.org(whichOnto="BP", feasibleGenes=NULL, mapping="org.Hs.eg.db", ID="symbol")
GOdata <- new("topGOdata",
ontology="BP",
allGenes=genes,
annot=annFUN.GO2genes,
GO2genes=allGO2genes,
geneSel=selection,
nodeSize=10)
results.ks <- runTest(GOdata, statistic="KS",algorithm="weight01")
goEnrichment <- GenTable(GOdata, KS=results.ks, orderBy="KS", topNodes=20)
goEnrichment$KS <- as.numeric(goEnrichment$KS)
goEnrichment <- goEnrichment[goEnrichment$KS<0.01,]
goEnrichment <- goEnrichment[goEnrichment$Annotated<200,]
goEnrichment <- goEnrichment %>% arrange(goEnrichment$KS)
goEnrichment <-goEnrichment[1:5,]
goEnrichment <- goEnrichment[,c("GO.ID","Term","KS")]
goEnrichment$Term <- gsub(" [a-z]*\\.\\.\\.$", "", goEnrichment$Term)
#goEnrichment$Term <- gsub("\\.\\.\\.$", "", goEnrichment$Term)
goEnrichment$Term <- paste(goEnrichment$GO.ID, goEnrichment$Term, sep=", ")
goEnrichment$celltype=tag
goEnrichment
}

TopDOgsea = function(df,n,my_entrez_gene_info) #eg n=0.1 for top 10%
{
tag=colnames(df)[2]
colnames(df)=c("gene","score")
df=df[order(-df$score), , drop = FALSE]
topn=round(dim(df)[1]*n)
df=df[1:topn,]
df=df[df$score > 0 ,]
indx=match(df$gene,my_entrez_gene_info$gene)
indx=indx[!is.na(indx)]
df=my_entrez_gene_info[indx,]
query=df$entrezID
x <- enrichDO(gene          = query,
              ont           = "DO",
              pvalueCutoff  = 0.1,
              pAdjustMethod = "BH",
              minGSSize     = 5,
              maxGSSize     = 500,
              qvalueCutoff  = 0.1,
              readable      = TRUE)
DO.df=head(x)

	if(nrow(DO.df)>0)
	{
		DO.df$ct=tag
		DO.df$qvalue=-1*log10(DO.df$qvalue)
		DO.df=DO.df[,c("Description","qvalue","ct")]
		return(DO.df)
	}
}

TopKEGGgsea = function(cent.mat,n,my_entrez_gene_info) #eg n=0.3 for 30%
{
	tag=colnames(df)[2]
	colnames(df)=c("gene","score")
	df=df[order(-df$score), , drop = FALSE]
	topn=round(dim(df)[1]*n)
	df=df[1:topn,]
	df=df[df$score > 0 ,]
	indx=match(df$gene,my_entrez_gene_info$gene)
	indx=indx[!is.na(indx)]
	df=my_entrez_gene_info[indx,]
	gene=df$entrezID
	x <- enrichKEGG(gene          = gene,
              organism     = 'hsa',
							keyType = "kegg",
              pvalueCutoff = 0.05,
							pAdjustMethod = "BH",
							minGSSize = 10,
  						maxGSSize = 500)
	KEGG.df=head(x)
	if(nrow(KEGG.df)>0)
	{
		KEGG.df$ct=tag
		KEGG.df$qvalue=-1*log10(KEGG.df$qvalue)
		KEGG.df=KEGG.df[,c("Description","qvalue","ct")]
		return(KEGG.df)
		}
}


#https://stackoverflow.com/questions/15624656/label-points-in-geom-point
plot_scatter=function(df,tag,th){

d=df
d$Name=rownames(d)
d$lfc=log2(d[,1]/d[,2])
p=ggplot(d, aes(x= AD, y = Ctrl, label = Name)) +
geom_point(color = dplyr::case_when(d$lfc > th ~ "green",
                                      d$lfc < -th ~ "red",
                                      TRUE ~ "#7570b3"),
             size = 3, alpha = 0.8) +
  geom_text_repel(data          = subset(d, lfc > th),
                  nudge_x       = 10 - subset(d, lfc > th)$lfc,
                  size          = 2,
                  box.padding   = 2.5,
                  point.padding = 0.5,
                  force         = 200,
                  segment.size  = 0.2,
                  segment.color = "grey50",
                  direction     = "x") +
  geom_label_repel(data         = subset(d, lfc < -th),
                  nudge_x       = 10 - subset(d, lfc < -th)$lfc,
                  size          = 2,
                  box.padding   = 1,
                  point.padding = 0.5,
                  force         = 200,
                  segment.size  = 0.2,
                  segment.color = "grey50",
                  direction     = "x") +
  theme_classic(base_size = 16)

#name=paste(tag,"scatter.pdf", sep=".")
#pdf(file=name)
p
#dev.off()
}
