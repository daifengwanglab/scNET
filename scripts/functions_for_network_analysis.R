

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


find_target_pairs_matrix=function(net) #network
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
	jacc
}


detect_modules = function(matrix, msize)
{

	dissMatrix = 1 - matrix
	# Call the hierarchical clustering function
	geneTree = flashClust(as.dist(dissMatrix),method="average");
	minModuleSize = msize;
	# Module identification using dynamic tree cut:
	dynamicMods = cutreeDynamic(dendro = geneTree,  method="tree", minClusterSize = minModuleSize)
	dynamicColors = labels2colors(dynamicMods)
	#restGenes= (dynamicColors != "grey")
	#diag(dissMatrix) = NA
	#TOMplot(dissMatrix, geneTree, as.character(dynamicColors))

	modules=cbind(as.data.frame(dynamicMods),rownames(matrix))
	modules=modules[,c(2,1)]
	colnames(modules)=c("gene","moduleID")
	modules
}
