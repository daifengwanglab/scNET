

scNET=function(inlist,tag,nedges)
{
  if(length(inlist)<2)
  {
    consensus=rank_network_edges(as.data.frame(list[1]),nedges,tag)
    consensus=consensus %>% separate(edge, c("TF","target"),sep="-")
  }
  else
  {
    consensus=ara(inlist,tag,nedges) #consensus network
  }
  pr=get_centrality(consensus,"pr",tag)
  betweenness=get_centrality(consensus,"betweenness",tag)
  hub_score=get_centrality(consensus,"hub_score",tag)
  closeness=get_centrality(consensus,"closeness",tag)
  degree=get_centrality(consensus,"degree",tag)
  concomp=get_centrality(consensus,"components",tag)
  target_pairs=find_target_pairs_matrix(consensus,0.5) #0.5 is arbitrary; test a few values
  modules=detect_modules(target_pairs)
  outlist=list(consensus=consensus,pagerank=pr,betweenness=betweenness,degree=degree,
  hub_score=hub_score,closeness=closeness,concomp=concomp,
  target_pairs=target_pairs,modules=modules)
  print(paste("creating output for: ",celltypes[i]))
  return(outlist)
}
