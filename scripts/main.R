source('~/work/scNET/scripts/functions_for_network_analysis.R')

ex.list=list(ex_genie,ex_grnbst,ex_pidc)

scNET=function(inlist,tag,nedges)
{
  consensus=ara(inlist,tag,nedges) #consensus network
  pr=get_centrality(consensus,"pr",tag)
  betweenness=get_centrality(consensus,"betweenness",tag)
  hub_score=get_centrality(consensus,"hub_score",tag)
  closeness=get_centrality(consensus,"closeness",tag)
  degree=get_centrality(consensus,"degree",tag)
  target_pairs=find_target_pairs_matrix(consensus,0.5) #0.5 is arbitrary; test a few values
  modules=detect_modules(target_pairs)
  outlist=list(consensus=consensus,pr=pr,bet=betweenness,deg=degree,hs=hub_score,cl=closeness,tp=target_pairs,mod=modules)
  return(outlist)
}


ex.target.graph=get_target_graph(ex.consensus,0.5)
ex.igraph.cluster=cluster_fast_greedy(ex.target.graph)
ex.modularity=modularity(ex.target.graph, membership(ex.igraph.cluster))
ex.density=edge_density(ex.target.graph)
