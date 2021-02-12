
rm(list=ls())

source('~/work/scNetAnalysis.beta/scripts/load_libraries.R')
source('~/work/scNetAnalysis.beta/scripts/read_data.R')
source('~/work/scNetAnalysis.beta/scripts/functions_for_network_analysis.R')


nedges=200000#no of edges to keep in the consensus

#create a list of networks for each cell type
#Celltype1:

ex.list=list(ex_genie,ex_grnbst,ex_pidc)
ex.consensus=ara(ex.list,"ex",nedges) #consensus network

ex.loregic.out=edgelist2triplets(ex.consensus[,1:2])
ex.pr=get_pageRank(ex.consensus,"ex")
ex.loregic.pr=calculate_triplet_hubScores(ex.loregic.out,ex.pr)

ex.target_pairs=find_target_pairs(ex.consensus,0.5) #0.5 is arbitrary; test a few values
ex.modules=detect_modules(ex.target_pairs)

###############################################

#Celltype2:

in.list=list(in_genie,in_grnbst,in_pidc)
in.consensus=ara(in.list,"in",nedges) #consensus network

in.loregic.out=edgelist2triplets(in.consensus[,1:2])
in.pr=get_pageRank(in.consensus,"in")
in.loregic.pr=calculate_triplet_hubScores(in.loregic.out,in.pr)

in.target_pairs=find_target_pairs(in.consensus,0.5)
in.modules=detect_modules(in.target_pairs)

############


oligo.list=list(oligo_genie,oligo_grnbst,oligo_pidc)
oligo.consensus=ara(oligo.list,"oligo",nedges) #consensus network

oligo.loregic.out=edgelist2triplets(oligo.consensus[,1:2])
oligo.pr=get_pageRank(oligo.consensus,"oligo")
oligo.loregic.pr=calculate_triplet_hubScores(oligo.loregic.out,oligo.pr)

oligo.target_pairs=find_target_pairs(oligo.consensus,0.5)
oligo.modules=detect_modules(oligo.target_pairs)

####################


micro.list=list(micro_genie,micro_grnbst,micro_pidc)
micro.consensus=ara(micro.list,"micro",nedges) #consensus network

micro.loregic.out=edgelist2triplets(micro.consensus[,1:2])
micro.pr=get_pageRank(micro.consensus,"micro")
micro.loregic.pr=calculate_triplet_hubScores(micro.loregic.out,micro.pr)

micro.target_pairs=find_target_pairs(micro.consensus,0.5)
micro.modules=detect_modules(micro.target_pairs)


###########################
