# scNET: single cell network biology for understanding cell-type functional genomics


The following outlines the pipeline implemented for post-inference analysis of brain cell type specific gene regulatory networks. We developed this pipeline to highlight (and assign a biological meaning to) similarities and differences in global and local network topologies in a bunch of single cell gene regulatory networks.

If you wish to use this pipeline, start with a GRN inferred using any method of choice as input. Going with the trend, we have also added a module that rank aggregates networks inferred using different algorithm into a single consensus network for each cell type.

**Citation**

If you use this pipeline, please cite: *(under development)*

## Download code
The code has been tested on R version 4.0.3 on Mac and Linux OS.

```r
git clone https://github.com/cngupta/scNET.git
cd scNET
mkdir results
```
___

## Example
The first step is to load external libraries and functions required to execute the pipeline:

```r
source('scripts/load_libraries.R')
source('scripts/functions_for_network_analysis.R')
```

##### Read demo data
```r
# Four brain cell type networks (these are random networks for demo purposes)
# Each inferred using three different algorithms (GENIE3, PIDC, GRNBoost2)

source('scripts/read_demo_data.R')
```

#### Consensus Network


The first step is to create a consensus network using the `ara` function that implements the average rank aggregation method on
networks inferred from different algorithms. The `ara` function expects as input a list of networks as data frames, the cell type tag, and the total number of edges to retain in the consensus.

This step is optional if you have inferred the network using only one algorithm. Although it is expected that some filtering step is taken by the user to reduce the computational burden on downstream analysis, especially if a large number of cell types are to be analyzed.    

```r
#no. of edges to keep in the consensus; test different values
nedges=1000 #just a random number of demo

#create a list with cell type networks as the elements
ex.list=list(ex_genie,ex_grnbst,ex_pidc)

#call the ara function with the list of data frames, the cell type tag, and the total number of edges as arguments.
ex.consensus=ara(ex.list,"ex",nedges)

```

#### Node Centrality
The `get_pageRank` function calculates node centrality, currently using the Page Rank algorithm.    
```r
#modify to let user choose option (pr, hubscore, degree, betweenness)
ex.pr=get_pageRank(ex.consensus,"ex")
```

#### Regulatory triplets  
The next step is to find triplets in the network. We use the `edgelist2triplets` function of the `Loregic` package to report network motifs comprising of two TFs and a target gene (TF can also be a target gene).
```r
ex.loregic.out=edgelist2triplets(ex.consensus[,1:2])
```

#### Triplet importance
Next, the `calculate_triplet_hubScores` function calculates the importance of every triplet by reporting the sum of page rank scores of each gene in the triplet.
```r
ex.loregic.PRscores=calculate_triplet_hubScores(ex.loregic.out,ex.pr)
```

#### Cluster target genes into modules
To cluster genes into functional modules, the first step is to connect gene-pairs if they have a high overlap between their predicted regulators. The `find_target_pairs` achieves this, based on a user defined threshold as one of the argument. Next, the output of `find_target_pairs` is passed into the function `detect_modules` to cluster genes into modules. The `detect_modules` function uses WGCNA's `TOMsimilarity` function and heirarchical clustering to report modules.
```r

#The threshold of 0.5 indicates 50% overlap between the predicted regulators of every target gene-pair. This threshold should ideally be tested for a range of values.

ex.target_pairs=find_target_pairs(ex.consensus,0.5)
ex.modules=detect_modules(ex.target_pairs)
```

#### Store results
```r

```

Repeat these steps for each cell type in your dataset. Then, follow the steps below
if you wish to compare networks and characterize each cell type.  
#### Compare local and global network topologies across cell types
*(under development)*
