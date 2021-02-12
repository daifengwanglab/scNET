# scNET: single cell network biology for understanding cell-type functional genomics


The following outlines the pipeline implemented for post-inference analysis of brain cell type specific gene regulatory networks. We developed this pipeline to highlight (and assign a biological meaning to) similarities and differences in global and local network topologies in a bunch of single cell gene regulatory networks.

If you wish to use this pipeline, start with a GRN inferred using any method of choice as input. Going with the trend, we have also added a module that rank aggregates networks inferred using different algorithm into a single consensus network for each cell type.

**Citation**

If you use this pipeline, please cite:

## Download code
The code has been tested on R version 4.0.3 on Mac and Linux OS.

```r
git clone https://github.com/cngupta/scNET.git
cd scNET
```
___

The main functions perform the following tasks:  


*Steps 1-4 must be applied to every cell type network before proceeding to step 5*

0. (*Optional*) Construct cell-type consensus GRN by **rank aggregating** networks inferred using different algorithms.
  -  Follows the wisdom of crowds approach to infer a consensus GRN


1. Calculate node **centrality**
  - Uses the **Page Rank algorithm** implemented in the igraph library to calculate node centralities.  


2. Identify **network motifs**
  * Uses the **Loregic** algorithm to identify coordinated TF activities as triplets (TF1-TF2-target)  


3. Estimate **Motif importance**   
  * Calculates the sum of node centrality within each triplet


4. **Clustering**
  * Clusters target genes into functional modules based on the similarities within their regulators.    


5. **Compare** cell type networks  
  * Finds common and unique network motifs and modules across cell type Networks

___

## Example
The first step is to load external libraries and functions required to execute the pipeline:

```r
source('scripts/load_libraries.R')
source('scripts/functions_for_network_analysis.R')
```

##### Read demo data
```r
# Four brain cell type networks (microglia, oligodendrocytes, inhibitory and excitatory neurons)

# Each inferred using three different algorithms (GENIE3, PIDC, GRNBoost2)

source('scripts/read_demo_data.R')
```

#### Consensus Network


The first step is to create a consensus network using the `ara` function that implements the average rank aggregation method on
networks inferred from different algorithms. The `ara` function expects as input a list of networks as data frames, the cell type tag, and the total number of edges to retain in the consensus.

This step is optional if you have inferred the network using only one algorithm. Although it is expected that some filtering step is taken by the user to reduce the computational burden on downstream analysis, especially if a large number of cell types are to be analyzed.    

```r
#no. of edges to keep in the consensus; test different values
nedges=200000 #top 10% of all expected edges

#create a list with cell type networks as the elements
ex.list=list(ex_genie,ex_grnbst,ex_pidc)

#call the ara function to crteate the consensus network.
ex.consensus=ara(ex.list,"ex",100000)

```

#### Node Centrality
```r
#modify to let user choose option (pr, hubscore, degree, betweenness)
ex.pr=get_pageRank(ex.consensus,"ex")
```
#### Regulatory triplets  
```r
#two TFs target the same gene
ex.loregic.out=edgelist2triplets(ex.consensus[,1:2])
```
#### Triplet centrality
```r
#Sum centrality score of each node in every triplet
ex.loregic.PRscores=calculate_triplet_hubScores(ex.loregic.out,ex.pr)
```

#### Cluster target genes into modules
```r
#Uses WGCNA's TOMsimilarity function and heirarchical clustering
#threshold of 0.5 is arbitrary in this example. Ideally one should
#test different values

#returns a data frame with three columsn (g1 g2 jaccard)
ex.target_pairs=find_target_pairs(ex.consensus,0.5) #50% overlap

#returns a df with first col as gene and second col as module ID
ex.modules=detect_modules(ex.target_pairs)
```

#### Compare local and global network topologies across cell types
