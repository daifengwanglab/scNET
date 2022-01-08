# Single-cell network biology characterizes cell type gene regulation for drug repurposing and phenotype prediction in Alzheimer’s disease

## Summary
Dysregulation of gene expression in Alzheimer’s disease (AD) remains elusive, especially at the cell type level. Gene regulatory network, a key molecular mechanism linking transcription factors (TFs) and regulatory elements to govern target gene expression, can change across cell types in the human brain and thus serve as a model for studying gene dysregulation in AD. However, it is still challenging to understand how cell type networks work abnormally under AD. To address this, we integrated single-cell multi-omics data and predicted the gene regulatory networks in AD and control for four major cell types, excitatory and inhibitory neurons, microglia and oligodendrocytes. Importantly, we applied network biology approaches to analyze the changes of network characteristics across these cell types, and between AD and control. For instance, many hub TFs target different genes between AD and control (rewiring). Also, these networks show strong hierarchical structures in which top TFs (master regulators) are largely common across cell types, whereas different TFs operate at the middle levels in some cell types (e.g., microglia). The regulatory logics of enriched network motifs (e.g., feed-forward loops) further uncover cell type-specific TF-TF cooperativities in gene regulation. The cell type networks are highly modular and several network modules with cell-type-specific expression changes in AD pathology are enriched with AD-risk genes and putative targets of approved and pending AD drugs, suggesting possible cell-type genomic medicine in AD. Finally, using the cell type gene regulatory networks, we developed machine learning models to classify and prioritize additional AD genes. We found that top prioritized genes predict clinical phenotypes (e.g., cognitive impairment) with reasonable accuracy. Overall, this single-cell network biology analysis provides a comprehensive map linking genes, regulatory networks, cell types and drug targets and reveals dysregulated cell type gene dysregulatory mechanisms in AD.

## Flow chart
![alt text](https://github.com/cngupta/scNET/blob/master/workflow.png)

## System requirements

The analysis is based on R version 4.0.3. For the gene regulatory network, a *Linux* system with 32 GB RAM and 32GB storage would be enough. For network analysis, a standard computer should be enough.

## Software installation

Packages needed for the whole analysis.

- `igraph` (R package, https://igraph.org/r/)
- `mfinder` (motif finding tool from https://www.weizmann.ac.il/mcb/UriAlon/download/network-motif-software)
- `Loregic` (R package for regulatory logics, https://github.com/gersteinlab/Loregic.git)
- `WGCNA` (R package for network module detection, https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/)
- `randomForest` (R package for binary classification)
- other supporting packages should install automatically. Check `load_libraries.R` file in the `scripts` directory

## Download code
The code has been tested on R version 4.0.3 on Linux and Mac OS.
```r
git clone https://github.com/daifengwanglab/scNET
cd scNET
```

## Download data
Download brain cell type gene regulatory networks for four cell types from Zenodo `(DOI:10.5281/zenodo.5829585)` and store them in the `data` directory. These networks link TFs to target genes based on open chromatin data, chromatin interaction maps, and gene expression correlations in control and AD cells. `scGRNom` pipeline was used to predict the networks.


```r
cd ADnets/
mkdir results
```
Run the following lines of code in the R console. The outputs will be stored in the `results` directory

### 1. Get gene centrality
```r
source('../scripts/get_centrality.R')
source('../scripts/get_h_metric.R')
```

### 2. Get network motifs.
This analysis is run outside of R using the `mfinder` tool. Please refer to mfinder manual on the link provided above.
Note that mfinder uses numeric gene IDs. The script `convert_symbols_to_entrezID.R` will assign numeric IDs to gene symbols and generate a `*.node_symbol-integar.key.txt` file for this mapping to be used later.
```r
Rscript ../scripts/convert_symbols_to_entrezID.R data/MIT.AD.Mic.grn.txt Mic.AD
/path/to/mfinder Mic.Ctrl.TF-TF.entrez.txt -s 3 -r 1000 -f Mic.AD -ospmem 38
mv Mic.AD_MEMBERS.txt results/
mv Mic.AD_MEMBERS.txt.node_symbol-integar.key.txt results/
rm Mic.AD.TF-TF.entrez.txt
```

### 3. Get regulatory logics
An example of how to get regulatory logics for a single cell type GRN using `Loregic`.
```r
Rscript run_loregic_FFL.R results/Mic.AD_MEMBERS.txt results/Mic.AD_MEMBERS.txt.node_symbol-integar.key.txt Mic AD 99
```

### 4. Get co-regulatory network modules
```r
source('../scripts/get_module_coregnet.R')
```

### 5. Prioritize network genes using random forest classifier
```r
source('../scripts/network_based_classifier.R')
```
