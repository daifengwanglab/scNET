# Single-cell network biology characterizes cell type gene regulation for drug repurposing and phenotype prediction in Alzheimer’s disease

## Summary
Dysregulation of gene expression in Alzheimer’s disease (AD) remains elusive, especially at the cell type level. Gene regulatory network, a key molecular mechanism linking transcription factors (TFs) and regulatory elements to govern target gene expression, can change across cell types in the human brain and thus serve as a model for studying gene dysregulation in AD. However, it is still challenging to understand how cell type networks work abnormally under AD. To address this, we integrated single-cell multi-omics data and predicted the gene regulatory networks in AD and control for four major cell types, excitatory and inhibitory neurons, microglia and oligodendrocytes. Importantly, we applied network biology approaches to analyze the changes of network characteristics across these cell types, and between AD and control. For instance, many hub TFs target different genes between AD and control (rewiring). Also, these networks show strong hierarchical structures in which top TFs (master regulators) are largely common across cell types, whereas different TFs operate at the middle levels in some cell types (e.g., microglia). The regulatory logics of enriched network motifs (e.g., feed-forward loops) further uncover cell type-specific TF-TF cooperativities in gene regulation. The cell type networks are highly modular and several network modules with cell-type-specific expression changes in AD pathology are enriched with AD-risk genes and putative targets of approved and pending AD drugs, suggesting possible cell-type genomic medicine in AD. Finally, using the cell type gene regulatory networks, we developed machine learning models to classify and prioritize additional AD genes. We found that top prioritized genes predict clinical phenotypes (e.g., cognitive impairment) with reasonable accuracy. Overall, this single-cell network biology analysis provides a comprehensive map linking genes, regulatory networks, cell types and drug targets and reveals dysregulated cell type gene dysregulatory mechanisms in AD.

## Flow chart
![alt text](https://github.com/daifengwanglab/scNET/workflow.png?raw=true)

## System requirements

The analysis is based on R version 4.0.3. For the gene regulatory network, a *Linux* system with 32 GB RAM and 32GB storage would be enough. For network analysis, a standard computer should be enough.

## Software installation guide

Packages needed for the whole analysis.

- `igraph` (R package, https://igraph.org/r/)
- `mfinder` (motif finding tool from https://www.weizmann.ac.il/mcb/UriAlon/download/network-motif-software)
- `Loregic` (R package for regulatory logics)


## Download code
The code has been tested on R version 4.0.3 on Mac and Linux OS.
```r
git clone https://github.com/daifengwanglab/scNET
cd scNET
mkdir results
```

## Demo for analysis of brain cell type network characteristics

A demo for aligning single-cell multi-modal data is available at `demo_brain_ct/`

```r
cd demo_brain_ct/
```
Store all networks with .txt extension in the 'demo_data' directory. The demo contains outputs from scGRNom ().

Run the following lines of code in

Get gene centrality
```r
source('../scripts/get_centrality.R')
source('../scripts/get_h_metric.R')
```

Get network motifs. This analysis is run outside of R using the mfinder tool. Please refer to mfinder manual on the link provided above.
```r
Rscript ../scripts/convert_symbols_to_entrezID.R demo_data/MIT.AD.Ex.grn.demo.txt Mic.AD
/path/to/mfinder Mic.Ctrl.TF-TF.entrez.txt -s 3 -r 1000 -f Mic.AD -ospmem 38
mv Mic.AD_MEMBERS.txt results/
mv Mic.AD.node_symbol-integar.key.txt results/
```

Get co-regulatory network modules
```r
source('../scripts/find_module_coregnet.R')
```

Prioritize network genes using random forest classifier
```r
source('../scripts/network_based_classifier.R')
```
