
# needs brush up

library(pheatmap)
library(RColorBrewer)
library(viridis)
library(clusterProfiler)
library(DOSE)
#library(ComplexHeatmap)
library("data.table")
library(ggplot2)
library(hypeR)
library(GSA)
library(mGSZ)
library(ggrepel)


list.of.packages <- c("dplyr", "tidyr",
"igraph", "Loregic",
"reshape2","WGCNA",
"flashClust","arsenal",
"piano","topGO","org.Hs.eg.db")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) BiocManager::install(new.packages)

lapply(list.of.packages, require, character.only = TRUE)
