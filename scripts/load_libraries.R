
# needs brush up




library(pheatmap)
library(RColorBrewer)
library(viridis)
library(clusterProfiler)
library(DOSE)
library(ComplexHeatmap)
library(ggplot2)
library(hypeR)
library(GSA)
library(mGSZ)
library(ggrepel)
library(ggsci)
library(data.table)
library(caret)
library(easyalluvial)
library('alluvial')
library(ggsankey)
library(dplyr)
library(wordcloud)
library(disgenet2r)
library(ggplot2)
library(ComplexUpset)
library(readxl)

npgcolors=pal_npg("nrc", alpha = 1)(10)
troncolors=pal_tron("legacy", alpha = 1)(7)

list.of.packages <- c("dplyr", "tidyr",
"igraph", "Loregic",
"reshape2","WGCNA",
"flashClust","arsenal",
"piano","topGO","org.Hs.eg.db")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) BiocManager::install(new.packages)

lapply(list.of.packages, require, character.only = TRUE)
