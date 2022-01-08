

list.of.packages <- c("dplyr", "tidyr","stringr",
"igraph", "Loregic","randomForest",
"reshape2","WGCNA","PPROC",
"flashClust","data.table","ArrayBin","tidyverse")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) BiocManager::install(new.packages)

lapply(list.of.packages, require, character.only = TRUE)
