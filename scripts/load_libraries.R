
# needs brush up


list.of.packages <- c("dplyr", "tidyr",
"igraph", "Loregic",
"reshape2","WGCNA",
"flashClust","arsenal",
"piano","topGO","org.Hs.eg.db")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) BiocManager::install(new.packages)

lapply(list.of.packages, require, character.only = TRUE)
