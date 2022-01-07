
library(data.table)
library(caret)
library(dplyr)
library(readxl)
library(caret)
library(PRROC)
library(e1071)
library(randomForest)
library(stringr)
library(tidyverse)
library(Loregic)
library(ArrayBin)
library(tidyr)


list.of.packages <- c("dplyr", "tidyr",
"igraph", "Loregic",
"reshape2","WGCNA",
"flashClust",)

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) BiocManager::install(new.packages)

lapply(list.of.packages, require, character.only = TRUE)
