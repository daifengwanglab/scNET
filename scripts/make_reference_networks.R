rm(list=ls())
source('/ua/cgupta8/tools/scGRNom/code/scGRNom_getTF.R')
source('/ua/cgupta8/tools/scGRNom/code/scGRNom_getNt.R')
source('/ua/cgupta8/tools/scGRNom/code/scGRNom_interaction.R')

#load('../data/gematrix/MIT_imputed_gexpr.rdata')

library(readxl)
library(data.table)

cells=c("Microglia","Neuronal","Oligo")

for (i in 1:length(cells))
{
  print (cells[i])
  sheet1=paste(cells[i],"interactome", sep=" ")
  sheet2=paste(cells[i],"enhancers", sep=" ")
  interactome_data = read_xlsx("/ua/cgupta8/tools/scGRNom/data/PLAC-seq promoter interactome map.xlsx",sheet = sheet1, skip = 2)[,1:6]
  enhancers = read_xlsx("/ua/cgupta8/tools/scGRNom/data/PLAC-seq promoter interactome map.xlsx",sheet = sheet2,skip = 2,col_names = c('chr','start','end'))
  df1 = scGRNom_interaction(interactome_data,enhancers)
  df2 = scGRNom_getTF(df1,num_cores = 24)
  name=paste("df2",cells[i],sep=".")
  assign(name, df2)
}

save.image("BrainCelltypeReferenceNetworks.RData")
