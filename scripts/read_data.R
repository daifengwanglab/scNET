# Regulatory network format;: 3 columns (TF target score)
rm(list=ls())
home_dir=getwd()
#network data directory; change accordingly

scGRNom_data_dir=c("/Users/chiraggupta/work/scNET_manuscript/AD_MIT/data")
#scGRNom_data_dir=c("~/work/scNET_manuscript/data/Ting_Jin/AD_networks/AD01103_GRN/")


setwd(scGRNom_data_dir)
files= list.files(path=".", pattern="*.txt", full.names=TRUE)
for(i in 1:length(files))
{
  name=gsub(".grn.txt","",gsub("./MIT.","",files[i]))
  name=paste(name,"network",sep=".")
  net=read.table(files[i], header=T, sep="\t")
  net=net[net$mse<0.1 & net$abs_coef > 0.01,]
  assign(name, net)
}

setwd(home_dir)
