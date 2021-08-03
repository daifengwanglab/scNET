# Regulatory network format;: 3 columns (TF target score)
rm(list=ls())
home_dir=getwd()
#network data directory; change accordingly

#scGRNom_data_dir=c("/Users/chiraggupta/work/scNET_manuscript/AD_MIT/data/LakeCtrl")

scGRNom_data_dir=c("/Users/chiraggupta/work/scNET_manuscript/AD_MIT/data")

##genes=c("LINGO1", "CNTNAP2","ERBB2IP","NEGR1","BEX1","NTNG1","SLC17A7")
setwd(scGRNom_data_dir)
files= list.files(path=".", pattern="*.txt", full.names=TRUE)
for(i in 1:length(files))
{
  name=gsub(".grn.txt","",gsub("./MIT.","",files[i]))
  name=paste(name,"network",sep=".")
  net=read.table(files[i], header=T, sep="\t")
  net=net[net$mse<0.1 & net$abs_coef > 0.01,]
  print(paste(name,nrow(net),sep=":"))
  assign(name, net)
}

setwd(home_dir)
