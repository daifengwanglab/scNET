# Regulatory network format;: 3 columns (TF target score)
rm(list=ls())
home_dir=getwd()

#network data directory; change accordingly

data_dir=c("./data/")

setwd(data_dir)
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
