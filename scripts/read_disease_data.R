# Regulatory network format;: 3 columns (TF target score)
home_dir=getwd()
#network data directory; change accordingly
scGRNom_data_dir=c("~/work/scNET_manuscript/data/Ting_Jin/AD_networks/AD01103_GRN/")


setwd(scGRNom_data_dir)
files= list.files(path=".", pattern="*.csv", full.names=TRUE)
for(i in 1:length(files))
{
  name=gsub(".csv","",gsub("./","",files[i]))
  name=paste(name,"network",sep=".")
  assign(name, read.table(files[i], header=T, sep=","))
}

setwd(home_dir)
