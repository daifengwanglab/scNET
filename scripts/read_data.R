# Regulatory network format;: 3 columns (TF target score)


#genome data directory;
genome_data="~/work/scNET/demo_data/genome"
TFs=read.table(paste(genome_data,"TF_names_v_1.01.txt",sep="/"))

home_dir=getwd()
#network data directory; change accordingly
genie_net_data="~/work/scNET_manuscript/data/Ting_Jin/GENIE3"
GRNboost2_net_data="~/work/scNET_manuscript/data/Ting_Jin/GRNboost2"
pidc_net_data="~/work/scNET_manuscript/data/Ting_Jin/PIDC"


setwd(genie_net_data)
files= list.files(path=".", pattern="*.txt", full.names=TRUE)
for(i in 1:length(files))
{
  name=gsub(".txt","",gsub("./","",files[i]))
  assign(name, read.table(files[i], header=T))
}

setwd(GRNboost2_net_data)
files= list.files(path=".", pattern="*.csv", full.names=TRUE)
for(i in 1:length(files))
{
  name=gsub(".csv","",gsub("./","",files[i]))
  assign(name, read.table(files[i], header=T, sep=","))
}

setwd(pidc_net_data)
files= list.files(path=".", pattern="*.txt", full.names=TRUE)
for(i in 1:length(files))
{
  name=gsub(".txt","",gsub("./","",files[i]))
  assign(name, read.table(files[i], header=F, sep="\t"))
}

setwd(home_dir)
