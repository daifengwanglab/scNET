source('~/work/scNET-devel/scripts/read_data.R')
source('~/work/scNET-devel/scripts/functions_for_network_analysis.R')
source('~/work/scNET-devel/scripts/load_libraries.R')

nets=ls(pattern="*\\.network")

filter nets to select only edges with BS in promoter
for(i in 1:length(nets))
{
  name=nets[i]
  tag=gsub(".network","",name)
  name=paste(tag,"filtered",sep=".")
  net=as.data.frame(lapply(nets[i],get))

  TF.net1=net[net$TFbs %in% "promoter", ]
  TF.net2=net[net$TFbs %in% "both", ]
  TF.net=rbind(TF.net1,TF.net2)
  TF.net=net[,c("TF","TG","mse")]

  TF.net=distinct(TF.net)
  assign(name, TF.net)

}


load('filtered_nets.RData')
source('HirNet_function.R')
data=AD.Mic.filtered[,1:2]
cal_hier_score26(data=data, kmax=10000, ptim=100, anneal.coeff=1e-6,  myoutf = "AD_MIC.result.txt")
