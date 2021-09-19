
# mfinder _OUT files with extra text removed

rm(list=ls())
source('~/work/scNET-devel/scripts/load_libraries.R')
source('~/work/scNET-devel/scripts/functions_for_network_analysis.R')
home_dir=getwd()


data_dir=c("/Users/chiraggupta/work/scNET_manuscript/data/motifs/1000_rands_with_100_sampling")

#data_dir=c("/Users/chiraggupta/work/scNET_manuscript/data/motifs/fullnets")

setwd(data_dir)


motif_count.tbl=data.frame(mID=NULL,nREAL=NULL,nRAND=NULL,zscore=NULL,pval=NULL,cREAL=NULL,uniq=NULL,ct=NULL,cond=NULL)

files= list.files(path=".", pattern="*.txt.for_plot", full.names=TRUE)
for(i in 1:length(files))
{
  name=gsub("_OUT.txt.for_plot","",files[i])
  name=gsub("\\./","",name)
  name=paste(name,"motif_counts",sep=".")
  df=read.table(files[i])
  colnames(df)=c("mID","nREAL","nRAND","zscore","pval","cREAL","uniq")
  df$ct=ifelse(name %like% "Ex","Ex",ifelse(name %like% "In","In",ifelse(name %like% "Mic","Mic","Oli")))
  df$cond=ifelse(name %like% "AD","AD","Ctrl")
  motif_count.tbl=rbind(motif_count.tbl,df)
}
setwd(home_dir)

#filter rubbish zscores
motif_count.tbl$zscore=ifelse(motif_count.tbl$zscore == 888888.00,0,motif_count.tbl$zscore)

p=ggplot(motif_count.tbl,aes(x=as.character(mID),y=zscore,fill=cond))+geom_bar(stat="identity",position="dodge")+facet_wrap(~ct,ncol=1)

ggsave(p,filename="Figures/p.motifs.fullnets.pdf", device="pdf",width=4,height=6,units="in")


#for full nets
ggsave(p,filename="Figures/p.motifs.fullnets.pdf", device="pdf",width=4,height=6,units="in")
