
rm(list=ls())
source('~/work/scNET-devel/scripts/load_libraries.R')
source('~/work/scNET-devel/scripts/read_data.R')
source('~/work/scNET-devel/scripts/functions_for_network_analysis.R')


#feature matrix from AD-gene classifier
fw=read.table("Mic.AD.Feature_weights.mat", header=T)
fw$average=abs(fw$average)
fw=fw[order(fw$average),]

positives=rownames(tail(fw,dim(fw)[1]*.1))


net=AD.Mic.network
tag="AD.Mic"
TF.net=net[,c("TF","TG","mse")]
TF.net=distinct(TF.net)
TF.net=TF.net[TF.net$TG %in% TF.net$TF,]
TF.net.igraph=graph_from_data_frame(TF.net, directed = TRUE, vertices = NULL)
out.deg=get_centrality(TF.net,"degree_out",tag)
in.deg=get_centrality(TF.net,"degree_in",tag)
df=merge(out.deg,in.deg,by="gene",all=TRUE)
colnames(df)=c("gene","out.d","in.d")
#df$o_minus_i=df$out.d-df$in.d
#df$o_plus_i=df$out.d+df$in.d
df$h=(df$out.d-df$in.d)/(df$out.d+df$in.d)
df=df[,c("gene","h")]



fw$gene=rownames(fw)
fw=fw[,c("gene","average")]
rownames(fw)=NULL
fw.df=left_join(fw,df)
fw.df=na.omit(fw.df)


fw.df$levels = cut(fw.df[,3],c(-Inf,-0.33,0.33,1),labels=c("Bottom","Middle","Top"))


p1=ggplot(fw.df, aes(x = factor(levels), y = average)) +
  geom_bar(stat = "summary", fun = "mean")+
  labs(y="Feature importance score in \n AD gene classifier ",x="Hierarchy level")+
  theme_bw(base_size=12)+theme(legend.position="top")

  ggsave(p1,filename="Figures/p.feat_imp.hierarchy.barplot.pdf", device="pdf",width=2,height=3,units="in")



positves.net=TF.net[TF.net$TF %in% positives & TF.net$TG %in% positives,]
write.table(positves.net,file="AD.mic.imp.feat.top10perc.dat",sep="\t",col.names=T,row.names=F,quote=F)


write.table(fw.df,file="D.mic.imp.feat.attributes.txt",sep="\t",col.names=T,row.names=F,quote=F)
