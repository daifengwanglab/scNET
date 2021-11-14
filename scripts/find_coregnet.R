rm(list=ls())
source('~/work/scNET-devel/scripts/load_libraries.R')
source('~/work/scNET-devel/scripts/read_data.R')
source('~/work/scNET-devel/scripts/functions_for_network_analysis.R')

npgcolors=pal_npg("nrc", alpha = 1)(10)

#fig 2a in https://www.nature.com/articles/s12276-020-00528-0#Sec2
#fig 1c in https://www.nature.com/articles/s41467-019-10591-5/figures/1

#Read GO data
data=GSA.read.gmt('/Users/chiraggupta/work/scNET_manuscript/genome/genesets/BaderLab/Human_GO_bp_with_GO_iea_symbol.gmt')
genesets=data$genesets
names(genesets)=data$geneset.descriptions

#number of random nets needed
nrandnets=1
set.seed(123)

nets=ls(pattern="*\\.network")

#th=c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)
#msize=c(10,20,30,40,50,60,70,80,90,100)


tm=c(0.2)
msize=c(30)

for(i in 1:length(nets))
{
  for (j in 1:length(th))
  {
    for (k in 1:length(msize))
    {
      name=nets[i]
      celltype=gsub(".network","",name)
      net=as.data.frame(lapply(nets[i],get))
      net=net[,c("TF","TG","abs_coef")]
      net=distinct(net)
      name=paste(name,"JI.coreg.mat",sep=".")
      assign(name, find_target_pairs_matrix(net))
      mat=get(name)
      mat[mat < th[j]] <-0 #remove edges with less than th% overlap
      name=gsub(".mat",".modules", name)
      tag=paste("JI",th[j],sep="_")
      name=paste(name,tag,sep=".")
      tag=paste("ModSize",msize[k],sep="_")
      name=paste(name,tag,sep=".")
      assign(name, detect_modules(mat,msize[k]))
    }
  }
}


#calculate number of genes and modules in each network
list=ls(pattern="*\\.modules")
module_gene_count.tbl=data.frame(celltype=NULL,Nmodules=NULL,Ngenes=NULL,cond=NULL,th=NULL,msize=NULL)
for (i in 1:length(list))
{
  df=get(list[i])
  df=df[as.character(df$moduleID) != "0",]
  name=list[i]
  name=gsub(".network.JI.coreg.modules.","",name)
  name=gsub("JI_"," ",gsub(".ModSize_"," ",name))
  th=sapply(strsplit(name, " "), head, 3)[2]
  msize=sapply(strsplit(name, " "), head, 3)[3]
  name=sapply(strsplit(name, " "), head, 3)[1]
  ct=ifelse(name %like% "Ex","Ex",ifelse(name %like% "In","In",ifelse(name %like% "Mic","Mic","Oli")))
  cond=ifelse(name %like% "AD","AD","Ctrl")
  tmp=data.frame(celltype=ct,Nmodules=length(unique(df$moduleID)),Ngenes=length(unique(df$gene)),cond=cond,th=th,msize=msize)
  module_gene_count.tbl=rbind(module_gene_count.tbl,tmp)
}


tmp=ggplot(module_gene_count.tbl,aes(x=celltype,y=Nmodules,fill=cond))+
geom_bar(stat="identity",position="dodge")+facet_grid(th~msize)
scale_fill_manual(values=c("AD"=npgcolors[1],"Ctrl"=npgcolors[2]))+
theme_bw(base_size=12) + labs(y="# of co-regulation modules",x="cell types")+
theme(legend.position = "top")



p1=ggplot(module_gene_count.tbl,aes(x=celltype,y=Nmodules,fill=cond))+
geom_bar(stat="identity",position="dodge")+facet_grid(th~as.numeric(msize))+scale_fill_manual(values=c("AD"=npgcolors[1],"Ctrl"=npgcolors[2]))+
theme_bw(base_size=12) + labs(y="# of co-regulation modules",x="cell types")+
theme(legend.position = "top")+theme(axis.text.x=element_text(angle=90))

 p2=ggplot(module_gene_count.tbl,aes(x=celltype,y=Ngenes))+
 geom_point(aes(colour=factor(cond),
      fill = factor(cond)), shape=21, size = 0.4)+
 scale_fill_manual(values=c("AD"=npgcolors[1],"Ctrl"=npgcolors[2]))+
 facet_grid(th~as.numeric(msize)) +
 theme_bw(base_size=12) + labs(y="# genes",x="cell types")+
 theme(legend.position = "top")+theme(axis.text.x=element_text(angle=90))

ggsave(p1,filename="Figures/p.no_mods_per_ct.pdf", device="pdf",width=6,height=6,units="in")
ggsave(p2,filename="Figures/p.no_genes_per_mod.pdf", device="pdf",width=6,height=6,units="in")


#select the th and msize parameters (here th = 0.2 and msize =30)
#select those modules
pattern="*.network.JI.coreg.modules.JI_0.2.ModSize_30"
list=ls(pattern=pattern)
module_gene_count.tbl=data.frame(celltype=NULL,Nmodules=NULL,Ngenes=NULL,cond=NULL,th=NULL,msize=NULL)
for (i in 1:length(list))
{
  df=get(list[i])
  df=df[as.character(df$moduleID) != "0",]
  name=list[i]
  name=gsub(".network.JI.coreg.modules.","",name)
  name=gsub("JI_"," ",gsub(".ModSize_"," ",name))
  th=sapply(strsplit(name, " "), head, 3)[2]
  msize=sapply(strsplit(name, " "), head, 3)[3]
  name=sapply(strsplit(name, " "), head, 3)[1]
  ct=ifelse(name %like% "Ex","Ex",ifelse(name %like% "In","In",ifelse(name %like% "Mic","Mic","Oli")))
  cond=ifelse(name %like% "AD","AD","Ctrl")
  tmp=data.frame(celltype=ct,Nmodules=length(unique(df$moduleID)),Ngenes=length(unique(df$gene)),cond=cond,th=th,msize=msize)
  module_gene_count.tbl=rbind(module_gene_count.tbl,tmp)
}

#check the end of this script to find where GO table came from
tmp=GO.diff.cent.enrich.tbl[ GO.diff.cent.enrich.tbl$th %in% "0.2" & GO.diff.cent.enrich.tbl$msize %in% "30" & GO.diff.cent.enrich.tbl$geneset < 500 & GO.diff.cent.enrich.tbl$geneset > 10 ,]
tmp=tmp[,c("label","cell")]
tmp.df=as.data.frame(table(tmp$cell))

module_gene_count.tbl$nGOBP=tmp.df$Freq
p3=ggplot(module_gene_count.tbl,aes(x=celltype,y=Nmodules,fill=cond))+
geom_bar(stat="identity",position="dodge")+scale_fill_manual(values=c("AD"=npgcolors[1],"Ctrl"=npgcolors[2]), name="")+
theme_bw(base_size=12) + labs(y="# of \n co-regulation modules",x=" ")+
theme(legend.position = "top")

p4=ggplot(module_gene_count.tbl,aes(x=celltype,y=nGOBP,fill=cond))+
geom_bar(stat="identity",position="dodge")+scale_fill_manual(values=c("AD"=npgcolors[1],"Ctrl"=npgcolors[2]), name="")+
theme_bw(base_size=12) + labs(y="# of GOBPs enriched in modules",x=" ")+
theme(legend.position = "top")

ggsave(p3,filename="../Figures/p.no_mods_per_ct_ji0.2_msize30.pdf", device="pdf",width=2,height=2,units="in")
ggsave(p4,filename="../Figures/p.no_GOBP_in_mods_per_ct_ji0.2_msize30.pdf", device="pdf",width=2,height=2,units="in")



#mic AD module expression fold change

all.deg.no=read_xlsx("~/work/scNET_manuscript/data/gematrix/Diff.Exp.Genes.DataS2.MIT.xlsx",sheet="Mic",skip=1)[,1:9]
all.deg.no=as.data.frame(all.deg.no[,c(1,5)])
colnames(all.deg.no)=c("gene","nopath")

all.deg.noearly=read_xlsx("~/work/scNET_manuscript/data/gematrix/Diff.Exp.Genes.DataS2.MIT.xlsx",sheet="Mic",skip=1)[,12:20]
all.deg.noearly=as.data.frame(all.deg.noearly[,c(1,5)])
colnames(all.deg.noearly)=c("gene","noveraly")

all.deg.earlylate=read_xlsx("~/work/scNET_manuscript/data/gematrix/Diff.Exp.Genes.DataS2.MIT.xlsx",sheet="Mic",skip=1)[,23:31]
all.deg.earlylate=as.data.frame(all.deg.earlylate[,c(1,5)])
colnames(all.deg.earlylate)=c("gene","earlyvlate")

#for enrichment analysis
tmp=left_join(all.deg.no,all.deg.noearly)
tmp=left_join(tmp,all.deg.earlylate)

write.table(tmp,file="all.deg.MIT.FC.mat", col.names=TRUE, row.names=FALSE,sep="\t",quote=F)

df=AD.Mic.network.JI.coreg.modules.JI_0.2.ModSize_30
write.table(df[as.character(df$moduleID) != "0",],file="AD.Mic.network.JI.coreg.modules.JI_0.4.ModSize_10.genesets", col.names=FALSE, row.names=FALSE,sep="\t",quote=F)

#read enrichment results
npgcolors=pal_npg("nrc", alpha = 1)(10)

mat=read.table("all.deg.MIT.FC.mat-AD.Mic.network.JI.coreg.modules.JI_0.2.ModSize_30.fdr-based.zscore.mat", header=T, row.names=1)
mat=mat[,-c(1:2)]
rownames(mat)=paste("M",rownames(mat),sep="")
colnames(mat)=c("no pathology","no vs. early","early vs. late")

paletteLength=100
#myColor = colorRampPalette(c(npgcolors[8], "white", npgcolors[4]))(paletteLength)
myColor = colorRampPalette(c("#B22222", "white", "#0A24CC"))(paletteLength)
myBreaks <- c(seq(-2, 0, length.out=ceiling(paletteLength/2) + 1),
+ seq(max(as.matrix(mat))/paletteLength, 2, length.out=floor(paletteLength/2)))

pdf(file="Figures/p.module_FC_AD_Path.pdf")
pheatmap(t(as.matrix(mat)),cellwidth=12,cellheight=8,show_rownames=TRUE,
 breaks=myBreaks,
 col=myColor,
 cluster_cols=F,
 cluster_rows=F,
 border_color = "black",
 angle_col=c("45")
 )
dev.off()


#select module 6 from heatmap and extract GO
#run  the script DO_module_enrichment.R to get "GO.diff.cent.enrich.tbl" object
GO.diff.cent.enrich.filt.tbl=GO.diff.cent.enrich.tbl[GO.diff.cent.enrich.tbl$geneset < 500 & GO.diff.cent.enrich.tbl$geneset > 10,]
write.table(GO.diff.cent.enrich.filt.tbl, file="GO.diff.cent.enrich.filt.tbl", sep="\t",quote=F, row.names=F, col.names=T)



#module 2 GO BP barplot
 mat=read.table("../data/mic.module2.GOBP.txt", sep="\t", header=T)
mat$logp=-1*log10(mat$pval)
 p=ggplot(mat,aes(x=reorder(label, logp),y=logp))+
 geom_bar( stat="identity", fill="grey47")+coord_flip()+
 theme_bw(base_size=12) + labs(y="-1*log(pvalue)",x="Biological process")

 ggsave(p,filename="Figures/p.mod2_GOBP_bar.pdf", device="pdf",width=2,height=2,units="in")
