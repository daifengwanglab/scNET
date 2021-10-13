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

th=c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)
msize=c(10,20,30,40,50,60,70,80,90,100)

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


#select the th and msize parameters (here th = 0.4 and msize =10)
#select those modules
pattern="*.network.JI.coreg.modules.JI_0.4.ModSize_10"
list=ls(pattern=pattern)
list=Filter(function(x) !any(grepl("ModSize_100", x)), list)
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


p3=ggplot(module_gene_count.tbl,aes(x=celltype,y=Nmodules,fill=cond))+
geom_bar(stat="identity",position="dodge")+scale_fill_manual(values=c("AD"=npgcolors[1],"Ctrl"=npgcolors[2]))+
theme_bw(base_size=12) + labs(y="# of co-regulation modules",x="cell types")+
theme(legend.position = "top")

ggsave(p3,filename="Figures/p.no_mods_per_ct_ji0.4_msize10.pdf", device="pdf",width=4,height=4,units="in")


#calculate module denisty
Modules.df=data.frame(Modulename=NULL,cell=NULL)

#calculate edge density per module for ctrl nets
list=Filter(function(x) !any(grepl("AD", x)), list)
for(i in 1:length(list))
{
  Module.density.tbl=data.frame("Modulename"=NULL,"Ctrl"=NULL,"AD"=NULL,lfc=NULL)
  cell=gsub(".network.JI.coreg.modules","",list[i])
  cell=gsub("Ctrl.","",cell)
  df=get(list[i])
  allmodules=unique(df$moduleID)
  #remove module0
  allmodules= Filter(function(x) !any(grepl("0", x)), allmodules)
  for (j in 1:length(allmodules))
  {
      module=df[df$moduleID %in% allmodules[j],]$gene
      matname=paste("Ctrl",cell,sep=".")
      matname=paste(matname,".network.JI.coreg.mat",sep="")
      mat=get(matname)
      indx.c=match(module,colnames(mat))
      indx.r=match(module,rownames(mat))
      indx.c=indx.c[!is.na(indx.c)]
      indx.r=indx.r[!is.na(indx.r)]
      module.mat=mat[indx.r,indx.c]
      g=graph_from_adjacency_matrix(module.mat,weighted=TRUE, diag=FALSE, mode='undirected')
      df.2= get.data.frame(igraph::simplify(g,remove.multiple = TRUE, remove.loops = TRUE))
      e.density.Ctrl=sum(df.2$weight)/length(unique(module))

      matname=paste("AD",cell,sep=".")
      matname=paste(matname,".network.JI.coreg.mat",sep="")
      mat=get(matname)
      indx.c=match(module,colnames(mat))
      indx.r=match(module,rownames(mat))
      indx.c=indx.c[!is.na(indx.c)]
      indx.r=indx.r[!is.na(indx.r)]
      module.mat=mat[indx.r,indx.c]
      g=graph_from_adjacency_matrix(module.mat,weighted=TRUE, diag=FALSE, mode='undirected')
      df.2= get.data.frame(igraph::simplify(g,remove.multiple = TRUE, remove.loops = TRUE))
      e.density.AD=sum(df.2$weight)/length(unique(module))

      df.3=data.frame("Modulename"=allmodules[j],"Ctrl"=e.density.Ctrl,"AD"=e.density.AD)
      df.3$lfc = log2(df.3[,colnames(df.3)%like% "AD"]/ df.3[,colnames(df.3)%like% "Ctrl"])
      Module.density.tbl=rbind(Module.density.tbl,df.3)
  }
  colnames(Module.density.tbl)=c("Modulename","Ctrl","AD","lfc")
  Module.density.tbl$Cell=cell
  name=paste(cell,"Module.density.tbl",sep=".")
  assign(name,Module.density.tbl)
}


#box plot
for_boxplot=rbind(Mic.Module.density.tbl,Oli.Module.density.tbl)
for_boxplot=rbind(for_boxplot,Ex.Module.density.tbl)
for_boxplot=rbind(for_boxplot,In.Module.density.tbl)

npgcolors=pal_npg("nrc", alpha = 1)(10)
p.lfc_density_module_boxplot=ggplot(for_boxplot,aes(x=Cell,y=lfc)) +
  geom_boxplot()+
  labs(y="change in module edge density",x="Cell types")+
 theme_bw(base_size=12)+theme(legend.position="top")
#ggsave(p.lfc_density_module_boxplot,filename="Figures/p.lfc_density_module_boxplot.pdf", device="pdf",width=3,height=3,units="in")


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
tmp=left_join(all.deg.earlylate,all.deg.no)
tmp=left_join(tmp,all.deg.noearly)
write.table(tmp,file="all.deg.MIT.FC.mat", col.names=TRUE, row.names=FALSE,sep="\t",quote=F)


df=AAD.Mic.network.JI.coreg.modules.JI_0.4.ModSize_10
allmodules=unique(df$moduleID)
#remove module0
allmodules= Filter(function(x) !any(grepl("0", x)), allmodules)
Module.fc.tbl=data.frame("gene"=NULL,"FC"=NULL,"module"=NULL,"Pathology"=NULL)
for (j in 1:length(allmodules))
{
    module=df[df$moduleID %in% allmodules[j],]$gene
    module.de=all.deg.no[all.deg.no$gene %in% module,]
    module.de$module=allmodules[j]
    module.de$Pathology="no-pathology vs pathology"
    Module.fc.tbl=rbind(Module.fc.tbl,module.de)

    module.de=all.deg.noearly[all.deg.noearly$gene %in% module,]
    module.de$module=allmodules[j]
    module.de$Pathology="no-pathology vs early-pathology"
    Module.fc.tbl=rbind(Module.fc.tbl,module.de)

    module.de=all.deg.earlylate[all.deg.earlylate$gene %in% module,]
    module.de$module=allmodules[j]
    module.de$Pathology=" early-pathology vs late-pathology"
    Module.fc.tbl=rbind(Module.fc.tbl,module.de)
}

p.module_fc_boxplot=ggplot(Module.fc.tbl,aes(x=as.factor(module),y=FC)) + geom_boxplot()+
  labs(y="change in GO edge density",x="Cell types")+facet_wrap(~Pathology, nrow=3)+
 theme_bw(base_size=12)+theme(legend.position="top")
#ggsave(p.lfc_density_GOBP_boxplot,filename="Figures/p.lfc_density_GOBP_boxplot.pdf", device="pdf",width=3,height=3,units="in")
