rm(list=ls())
#heirarchy

source('~/work/scNET-devel/scripts/load_libraries.R')

celltypes=c("MIC","Oli","Ex","In")

##https://stackoverflow.com/questions/51991825/r-alluvial-plots-for-combinations

#https://github.com/davidsjoberg/ggsankey

mic.ad=read.table('heirarchy/AD_MIC.result.Lev6.txt', header=T)
colnames(mic.ad)=gsub("X","Lev",colnames(mic.ad))
mic.ad$Gene=rownames(mic.ad)
rownames(mic.ad)=NULL
mic.ad=as.data.table(mic.ad)
mic.ad[,maximum_column :=  names(.SD)[max.col(.SD)], .SDcols = 1:6]
mic.ad.gene.lev=mic.ad[,c("Gene","maximum_column")]
colnames(mic.ad.gene.lev)=c("Gene","Mic.AD")
mic.ctrl=read.table('heirarchy/Ctrl_MIC.result.Lev6.txt', header=T)
colnames(mic.ctrl)=gsub("X","Lev",colnames(mic.ctrl))
mic.ctrl$Gene=rownames(mic.ctrl)
rownames(mic.ctrl)=NULL
mic.ctrl=as.data.table(mic.ctrl)
mic.ctrl[,maximum_column :=  names(.SD)[max.col(.SD)], .SDcols = 1:6]
mic.ctrl.gene.lev=mic.ctrl[,c("Gene","maximum_column")]
colnames(mic.ctrl.gene.lev)=c("Gene","Mic.Ctrl")
Mic=merge(mic.ctrl.gene.lev, mic.ad.gene.lev)
df = Mic %>% make_long(colnames(Mic[,2:3]))
p.mic=ggplot(df, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = factor(node), label = node)) +
  geom_sankey(flow.alpha = .6,
              node.color = "gray30") +
  scale_fill_npg() +
  theme_sankey(base_size = 22) +
  labs(x = NULL) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = .5))

pdf(file="Figures/mic.sankey.pdf")
p.mic
dev.off()

Oli.ad=read.table('heirarchy/AD_Oli.result.Lev6.txt', header=T)
colnames(Oli.ad)=gsub("X","Lev",colnames(Oli.ad))
Oli.ad$Gene=rownames(Oli.ad)
rownames(Oli.ad)=NULL
Oli.ad=as.data.table(Oli.ad)
Oli.ad[,maximum_column :=  names(.SD)[max.col(.SD)], .SDcols = 1:6]
Oli.ad.gene.lev=Oli.ad[,c("Gene","maximum_column")]
colnames(Oli.ad.gene.lev)=c("Gene","Oli.AD")
Oli.ctrl=read.table('heirarchy/Ctrl_Oli.result.Lev6.txt', header=T)
colnames(Oli.ctrl)=gsub("X","Lev",colnames(Oli.ctrl))
Oli.ctrl$Gene=rownames(Oli.ctrl)
rownames(Oli.ctrl)=NULL
Oli.ctrl=as.data.table(Oli.ctrl)
Oli.ctrl[,maximum_column :=  names(.SD)[max.col(.SD)], .SDcols = 1:6]
Oli.ctrl.gene.lev=Oli.ctrl[,c("Gene","maximum_column")]
colnames(Oli.ctrl.gene.lev)=c("Gene","Oli.Ctrl")
Oli=merge(Oli.ctrl.gene.lev, Oli.ad.gene.lev)
df = Oli %>% make_long(colnames(Oli[,2:3]))
p.Oli=ggplot(df, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = factor(node), label = node)) +
  geom_sankey(flow.alpha = .6,
              node.color = "gray30") +
  scale_fill_npg() +
  theme_sankey(base_size = 22) +
  labs(x = NULL) +
  theme(legend.position = "none", plot.title = element_text(hjust = .5),)

  pdf(file="Figures/oli.sankey.pdf")
  p.Oli
  dev.off()

Ex.ad=read.table('heirarchy/AD_Ex.result.Lev6.txt', header=T)
colnames(Ex.ad)=gsub("X","Lev",colnames(Ex.ad))
Ex.ad$Gene=rownames(Ex.ad)
rownames(Ex.ad)=NULL
Ex.ad=as.data.table(Ex.ad)
Ex.ad[,maximum_column :=  names(.SD)[max.col(.SD)], .SDcols = 1:6]
Ex.ad.gene.lev=Ex.ad[,c("Gene","maximum_column")]
colnames(Ex.ad.gene.lev)=c("Gene","Ex.AD")
Ex.ctrl=read.table('heirarchy/Ctrl_Ex.result.Lev6.txt', header=T)
colnames(Ex.ctrl)=gsub("X","Lev",colnames(Ex.ctrl))
Ex.ctrl$Gene=rownames(Ex.ctrl)
rownames(Ex.ctrl)=NULL
Ex.ctrl=as.data.table(Ex.ctrl)
Ex.ctrl[,maximum_column :=  names(.SD)[max.col(.SD)], .SDcols = 1:6]
Ex.ctrl.gene.lev=Ex.ctrl[,c("Gene","maximum_column")]
colnames(Ex.ctrl.gene.lev)=c("Gene","Ex.Ctrl")
Ex=merge(Ex.ctrl.gene.lev, Ex.ad.gene.lev)
df = Ex %>% make_long(colnames(Ex[,2:3]))
p.Ex=ggplot(df, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = factor(node), label = node)) +
  geom_sankey(flow.alpha = .6,
              node.color = "gray30") +
    scale_fill_npg() +
  theme_sankey(base_size = 22) +
  labs(x = NULL) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = .5))

        pdf(file="Figures/Ex.sankey.pdf")
        p.Ex
        dev.off()

In.ad=read.table('heirarchy/AD_In.result.Lev6.txt', header=T)
colnames(In.ad)=gsub("X","Lev",colnames(In.ad))
In.ad$Gene=rownames(In.ad)
rownames(In.ad)=NULL
In.ad=as.data.table(In.ad)
In.ad[,maximum_column :=  names(.SD)[max.col(.SD)], .SDcols = 1:6]
In.ad.gene.lev=In.ad[,c("Gene","maximum_column")]
colnames(In.ad.gene.lev)=c("Gene","In.AD")
In.ctrl=read.table('heirarchy/Ctrl_In.result.Lev6.txt', header=T)
colnames(In.ctrl)=gsub("X","Lev",colnames(In.ctrl))
In.ctrl$Gene=rownames(In.ctrl)
rownames(In.ctrl)=NULL
In.ctrl=as.data.table(In.ctrl)
In.ctrl[,maximum_column :=  names(.SD)[max.col(.SD)], .SDcols = 1:6]
In.ctrl.gene.lev=In.ctrl[,c("Gene","maximum_column")]
colnames(In.ctrl.gene.lev)=c("Gene","In.Ctrl")
In=merge(In.ctrl.gene.lev, In.ad.gene.lev)
df = In %>% make_long(colnames(In[,2:3]))
p.In=ggplot(df, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = factor(node), label = node)) +
  geom_sankey(flow.alpha = .6,
              node.color = "gray30") +
    scale_fill_npg() +
  theme_sankey(base_size = 18) +
  labs(x = NULL) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = .5))

        pdf(file="Figures/In.sankey.pdf")
        p.In
        dev.off()

#combine all and find enrichment

list=ls(pattern=".gene.lev")
data = Reduce(function(x, y) merge(x, y, all=T), lapply(list,get), accumulate=F)

data[is.na(data)]="None"
data[data == "Lev1"] = "Bottom"
data[data == "Lev2" | data == "Lev3" | data == "Lev4" | data == "Lev5" ] = "Middle"
data[data == "Lev6"] = "Top"

write.table(data, file="~/work/scNET_manuscript/AD_MIT/heirarchy/all_genes_heirarchy_levs.txt",sep="\t",row.names=F, col.names=T, quote=F)

#Level enriehment


#GO enrichment analysis
GOdata=GSA.read.gmt('~/work/scNET_manuscript/genome/genesets/GO_annotations-9606-inferred-allev.gmt')
genesets=GOdata$genesets
names(genesets)=GOdata$geneset.descriptions

levels=c("Top","Bottom","Middle")
celltypes=c("Mic","Oli","Ex","In")

df=as.data.frame(data)
rownames(df)=df$Gene
df$Gene=NULL


heir.GO.enrich.tbl=data.frame("label"=NULL,"pval"=NULL,"fdr"=NULL,"signature"=NULL,"geneset"=NULL,
"overlap"=NULL,"background"=NULL,"hits"=NULL,"cell"=NULL,"Level"=NULL)

for(i in 1:length(celltypes))
{
    for (j in 1:length(levels))
    {
      cell.df=df[,colnames(df) %like% celltypes[i]]
      colnames(cell.df)=gsub(celltypes[i],"",colnames(cell.df))
      colnames(cell.df)=gsub("\\.","",colnames(cell.df))
      cell.df=cell.df %>% dplyr::filter(AD==levels[j] & Ctrl == levels[j])
      signature=rownames(cell.df)
      hyp_obj = hypeR(signature, genesets ,fdr=0.001)$data
      if(nrow(hyp_obj) > 0)
      {
        hyp_obj$cell=celltypes[i]
        hyp_obj$Level=levels[j]
        hyp_obj=hyp_obj[1:5,]
      }
      heir.GO.enrich.tbl=rbind(heir.GO.enrich.tbl,hyp_obj)
    }
}

tmp=heir.GO.enrich.tbl[,c("label","fdr","cell","Level")]
tmp=tmp[!is.na(tmp$fdr),]
tmp$logFDR=-log10(tmp$fdr)
tmp$label=gsub("_"," ",tmp$label)
p.heir.lvlGO=ggplot(tmp, aes(x=cell, y=label))+geom_point(aes(size=logFDR))+
facet_wrap(~Level, nrow=1) +
theme_bw(base_size=12)+theme(text = element_text(size = 12),axis.text.x=element_text(angle=90))+
labs(x="Levels",y="Biological process")+scale_y_discrete(label = function(x) stringr::str_trunc(x, 50))

pdf(file="Figures/heir.lev.GOBP.fullnet.pdf")
p.heir.lvlGO
dev.off()


#genesets <- msigdb_gsets("Homo sapiens", "C2", "CP:KEGG", clean=TRUE)

#difference enrichment
Ex=data[data$Ex.AD!=data$Ex.Ctrl,]$Gene
In=data[data$In.AD!=data$In.Ctrl,]$Gene
Mic=data[data$Mic.AD!=data$Mic.Ctrl,]$Gene
Oli=data[data$Oli.AD!=data$Oli.Ctrl,]$Gene

Ex.hyp = hypeR(Ex, genesets,fdr=0.05)$data
Ex.hyp$cell="Ex"
Mic.hyp = hypeR(Mic, genesets,fdr=0.05)$data
Mic.hyp$cell="Mic"
Oli.hyp = hypeR(Oli, genesets,fdr=0.05)$data
Oli.hyp$cell="Oli"
In.hyp = hypeR(In, genesets,fdr=0.05)$data
In.hyp$cell="In"
tbl=rbind(Ex.hyp,Mic.hyp,Oli.hyp,In.hyp)
p.heir.GOenrich=ggplot(tbl, aes(y=label, x=cell)) + geom_point(aes(size=-log10(fdr))) +
theme_minimal()


source('~/work/scNET_manuscript/get_api.R')
#disease ontology enrichment
#res_enrich = disease_enrichment( entities =Ex, vocabulary = "HGNC", database = "ALL")
#Ex.table = res_enrich@qresult[res_enrich@qresult$FDR<0.01, c("Description", "FDR", "Ratio",  "BgRatio")]
#Ex.table$cell="Ex"

#res_enrich = disease_enrichment( entities =In, vocabulary = "HGNC", database = "ALL")
#In.table = res_enrich@qresult[res_enrich@qresult$FDR<0.01, c("Description", "FDR", "Ratio",  "BgRatio")]
#In.table$cell="In"

res_enrich = disease_enrichment( entities =Mic, vocabulary = "HGNC", database = "ALL", universe=data$Gene)
Mic.table = res_enrich@qresult[res_enrich@qresult$FDR<0.01, c("Description", "FDR", "Ratio",  "BgRatio")]
Mic.table$cell="Mic"

res_enrich = disease_enrichment( entities =Oli, vocabulary = "HGNC", database = "ALL",universe=data$Gene)
Oli.table = res_enrich@qresult[res_enrich@qresult$FDR<0.01, c("Description", "FDR", "Ratio",  "BgRatio")]
Oli.table$cell="Oli"

DO.tbl=rbind(Mic.table,Oli.table)

p.heir.DOenrich=ggplot(DO.tbl, aes(y=Description, x=cell)) + geom_point(aes(size=-log10(FDR))) +
theme_minimal()
