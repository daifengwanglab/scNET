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
  theme_sankey(base_size = 18) +
  labs(x = NULL) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = .5))


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
  theme_sankey(base_size = 18) +
  labs(x = NULL) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = .5))


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
  theme_sankey(base_size = 18) +
  labs(x = NULL) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = .5))



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


pdf(file="Figures/Hei.In.pdf")
p.In
dev.off()

pdf(file="Figures/Hei.Ex.pdf")
p.Ex
dev.off()

pdf(file="Figures/Hei.Mic.pdf")
p.mic
dev.off()

pdf(file="Figures/Hei.Oli.pdf")
p.Oli
dev.off()


#combine all and find enrichment

list=ls(pattern=".gene.lev")
data = Reduce(function(x, y) merge(x, y, all=T), lapply(list,get), accumulate=F)

 #df = data %>% make_long(colnames(data[,2:ncol(data)]))
 p.data=ggplot(df, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = factor(node), label = node)) +
   geom_sankey(flow.alpha = .6,
               node.color = "gray30") +
     scale_fill_npg() +
   theme_sankey(base_size = 18) +
   labs(x = NULL) +
   theme(legend.position = "none",
         plot.title = element_text(hjust = .5))

data[is.na(data)]="None"
write.table(data, file="~/work/scNET_manuscript/AD_MIT/heirarchy/all_genes_heirarchy_levs.txt",sep="\t",row.names=F, col.names=T, quote=F)

#enrichment
Ex=data[data$Ex.AD!=data$Ex.Ctrl,]$Gene
In=data[data$In.AD!=data$In.Ctrl,]$Gene
Mic=data[data$Mic.AD!=data$Mic.Ctrl,]$Gene
Oli=data[data$Oli.AD!=data$Oli.Ctrl,]$Gene

#GO enrichment analysis
GOdata=GSA.read.gmt('~/work/scNET_manuscript/genome/genesets/GO_annotations-9606-inferred-allev.gmt')
genesets=GOdata$genesets
names(genesets)=GOdata$geneset.descriptions
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


#disease ontology enrichment
res_enrich = disease_enrichment( entities =Ex, vocabulary = "HGNC", database = "ALL")
Ex.table = res_enrich@qresult[res_enrich@qresult$FDR<0.01, c("Description", "FDR", "Ratio",  "BgRatio")]
Ex.table$cell="Ex"

res_enrich = disease_enrichment( entities =In, vocabulary = "HGNC", database = "ALL")
In.table = res_enrich@qresult[res_enrich@qresult$FDR<0.01, c("Description", "FDR", "Ratio",  "BgRatio")]
In.table$cell="In"

res_enrich = disease_enrichment( entities =Mic, vocabulary = "HGNC", database = "ALL")
Mic.table = res_enrich@qresult[res_enrich@qresult$FDR<0.01, c("Description", "FDR", "Ratio",  "BgRatio")]
Mic.table$cell="Mic"

res_enrich = disease_enrichment( entities =Oli, vocabulary = "HGNC", database = "ALL")
Oli.table = res_enrich@qresult[res_enrich@qresult$FDR<0.01, c("Description", "FDR", "Ratio",  "BgRatio")]
Oli.table$cell="Oli"

DO.tbl=rbind(Ex.table,In.table,Mic.table,Oli.table)

p.heir.DOenrich=ggplot(DO.tbl, aes(y=Description, x=cell)) + geom_point(aes(size=-log10(FDR))) +
theme_minimal()


Ex.data2 <- gene2disease(gene = Ex,score =c(0.2, 1),api_key=api)
pdf(file="Figures/Ex.heir.diseaseclass.heatmap.pdf")
plot(Ex.data2,class="DiseaseClass", nchars=30)
dev.off()


In.data2 <- gene2disease(gene = In,score =c(0.2, 1),api_key=api)
pdf(file="Figures/In.heir.diseaseclass.heatmap.pdf")
plot(In.data2,class="DiseaseClass", nchars=30)
dev.off()


Mic.data2 <- gene2disease(gene = Mic,score =c(0.2, 1),api_key=api)
pdf(file="Figures/Mic.heir.diseaseclass.heatmap.pdf")
plot(Mic.data2,class="DiseaseClass", nchars=30)
dev.off()


Oli.data2 <- gene2disease(gene = Oli,score =c(0.2, 1),api_key=api)
pdf(file="Figures/Oli.heir.diseaseclass.heatmap.pdf")
plot(Oli.data2,class="DiseaseClass", nchars=30)
dev.off()
