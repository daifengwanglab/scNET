#heirarchy
library(data.table)
library(caret)
library(easyalluvial)
library('alluvial')


celltypes=c("MIC","Oli","Ex","In")

##https://stackoverflow.com/questions/51991825/r-alluvial-plots-for-combinations


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
p.mic=easyalluvial::alluvial_wide( Mic, id = Gene, fill_by = 'first_variable')


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
p.Oli=easyalluvial::alluvial_wide( Oli, id = Gene, fill_by = 'first_variable')

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
p.Ex=easyalluvial::alluvial_wide( Ex, id = Gene, fill_by = 'first_variable')


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
p.In=easyalluvial::alluvial_wide( In, id = Gene, fill_by = 'first_variable')
