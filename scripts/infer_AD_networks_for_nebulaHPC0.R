rm(list=ls())

source('~/tools/scGRNom/code/scGRNom_getNt.R')
load('~/scGRN/data/gematrix/MIT_imputed_gexpr.rdata')
load('~/BrainCelltypeReferenceNetworks.RData')


library(readxl)
library(data.table)

chromatin_access_regions <- read_xlsx("~/tools/scGRNom/data/open_chromatin_regions.xlsx", sheet = 'Feature Binarization Peaks',skip = 16)
open_chrom_regions <- chromatin_access_regions[which(chromatin_access_regions$Microglia == 1),]
open_chrom_regions <- data.frame(na.omit(open_chrom_regions[,c('hg38_Chromosome', 'hg38_Start', 'hg38_Stop')]))

#AD Mic
AD_cells.tmp=AD_cells[,c("TAG","broad.cell.type")]
AD_cells.tmp$new.broad.cell.type = sub('[.]', '.', make.names(AD_cells.tmp$broad.cell.type, unique=TRUE))
gexpr_AD.tmp=gexpr_AD
gexpr_AD.tmp$newrows=rownames(gexpr_AD.tmp)
table=setDT(gexpr_AD.tmp)
lookup=setDT(AD_cells.tmp)
setkey(gexpr_AD.tmp, "newrows")
setkey(AD_cells.tmp, "TAG")
joined=gexpr_AD.tmp[AD_cells.tmp]
joined.df=as.data.frame(joined)
rownames(joined.df)=joined.df$new.broad.cell.type
joined.df=joined.df[,1:ncol(gexpr_AD)]
Mic.joined.df=joined.df[rownames(joined.df) %like% "Mic", ]
gexpr <- t(Mic.joined.df[,log10(colSums(Mic.joined.df)+1)> 1])

df3 <- scGRNom_getNt(df = df2.Microglia, gexpr = gexpr, open_chrom = open_chrom_regions, extension_bps = 2000,num_cores = 24)
df3$abs_coef=abs(df3$coef)
write.table(df3,file='MIT.AD.Mic.grn.txt',sep="\t",col.names=TRUE, row.names=TRUE, quote=F)

rm(gexpr, df3)

##Ctrl Mic
CTL_cells.tmp=CTL_cells[,c("TAG","broad.cell.type")]
CTL_cells.tmp$new.broad.cell.type = sub('[.]', '.', make.names(CTL_cells.tmp$broad.cell.type, unique=TRUE))
gexpr_CTL.tmp=gexpr_CTL
gexpr_CTL.tmp$newrows=rownames(gexpr_CTL.tmp)
table=setDT(gexpr_CTL.tmp)
lookup=setDT(CTL_cells.tmp)
setkey(gexpr_CTL.tmp, "newrows")
setkey(CTL_cells.tmp, "TAG")
joined=gexpr_CTL.tmp[CTL_cells.tmp]
joined.df=as.data.frame(joined)
rownames(joined.df)=joined.df$new.broad.cell.type
joined.df=joined.df[,1:ncol(gexpr_CTL)]
Mic.joined.df=joined.df[rownames(joined.df) %like% "Mic", ]
gexpr <- t(Mic.joined.df[,log10(colSums(Mic.joined.df)+1)> 1])

df3 <- scGRNom_getNt(df = df2.Microglia, gexpr = gexpr, open_chrom = open_chrom_regions, extension_bps = 2000,num_cores = 24)
df3$abs_coef=abs(df3$coef)
write.table(df3,file='MIT.Ctrl.Mic.grn.txt',sep="\t",col.names=TRUE, row.names=TRUE, quote=F)
