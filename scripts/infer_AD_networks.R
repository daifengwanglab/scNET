source('~/work/tools/scGRNom/code/scGRNom_getTF.R')
source('~/work/tools/scGRNom/code/scGRNom_getNt.R')
source('~/work/tools/scGRNom/code/scGRNom_interaction.R')

load('~/work/scNET_manuscript/data/Ting_Jin/AD_MIT/MIT_imputed_gexpr.rdata')

library(readxl)
interactome_data = read_xlsx("~/work/tools/scGRNom/data/PLAC-seq promoter interactome map.xlsx",sheet = 'Microglia interactome',skip = 2)[,1:6]
enhancers = read_xlsx("~/work/tools/scGRNom/data/PLAC-seq promoter interactome map.xlsx",sheet = 'Microglia enhancers',skip = 2,col_names = c('chr','start','end'))


df1 = scGRNom_interaction(interactome_data,enhancers)
df2 = scGRNom_getTF(df1)

#assign ct names to gexpr mat
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
