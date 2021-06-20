rm(list=ls())
source('~/work/scNET-devel/scripts/load_libraries.R')
source('~/work/scNET-devel/scripts/read_data.R')
source('~/work/scNET-devel/scripts/functions_for_network_analysis.R')
source('~/work/scNET-devel/scripts/find_coregnet.R')


diseasegenes=read.csv("~/Desktop/all_gene_disease_associations.tsv", header=T, sep="\t")

nets=ls(pattern="*.coregnet.igraph")

random=sample(unique(diseasegenes$geneSymbol), 4000)

AD=diseasegenes[diseasegenes$diseaseName %like% "Alzheimer's Disease",]
AD=unique(AD$geneSymbol)

BD=diseasegenes[diseasegenes$diseaseName %like% "Bipolar Disorder",]
BD=unique(BD$geneSymbol)

Epilepsy=diseasegenes[diseasegenes$diseaseName %like% "Epilepsy",]
Epilepsy=unique(Epilepsy$geneSymbol)

Park=diseasegenes[diseasegenes$diseaseName %like% "Parkinson",]
Park=unique(Park$geneSymbol)

HD=diseasegenes[diseasegenes$diseaseName %like% "Huntington",]
HD=unique(HD$geneSymbol)

ALS=diseasegenes[diseasegenes$diseaseName %like% "sclerosis",]
ALS=unique(ALS$geneSymbol)


Skin=diseasegenes[diseasegenes$diseaseName %like% "Skin lesion",]
Skin=unique(Skin$geneSymbol)

DM=diseasegenes[diseasegenes$diseaseName %like% "Diabetes Mellitus",]
DM=unique(DM$geneSymbol)


density.tbl=data.frame(cell=NULL,density=NULL, disease=NULL)

for(i in 1:length(nets))
{
  name=nets[i]
  tag=gsub(".igraph","",name)
  name=paste(tag,"density",sep=".")
  allNodes=as_ids(V(get(nets[i])))
  geneset=intersect(random,allNodes)
  assign(name, calculate_geneset_density(get(nets[i]),geneset))
  newtbl=data.frame(cell=tag,density=get(name),disease="Random4")
  density.tbl=rbind(density.tbl,newtbl)
}


neurodegen=c("AD","ALS","BD","Epilepsy","HD","Parkinson's","Random")
density.tbl=density.tbl %>% filter(disease %in% neurodegen)

#https://nitinahuja.github.io/2017/heatmaps-in-r/
ggplot(density.tbl , aes(y=cell, x=disease, fill = density)) +
    geom_tile(color = "white", size = 0.1) +
    scale_x_discrete(expand=c(0,0)) +
    scale_y_discrete(expand=c(0,0)) +
    scale_fill_viridis(name="Density", option = "plasma") +
    coord_equal() +
    theme_tufte(base_family="Helvetica") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
