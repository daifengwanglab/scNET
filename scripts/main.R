rm(list=ls())

#

#step0
source('~/work/scNET-devel/scripts/functions_for_network_analysis.R')
source('~/work/scNET-devel/scripts/load_libraries.R')
source('~/work/scNET-devel/scripts/read_data.R')

out_dir="./Figures"

#step1
#source('~/work/scNET-devel/scripts/get_topological_properties.R')
#ggsave(p.topology,file=paste(out_dir,"topology.pdf",sep="/"),width = 4, height = 3, dpi = 300, units = "in", device='pdf')


#step2
source('~/work/scNET-devel/scripts/get_centrality.R')
source('~/work/scNET-devel/scripts/plot_upset_cent.R')
pdf(file="Figures/betweenness.upset.pdf")
p.bet
dev.off()
pdf(file="Figures/indegree.upset.pdf")
p.in
dev.off()
pdf(file="Figures/outdegree.upset.pdf")
p.out
dev.off()

source('~/work/scNET-devel/scripts/plot_centrality.R')


#step3
#get h metric and DE
#get rewiring



#step4
source('~/work/scNET-devel/scripts/plot_heirarchy_scores.R')
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
pdf(file="Figures/Heir.lvl.GO.pdf")
p.heir.lvlGO
dev.off()


#step4
source('~/work/scNET-devel/scripts/find_coregnet.R')


#step5
source('~/work/scNET-devel/scripts/get_motifs.R')

#step6
#loregic
