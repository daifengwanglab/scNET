rm(list=ls())

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
pdf(file=paste(out_dir,"feature_heatmap.pdf",sep="/"))
draw(p.heatmap, heatmap_legend_side = "left", annotation_legend_side = "bottom")
dev.off()
pdf(file=paste(out_dir,"Centrality.GOgsea.pdf",sep="/"))
p.GOgsea
dev.off()
pdf(file=paste(out_dir,"Centrality.Kegg.pdf",sep="/"))
p.Kegg
dev.off()

#step3
source('~/work/scNET-devel/scripts/get_heirarchy_scores.R')

#step4
source('~/work/scNET-devel/scripts/find_coregnet.R')


#step5
source('~/work/scNET-devel/scripts/get_motifs.R')

#step6
#loregic
