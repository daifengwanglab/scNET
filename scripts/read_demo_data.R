# Regulatory network format;: 3 columns (TF target score)


#network data directory; change accordingly
net_data="./demo_data"

#genome data directory; change accordingly
genome_data="./demo_data/genome"


#Read genome data
TFs=read.table(paste(genome_data,"TF_names_v_1.01.txt",sep="/"))
GOBP=read.table(paste(genome_data,"Human_GO_bp_with_GO_iea_symbol.mt3lt500.genesets",sep="/"), quote="")
GOBP.gs=loadGSC(GOBP)

#Read network data; careful of headers
ex_genie=read.table(paste(net_data,"ex_genie.demo.txt",sep="/"))
ex_grnbst=read.table(paste(net_data,"ex_grnbst.demo.txt",sep="/"))
ex_pidc=read.table(paste(net_data,"ex_pidc.demo.txt", sep="/"))

in_genie=read.table(paste(net_data,"in_genie.demo.txt",sep="/"))
in_grnbst=read.table(paste(net_data,"in_grnbst.demo.txt",sep="/"))
in_pidc=read.table(paste(net_data,"in_pidc.demo.txt", sep="/"))

oligo_genie=read.table(paste(net_data,"oligo_genie.demo.txt",sep="/"))
oligo_grnbst=read.table(paste(net_data,"oligo_grnbst.demo.txt",sep="/"))
oligo_pidc=read.table(paste(net_data,"oligo_pidc.demo.txt", sep="/"))

micro_genie=read.table(paste(net_data,"micro_genie.demo.txt",sep="/"))
micro_grnbst=read.table(paste(net_data,"micro_grnbst.demo.txt",sep="/"))
micro_pidc=read.table(paste(net_data,"micro_pidc.demo.txt", sep="/"))
