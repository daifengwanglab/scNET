# Regulatory network format;: 3 columns (TF target score)


#network data directory; change accordingly
net_data="~/work/scNET_manuscript/data/Brain_GRNs/mufang_yiang_cell_type_networks"

#genome data directory;
genome_data="~/work/scNET/demo_data/genome"


#Read genome data
TFs=read.table(paste(genome_data,"TF_names_v_1.01.txt",sep="/"))

#Read network data; careful of headers
ex_genie=read.table(paste(net_data,"ex/ex_GENIE3.txt",sep="/"), header=T)
ex_grnbst=read.table(paste(net_data,"ex/ex_GRNboost2_GRN.csv",sep="/"),header=T,sep=",")
ex_pidc=read.table(paste(net_data,"ex/ex_PIDC_GRN.txt", sep="/"))

in_genie=read.table(paste(net_data,"in/in_GENIE3.txt",sep="/"), header=T)
in_grnbst=read.table(paste(net_data,"in/in_GRNboost2_GRN.csv",sep="/"),header=T,sep=",")
in_pidc=read.table(paste(net_data,"in/in_PIDC_GRN.txt", sep="/"))

oligo_genie=read.table(paste(net_data,"oligo/oligo_GENIE3.txt",sep="/"), header=T)
oligo_grnbst=read.table(paste(net_data,"oligo/oligo_GRNboost2_GRN.csv",sep="/"),header=T,sep=",")
oligo_pidc=read.table(paste(net_data,"oligo/oligo_PIDC_GRN.txt", sep="/"))

micro_genie=read.table(paste(net_data,"micro/micro_GENIE3.txt",sep="/"), header=T)
micro_grnbst=read.table(paste(net_data,"micro/micro_GRNboost2_GRN.csv",sep="/"),header=T,sep=",")
micro_pidc=read.table(paste(net_data,"micro/micro_PIDC_GRN.txt", sep="/"))

#read genesets
GOBP=read.table("~/work/scNET_manuscript/genome/genesets/BaderLab/Human_GO_bp_with_GO_iea_symbol.mt3lt500.genesets", quote="")
GOBP$V2=gsub("%.*","",GOBP$V2)
gsc=loadGSC(GOBP)
