rm(list=ls())

#predictions
fw=read.table("~/work/scNET_manuscript/AD_MIT/supp_data/AD_gene_prediction.mic.txt", header=T)
fw$gene=rownames(fw)
rownames(fw)=NULL
fw=fw[order(-fw[,2]),]
top=fw[1:round(nrow(fw)*0.2),]$gene

#Read GO data
data=GSA.read.gmt('/Users/chiraggupta/work/scNET_manuscript/genome/genesets/BaderLab/Human_GO_bp_with_GO_iea_symbol.gmt')
genesets=data$genesets
names(genesets)=data$geneset.descriptions



hyp_obj = hypeR(top, genesets,fdr=0.001, background=fw$gene)
hyp_df =  hyp_obj$data
hyp_df = hyp_df[hyp_df$geneset < 200,]

hyp_df$logFDR=-1*log10(hyp_df$fdr)


#select top 20 terms

hyp_df = hyp_df[,c("label","fdr","logFDR")]
hyp_df = hyp_df[order(-hyp_df$logFDR),]
remove_terms=c("hematopoietic or lymphoid organ development","cellular response to lipid")
hyp_df=hyp_df[!(hyp_df$label %in% remove_terms) ,]

hyp_df=hyp_df[!(hyp_df$label %like% "regulation of") & !(hyp_df$label %like% "biosynthetic process") & !(hyp_df$label %like% "transcription"),]
write.table(hyp_df,file="~/work/scNET_manuscript/AD_MIT/supp_data/top_genes_bp.txt",sep="\t",col.names=T, row.names=F, quote=F)
hyp_df=head(hyp_df,20)


p.topgene.enrich=ggplot(hyp_df,aes(y=reorder(label,logFDR),x=logFDR, fill=reorder(label,-logFDR))) +
geom_bar(stat = "summary") + scale_fill_grey() +
labs(y="Biological process",x="-1*log(FDR)")+
theme_bw(base_size=12)+theme(legend.position="none")
ggsave(p.topgene.enrich,filename="Figures/p.topgene.enrich.pdf", device="pdf",width=3,height=3,units="in")
