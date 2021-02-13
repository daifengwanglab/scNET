source('~/work/scNetAnalysis.beta/scripts/load_libraries.R')
source('~/work/scNetAnalysis.beta/scripts/read_data.R')

create_demo_networks = function(net) #network and threshold for # of edges to report
{
	#rank edges of each network from 1 to n
	colnames(net)=c("regulatoryGene","targetGene","weight")
	net= net %>% filter(regulatoryGene %in% TFs$V1)
	net
}
write.table(head(create_demo_networks(ex_genie),100000), file="ex_genie.demo.txt", sep="\t", quote=F, row.names=FALSE, col.names=FALSE)
write.table(head(create_demo_networks(in_genie),100000), file="in_genie.demo.txt", sep="\t", quote=F, row.names=FALSE, col.names=FALSE)
write.table(head(create_demo_networks(oligo_genie),100000), file="oligo_genie.demo.txt", sep="\t", quote=F, row.names=FALSE, col.names=FALSE)
write.table(head(create_demo_networks(micro_genie),100000), file="micro_genie.demo.txt", sep="\t", quote=F, row.names=FALSE, col.names=FALSE)

write.table(head(create_demo_networks(ex_grnbst),100000), file="ex_grnbst.demo.txt", sep="\t", quote=F, row.names=FALSE, col.names=FALSE)
write.table(head(create_demo_networks(in_grnbst),100000), file="in_grnbst.demo.txt", sep="\t", quote=F, row.names=FALSE, col.names=FALSE)
write.table(head(create_demo_networks(oligo_grnbst),100000), file="oligo_grnbst.demo.txt", sep="\t", quote=F, row.names=FALSE, col.names=FALSE)
write.table(head(create_demo_networks(micro_grnbst),100000), file="micro_grnbst.demo.txt", sep="\t", quote=F, row.names=FALSE, col.names=FALSE)


write.table(head(create_demo_networks(ex_pidc),100000), file="ex_pidc.demo.txt", sep="\t", quote=F, row.names=FALSE, col.names=FALSE)
write.table(head(create_demo_networks(in_pidc),100000), file="in_pidc.demo.txt", sep="\t", quote=F, row.names=FALSE, col.names=FALSE)
write.table(head(create_demo_networks(oligo_pidc),100000), file="oligo_pidc.demo.txt", sep="\t", quote=F, row.names=FALSE, col.names=FALSE)
write.table(head(create_demo_networks(micro_pidc),100000), file="micro_pidc.demo.txt", sep="\t", quote=F, row.names=FALSE, col.names=FALSE)
