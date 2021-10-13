
rm(list=ls())
source('~/work/scNET-devel/scripts/functions_for_network_analysis.R')
source('~/work/scNET-devel/scripts/load_libraries.R')


###plotting overlaps

cent.mat=read.table('centralities/centrality_matrix.all-cellTypes.txt', header=T) #first col is gene name
cent.mat[is.na(cent.mat)]=0


for (i in 2:ncol(cent.mat))
{
  name=colnames(cent.mat)[i]
  list = GetTopCentralGenes(cent.mat[,c("gene",name)],0.1)
  name=paste(name,"topCent",sep=".")
  assign(name, list)
}
list=ls(pattern="topCent")
newlist=lapply(list,function(x) get(x))
names(newlist)=list
mat=list_to_matrix(newlist)

#betweenness
c=mat[,colnames(mat) %like% "betweenness"]
c= c[as.logical(rowSums(c != 0)), ]
colnames(c)=gsub(".betweenness.topCent","",colnames(c))
colnames(c)=gsub("\\."," ",colnames(c))
p.bet=upset(
            as.data.frame(c), colnames(c), width_ratio=0.1,name='Betweenness', min_size=4,
            sort_sets=FALSE,
            stripes=upset_stripes(
              geom=geom_segment(size=5),
              colors=c(npgcolors[1],npgcolors[1],npgcolors[1],npgcolors[1],
              npgcolors[2],npgcolors[2],npgcolors[2],npgcolors[2])
            ),
            themes=upset_modify_themes(
              list(
                'intersections_matrix'=theme(text=element_text(size=10)),
                'overall_sizes'=theme(text=element_text(size=10),axis.text.x=element_text(angle=90))
              )
            )
      )

pdf(file="Figures/p.betweenness.upset.pdf", width=4,height=4)
p.bet
dev.off()

#make network fig for common Ex betweenness genes
m=make_comb_mat(c)
ex.bet.Ad.ctrl.common=extract_comb(m,"10001000")
ex.net.AD=read.table("/Users/chiraggupta/work/scNET_manuscript/AD_MIT/data/MIT.AD.Ex.grn.txt",header=T, sep="\t")
ex.net.AD=ex.net.AD[ex.net.AD$mse<0.1 & ex.net.AD$abs_coef > 0.5,]
ex.net.AD.between=ex.net.AD[ex.net.AD$TF %in% ex.bet.Ad.ctrl.common,]
ex.net.AD.between=distinct(ex.net.AD.between[,c("TF","TG","abs_coef")])
write.table(ex.net.AD.between,file="ex.net.AD.between.net.txt",sep="\t",quote=F,col.names=T,row.names=F)

ex.net.Ctrl=read.table("/Users/chiraggupta/work/scNET_manuscript/AD_MIT/data/MIT.Ctrl.Ex.grn.txt",header=T, sep="\t")
ex.net.Ctrl=ex.net.Ctrl[ex.net.Ctrl$mse<0.1 & ex.net.Ctrl$abs_coef > 0.5,]
ex.net.Ctrl.between=ex.net.Ctrl[ex.net.Ctrl$TF %in% ex.bet.Ad.ctrl.common,]
ex.net.Ctrl.between=distinct(ex.net.Ctrl.between[,c("TF","TG","abs_coef")])
write.table(ex.net.Ctrl.between,file="ex.net.Ctrl.between.net.txt",sep="\t",quote=F,col.names=T,row.names=F)


#indegree
c=mat[,colnames(mat) %like% "degree_in"]
c= c[as.logical(rowSums(c != 0)), ]
colnames(c)=gsub(".degree_in.topCent","",colnames(c))
colnames(c)=gsub("\\."," ",colnames(c))
p.in=upset(
            as.data.frame(c), colnames(c), width_ratio=0.1,name='In-degree',min_size=5,
            sort_sets=FALSE,
            stripes=upset_stripes(
              geom=geom_segment(size=5),
              colors=c(npgcolors[1],npgcolors[1],npgcolors[1],npgcolors[1],
              npgcolors[2],npgcolors[2],npgcolors[2],npgcolors[2])
            ),
            themes=upset_modify_themes(
              list(
                'intersections_matrix'=theme(text=element_text(size=10)),
                'overall_sizes'=theme(text=element_text(size=10),axis.text.x=element_text(angle=90))
              )
            )
      )

pdf(file="Figures/p.indegree.upset.pdf", width=4,height=4)
p.in
dev.off()

c=mat[,colnames(mat) %like% "degree_out"]
c= c[as.logical(rowSums(c != 0)), ]
colnames(c)=gsub(".degree_out.topCent","",colnames(c))
colnames(c)=gsub("\\."," ",colnames(c))
p.out=upset(
            as.data.frame(c), colnames(c), width_ratio=0.1,name='Out-degree',min_size=5,
            sort_sets=FALSE,
            stripes=upset_stripes(
              geom=geom_segment(size=5),
              colors=c(npgcolors[1],npgcolors[1],npgcolors[1],npgcolors[1],
              npgcolors[2],npgcolors[2],npgcolors[2],npgcolors[2])
            ),
            themes=upset_modify_themes(
              list(
                'intersections_matrix'=theme(text=element_text(size=10)),
                'overall_sizes'=theme(text=element_text(size=10),axis.text.x=element_text(angle=90))
              )
            )
      )

pdf(file="Figures/p.outdegree.upset.pdf", width=4,height=4)
p.out
dev.off()
