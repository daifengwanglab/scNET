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
            as.data.frame(c), colnames(c), width_ratio=0.1,name='Betweenness',min_size=3,
            sort_sets=FALSE,
            stripes=upset_stripes(
              geom=geom_segment(size=5),
              colors=c(npgcolors[1],npgcolors[1],npgcolors[1],npgcolors[1],
              npgcolors[2],npgcolors[2],npgcolors[2],npgcolors[2])
            ),
            themes=upset_modify_themes(
              list(
                'intersections_matrix'=theme(text=element_text(size=20)),
                'overall_sizes'=theme(text=element_text(size=10),axis.text.x=element_text(angle=90))
              )
            )
      )


#indegree
c=mat[,colnames(mat) %like% "degree_in"]
c= c[as.logical(rowSums(c != 0)), ]
colnames(c)=gsub(".degree_in.topCent","",colnames(c))
colnames(c)=gsub("\\."," ",colnames(c))
p.in=upset(
            as.data.frame(c), colnames(c), width_ratio=0.1,name='In-degree',min_size=3,
            sort_sets=FALSE,
            stripes=upset_stripes(
              geom=geom_segment(size=5),
              colors=c(npgcolors[1],npgcolors[1],npgcolors[1],npgcolors[1],
              npgcolors[2],npgcolors[2],npgcolors[2],npgcolors[2])
            ),
            themes=upset_modify_themes(
              list(
                'intersections_matrix'=theme(text=element_text(size=20)),
                'overall_sizes'=theme(text=element_text(size=10),axis.text.x=element_text(angle=90))
              )
            )
      )



c=mat[,colnames(mat) %like% "degree_out"]
c= c[as.logical(rowSums(c != 0)), ]
colnames(c)=gsub(".degree_out.topCent","",colnames(c))
colnames(c)=gsub("\\."," ",colnames(c))
p.out=upset(
            as.data.frame(c), colnames(c), width_ratio=0.1,name='Out-degree',min_size=3,
            sort_sets=FALSE,
            stripes=upset_stripes(
              geom=geom_segment(size=5),
              colors=c(npgcolors[1],npgcolors[1],npgcolors[1],npgcolors[1],
              npgcolors[2],npgcolors[2],npgcolors[2],npgcolors[2])
            ),
            themes=upset_modify_themes(
              list(
                'intersections_matrix'=theme(text=element_text(size=20)),
                'overall_sizes'=theme(text=element_text(size=10),axis.text.x=element_text(angle=90))
              )
            )
      )
