#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE) #loregic_out_file binarized_gexprmat celltype (e.g: Ex1)

library(tidyverse)
library(Loregic)

loregicOut=read.table(args[1], header=T, sep="\t")
bin.mat=read.table(args[2], header=T, row.names=1, sep="\t")

tmp=loregicOut %>% unite(x, c(RF1,RF2,target), sep="_")
rownames(tmp)=tmp$x
tmp=tmp[,-1]
mat=tmp[apply(tmp[,], 1, function(x) !all(x==0)),] # remove gates with no scores
colnames(mat)=c("T=0","AND","RF1*~RF2","RF1","~RF1*RF2","RF2","XOR","OR","NOR","XNOR","~RF2","RF1+~RF2","~RF1","~RF1+RF2","NAND","T=1")

l=apply(mat,1,function(x) which(x==max(x))) #store max of each row in a list

DF=data.frame(matrix(NA, nrow = 1, ncol = 16))
colnames(DF)=colnames(mat)

#for each list, check of the number of items in more than 2; if yes then gate-inconsistent
for(i in 1:length(l))
{
    nGates=length(l[[i]])
    if(nGates<2)
    {
      DF=rbind(DF,mat[i,]) #store all gate consistent in a DF df
    }
}

DF=DF[-1,] #stores all gate consistent triplets
tmp=DF
tmp$triplet=rownames(tmp)
tmp=tmp %>% separate(triplet, c("RF1","RF2","target"),sep="_")
rownames(tmp)=c()
gate_consistent_trips=tmp

targets=unique(gate_consistent_trips$target)
pvalues=data.frame(matrix(ncol=2,nrow=1))
colnames(pvalues)=c("pvalue","qvalue")



for(i in 1:nrow(gate_consistent_trips))
{
    print (paste("triplet:",i))
    count=0
    trip.gate=colnames(gate_consistent_trips[,1:14])[max.col(gate_consistent_trips[i,1:14])]
    trip=gate_consistent_trips[i,c("RF1","RF2","target")]
    for(j in 1:10)
    {
      rand.target=targets[[sample(1:length(targets), 1)]]
      rand.trip=trip
      rand.trip$target=rand.target
      rand.loregic.out=loregic(trip, bin.mat)
      rand.gate=colnames(rand.loregic.out[,4:19])[max.col(rand.loregic.out[1,4:19])]
      compare=identical(rand.gate,trip.gate)
      if(compare=='TRUE')
      {
        count=count+1
      }
    }
    pval=(1+count)/10
    pvalues=rbind(pvalues,pval)
}

pvalues=pvalues[-1,]
gate_consistent_trips.pvals=merge(gate_consistent_trips,pvalues)

write.table(gate_consistent_trips.pvals, file=paste(args[2],"gate_consistent.txt",sep="."), sep="\t",col.names=TRUE, row.names=FALSE)

count occurences per gate
df=data.frame(matrix(NA, nrow = 1, ncol = 1))
colnames(df)=c("Count")
for (i in 1:ncol(DF))
{
  count=as.data.frame(dim(DF[which(DF[,i]>0),])[1])
  rownames(count)=colnames(DF)[i]
  colnames(count)="Count"
  df=rbind(df,count)
}
write.table(df, file=paste(args[2],"gate_consistent_counts.txt",sep="."), sep="\t",col.names=TRUE, row.names=FALSE)
