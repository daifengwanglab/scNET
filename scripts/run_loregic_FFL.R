library(dplyr)
library(Loregic)
library(dplyr)
library(ArrayBin)
library(tidyr)
library(data.table)

args = commandArgs(trailingOnly=TRUE) #FFL_file #idmap_file celltype condition #permut_for_loregic

data=read.table(args[1],skip=5, header=F)
ids=read.table(args[2],header=T)
ct=args[3]
cond=args[4]
gexpr=args[5]
npermut=as.numeric(args[6])


colnames(data)=c("RF1","RF2","target")

#map ids to symbols
new=data
new[] = lapply(data, function(x) ids$symbol[match(x, ids$id)])
new=distinct(new)

#load expression matrix
ct.gexpr.bin=read.table(gexpr, header=T, row.names=1)

nodes=append(unique(as.character(new$RF1)),unique(as.character(new$RF2)))
nodes=append(unique(nodes),unique(as.character(new$target)))
nodes=unique(nodes)
indx=match(nodes,rownames(gexpr))
indx=indx[!is.na(indx)]
nodes.ct.gexpr=gexpr[indx,]
nodes.ct.gexpr.bin=binarize.array(nodes.ct.gexpr[,1:100])
ct.loregic.out=loregic(new,nodes.ct.gexpr.bin)

filename=paste(ct,cond,sep=".")
filename=paste(filename,"ffl_logics.txt",sep=".")
write.table(ct.loregic.out,file=filename, col.names=TRUE, row.names=FALSE, sep="\t", quote=F)

loregicOut = ct.loregic.out

#find gate consistent trips and calculate pvalues
rownames(loregicOut)=paste(loregicOut$RF1,loregicOut$RF2,loregicOut$target,sep="_")
loregicOut=loregicOut[,-c(1:3)]
mat=loregicOut[apply(loregicOut[,], 1, function(x) !all(x==0)),] # remove gates with no scores
#colnames(mat)=c("T=0","AND","RF1*~RF2","RF1","~RF1*RF2","RF2","XOR","OR","NOR","XNOR","~RF2","RF1+~RF2","~RF1","~RF1+RF2","NAND","T=1")

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

#count occurences per gate
#df=data.frame(matrix(NA, nrow = 1, ncol = 1))
#colnames(df)=c("Count")
df=data.frame(Count=NULL)
for (i in 1:ncol(DF))
{
  count=as.data.frame(dim(DF[which(DF[,i]>0),])[1])
  rownames(count)=colnames(DF)[i]
  colnames(count)="Count"
  df=rbind(df,count)
}

df$cell=ct

targets=unique(gate_consistent_trips$target)
pvalues=data.frame(matrix(ncol=2,nrow=1))
colnames(pvalues)=c("pvalue","qvalue")


print("nrow(gate_consistent_trips)")
for(i in 1:nrow(gate_consistent_trips))
{
    print (paste("triplet:",i))
    count=1
    trip.gate=colnames(gate_consistent_trips[,1:14])[max.col(gate_consistent_trips[i,1:14])]
    trip=gate_consistent_trips[i,c("RF1","RF2","target")]
    for(j in 1:npermut)
    {
      rand.target=targets[[sample(1:length(targets), 1)]]
      rand.trip=trip
      rand.trip$target=rand.target
      rand.loregic.out=loregic(trip, nodes.ct.gexpr.bin)
      rand.gate=colnames(rand.loregic.out[,4:19])[max.col(rand.loregic.out[1,4:19])]
      compare=identical(rand.gate,trip.gate)
      if(compare=='TRUE')
      {
        count=count+1
      }
    }

    print (count)
    pval=count/(npermut+1)
    pvalues=rbind(pvalues,pval)
}

pvalues=pvalues[-1,]
gate_consistent_trips$pvalue=pvalues$pvalue

filename=paste(ct,cond,sep=".")
filename=paste(filename,"gate_consistent_FLL_trips.txt",sep=".")
filename=paste("results",filename,sep="/")
write.table(gate_consistent_trips, file=filename, sep="\t",col.names=TRUE, row.names=FALSE)
