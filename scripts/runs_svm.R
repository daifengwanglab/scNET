#source code modified from: https://rpubs.com/markloessi/506713
#https://johanndejong.wordpress.com/2018/04/04/nested-cross-validation/
rm(list=ls())
library(kernlab)
library(caret)
library(mlbench)
library(PRROC)

load("data/gematrix/gexpr_All.rdata")
gexpr=t(gexpr)

cent.mat=read.table("results/centralities/centrality_matrix.all-cellTypes.txt", header=T)


markers=read.table("~/work/scNET_manuscript/data/marker_genes.txt", header=T)
markers$Cluster=gsub("Oligo","Oli",markers$Cluster)
markers$Cluster=gsub("Microglia","Mic",markers$Cluster)



get_top_cent.genes = function(cent.df,n) #eg n=0.3 for 30%
{
  name=colnames(cent.df)[2]
  tag=strsplit(name,"[.]")
  colnames(cent.df)=c("gene","score")
  cent.df=cent.df[order(-cent.df$score), , drop = FALSE]
  topn=round(dim(cent.df)[1]*n)
  cent.df=cent.df[1:topn,]
  cent.df=cent.df[cent.df$score > 0 ,]
  cent.df$score=name
  colnames(cent.df)=c("Gene","Cluster")
  #ct.markers=markers[markers$Cluster==tag[[1]][1],]
  #tbl=rbind(cent.df,ct.markers)
  #tbl
  cent.df
}

cells=c("gene","Mic.betweenness","Oli.betweenness","Ex1.betweenness","In1a.betweenness")
betweeness.mat=dplyr::select(as.data.frame(cent.mat),contains(cells))
colnames(betweeness.mat)=gsub(".betweenness","",colnames(betweeness.mat))

for (i in 2:ncol(betweeness.mat))
{
  cl=i-1 #label
  ct=paste(colnames(betweeness.mat)[i])
  df=betweeness.mat[,c(1,i)]
  assign(ct, get_top_cent.genes(df,0.1))
  tbl=get(ct)
  indx=match(tbl$Gene, rownames(gexpr))
  indx=indx[!is.na(indx)]
  tmp=gexpr[indx,]
  tmp=dplyr::select(as.data.frame(tmp),contains(ct))
  cent.expr.mat=tmp[,1:100]#biased by # of samples
  name=paste(ct,"cent.list",sep=".")
  assign(name, rownames(cent.expr.mat))
  cent.expr.mat=as.data.frame(t(cent.expr.mat))
  cent.expr.mat$Class=cl
  name=paste(ct,"tbl",sep=".")
  assign(name,cent.expr.mat)
}
list=list(Ex1=Ex1.cent.list,In1a=In1a.cent.list,Mic=Mic.cent.list,Oli=Oli.cent.list)



df1=gexpr[indx,]
df1=dplyr::select(as.data.frame(df1),contains("In1a"))
#df=df[,1:200]
df1=as.data.frame(t(df1))
df1$Class=4


#data=rbind(cent.expr.mat,df)
data=rbind(data,df)


df=gexpr[indx,]
df=dplyr::select(as.data.frame(df),-contains("Ex1"))
df=df[,1:200]
df=as.data.frame(t(df))
df$Class=0

data=rbind(Oli.bet.expr.mat,df)
#data

# Set the random seed for reproducibility
set.seed(111)

y <- data$Class
X <- data
X$Class <- NULL

train_and_validate <- function(data,fold,C){
  # Train an SVM, excluding the fold
  fit <- randomForest(Class ~ .,data = data[-fold,],ntrees = C)
  # Predict the fold
  yh <- predict(fit, newdata = data[fold,])
  # Compare the predictions to the labels
  posneg <- split(yh[,1], data$Class[fold])
  # Return the AUC under the ROC
  roc.curve(posneg[[1]], posneg[[2]])$auc
}

# Function for doing a k-fold cross-validation for each C in CC
cv <- function(data,k,CC,seed = NULL){

  # For each value of the hyperparameter C ...
  auc <- lapply(CC, function(C) {
    folds <- createFolds(data$Class, k = k)
    # For each fold ...
    sapply(folds, function(fold) {
      # Train an SVM, and validate on the fold
      train_and_validate(data,fold,C)
    })
  })
  auc
}

ncv <- function(data,k,CC) {
  folds <- createFolds(data$Class, k = k)
  # For each fold ...
  auc <- sapply(folds, function(fold) {
    # Do a cross-validation for each C
    auc <- cv(data[-fold,],k,CC)
    # Select the C with the highest AUC
    C <- CC[which.max(sapply(auc, mean))]
    # Test this C on the test data
    train_and_validate(data,fold = fold,C = C)
  })
  auc
}

auc <- ncv(
  data = data,
  k = 5,
  CC = 100
  )
