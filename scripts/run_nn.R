#source code modified from: https://www.r-bloggers.com/2017/02/multilabel-classification-with-neuralnet-package/
#https://johanndejong.wordpress.com/2018/04/04/nested-cross-validation/
rm(list=ls())

require(neuralnet)
require(nnet)
require(ggplot2)
set.seed(10)


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

i=59
ct=paste(colnames(cent.mat)[i])
ct=paste(ct,"tbl",sep=".")
df=cent.mat[,c(1,i)]
assign(ct, get_top_cent.genes(df,0.1))
tbl=get(ct)
indx=match(tbl$Gene, rownames(gexpr))
indx=indx[!is.na(indx)]
tmp=gexpr[indx,]
tmp=dplyr::select(as.data.frame(tmp),contains("Mic"))
Oli.bet.expr.mat=tmp[,1:200]
Oli.bet.expr.mat=as.data.frame(t(Oli.bet.expr.mat))
Oli.bet.expr.mat$Class=1

df=gexpr[indx,]
df=dplyr::select(as.data.frame(df),-contains("Mic"))
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
  fit <- ksvm(
    Class ~ .,
    data = data[-fold,],
    kernel = "vanilladot",
    kpar = list(),
    C = C,
    type="C-svc",
    prob.model = TRUE
  )
  # Predict the fold
  yh <- predict(fit, newdata = data[fold,], type = "probabilities")
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
