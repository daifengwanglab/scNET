
rm(list=ls())
rm(list=ls())

source('../scripts/read_data.R')
source('../scripts/functions_for_network_analysis.R')


celltypes=c("Mic","Oli","Ex","In")


#disease gene database
alz=read.table("../alz.txt")

train_and_validate = function( data, fold, C)
{
  fit = randomForest(formula= as.factor(Class) ~ ., data = data[-fold,],importance=TRUE)

  # Predict the fold
  yh = predict(fit, newdata = data[fold,])

  # Return the balanced accuracy
  conf.mat=caret::confusionMatrix(yh, as.factor(data[fold,]$Class))
  acc=conf.mat$byClass['Balanced Accuracy']
  acc
}

# Function for doing a k-fold cross-validation for each mtry in mtry
cv = function(data, k, mtry)
{
  # For each value of mtry ...
  auc = lapply(mtry, function(C)
  {
    folds = createFolds(data$Class, k = k)
    sapply(folds, function(fold)
    {
      train_and_validate(data,fold,C)
    })
  })
  auc
}

cv_feat_imp = function(data, k)
{
    folds <- createFolds(data$Class, k = k)

    # For each fold ...
    imp=sapply(folds, function(fold)
    {
      get_weights(data,fold)
    })
  imp
}


nets=ls(pattern="*\\.network")

#label prediction for all cell types

for(i in 1:length(nets))
{

name=nets[i]
celltype=gsub(".network","",name)
net=as.data.frame(lapply(nets[i],get))

net=net[,c("TF","TG","abs_coef")]
net=distinct(net)

mat=tidyr::pivot_wider(net, names_from = TF, values_from = abs_coef)
#mat[is.na(mat)]=0
mat[is.na(mat)]=min(net$abs_coef)*0.01
mat=as.data.frame(mat)


#add TFnodes with no indegrees
TF.nodes=unique(net$TF)
df = data.frame(matrix(ncol = ncol(mat), nrow = length(TF.nodes)))
colnames(df)=colnames(mat)
df$TG=TF.nodes
df[is.na(df)]=min(net$abs_coef)*0.01
#df[is.na(df)]=0
df=df[!(df$TG %in% intersect(TF.nodes,mat$TG)),]
mat=rbind(mat,df)

#label alz genes as positives in the network matrix
mat$Class=ifelse(mat$TG %in% alz,1,-1)

mat.positive=mat[mat$Class == 1,]

print(dim(mat.positive))

##create negative labels as those that are not positives and also not related to mental disorders
#mat.negative=mat[mat$Class == -1 & (!mat$TG %in% uncorrDis.genes),]
mat.negative=mat[mat$Class == -1 ,]


set.seed(123)
#randomly select rows equal to the number of positives
mat.tmp=sample_n(mat.negative,dim(mat.positive)[1])

#train matrix
mat.feat=rbind(mat.positive,mat.tmp)

#test matrix
#final_test=mat[!(mat$TG %in% mat.feat$TG),]
final_test=mat
colnames(final_test)=gsub("-","_",colnames(final_test))
rownames(final_test)=final_test$TG
final_test$TG=NULL

# Set the random seed for reproducibility
data=mat.feat[,-1]
colnames(data)=gsub("-","_",colnames(data))

fit = randomForest(formula= as.factor(Class) ~ ., data = data,importance=TRUE)

# Predict the fold
yh = as.data.frame(predict(fit, newdata = final_test, type="prob"))
yh=yh[order(yh[,1]),]
yh$prior=ifelse(rownames(yh) %in% alz.dis,1,0)


name=paste(celltype,"full_AD_gene_prediction.txt",sep=".")
filename=paste("results",name,sep="")
write.table(yh,file=filename,row.names=T,col.names=T,sep="\t",quote=F)

}#end of all ct label pred for loop
