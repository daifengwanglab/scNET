#source code modified from: https://rpubs.com/markloessi/506713
#https://johanndejong.wordpress.com/2018/04/04/nested-cross-validation/
#https://stackoverflow.com/questions/48379502/generate-a-confusion-matrix-for-svm-in-e1071-for-cv-results


rm(list=ls())
source('~/work/scNET-devel/scripts/load_libraries.R')
#source('~/work/scNET-devel/scripts/read_data.R')
source('~/work/scNET-devel/scripts/functions_for_network_analysis.R')


#feature matrix from AD-gene classifier
fw=read.table("Mic.AD.Feature_weights.mat", header=T)
fw$absaverage=abs(fw$average)
fw=fw[order(fw$average),]
positives=rownames(tail(fw,dim(fw)[1]*.1))
negatives=rownames(head(fw,dim(fw)[1]*.1))
labels=union(positives,negatives)



#read rosmap data
rosmap.ge=read.table("~/Desktop/finalAlzheimersDLPFCDataSetForSaniyaToUse.csv", sep=",", header=T)
rosmap.ge$entrezID=NULL
colnames(rosmap.ge)=gsub("X","",colnames(rosmap.ge))

rosmap.pheno=read.table("~/Desktop/AlzheimersDLPFCPhenotypesUpdated.csv", sep=",", header=T, )[-c(1:5)]


#cogdx
CERADScore.patients=rosmap.pheno[rosmap.pheno$cogdx==4 | rosmap.pheno$cogdx==5,]

CERADScore.patients.ge=rosmap.ge[,colnames(rosmap.ge) %in% c(CERADScore.patients$Patient,"geneName")]

CERADScore.patients.ge.labels=CERADScore.patients.ge[CERADScore.patients.ge$geneName %in% labels,]
tmp=scale(CERADScore.patients.ge.labels[,-1])
rownames(tmp)=CERADScore.patients.ge.labels$geneName
CERADScore.patients.ge.labels=as.data.frame(tmp)

colnames(CERADScore.patients.ge.labels)=paste("X",colnames(CERADScore.patients.ge.labels),sep="")


CERADScore.patients.ge.labels$Class=ifelse(rownames(CERADScore.patients.ge.labels)%in% positives,1,-1)

data.cogdx=CERADScore.patients.ge.labels

#cerad
CERADScore.patients=rosmap.pheno[rosmap.pheno$CERADScore==1,]

CERADScore.patients.ge=rosmap.ge[,colnames(rosmap.ge) %in% c(CERADScore.patients$Patient,"geneName")]

CERADScore.patients.ge.labels=CERADScore.patients.ge[CERADScore.patients.ge$geneName %in% labels,]
tmp=scale(CERADScore.patients.ge.labels[,-1])
rownames(tmp)=CERADScore.patients.ge.labels$geneName
CERADScore.patients.ge.labels=as.data.frame(tmp)

colnames(CERADScore.patients.ge.labels)=paste("X",colnames(CERADScore.patients.ge.labels),sep="")


CERADScore.patients.ge.labels$Class=ifelse(rownames(CERADScore.patients.ge.labels)%in% positives,1,-1)

data.cerad=CERADScore.patients.ge.labels




train_and_validate = function( data, fold, C) #returns average balanced acc and feat, imp. scores for each fold
{
  fit = randomForest(formula= as.factor(Class) ~ ., data = data[-fold,],importance=TRUE)
  # Predict the fold
  yh = predict(fit, newdata = data[fold,])

  # Compare the predictions to the labels
  #posneg = split(yh, data$Class[fold])

  # Return the AUC under the ROC
  #roc.curve(posneg[[1]], posneg[[2]])$auc
  conf.mat=caret::confusionMatrix(yh, as.factor(data[fold,]$Class))
  acc=conf.mat$byClass['Balanced Accuracy']
  acc
}

# Function for doing a k-fold cross-validation for each C in CC
cv = function(data, k, CC, seed = 123)
{
  # For each value of the hyperparameter C ...
  auc = lapply(CC, function(C)
  {
    folds <- createFolds(data$Class, k = k)

    # For each fold ...
    sapply(folds, function(fold)
    {
      # Train an SVM, and validate on the fold
      train_and_validate(data,fold,C)
    })#end sapply
  })#end lapply
  auc
}

k=5
  auc.cerad = cv(
    data = data.cerad,
    k = k,
    CC=10,
    seed = 123
  )

  auc.cogdx = cv(
    data = data.cogdx,
    k = k,
    CC=10,
    seed = 123
  )

  auc.rand = cv(
    data = transform( data.cogdx, Class = sample(Class)),
    k = k,
    CC=10,
    seed = 123
  )


  auc=append(auc.cerad,auc.cogdx)
  auc=append(auc,auc.rand)
  names(auc)=c("cerad","cogdx","random")
  df=melt(do.call(rbind,auc))

  p.pheno.acc.boxplot=ggplot(df,aes(x=Var1,y=value)) +
    geom_boxplot(notch=FALSE)+
    labs(y="Balanced accuracy in predicting \n known AD phenotyes",x="Phenotypes")+
   theme_bw(base_size=12)+theme(legend.position="top")
  ggsave(p.pheno.acc.boxplot,filename="Figures/p.pheno.acc.boxplot.pdf", device="pdf",width=3,height=3,units="in")
