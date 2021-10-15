#source code modified from: https://rpubs.com/markloessi/506713
#https://johanndejong.wordpress.com/2018/04/04/nested-cross-validation/
#https://stackoverflow.com/questions/48379502/generate-a-confusion-matrix-for-svm-in-e1071-for-cv-results

set.seed(123)

rm(list=ls())
source('~/work/scNET-devel/scripts/load_libraries.R')
#source('~/work/scNET-devel/scripts/read_data.R')
source('~/work/scNET-devel/scripts/functions_for_network_analysis.R')


##create a feature matrix using output (feat importance scores) from the AD-gene classifier
fw=read.table("Mic.AD.Feature_weights.mat", header=T)
fw=fw[order(fw$average),]
features.tfs=rownames(tail(fw,dim(fw)[1]*.1))



#read rosmap data
rosmap.ge=read.table("~/Desktop/finalAlzheimersDLPFCDataSetForSaniyaToUse.csv", sep=",", header=T)
rosmap.ge$entrezID=NULL
colnames(rosmap.ge)=gsub("X","",colnames(rosmap.ge))

features=rosmap.ge[rosmap.ge$geneName %in% features.tfs, ]
rownames(features)=features$geneName
features$geneName=NULL
features=as.data.frame(scale(t(features)))
#features=scale(features)

#phenotypes
rosmap.pheno=read.table("~/Desktop/AlzheimersDLPFCPhenotypesUpdatedForSaniya.csv", sep=",", header=T, )

#cerad
patients.cerad.pos=rosmap.pheno[rosmap.pheno$CERADScore==1,]$Patient
patients.cerad.neg=rosmap.pheno[rosmap.pheno$CERADScore==4 ,]$Patient

#add labels for cerad
cerad.train=features
cerad.train=features[rownames(cerad.train) %in% patients.cerad.pos |rownames(cerad.train) %in% patients.cerad.neg, ]
cerad.train$Class=ifelse(rownames(cerad.train) %in% patients.cerad.pos,"yes","no")


#cogdx
patients.cogdx.pos=rosmap.pheno[rosmap.pheno$cogdx==4 | rosmap.pheno$cogdx==5 | rosmap.pheno$cogdx==6,]$Patient
patients.cogdx.neg=rosmap.pheno[rosmap.pheno$cogdx==1 ,]$Patient # | rosmap.pheno$cogdx==2 | rosmap.pheno$cogdx==3,

#add labels for cogdx
cogdx.train=features
cogdx.train=features[rownames(cogdx.train) %in% patients.cogdx.pos |rownames(cogdx.train) %in% patients.cogdx.neg, ]
cogdx.train$Class=ifelse(rownames(cogdx.train) %in% patients.cogdx.pos,"yes","no")


#dcfdx Clinical cognitive diagnosis summary
patients.dcfdx.pos=rosmap.pheno[rosmap.pheno$dcfdx==1 | rosmap.pheno$dcfdx==2 | rosmap.pheno$dcfdx==3,]$Patient
patients.dcfdx.neg=rosmap.pheno[rosmap.pheno$dcfdx==4 | rosmap.pheno$dcfdx==5 |  rosmap.pheno$dcfdx==6  ,]$Patient  #| rosmap.pheno$dcfdx==2 | rosmap.pheno$dcfdx==3,

#add labels for dcfdx
dcfdx.train=features
dcfdx.train=features[rownames(dcfdx.train) %in% patients.dcfdx.pos |rownames(dcfdx.train) %in% patients.dcfdx.neg, ]
dcfdx.train$Class=ifelse(rownames(dcfdx.train) %in% patients.dcfdx.pos,"yes","no")




#braak
patients.braak.pos=rosmap.pheno[rosmap.pheno$Braak.Progression==0 | rosmap.pheno$Braak.Progression==1 | rosmap.pheno$Braak.Progression==2,]$Patient
patients.braak.neg=rosmap.pheno[rosmap.pheno$Braak.Progression==5 | rosmap.pheno$Braak.Progression==6 ,]$Patient

#add labels for braak
braak.train=features
braak.train=features[rownames(braak.train) %in% patients.braak.pos |rownames(braak.train) %in% patients.braak.neg, ]
braak.train$Class=ifelse(rownames(braak.train) %in% patients.braak.pos,"yes","no")




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
cv = function(data, k, CC)
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
mtry=1

for(j in 1:10)
{
  #auc.cerad = cv(data = cerad.train, k = 5,C=mtry)
  #name=paste("auc.cerad",j,sep="_")
  #assign(name,auc.cerad)

  #auc.cogdx = cv(data = cogdx.train,k = k,C=mtry)
  #name=paste("auc.cogdx",j,sep="_")
  #assign(name,auc.cogdx)

  auc.braak = cv(data = braak.train,k = k,C=mtry)
  name=paste("auc.braak",j,sep="_")
  assign(name,auc.braak)

  #auc.dcfdx = cv(data = dcfdx.train,k = k,C=mtry)
  #name=paste("auc.dcfdx",j,sep="_")
  #assign(name,auc.dcfdx)

  auc.random = cv(data = transform( braak.train, Class = sample(Class)),k = k,C=mtry)
  name=paste("auc.random",j,sep="_")
  assign(name,auc.random)

}

#phenotypes=c("auc.cerad","auc.cogdx","auc.braak","auc.random","auc.dcfdx")
phenotypes=c("auc.random","auc.dcfdx")

acc.tbl=data.frame(acc=NULL,pheno=NULL,cond=NULL)
for (i in 1:length(phenotypes))
{
    pattern=paste(phenotypes[i],"_*",sep="")
    list=ls(pattern=pattern)
    for(j in 1:length(list))
    {
      df=as.data.frame(unlist(get(list[j])))
      colnames(df)="acc"
      df$pheno=phenotypes[i]
      acc.tbl=rbind(acc.tbl,df)
    }
}


acc.tbl$pheno=gsub("auc.","",acc.tbl$pheno)
  p.pheno.acc.boxplot=ggplot(acc.tbl,aes(x=pheno,y=acc)) +
    geom_boxplot(notch=TRUE)+
    labs(y="Balanced accuracy in predicting \n AD phenotyes",x="Phenotypes")+
   theme_bw(base_size=12)+theme(legend.position="top")
  ggsave(p.pheno.acc.boxplot,filename="Figures/p.pheno.acc.boxplot.pdf", device="pdf",width=3,height=3,units="in")
