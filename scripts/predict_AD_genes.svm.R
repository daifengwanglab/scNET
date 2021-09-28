#source code modified from: https://rpubs.com/markloessi/506713
#https://johanndejong.wordpress.com/2018/04/04/nested-cross-validation/
#https://stackoverflow.com/questions/48379502/generate-a-confusion-matrix-for-svm-in-e1071-for-cv-results


rm(list=ls())
source('~/work/scNET-devel/scripts/load_libraries.R')
source('~/work/scNET-devel/scripts/read_data.R')
source('~/work/scNET-devel/scripts/functions_for_network_analysis.R')

library(randomForest)

#disease gene database
dis=read.table("/Users/chiraggupta/work/scNET_manuscript/genome/genesets/Disgenet/all_gene_disease_associations.tsv",header=T,sep="\t",quote="")
#select alz genes with DGA score > 0 (at least some evidence)
alz.dis=unique(dis[dis$diseaseName %like% "Alzheimer" & dis$score >= 0.1,]$geneSymbol)




###identify negative labels based on Dis-Dis-assoc
source('~/work/scNET_manuscript/get_api.R')
Sys.setenv(DISGENET_API_KEY= disgenet_api_key)

#find Dis-Dis assoc
neg = disease2disease_by_gene(disease  = "C0002395",database = "CURATED",ndiseases = 100)
qr = disgenet2r::extract(neg)
uncorrDis=qr[qr$jaccard_genes < 0.09,]$disease2_name
uncorrDis.genes=unique(dis[dis$diseaseName %in% uncorrDis ,]$geneSymbol)


alz=alz.dis

#function to get feature scores
get_weights = function( data, fold)
{
  fit = randomForest(formula= as.factor(Class) ~ ., data = data[-fold,],importance=TRUE)
  imp=round(importance(fit), 3)
  imp=imp[,3]
  imp
}


train_and_validate = function( data, fold, C) #returns average balanced acc and feat, imp. scores for each fold
{

  #fit = svm(as.factor(Class) ~ ., data = as.matrix(data[-fold,]), kernel = "linear", cost = C, scale = FALSE,type="C-classification")

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

cv_feat_imp = function(data, k, seed = 123)
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

#for each ct network

for(i in 1:length(nets))
{
  name=nets[i]
  celltype=gsub(".network","",name)
  net=as.data.frame(lapply(nets[i],get))
  net=net[,c("TF","TG","abs_coef")]
  net=distinct(net)
  mat=tidyr::pivot_wider(net, names_from = TF, values_from = abs_coef)
  mat[is.na(mat)]=0
  mat=as.data.frame(mat)
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

  #feature matrix
  mat.feat=rbind(mat.positive,mat.tmp)

  # Set the random seed for reproducibility
  data=mat.feat[,-1]
  colnames(data)=gsub("-","_",colnames(data))
#  y = as.factor(data$Class)
#  X = data
#  X$Class = NULL

  # Do the cross-validation for each C in CC
  auc = cv(
    data = data,
    k = 5,
    CC=10,
    seed = 123
  )
  name=paste(celltype,"auc",sep=".")

  assign(name,auc)

  f.imp=cv_feat_imp(data,5)
  name=paste(celltype,"imp",sep=".")

  assign(name,as.data.frame(f.imp))


}#end of ct for loop


#make boxplots

celltypes=c("Mic.AD","Oli.AD","Ex.AD","In.AD","Mic.Ctrl","Oli.Ctrl","Ex.Ctrl","In.Ctrl")
list=ls(pattern="*.auc")
acc.tbl=data.frame(acc=NULL,ct=NULL,cond=NULL)
for (i in 1:length(list))
{
    df=as.data.frame(unlist(get(list[i])))
    colnames(df)="acc"
    df$ct=list[i]
    df$cond=ifelse(df$ct %like% "AD","AD","Ctrl")
    acc.tbl=rbind(acc.tbl,df)
}
acc.tbl$ct=gsub(".auc","",acc.tbl$ct)
acc.tbl$ct=gsub("AD.","",gsub("Ctrl.","",acc.tbl$ct))

p.adgenes.acc.boxplot=ggplot(acc.tbl,aes(x=ct,y=acc,fill=cond)) +
  geom_boxplot(notch=FALSE)+
  labs(y="Balanced accuracy in predicting \n known AD genes",x="Cell type networks")+
 theme_bw(base_size=12)+theme(legend.position="top")
#ggsave(p.adgenes.acc.boxplot,filename="Figures/p.adgenes.acc.boxplot.pdf", device="pdf",width=3,height=3,units="in")


#select mic feature scores
AD.Mic.imp$average=rowMeans(AD.Mic.imp)
AD.Mic.imp=AD.Mic.imp[order(AD.Mic.imp$average),]
#write.table(AD.Mic.imp,file="Mic.AD.Feature_weights.mat",row.names=T,col.names=T,sep="\t",quote=F)
