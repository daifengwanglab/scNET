#source code modified from: https://rpubs.com/markloessi/506713
#https://johanndejong.wordpress.com/2018/04/04/nested-cross-validation/
#https://stackoverflow.com/questions/48379502/generate-a-confusion-matrix-for-svm-in-e1071-for-cv-results


rm(list=ls())
source('~/work/scNET-devel/scripts/load_libraries.R')
source('~/work/scNET-devel/scripts/read_data.R')
source('~/work/scNET-devel/scripts/functions_for_network_analysis.R')


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
  imp=imp[,4]
  imp
}


##returns average balanced acc and feat, imp. scores for each fold
#uses default mtry param for now: C is neglected

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

  #mtry=c("7","10","20","50","100")
  for(j in 1:10)
  {
    auc = cv(data = data, k = 5,mtry=1)
    name=paste(celltype,"auc",sep=".")
    name=paste(name,j,sep="_")
    assign(name,auc)
  }

  f.imp=cv_feat_imp(data,5)
  name=paste(celltype,"imp",sep=".")

  assign(name,as.data.frame(f.imp))


}#end of ct for loop


#make boxplots

celltypes=c("AD.Mic","AD.Oli","AD.Ex","AD.In","Ctrl.Mic","Ctrl.Oli","Ctrl.Ex","Ctrl.In")
acc.tbl=data.frame(acc=NULL,ct=NULL,cond=NULL)
for (i in 1:length(celltypes))
{
    pattern=paste(celltypes[i],"auc_*",sep=".")
    list=ls(pattern=pattern)
    for(j in 1:length(list))
    {
      df=as.data.frame(unlist(get(list[j])))
      colnames(df)="acc"
      df$ct=celltypes[i]
      df$cond=ifelse(df$ct %like% "AD","AD","Ctrl")
      acc.tbl=rbind(acc.tbl,df)
    }
}
acc.tbl$ct=gsub(".auc","",acc.tbl$ct)
acc.tbl$ct=gsub("AD.","",gsub("Ctrl.","",acc.tbl$ct))

npgcolors=pal_npg("nrc", alpha = 1)(10)
p.adgenes.acc.boxplot=ggplot(acc.tbl,aes(x=ct,y=acc,fill=cond)) +
  geom_boxplot(notch=TRUE)+
  labs(y="Balanced accuracy in predicting \n known AD genes \n Repeated 5 fold CV",x="Cell type networks")+
 theme_bw(base_size=12)+theme(legend.position="top")+scale_fill_manual(values=c("AD"=npgcolors[1],"Ctrl"=npgcolors[2]))
#ggsave(p.adgenes.acc.boxplot,filename="Figures/p.adgenes.acc.boxplot.pdf", device="pdf",width=3,height=3,units="in")


#select mic feature scores
AD.Mic.imp$average=rowMeans(AD.Mic.imp)
AD.Mic.imp=AD.Mic.imp[order(AD.Mic.imp$average),]
#write.table(AD.Mic.imp,file="Mic.AD.Feature_weights.mat",row.names=T,col.names=T,sep="\t",quote=F)


#concatenate networks as integrated features
nets=ls(pattern="*\\.network")

list.ctrl=Filter(function(x) !any(grepl("AD", x)), nets)
list.AD=Filter(function(x) !any(grepl("Ctrl", x)), nets)

nets=list.ctrl
i=1
name=nets[i]
celltype=gsub(".network","",name)
net=as.data.frame(lapply(nets[i],get))
net=net[,c("TF","TG","abs_coef")]
net=distinct(net)
net$TF=paste(net$TF,celltype,sep="_")
integrated.net=net

for(i in 2:length(nets))
{
  name=nets[i]
  celltype=gsub(".network","",name)
  net=as.data.frame(lapply(nets[i],get))
  net=net[,c("TF","TG","abs_coef")]
  net=distinct(net)
  net$TF=paste(net$TF,celltype,sep="_")
#  net$TG=paste(net$TG,celltype,sep="_")
  integrated.net=rbind(integrated.net,net)
}

mat=tidyr::pivot_wider(integrated.net, names_from = TF, values_from = abs_coef)
mat[is.na(mat)]=0
mat=as.data.frame(mat)

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

#mtry=c("7","10","20","50","100")
for(i in 1:2)
{
  auc = cv(data = data, k = 5,mtry=1)
  name=paste("integrated","auc",sep=".")
  name=paste(name,i,sep="_")
  assign(name,auc)
}

#f.imp=cv_feat_imp(data,5)
#f.imp.integrated=f.imp
acc.tbl=data.frame(acc=NULL,ct=NULL,cond=NULL)

pattern=paste("integrated","auc_*",sep=".")
list=ls(pattern=pattern)
for(j in 1:length(list))
{
  df=as.data.frame(unlist(get(list[j])))
  colnames(df)="acc"
  acc.tbl=rbind(acc.tbl,df)
}
