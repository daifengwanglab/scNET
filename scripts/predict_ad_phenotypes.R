#source code modified from: https://rpubs.com/markloessi/506713
#https://johanndejong.wordpress.com/2018/04/04/nested-cross-validation/
#https://stackoverflow.com/questions/48379502/generate-a-confusion-matrix-for-svm-in-e1071-for-cv-results

set.seed(123)

#rm(list=ls())
source('~/work/scNET-devel/scripts/load_libraries.R')
#source('~/work/scNET-devel/scripts/read_data.R')
source('~/work/scNET-devel/scripts/functions_for_network_analysis.R')


##create a feature matrix using output (feat importance scores) from the AD-gene classifier

#fw=read.table("Mic.AD.Feature_weights.mat", header=T)
#fw=fw[order(fw$average),]
#features.tfs=rownames(tail(fw,dim(fw)[1]*.1))

##create a feature matrix using output (TG classification) from the AD-gene classifier

fw=read.table("~/work/scNET_manuscript/AD_MIT/supp_data/AD_gene_prediction.mic.txt", header=T)
features.tfs=rownames(fw[fw[,2]>0.60,])




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
patients.cerad.pos=rosmap.pheno[rosmap.pheno$CERADScore==1 ,]$Patient

#patients.cerad.neg=rosmap.pheno[rosmap.pheno$CERADScore==4 ,]$Patient
patients.cerad.neg=rosmap.pheno[rosmap.pheno$CERADScore==4 | rosmap.pheno$CERADScore==3 ,]$Patient


#add labels for cerad
cerad.train=features
cerad.train=features[rownames(cerad.train) %in% patients.cerad.pos |rownames(cerad.train) %in% patients.cerad.neg, ]
cerad.train$Class=ifelse(rownames(cerad.train) %in% patients.cerad.pos,"yes","no")
colnames(cerad.train)=gsub("-","_",colnames(cerad.train))

#cogdx
#patients.cogdx.pos=rosmap.pheno[rosmap.pheno$cogdx==4 | rosmap.pheno$cogdx==5 | rosmap.pheno$cogdx==6,]$Patient
patients.cogdx.pos=rosmap.pheno[rosmap.pheno$cogdx==4 | rosmap.pheno$cogdx==5 ,]$Patient
patients.cogdx.neg=rosmap.pheno[rosmap.pheno$cogdx==1 ,]$Patient # | rosmap.pheno$cogdx==2 | rosmap.pheno$cogdx==3,

#add labels for cogdx
cogdx.train=features
cogdx.train=features[rownames(cogdx.train) %in% patients.cogdx.pos |rownames(cogdx.train) %in% patients.cogdx.neg, ]
cogdx.train$Class=ifelse(rownames(cogdx.train) %in% patients.cogdx.pos,"yes","no")
colnames(cogdx.train)=gsub("-","_",colnames(cogdx.train))


#dcfdx Clinical cognitive diagnosis summary
#patients.dcfdx.pos=rosmap.pheno[rosmap.pheno$dcfdx==1 | rosmap.pheno$dcfdx==2 | rosmap.pheno$dcfdx==3,]$Patient
patients.dcfdx.pos=rosmap.pheno[rosmap.pheno$dcfdx==1,]$Patient
#patients.dcfdx.neg=rosmap.pheno[rosmap.pheno$dcfdx==4 | rosmap.pheno$dcfdx==5 |  rosmap.pheno$dcfdx==6  ,]$Patient  #| rosmap.pheno$dcfdx==2 | rosmap.pheno$dcfdx==3,
patients.dcfdx.neg=rosmap.pheno[rosmap.pheno$dcfdx==4 | rosmap.pheno$dcfdx==5  ,]$Patient

#add labels for dcfdx
dcfdx.train=features
dcfdx.train=features[rownames(dcfdx.train) %in% patients.dcfdx.pos |rownames(dcfdx.train) %in% patients.dcfdx.neg, ]
dcfdx.train$Class=ifelse(rownames(dcfdx.train) %in% patients.dcfdx.pos,"yes","no")
colnames(dcfdx.train)=gsub("-","_",colnames(dcfdx.train))




#braak
patients.braak.pos=rosmap.pheno[rosmap.pheno$Braak.Progression==0 | rosmap.pheno$Braak.Progression==1 | rosmap.pheno$Braak.Progression==2 ,]$Patient
patients.braak.neg=rosmap.pheno[rosmap.pheno$Braak.Progression==5 | rosmap.pheno$Braak.Progression==6 ,]$Patient

#add labels for braak
braak.train=features
braak.train=features[rownames(braak.train) %in% patients.braak.pos |rownames(braak.train) %in% patients.braak.neg, ]
braak.train$Class=ifelse(rownames(braak.train) %in% patients.braak.pos,"yes","no")
colnames(braak.train)=gsub("-","_",colnames(braak.train))


#random
patients.rand.pos=rosmap.pheno[sample(nrow(rosmap.pheno), 100), ]$Patient
tmp=rosmap.pheno[!(rosmap.pheno$Patient %in% patients.rand.pos),]
patients.rand.neg=tmp[sample(nrow(tmp),100),]$Patient

#add labels for braak
rand.train=features
rand.train=features[rownames(rand.train) %in% patients.rand.pos |rownames(rand.train) %in% patients.rand.neg, ]
rand.train$Class=ifelse(rownames(rand.train) %in% patients.rand.pos,"yes","no")
colnames(rand.train)=gsub("-","_",colnames(rand.train))



k=5
mtry=1

for(j in 1:10)
{
  auc.cerad = cv(data = cerad.train, k = 5,C=mtry)
  name=paste("auc.cerad",j,sep="_")
  assign(name,auc.cerad)

  auc.cogdx = cv(data = cogdx.train,k = k,C=mtry)
  name=paste("auc.cogdx",j,sep="_")
  assign(name,auc.cogdx)

  auc.braak = cv(data = braak.train,k = k,C=mtry)
  name=paste("auc.braak",j,sep="_")
  assign(name,auc.braak)

  auc.dcfdx = cv(data = dcfdx.train,k = k,C=mtry)
  name=paste("auc.dcfdx",j,sep="_")
  assign(name,auc.dcfdx)

  auc.random = cv(data =  rand.train,k = k,C=mtry)
  name=paste("auc.random",j,sep="_")
  assign(name,auc.random)

}

phenotypes=c("auc.cerad","auc.cogdx","auc.braak","auc.random","auc.dcfdx")
#phenotypes=c("auc.random","auc.braak")

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
  p.pheno.acc.boxplot=ggplot(acc.tbl,aes(x=reorder(pheno,-acc),y=acc, fill=reorder(pheno,-acc))) +
    geom_boxplot(notch=TRUE)+scale_fill_grey() +
    labs(y="Balanced accuracy in predicting \n AD phenotyes",x="Phenotypes")+
   theme_bw(base_size=12)+theme(legend.position="none")
ggsave(p.pheno.acc.boxplot,filename="Figures/p.pheno.acc.boxplot.pdf", device="pdf",width=3,height=3,units="in")
