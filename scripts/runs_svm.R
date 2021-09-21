#source code modified from: https://rpubs.com/markloessi/506713
#https://johanndejong.wordpress.com/2018/04/04/nested-cross-validation/


rm(list=ls())
source('~/work/scNET-devel/scripts/load_libraries.R')
source('~/work/scNET-devel/scripts/read_data.R')
source('~/work/scNET-devel/scripts/functions_for_network_analysis.R')


#gwas catalog

gwas=read.table("/Users/chiraggupta/work/scNET_manuscript/genome/genesets/GWAS_catalog/full",header=T,sep="\t",quote="")
alz.gwas=gwas[gwas$DISEASE %like% "Alzheimer",c("MAPPED_GENE")]

alz.gwas=unlist(strsplit(alz.gwas, ","))
alz.gwas=unlist(strsplit(alz.gwas, ";"))
alz.gwas=unlist(strsplit(alz.gwas, "-"))
alz.gwas=gsub(" ","",alz.gwas)
alz.gwas=gsub("\"","",alz.gwas)
alz.gwas=unique(alz.gwas)

#disease gene database
dis=read.table("/Users/chiraggupta/work/scNET_manuscript/genome/genesets/Disgenet/all_gene_disease_associations.tsv",header=T,sep="\t",quote="")
#select alz genes with DGA score > 0 (at least some evidence)
alz.dis=unique(dis[dis$diseaseName %like% "Alzheimer" & dis$score > 0,]$geneSymbol)


alz=unique(append(alz.gwas,alz.dis))


###identify negative labels based on Dis-Dis-assoc
source('~/work/scNET_manuscript/get_api.R')
Sys.setenv(DISGENET_API_KEY= disgenet_api_key)

#find DD assoc
neg = disease2disease_by_gene(disease  = "C0002395",database = "CURATED",ndiseases = 100)
qr = disgenet2r::extract(neg)
uncorrDis=qr[qr$jaccard_genes < 0.09,]$disease2_name
uncorrDis.genes=unique(dis[dis$diseaseName %in% uncorrDis & dis$score > 0 ,]$geneSymbol)

#other diseases with overlaps with AD
corrDis=qr[qr$jaccard_genes >0.8,]$disease2_name
corrDis.genes=unique(dis[dis$diseaseName %in% corrDis & dis$score >= 0 ,]$geneSymbol)

#combine alz genes with genes associated with correlated diseases to AD
alz=unique(append(alz,corrDis.genes))

train_and_validate = function( data, fold, C)
{
  fit = ksvm(
    Class ~ .,
    data = data[-fold,],
    kernel = "vanilladot",
    C = C,
    scale=TRUE,
    prob.model = TRUE
  )

  # Predict the fold
  yh = predict(fit, newdata = data[fold,])

  # Compare the predictions to the labels
  posneg <- split(yh[,1], data$Class[fold])

  # Return the AUC under the ROC
  roc.curve(posneg[[1]], posneg[[2]])$auc
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

  #label alz genes as positives in the network matrix
  mat$Class=ifelse(mat$TG %in% alz,1,0)
  mat.positive=mat[mat$Class == 1,]


  ##create negative labels as those that are not positives and also not related to mental disorders
  mat.negative=mat[mat$Class == 0 & (!mat$TG %in% uncorrDis.genes),]

  set.seed(100)
  #randomly select rows equal to the number of positives
  mat.tmp=sample_n(mat.negative,dim(mat.positive)[1])

  #feature matrix
  mat.feat=rbind(mat.positive,mat.tmp)

  # Set the random seed for reproducibility
  data=as.data.frame(mat.feat[,-1])
  y = data$Class
  X = data
  X$Class = NULL

  # Do the cross-validation for each C in CC
  auc = cv(
    data = data,
    k = 5,
    #CC = 2^seq(log2(.01), log2(10), length.out = 21),
    #CC = 2^seq(log2(0.01), log2(10), length.out = 20),
    #CC=2^seq(log2(1), log2(100), length.out = 5),
    CC=c(0.1,1,10,100),
    seed = 111
  )
name=paste(celltype,"auc",sep=".")
assign(name,auc)
}#end of ct for loop
