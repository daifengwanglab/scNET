#source code modified from: https://rpubs.com/markloessi/506713
#https://johanndejong.wordpress.com/2018/04/04/nested-cross-validation/


rm(list=ls())
source('~/work/scNET-devel/scripts/load_libraries.R')
source('~/work/scNET-devel/scripts/read_data.R')
source('~/work/scNET-devel/scripts/functions_for_network_analysis.R')


mat=tidyr::pivot_wider(net, names_from = TF, values_from = abs_coef)
mat[is.na(mat)]=0

#disease gene database
dis=read.table("/Users/chiraggupta/work/scNET_manuscript/genome/genesets/Disgenet/all_gene_disease_associations.tsv",header=T,sep="\t",quote="")

#select alz genes with DGA score > 0 (at least some evidence)
alz=unique(dis[dis$diseaseName %like% "Alzheimer" & dis$score >= 0,]$geneSymbol)

#label alz genes as positives in the network matrix
mat$Class=ifelse(mat$TG %in% alz,1,0)
mat.positive=mat[mat$Class == 1,]

#identify mental disorders genes
mat.MESHF03=unique(dis[dis$diseaseClass %like% "F03" &dis$score >= 0 ,]$geneSymbol)

##create negative labels as those that are not positives and also not related to mental disorders
mat.negative=mat[mat$Class == 0 & (!mat$TG %in% mat.MESHF03),]

set.seed(100)

#randomly select rows equal to the number of positives
mat.tmp=sample_n(mat.negative,dim(mat.positive)[1])

#feature matrix
mat.feat=rbind(mat.positive,mat.tmp)


# Set the random seed for reproducibility
data=as.data.frame(mat.feat[,-1])
y <- data$Class
X <- data
#X$Class <- NULL


train_and_validate = function( data, fold, C)
{
  fit = ksvm(
    Class ~ .,
    data = data[-fold,],
    kernel = "vanilladot",
    C = C,
    scale=FALSE,
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

# Do the cross-validation for each C in CC
auc = cv(
  data = data,
  k = 5,
  #CC = 2^seq(log2(.01), log2(10), length.out = 21),
CC = 2^seq(log2(.1), log2(10), length.out = 10),
  seed = 111
)






#####
train_and_validate <- function(data,fold,C)
{
  # Train an SVM, excluding the fold
  fit = ksvm(Class ~ .,
    data = data[-fold,],
    kernel = "vanilladot",
    kpar = list(),
    C = 1,
    prob.model = TRUE,
    Class.weights = 1 / table(data$Class[-fold])
    )
  # Predict the fold
  yh = predict(fit, newdata = data[fold,])
  # Compare the predictions to the labels
  posneg = split(yh[,1], data$Class[fold])
  # Return the AUC under the ROC
  roc.curve(posneg[[1]], posneg[[2]])$auc
}

# Function for doing a k-fold cross-validation for each C in CC
cv = function(data,k,CC,seed = NULL)
{
  # For each value of the hyperparameter C ...
  auc = lapply(CC, function(C) {
    folds = createFolds(data$Class, k = k)
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
