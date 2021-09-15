#source code modified from: https://rpubs.com/markloessi/506713
#https://johanndejong.wordpress.com/2018/04/04/nested-cross-validation/
rm(list=ls())
library(kernlab)
library(caret)
library(mlbench)
library(PRROC)



# Set the random seed for reproducibility
set.seed(123)

y <- data$Class
X <- data
X$Class <- NULL

train_and_validate <- function(data,fold,C)
{
  # Train an SVM, excluding the fold
  fit = randomForest(Class ~ .,data = data[-fold,],ntrees = C)
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
