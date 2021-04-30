#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(tidyverse)
library(knockoff)
library(rknn)
library(e1071)    
library(parallel)
cl <- makeCluster(10)

## Input parameters
labelset_ind = 3                ## 1:Y_species, 2:Y_empiric_treatments, 3:Y_MRSA_MSSA, 4:Y_gram_positive_negative, 5:Y_gram_positive, 6:Y_gram_negative
DWT_or_raw = 0                  ## 0:DWT, 1: raw
fold_id = 1                     ## 1, 2, 3, 4, or 5

## Input parameters
labelset_ind = as.integer(args[1])
DWT_or_raw = as.integer(args[2])
fold_id = as.integer(args[3])

ungrouped_or_grouped_lasso = 1 # Use grouped lasso for multinomial regression
Y_labelsets = c("Y_species", "Y_empiric_treatments", "Y_MRSA_MSSA", "Y_gram_positive_negative", "Y_gram_positive", "Y_gram_negative")
classes = c()              ## INPUT SPECIES CLASS LABELS TO INCLUDE. If selecting all data, set as *c()*

#### Register parallel
num.cores <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK")) - 1
cat(sprintf("num cores = %d\n", num.cores))
cl <- makeCluster(10)

## load X and Xknockoffs datasets
if (DWT_or_raw == 0) {
  X_dwt <- read_csv("../datasets/X_dwt_coif4lvl5_pywt.csv", col_names=FALSE, col_types=cols())
  Xk_dwt <- read_csv("../datasets/Xk-m_dwt_coif4lvl5_pywt.csv", col_names=FALSE, col_types=cols())
} else {
  X_raw <- read_csv("../datasets/X_raw.csv", col_names=FALSE, col_types=cols())
}
Y_all <- read_csv("../datasets/Y.csv", col_names = TRUE, col_types=cols())

## build datasets ----------------------------------------------------------
samples_per_class = 2000
if (length(classes)==0) {classes = seq(0,29)} ## select all data if no classes input
idx_list = c()
for (class in classes) {
  idx_start = class*samples_per_class + 1
  idx_end = (class+1)*samples_per_class
  idx_list <- append(idx_list, c(idx_start:idx_end))
}

if (DWT_or_raw == 0) {
    datasetX <- as.matrix(X_dwt[idx_list,])
    datasetXk <- as.matrix(Xk_dwt[idx_list,])
} else {
    datasetX <- as.matrix(X_raw[idx_list,])
}
datasetY <- as.numeric(unlist(Y_all[idx_list,Y_labelsets[labelset_ind]]))
## filter datasets to remove rows with NaN labels
nnan_idx = which(is.nan(datasetY) == 0)
datasetX <- datasetX[nnan_idx,]
datasetXk <- datasetXk[nnan_idx,]
colnames(datasetX) <- as.character(parse_number(colnames(datasetX)))
colnames(datasetXk) <- as.character(parse_number(colnames(datasetXk))+ncol(datasetXk))
datasetY <- datasetY[nnan_idx]

## Determine model and file names
nclasses = length(unique(datasetY))
out_filename = sprintf("%dclass_dwt%d", nclasses, DWT_or_raw)
if(nclasses==2) {
    binomial_or_multinomial = 0
} else {
    binomial_or_multinomial = 1
}

## concatenate X and Xk
datasetXXk <- cbind(datasetX, datasetXk)

## split into train/test ---------------------------------------------------
set.seed(123)
folds <- sample(cut(seq(1,nrow(datasetXXk)),breaks=5,labels=FALSE))
train_ind = which(folds!=fold_id)

## Split the data
X_train = datasetXXk[train_ind,]
X_test = datasetXXk[-train_ind,]
Y_train = datasetY[train_ind]
Y_test = datasetY[-train_ind]

## List of FDR levels
fdr_list = c(0.01, 0.02, seq(0.05, 1, 0.05))

##rows = sample(1:nrow(X_train), 1000)
rows = 1:nrow(X_train)

#################
## Naive Bayes ##
#################

cat(sprintf("Training Naive Bayes model on knockoff-augmented data...\n"))
nbc.fit = naiveBayes(scale(X_train)[rows,], factor(Y_train[rows]))
beta = sapply(nbc.fit$tables, function(tab) abs(tab[1,1]-tab[2,1]))
p = ncol(X_train)/2
W = beta[1:p] - beta[(p+1):(2*p)]

## Apply knockoff filter
err_list = c()
nvar_out_list = c()
nvar_in_list = c()
for(fdr_thresh in fdr_list) {
    tau = knockoff.threshold(W, fdr=fdr_thresh, offset=0)
    selected = which(W>=tau)
    cat(sprintf("FDR level %.2f: selected %d variables.\n", fdr_thresh, length(selected)))
    ## Save list of selected features
    out.file <- sprintf("selected_features/%s_nb_%s_fold%d.txt", out_filename, fdr_thresh, fold_id)
    tibble(Variable = selected) %>% write_delim(out.file, delim=" ")
    ## Estimate prediction error
    nbc.fit = naiveBayes(X_train[rows,selected], factor(Y_train)[rows])
    y.hat = predict(nbc.fit, X_test[,selected])
    err = mean(y.hat!=Y_test)
    err_list = c(err_list, err)
    nvar_out_list = c(nvar_out_list, length(selected))
    nvar_in_list = c(nvar_in_list, length(selected))
}
df.1 = tibble(Method="Naive Bayes", FDR=fdr_list, NvarIn=nvar_in_list, NvarOut=nvar_out_list, Test=err_list)

## Apply without knockoffs
cat(sprintf("Training Naive Bayes model without knockoffs...\n"))
nbc.fit = naiveBayes(X_train[rows,1:p], factor(Y_train[rows]))
y.hat = predict(nbc.fit, X_test[,selected])
err = mean(y.hat!=Y_test)
df.2 = tibble(Method="Naive Bayes", FDR=NA, NvarIn=p, NvarOut=p, Test=err)


#########################
## K-Nearest neighbors ##
#########################

cat(sprintf("Training KNN model on knockoff-augmented data...\n"))
knn.fit = rknnBeg(X_train[rows,], factor(Y_train)[rows], pk=0.5, stopat=10, k=1, cluster=cl, r=100)
beta.vals = knn.fit$vars[[which.max(knn.fit$mean_accuracy)]]
beta.idx = as.integer(names(beta.vals))
beta = rep(0, 2*p)
beta[beta.idx] = beta.vals
p = ncol(X_train)/2
W = beta[1:p] - beta[(p+1):(2*p)]

## Apply knockoff filter
err_list = c()
nvar_out_list = c()
nvar_in_list = c()
for(fdr_thresh in fdr_list) {
    tau = knockoff.threshold(W, fdr=fdr_thresh, offset=0)
    selected = which(W>=tau)
    cat(sprintf("FDR level %.2f: selected %d variables.\n", fdr_thresh, length(selected)))
    ## Save list of selected features
    out.file <- sprintf("selected_features/%s_knn_%s_fold%d.txt", out_filename, fdr_thresh, fold_id)
    tibble(Variable = selected) %>% write_delim(out.file, delim=" ")
    ## Estimate prediction error
    knn.fit = rknn(X_train[rows,selected], X_test[,selected], factor(Y_train)[rows])
    y.hat = knn.fit$pred
    err = mean(y.hat!=Y_test)
    err_list = c(err_list, err)
    nvar_out_list = c(nvar_out_list, length(selected))
    nvar_in_list = c(nvar_in_list, length(selected))
}
df.3 = tibble(Method="KNN", FDR=fdr_list, NvarIn=nvar_in_list, NvarOut=nvar_out_list, Test=err_list)

## Apply without knockoffs
cat(sprintf("Training KNN model without knockoffs...\n"))
knn.fit = rknn(X_train[rows,1:p], X_test[,1:p], factor(Y_train)[rows])
y.hat = knn.fit$pred
err = mean(y.hat!=Y_test)
df.4 = tibble(Method="KNN", FDR=NA, NvarIn=p, NvarOut=p, Test=err)

#####################
## Combine results ##
#####################
df = rbind(df.1, df.2, df.3, df.4)
out.file <- sprintf("test_errors/%s_others_fold%d.txt", out_filename, fold_id)
df %>% write_delim(out.file, delim=" ")
