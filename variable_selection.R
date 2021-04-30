#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(glmnet)
library(doParallel)
library(tidyverse)
library(knockoff)

## Input parameters
labelset_ind = 3                ## 1:Y_species, 2:Y_empiric_treatments, 3:Y_MRSA_MSSA, 4:Y_gram_positive_negative, 5:Y_gram_positive, 6:Y_gram_negative
DWT_or_raw = 0                  ## 0:DWT, 1: raw
ungrouped_or_grouped_lasso = 0 # Use grouped lasso for multinomial regression
fold_id = 1                     ## 1, 2, 3, 4, or 5

## Input parameters
labelset_ind = as.integer(args[1])
DWT_or_raw = as.integer(args[2])
ungrouped_or_grouped_lasso = as.integer(args[3]) # Use grouped lasso for multinomial regression
fold_id = as.integer(args[4])

Y_labelsets = c("Y_species", "Y_empiric_treatments", "Y_MRSA_MSSA", "Y_gram_positive_negative", "Y_gram_positive", "Y_gram_negative")
classes = c()              ## INPUT SPECIES CLASS LABELS TO INCLUDE. If selecting all data, set as *c()*

#### Register parallel
num.cores <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK")) - 1
cat(sprintf("num cores = %d\n", num.cores))
library(doParallel)
registerDoParallel(5)

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
datasetY <- datasetY[nnan_idx]

## Determine model and file names
nclasses = length(unique(datasetY))
out_filename = sprintf("%dclass_dwt%d_grp%d", nclasses, DWT_or_raw, ungrouped_or_grouped_lasso)
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

## perform lasso regression on X_Xk
parallel = TRUE
## model training ----------------------------------------------------------
cat(sprintf("Training model on %d x %d data...\n", nrow(X_train), ncol(X_train)))
if (binomial_or_multinomial == 0) {
    cvfit = cv.glmnet(X_train, Y_train, family = "binomial", type.measure = "class", parallel=parallel)
} else {
    if (ungrouped_or_grouped_lasso == 1) {
        cvfit = cv.glmnet(X_train, Y_train, family = "multinomial", type.multinomial = "grouped", parallel=parallel)
    } else {
        cvfit = cv.glmnet(X_train, Y_train, family = "multinomial", type.multinomial = "ungrouped", parallel=parallel)
    }
}
cat(sprintf("Model training completed.\n"))

lambda_min <- cvfit[["lambda.min"]]
coeffs <- coef(cvfit, s = "lambda.min")
mse_min <- cvfit[["cvm"]][which(cvfit[["lambda"]] == lambda_min)]
num_nonzero <- as.numeric(cvfit[["nzero"]][which(cvfit[["lambda"]] == lambda_min)])

## save model & coefficients
save(cvfit, coeffs, file = paste("./tmp/cvfit/", out_filename, ".RData", sep = "")) ## save model
if (binomial_or_multinomial==0) {
    coeffs_matrix <- as.matrix(coeffs)
} else {
    coeffs_matrix = c()
    for (i in seq_along(coeffs)) {
        c = as.matrix(coeffs[[i]])
        coeffs_matrix <- cbind(coeffs_matrix, c)
    }
}
write.csv(coeffs_matrix, file = paste("./tmp/fitted_coefficients/", out_filename, "_coeffs.csv", sep=""))

## model testing ---------------------------------------------------------
Ypred_train <- as.numeric(predict(cvfit, X_train, s = "lambda.min", type = "class"))
trainerror <- length(which(Y_train != Ypred_train)) / length(Y_train)
Ypred_test <- as.numeric(predict(cvfit, X_test, s = "lambda.min", type = "class"))
testerror <- length(which(Y_test != Ypred_test)) / length(Y_test)

cat(sprintf("%s %s %s %s %s\n", trainerror, testerror, mse_min, lambda_min, num_nonzero))


##############################################################################################

## load cvfit model
load(paste("./tmp/cvfit/", out_filename, ".RData", sep = ""))

## Get beta statistic & apply knockoff filter ------------------------------------------------------
coeffs <- coef(cvfit, s = "lambda.min")
nvar = dim(datasetXXk)[2]/2
##coeff_entry_pt = matrix(, nrow = nclasses, ncol = nvar)
std_devs <- apply(datasetX, 2, sd)

if (nclasses == 2) {
    beta = c()
    beta <- coeffs[2:dim(coeffs)[1],]
    beta_std <- beta * std_devs ## standardize coefficients so they can be compared against each other
    W = abs(beta_std[1:nvar]) - abs(beta_std[(nvar+1):(2*nvar)])
} else {
    beta <- abs(coeffs[[1]][2:(2*nvar+1)])
    for(b in 2:length(coeffs)) {
        beta <- beta + abs(coeffs[[b]][2:(2*nvar+1)])
    }
    beta_std <- beta * std_devs ## standardize coefficients so they can be compared against each other
    W = abs(beta_std[1:nvar]) - abs(beta_std[(nvar+1):(2*nvar)])
}

fdr_list = c(0.01, 0.02, seq(0.05, 1, 0.05))
## Apply knockoff filter ---------------------------------------------------
for(fdr_thresh in fdr_list) {
    tau = knockoff.threshold(W, fdr=fdr_thresh, offset=0)
    all_selected_features = which(W>=tau)
    cat(sprintf("FDR level %.2f: selected %d variables.\n", fdr_thresh, length(all_selected_features)))
    ## Save list of selected features
    out.file <- sprintf("selected_features/%s_knockoffs_%s_fold%d.txt", out_filename, fdr_thresh, fold_id)
    tibble(Variable = all_selected_features) %>% write_delim(out.file, delim=" ")
}

######################################
## Feature selection with the lasso ##
######################################

if(FALSE) {
    ## perform lasso regression on X_Xk
    ## model training ----------------------------------------------------------
    cat(sprintf("Training model on %d x %d data...\n", length(train_ind), ncol(datasetX)))
    if (binomial_or_multinomial == 0) {
        cvfit = cv.glmnet(datasetX[train_ind,], Y_train, family = "binomial", type.measure = "class", parallel=parallel)
    } else {
        if (ungrouped_or_grouped_lasso == 1) {
            cvfit = cv.glmnet(datasetX[train_ind,], Y_train, family = "multinomial", type.multinomial = "grouped", parallel=parallel)
        } else {
            cvfit = cv.glmnet(datasetX[train_ind,], Y_train, family = "multinomial", type.multinomial = "ungrouped", parallel=parallel)
        }
    }
    cat(sprintf("Model training completed.\n"))

    selected_lasso = which(beta!=0)
    cat(sprintf("The lasso selected %d variables.\n", length(selected_lasso)))
    ## Save list of selected features
    out.file <- sprintf("selected_features/%s_lasso_fold%d.txt", out_filename, fold_id)
    tibble(Variable = selected_lasso) %>% write_delim(out.file, delim=" ")
}
