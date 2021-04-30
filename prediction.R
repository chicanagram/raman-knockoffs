#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(glmnet)
library(doParallel)
library(tidyverse)

## Input parameters
labelset_ind = 2                ## 1:Y_species, 2:Y_empiric_treatments, 3:Y_MRSA_MSSA, 4:Y_gram_positive_negative, 5:Y_gram_positive, 6:Y_gram_negative
DWT_or_raw = 0                  ## 0:DWT, 1: raw
ungrouped_or_grouped_lasso = 0 # Use grouped lasso for multinomial regression
fold_id = 1                     ## 1, 2, 3, 4, or 5
PCA_YN = 0                      #0: no PCA, 1: PCA

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
registerDoParallel(10)

## load X dataset
if (DWT_or_raw == 0) {
  X_dwt <- read_csv("../datasets/X_dwt_coif4lvl5_pywt.csv", col_names=FALSE, col_types=cols())
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
} else {
    datasetX <- as.matrix(X_raw[idx_list,])
}
datasetY <- as.numeric(unlist(Y_all[idx_list,Y_labelsets[labelset_ind]]))
## filter datasets to remove rows with NaN labels
nnan_idx = which(is.nan(datasetY) == 0)
datasetX <- datasetX[nnan_idx,]
datasetY <- datasetY[nnan_idx]

# get PCA components of raw data ------------------------------------------
if (PCA_YN == 1) {
  res.pca <- prcomp(datasetX, scale = FALSE)
  #fviz_eig(res.pca)
  datasetX.ind <- get_pca_ind(res.pca)
  datasetX <- datasetX.ind[["coord"]][,1:20] # extract top 20 components
}

## Determine model and file names
nclasses = length(unique(datasetY))
out_filename = sprintf("%dclass_dwt%d_grp%d", nclasses, DWT_or_raw, ungrouped_or_grouped_lasso)
if(nclasses==2) {
    binomial_or_multinomial = 0
} else {
    binomial_or_multinomial = 1
}

## split into train/test ---------------------------------------------------
set.seed(123)
folds <- sample(cut(seq(1,nrow(datasetX)),breaks=5,labels=FALSE))
train_ind = which(folds!=fold_id)

## Split the data
X_train = datasetX[train_ind,]
X_test = datasetX[-train_ind,]
Y_train = datasetY[train_ind]
Y_test = datasetY[-train_ind]

prediction_error = function(selected_features) {
    ## Training indices
    idx <- sample(1:nrow(X_train), 1000)
    
    if(length(selected_features)>0) {
        ## perform lasso regression on X
        ## model training ----------------------------------------------------------
        cat(sprintf("Training model on %d x %d data...\n", nrow(X_train), length(selected_features)))
        dfmax = 1000
        parallel = TRUE
        if (binomial_or_multinomial == 0) {
            cvfit = cv.glmnet(X_train[,selected_features], Y_train, family = "binomial", type.measure = "class", parallel=parallel, dfmax=dfmax)
            nvar.out = sum(coef(cvfit, s="lambda.min")[-1]!=0)
        } else {
            if (ungrouped_or_grouped_lasso == 1) {
                cvfit = cv.glmnet(X_train[idx,selected_features], Y_train[idx], family="multinomial", type.multinomial="grouped",
                                  parallel=parallel, dfmax=dfmax)
            } else {
                cvfit = cv.glmnet(X_train[,selected_features], Y_train, family = "multinomial", type.multinomial = "ungrouped",
                                  parallel=parallel, dfmax=dfmax)
            }
            beta.list = coef(cvfit, s="lambda.min")
            nonzero = unique(do.call("c", lapply(beta.list, function(beta) which(beta[-1]!=0))))
            nvar.out = length(nonzero)
        }
        cat(sprintf("Model training completed.\n"))
        ## model testing ---------------------------------------------------------
        Ypred_test = as.numeric(predict(cvfit, X_test[,selected_features], s = "lambda.min", type = "class"))
    } else {
        Ypred_test = median(Y_train)
        nvar.out = 0
    }
    out <- c()
    out$err <- mean(Y_test != Ypred_test)
    out$nvar.out = nvar.out
    out$nvar.in = length(selected_features)
    return(out)
}

################################
## Predictions with knockoffs ##
################################
if(TRUE) {
##    fdr_list = c(0.01, 0.02, seq(0.05, 0.5, 0.05))
    fdr_list = c(0.1)
    err_list = c()
    nvar_out_list = c()
    nvar_in_list = c()
    for(fdr_thresh in fdr_list) {
        in.file = sprintf("selected_features/%s_fold%d_knockoffs_%s.txt", out_filename, fold_id, fdr_thresh)
        df = read_delim(in.file, delim=" ", col_types=cols())
        selected_features = df$Variable
        out = prediction_error(selected_features)
        cat(sprintf("FDR level %.2f: input %d variables, selected %d variables, error: %.3f.\n",
                    fdr_thresh, length(selected_features), out$nvar, out$err))
        err_list = c(err_list, out$err)
        nvar_out_list = c(nvar_out_list, out$nvar.out)
        nvar_in_list = c(nvar_in_list, out$nvar.in)
    }
    df.1 = tibble(Method="Knockoffs", FDR=fdr_list, NvarIn=nvar_in_list, NvarOut=nvar_out_list, Test=err_list)
} else {
    df.1 = tibble()
}

############################################
## Predictions without variable selection ##
############################################
if(FALSE) {
    selected_features = 1:ncol(X_test)
    res_lasso = prediction_error(selected_features)
    df.2 = tibble(Method="Lasso", FDR=NA, NvarIn=length(selected_features), NvarOut=res_lasso$nvar.out, Test=res_lasso$err)
    print(df.2)
} else {
    df.2 <- tibble()
}

#####################
## Combine results ##
#####################
df = rbind(df.1, df.2)
out.file <- sprintf("test_errors/%s_lasso_fold%d.txt", out_filename, fold_id)
df %>% write_delim(out.file, delim=" ")
