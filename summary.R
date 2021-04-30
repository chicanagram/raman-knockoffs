library(tidyverse)
library(kableExtra)
library(gridExtra)

## Input parameters
PCA_YN = 0                      #0: no PCA, 1: PCA


nclasses <- 2

DWT_or_raw = 0                  ## 0:DWT, 1: raw
df.lasso <- do.call("rbind", (lapply(c(30,8), function(nclasses) {
    in_filename = sprintf("%dclass_dwt%d_grp1", nclasses, DWT_or_raw)
    do.call("rbind", (lapply(1:5, function(fold_id) {
        #in.file <- sprintf("test_errors/%s_fold%d_lasso.txt", in_filename, fold_id)    
        in.file <- sprintf("test_errors/%s_lasso_fold%d.txt", in_filename, fold_id)    
        read_delim(in.file, delim=" ", col_types=cols()) %>%
            mutate(Fold=fold_id)
    }))) %>%
        mutate(Classes=nclasses, Raw=FALSE)
})))

DWT_or_raw = 1                  ## 0:DWT, 1: raw
df.lasso.raw <- do.call("rbind", (lapply(c(30,8,2), function(nclasses) {
    in_filename = sprintf("%dclass_dwt%d_grp0", nclasses, DWT_or_raw)
    do.call("rbind", (lapply(1:5, function(fold_id) {
        in.file <- sprintf("test_errors/%s_fold%d_lasso.txt", in_filename, fold_id)    
        read_delim(in.file, delim=" ", col_types=cols()) %>%
            mutate(Fold=fold_id)
    }))) %>%
        mutate(Classes=nclasses, Raw=TRUE)
})))

DWT_or_raw = 0                  ## 0:DWT, 1: raw
df.knockoffs <- do.call("rbind", (lapply(c(30,8,2), function(nclasses) {
    in_filename = sprintf("%dclass_dwt%d_grp0", nclasses, DWT_or_raw)
    do.call("rbind", (lapply(1:5, function(fold_id) {
        in.file <- sprintf("test_errors/%s_lasso_knockoffs_fold%d.txt", in_filename, fold_id)    
        read_delim(in.file, delim=" ", col_types=cols()) %>%
            mutate(Fold=fold_id)
    }))) %>%
        mutate(Classes=nclasses, Raw=FALSE)
})))

DWT_or_raw = 1                  ## 0:DWT, 1: raw
df.knockoffs.raw <- do.call("rbind", (lapply(c(30), function(nclasses) {
    in_filename = sprintf("%dclass_dwt%d_grp0", nclasses, DWT_or_raw)
    do.call("rbind", (lapply(1:5, function(fold_id) {
        in.file <- sprintf("test_errors/%s_lasso_knockoffs_fold%d.txt", in_filename, fold_id)    
        read_delim(in.file, delim=" ", col_types=cols()) %>%
            mutate(Fold=fold_id)
    }))) %>%
        mutate(Classes=nclasses, Raw=TRUE)
})))

df <- rbind(df.lasso, df.knockoffs, df.lasso.raw, df.knockoffs.raw) %>%
    mutate(Method = sprintf("%s-%s", Method, Raw)) %>% select(-Raw)

tb <- df %>%
    filter(FDR %in% c(0.1, NA)) %>%
    gather(NvarIn, NvarOut, Test, key="key", value="value") %>%
    group_by(Classes, Method, key) %>%
    summarise(Mean = mean(value), Sd = sd(value)) %>%
    pivot_wider(names_from=c("key"), values_from=c("Mean", "Sd")) %>%
    mutate(Input = factor(Method, c("Lasso-TRUE", "Knockoffs-TRUE", "Lasso-FALSE", "Knockoffs-FALSE"),
                          c("$X$", "$\\{X_j\\}_{j \\in \\hat{S}}$", "$X'$", "$\\{X'_j\\}_{j \\in \\hat{S}'}$"))) %>%
    mutate(Classes = factor(Classes, c(30,8,2), c("30 classes", "8 classes", "2 classes"))) %>%
    arrange(Classes, Input) %>%
    mutate(`Input features` = sprintf("%.0f (%.0f)", Mean_NvarIn, Sd_NvarIn),
           `Nonzero coefficients` = sprintf("%.0f (%.0f)", Mean_NvarOut, Sd_NvarOut),
           `Test error (\\%)` = sprintf("%.1f (%.1f)", Mean_Test*100, Sd_Test*100)) %>%
    ungroup()
tb

tb %>%
    select(Input, `Input features`, `Nonzero coefficients`, `Test error (\\%)`) %>%
    kable("latex", booktabs=TRUE, align="c", escape=FALSE) %>%
    pack_rows(index = table(tb$Classes))

###################################
## Results with other FDR values ##
###################################

fdr.values <- unique(df$FDR, na.rm=TRUE)

tb <- df %>%
    mutate(FDR = round(FDR,2)) %>%
    filter(Method %in% c("Lasso-FALSE", "Knockoffs-FALSE"), FDR %in% c(NA, fdr.values)) %>%
    gather(NvarIn, NvarOut, Test, key="key", value="value") %>%
    group_by(Classes, Method, FDR, key) %>%
    summarise(Mean = mean(value), Sd = sd(value)) %>%
    pivot_wider(names_from=c("key"), values_from=c("Mean", "Sd")) %>%
    mutate(Input = factor(Method, c("Lasso-FALSE", "Knockoffs-FALSE"),
                          c("$X'$", "$\\{X'_j\\}_{j \\in \\hat{S}'}$"))) %>%
    mutate(Classes = factor(Classes, c(30,8,2), c("30 classes", "8 classes", "2 classes"))) %>%
    arrange(Classes, Input, FDR) %>%
    mutate(`Input features` = sprintf("%.0f (%.0f)", Mean_NvarIn, Sd_NvarIn),
           `Nonzero coefficients` = sprintf("%.0f (%.0f)", Mean_NvarOut, Sd_NvarOut),
           `Test error (\\%)` = sprintf("%.1f (%.1f)", Mean_Test*100, Sd_Test*100)) %>%
    ungroup() %>%
    mutate(FDR = ifelse(is.na(FDR), "", FDR))    
tb

tb %>%
    select(Input, FDR, `Input features`, `Nonzero coefficients`, `Test error (\\%)`) %>%
    kable("latex", booktabs=TRUE, align="c", escape=FALSE) %>%
    pack_rows(index = table(tb$Classes))

tb.lasso <- tb %>% filter(FDR == "")

pp1 <- tb %>%
    filter(FDR != "") %>%
    mutate(FDR = parse_number(FDR)) %>%
    ggplot(aes(x=FDR, y=`Mean_Test`*100)) +
    geom_point() +
    geom_line() +
    geom_hline(data=tb.lasso, aes(yintercept=`Mean_Test`*100), linetype=2, color="blue") +
    coord_cartesian(ylim=c(4,12)) +
    xlab("FDR level") +
    ylab("Test error (%)") +
    facet_grid(.~Classes, scales="free_y") +
    theme_bw()
pp2 <- tb %>%
    filter(FDR != "") %>%
    mutate(FDR = parse_number(FDR)) %>%
    ggplot(aes(x=FDR, y=`Mean_NvarIn`)) +
    geom_point() +
    geom_line() +    
    geom_hline(data=tb.lasso, aes(yintercept=`Mean_NvarIn`), linetype=2, color="blue") +
    ylim(0, 1200) +
    xlab("FDR level") +
    ylab("Number of features") +
    facet_grid(.~Classes, scales="free_y") +
    theme_bw()
pp <- grid.arrange(pp1, pp2, ncol=1)
pp %>% ggsave(file=sprintf("../../manuscript/plot_fdr.png"), height=4, width=5, units = "in")


###################################
## Results with other algorithms ##
###################################

DWT_or_raw = 0                  ## 0:DWT, 1: raw
df.others <- do.call("rbind", (lapply(c(30,8,2), function(nclasses) {
    in_filename = sprintf("%dclass_dwt%d_others", nclasses, DWT_or_raw)
    do.call("rbind", (lapply(1:5, function(fold_id) {
        in.file <- sprintf("test_errors/%s_fold%d.txt", in_filename, fold_id)
        if(file.exists(in.file)) {
            tmp = read_delim(in.file, delim=" ", col_types=cols()) %>%
                mutate(Fold=fold_id)
        } else {
            tmp = tibble()
        }
        return(tmp)
    }))) %>%
        mutate(Classes=nclasses, Raw=FALSE)
}))) %>%
    mutate(Input = ifelse(is.na(FDR), "na", sprintf("%s",FDR)))
    

tb <- df.others %>%
    filter(FDR %in% c(0.1, NA)) %>%
    gather(NvarIn, NvarOut, Test, key="key", value="value") %>%
    group_by(Classes, Method, Input, key) %>%
    summarise(Mean = mean(value), Sd = sd(value)) %>%
    pivot_wider(names_from=c("key"), values_from=c("Mean", "Sd")) %>%
    mutate(Input = factor(Input, c("na", 0.1), c("$X'$", "$\\{X'_j\\}_{j \\in \\hat{S}'}$"))) %>%
    mutate(Classes = factor(Classes, c(30,8,2), c("30 classes", "8 classes", "2 classes"))) %>%
    arrange(Classes, Method, Input) %>%
    mutate(`Input features` = sprintf("%.0f (%.0f)", Mean_NvarIn, Sd_NvarIn),
           `Nonzero coefficients` = sprintf("%.0f (%.0f)", Mean_NvarOut, Sd_NvarOut),
           `Test error (\\%)` = sprintf("%.1f (%.1f)", Mean_Test*100, Sd_Test*100)) %>%
    ungroup()
tb

tb %>%
    select(Method, Input, `Input features`, `Test error (\\%)`) %>%
    kable("latex", booktabs=TRUE, align="c", escape=FALSE) %>%
    pack_rows(index = table(tb$Classes))

###################################
## Results with PCA ##
###################################

DWT_or_raw = 0                  ## 0:DWT, 1: raw
df.pca <- do.call("rbind", (lapply(c(30,8,2), function(nclasses) {
    in_filename = sprintf("%dclass_dwt%d_grp0", nclasses, DWT_or_raw)
    do.call("rbind", (lapply(1:5, function(fold_id) {
        in.file <- sprintf("test_errors/%s_fold%d_lasso_pca.txt", in_filename, fold_id)
        if(file.exists(in.file)) {
            tmp = read_delim(in.file, delim=" ", col_types=cols()) %>%
                mutate(Fold=fold_id)
        } else {
            tmp = tibble()
        }
        return(tmp)
    }))) %>%
        mutate(Classes=nclasses, Raw=FALSE)
}))) %>% mutate(Method = "PCA")

tb <- rbind(df.pca,df.knockoffs) %>%
    filter(FDR %in% c(0.1, NA)) %>%
    gather(NvarIn, NvarOut, Test, key="key", value="value") %>%
    group_by(Classes, Method, key) %>%
    summarise(Mean = mean(value), Sd = sd(value)) %>%
    pivot_wider(names_from=c("key"), values_from=c("Mean", "Sd")) %>%
    mutate(Input = factor(Method, c("PCA", "Knockoffs"), c("$X_{\\text{PCA}}$", "$\\{X'_j\\}_{j \\in \\hat{S}'}$"))) %>%
    mutate(Classes = factor(Classes, c(30,8,2), c("30 classes", "8 classes", "2 classes"))) %>%
    arrange(Classes, Input) %>%
    mutate(`Input features` = sprintf("%.0f (%.0f)", Mean_NvarIn, Sd_NvarIn),
           `Nonzero coefficients` = sprintf("%.0f (%.0f)", Mean_NvarOut, Sd_NvarOut),
           `Test error (\\%)` = sprintf("%.1f (%.1f)", Mean_Test*100, Sd_Test*100)) %>%
    ungroup()
tb

tb %>%
    select(Input, `Input features`, `Nonzero coefficients`, `Test error (\\%)`) %>%
    kable("latex", booktabs=TRUE, align="c", escape=FALSE) %>%
    pack_rows(index = table(tb$Classes))
