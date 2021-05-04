## Interpretable Classification of Bacterial Raman Spectra with Knockoff Wavelets

This repository contains Python and R code to reproduce the data analysis described in the [accompanying paper](https://arxiv.org/abs/2006.04937):
```
Interpretable Classification of Bacterial Raman Spectra with Knockoff Wavelets
Charmaine Chia, Matteo Sesia, Chi-Sing Ho, Stefanie S. Jeffrey, Jennifer Dionne, Emmanuel Candès, and Roger Howe
arXiv preprint (2021) https://arxiv.org/abs/2006.04937
```

Contents:
 - "generate_knockoffs.ipynb": wavelet feature construction and knockoff generation steps of the data analysis;
 - "variable_selection.R": feature selection with the knockoff filter, using wavelet features and logistic regression importance measures;
 - "prediction.R": fit a logistic regression model on the selected features to classify different bacterial species;
 - "prediction_pca.R": fit a logistic regression model using PCA instead of wavelet features;
 - "variable_selection_others.R": feature selection and prediction using wavelet features and other types of importance measures (Naive Bayes and Nearest neighbors);

Due to the large size of these data sets, fully reproducing the analysis described in the paper would be easiest with a computing cluster. Bash scripts are included to run the analysis on a cluster with a Slurm scheduler. 
