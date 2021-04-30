## Interpretable Classification of Bacterial Raman Spectra with Knockoff Wavelets

This repository contains Python and R code to reproduce the data analysis described in the accompanying paper:
```
Interpretable Classification of Bacterial Raman Spectra with Knockoff Wavelets
Charmaine Chia, Matteo Sesia, Chi-Sing Ho, Stefanie S. Jeffrey, Jennifer Dionne, Emmanuel Candès, and Roger Howe
arXiv preprint (2021)
```

Contents:
 - "generate_knockoffs.ipynb": wavelet feature construction and knockoff generation steps of the data analysis;
 - "variable_selection.R": feature selection with the knockoff filter, using wavelet features and logistic regression importance measures;
 - "prediction.R": fit a logistic regression model on the selected features to classify different bacterial species;
 - "prediction_pca.R": fit a logistic regression model using PCA instead of wavelet features;
 - "variable_selection_others.R": feature selection and prediction using wavelet features and other types of importance measures (Naive Bayes and Nearest neighbors);