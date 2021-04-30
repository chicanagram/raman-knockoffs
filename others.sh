#!/bin/bash

IND=$1
DWT=$2
FOLD=$3

ml gcc/8.3.0
module load openblas/0.3.8
ml r/4.0.0

export OMP_NUM_THREADS=1
Rscript --vanilla variable_selection_others.R $IND $DWT $FOLD
