#!/bin/bash

IND=$1
DWT=$2
GRP=$3
FOLD=$4

ml gcc/8.3.0
module load openblas/0.3.8
ml r/4.0.0

export OMP_NUM_THREADS=1
Rscript --vanilla prediction.R $IND $DWT $GRP $FOLD
