#!/bin/bash

# Parameters
CORES=6

# Slurm parameters
MEMO=50G                             # Memory required
TIME=23:00:00                       # Time required

# Assemble order prefix
ORDP="sbatch --mem="$MEMO" --nodes=1 --ntasks=1 --cpus-per-task="$CORES" --time="$TIME

# Create directory for log files
LOGS="logs"
mkdir -p $LOGS

mkdir -p "tmp"
mkdir -p "tmp/fitted_coefficients"
mkdir -p "tmp/cvfit"

OUT_DIR="selected_features"
mkdir -p $OUT_DIR

IND_LIST=(1 2 3)
DWT_LIST=(0 1)
GRP_LIST=(0)
FOLD_LIST=(1 2 3 4 5)

for IND in "${IND_LIST[@]}"; do
  for DWT in "${DWT_LIST[@]}"; do
    for GRP in "${GRP_LIST[@]}"; do
      for FOLD in "${FOLD_LIST[@]}"; do

        if [[ $IND == 1 ]]; then
          OUT_FILE=$OUT_DIR/"30class_dwt"$DWT"_grp"$GRP"_lasso_fold"$FOLD".txt"
        elif [[ $IND == 2 ]]; then
          OUT_FILE=$OUT_DIR/"8class_dwt"$DWT"_grp"$GRP"_lasso_fold"$FOLD".txt"
        else
          OUT_FILE=$OUT_DIR/"2class_dwt"$DWT"_grp"$GRP"_lasso_fold"$FOLD".txt"
        fi
        
        if [[ -f $OUT_FILE ]]; then
          echo "Found "$OUT_FILE
          RUN=0
        else
          RUN=1
          #echo $OUT_FILE
        fi

        if [[ $RUN == 1 ]]; then
          # Script to be run
          SCRIPT="variable_selection.sh $IND $DWT $GRP $FOLD"
          # Define job name
          JOBN=$IND"_"$DWT"_"$GRP"_"$FOLD
          OUTF=$LOGS"/"$JOBN".out"
          ERRF=$LOGS"/"$JOBN".err"
          # Assemble slurm order for this job
          ORD=$ORDP" -J "$JOBN" -o "$OUTF" -e "$ERRF" "$SCRIPT
          # Print order
          echo $ORD
          # Submit order
          $ORD
          # Run command now
          #./$SCRIPT

        fi
        
      done
    done
  done
done
