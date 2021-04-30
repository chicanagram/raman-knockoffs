#!/bin/bash

# Parameters
CORES=11

# Slurm parameters
MEMO=50G                             # Memory required
TIME=18:00:00                       # Time required

# Assemble order prefix
ORDP="sbatch --mem="$MEMO" --nodes=1 --ntasks=1 --cpus-per-task="$CORES" --time="$TIME

# Create directory for log files
LOGS="logs"
mkdir -p $LOGS

OUT_DIR="test_errors"
mkdir -p $OUT_DIR

IND_LIST=(1 2 3)
DWT_LIST=(0)
GRP_LIST=(0)
FOLD_LIST=(1 2 3 4 5)

for IND in "${IND_LIST[@]}"; do
  for DWT in "${DWT_LIST[@]}"; do
    for GRP in "${GRP_LIST[@]}"; do
    for FOLD in "${FOLD_LIST[@]}"; do

      RUN=1

      if [[ $RUN == 1 ]]; then
        # Script to be run
        SCRIPT="prediction_pca.sh $IND $DWT $GRP $FOLD"
        # Define job name
        JOBN="pca_"$IND"_"$DWT"_"$GRP"_"$FOLD
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
