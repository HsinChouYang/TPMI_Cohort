#!/bin/bash
# This script runs admixture analysis for different K values with cross-validation.
# Usage: ./run_admixture.sh

# Define variables
INPUT_FILE="ForAdmixture.bed"
THREADS=10
START_K=5
END_K=18

# Loop through each K value from START_K to END_K and run admixture
for (( K=$START_K; K<=$END_K; K++ ))
do
    nohup ./admixture --cv ${INPUT_FILE} ${K} -j${THREADS} > K${K}.log &
done


