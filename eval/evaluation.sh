#!/bin/bash

# Script Name: evaluation.sh
# Description: This script performs evaluation for schemes (CS*) and attributes (CA*).
# Author: Haomiao Luo
# Date: 2025-03-21

# Usage: ./evaluation.sh

# Specify the directory where the script is located
EVALUATION_SCHEMES_DIR="evaluation_schemes.py"
SUMMARIZATION_DIR="summarization.py"
LOG_FILE="$SCRIPT_DIR/evaluation.log"

# Initialize the scenarios 
SCHEME_LIST=("v1" "v2" "v3")
ATTRIBUTES_LIST=("Mechanism" "Cell Line" "Cell Type" "Tissue Type" "Organism")
READER_LIST=("REACH" "INDRA" "GPT" "LLAMA")
FILTERED_READER_LIST=("FLUTE/REACH" "FLUTE/INDRA" "FLUTE/GPT" "FLUTE/LLAMA")


# Create/clear log file 
> '$LOG_FILE'

# Logging function
log() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $*" | tee -a "$LOG_FILE"
}

log "Evaluation script started."

# Check if Python environment is active
if [ -z "${VIRTUAL_ENV}" ]; then
    log "ERROR: No virtual environment is activated."
    source ../../violin_test_20250321/bin/activate

else
    log "Environment activated: ${VIRTUAL_ENV}"
fi


# for s in "${SCHEME_LIST[@]}"
# do 
#     for r in "${READER_LIST[@]}"
#     do
#         python "$EVALUATION_SCHEMES_DIR" \
#         --scheme $s \
#         --reader FLUTE/$r \
#         --output ./output \
#         --attributes 'Regulator Compartment ID' 'Regulated Compartment ID' \
#         --log_file true

#         python "$SUMMARIZATION_DIR" \
#         --result_dir $s \
#         --reader_name FLUTE/$r/$s \
#         --output_dir ./FLUTE/$r/$s \
#         --scheme $s
#     done
# done

# for ca in "${ATTRIBUTES_LIST[@]}"
# do 
#     for r in "${READER_LIST[@]}"
#     do
#         python "$EVALUATION_SCHEMES_DIR" \
#         --scheme v1 \
#         --reader FLUTE/$r \
#         --output ./output \
#         --attributes 'Regulator Compartment ID' 'Regulated Compartment ID' $ca \
#         --log_file true

#         python "$SUMMARIZATION_DIR" \
#         --result_dir v1 \
#         --reader_name FLUTE/$r/v1 \
#         --output_dir ./FLUTE/$r/v1 \
#         --scheme v1
#     done

# Run the tests for evaluating data filtering (FLUTE tabs)
for s in "${SCHEME_LIST[@]}"
do
    for r in "${FILTERED_READER_LIST[@]}"
    do
        python "$EVALUATION_SCHEMES_DIR" \
        --scheme $s \
        --reader $r \
        --output ./output \
        --attributes 'Regulator Compartment ID' 'Regulated Compartment ID' \
        --log_file true
        
        python "$SUMMARIZATION_DIR" \
        --result_dir ./output/$r/$s \
        --reader_name $r/$s \
        --out_dir $r/$s \
        --scheme $s
    done
done

