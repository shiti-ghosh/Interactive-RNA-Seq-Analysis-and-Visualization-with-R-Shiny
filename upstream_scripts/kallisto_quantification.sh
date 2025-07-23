#!/bin/bash

# -----------------------------------------------------------------------------
# Script: kallisto_quant.sh
# Description: Loops through input FASTQ files and runs Kallisto quantification
# Author: Shiti Ghosh
# -----------------------------------------------------------------------------

# Step 1: Set input directory containing FASTQ files
# Please check your path
INPUT_DIR="/data/fastq_files"

# Step 2: Load the index file

INDEX="/reference/Mus_musculus.idx"


# Step 3: Output directory root
OUTPUT_ROOT="/data/output/kallisto_quant"
mkdir -p "$OUTPUT_ROOT"

# Step 4: Loop through each FASTQ file and run kallisto quant
for FILE in "$INPUT_DIR"/*; do
  # extract basename, sample name
  BASENAME=$(basename "$FILE")
  SAMPLE_NAME="${BASENAME%.*}"
  OUTPUT_DIR="${OUTPUT_ROOT}/${SAMPLE_NAME}"

  echo "Processing sample: $SAMPLE_NAME"

  mkdir -p "$OUTPUT_DIR"

  # Run Kallisto quantification
  kallisto quant -i "$INDEX" -l 200 -s 20 -o "$OUTPUT_DIR" --single "$FILE"

  # Rename output files for clarity
  mv "$OUTPUT_DIR/abundance.tsv" "$OUTPUT_DIR/${SAMPLE_NAME}.tsv"
  mv "$OUTPUT_DIR/abundance.h5" "$OUTPUT_DIR/${SAMPLE_NAME}.h5"
done

echo "All samples processed. Kallisto outputs saved in '$OUTPUT_ROOT/'"
