#!/bin/bash

# Standard CTC data preparation with 10,000 sample chunks
# Processes training POD5 files for all datasets

set -e

usage() {
    echo "Usage: $0 --input_pod5 <file> --reference <file> --output_dir <dir> --dataset_name <name>"
    echo ""
    echo "Required arguments:"
    echo "  --input_pod5     Training POD5 file (e.g. C9_128_train.pod5)"
    echo "  --reference      Reference FASTA file"
    echo "  --output_dir     Output directory for CTC files"
    echo "  --dataset_name   Dataset name (e.g. C9_128, FX_160)"
    echo ""
    echo "Example:"
    echo "  $0 --input_pod5 data/C9_128_train.pod5 --reference ref/C9_128.fa --output_dir ctc_standard/ --dataset_name C9_128"
    exit 1
}

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --input_pod5)
            INPUT_POD5="$2"
            shift 2
            ;;
        --reference)
            REFERENCE="$2"
            shift 2
            ;;
        --output_dir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        --dataset_name)
            DATASET_NAME="$2"
            shift 2
            ;;
        -h|--help)
            usage
            ;;
        *)
            echo "Unknown option $1"
            usage
            ;;
    esac
done

# Check required arguments
if [[ -z "$INPUT_POD5" || -z "$REFERENCE" || -z "$OUTPUT_DIR" || -z "$DATASET_NAME" ]]; then
    echo "Error: Missing required arguments"
    usage
fi

# Check files exist
if [[ ! -f "$INPUT_POD5" ]]; then
    echo "Error: Input POD5 file not found: $INPUT_POD5"
    exit 1
fi

if [[ ! -f "$REFERENCE" ]]; then
    echo "Error: Reference file not found: $REFERENCE"
    exit 1
fi

echo "Standard CTC Data Preparation"
echo "============================="
echo "Dataset: $DATASET_NAME"
echo "Input POD5: $INPUT_POD5"
echo "Reference: $REFERENCE"
echo "Output directory: $OUTPUT_DIR"
echo "Chunk size: 10,000 samples"
echo "Min accuracy: 0.8"
echo ""

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Create temporary directory for this dataset
TEMP_DIR="$OUTPUT_DIR/${DATASET_NAME}_temp"
mkdir -p "$TEMP_DIR"

# Copy POD5 file to temp directory
cp "$INPUT_POD5" "$TEMP_DIR/training_data.pod5"

echo "Generating CTC data..."

# Run bonito basecaller to generate CTC data
bonito basecaller dna_r10.4.1_e8.2_400bps_hac@v4.3.0 \
    "$TEMP_DIR" \
    --device cuda \
    --reference "$REFERENCE" \
    --min-accuracy 0.8 \
    --save-ctc \
    --chunksize 10000 \
    > "$OUTPUT_DIR/${DATASET_NAME}_ctc_log.txt" 2>&1

echo "CTC generation complete"

# Check if CTC files were generated
if ls "$TEMP_DIR"/*.npy 1> /dev/null 2>&1; then
    echo "Moving CTC files to output directory..."
    mv "$TEMP_DIR"/*.npy "$OUTPUT_DIR/"
    
    echo "CTC files generated:"
    ls -lh "$OUTPUT_DIR"/*.npy
    
    # Show dataset info
    python3 -c "
import numpy as np
import os
chunks = np.load('$OUTPUT_DIR/chunks.npy')
refs = np.load('$OUTPUT_DIR/references.npy')
ref_lengths = np.load('$OUTPUT_DIR/reference_lengths.npy')
print(f'CTC dataset for $DATASET_NAME:')
print(f'  Chunks: {chunks.shape}')
print(f'  References: {refs.shape}')
print(f'  Reference lengths: {ref_lengths.shape}')
print(f'  Chunk size: {chunks.shape[1]} (expected: 10,000)')
"
else
    echo "Error: No CTC files generated"
    echo "Check log file: $OUTPUT_DIR/${DATASET_NAME}_ctc_log.txt"
    exit 1
fi

# Clean up
rm -rf "$TEMP_DIR"

echo ""
echo "Standard CTC preparation complete for $DATASET_NAME"
echo "Files saved in: $OUTPUT_DIR"