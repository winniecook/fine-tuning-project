#!/bin/bash

# Custom CTC data preparation with dataset-specific chunk sizes
# Uses optimal chunk sizes from thesis: 35k, 55k, or 75k samples

set -e

usage() {
    echo "Usage: $0 --input_pod5 <file> --reference <file> --output_dir <dir> --dataset_name <name>"
    echo ""
    echo "Required arguments:"
    echo "  --input_pod5     Training POD5 file (e.g. C9_128_train.pod5)"
    echo "  --reference      Reference FASTA file"
    echo "  --output_dir     Output directory for CTC files"
    echo "  --dataset_name   Dataset name (determines chunk size)"
    echo ""
    echo "Supported datasets and chunk sizes:"
    echo "  FX_160, C9_128, C9_256  -> 35,000 samples"
    echo "  FX_640, C9_512          -> 55,000 samples" 
    echo "  C9_1024                 -> 75,000 samples"
    echo ""
    echo "Example:"
    echo "  $0 --input_pod5 data/C9_512_train.pod5 --reference ref/C9_512.fa --output_dir ctc_custom/ --dataset_name C9_512"
    exit 1
}

get_chunk_size() {
    case $1 in
        FX_160|C9_128|C9_256)
            echo "35000"
            ;;
        FX_640|C9_512)
            echo "55000"
            ;;
        C9_1024)
            echo "75000"
            ;;
        *)
            echo "Error: Unknown dataset $1"
            echo "Supported: FX_160, C9_128, C9_256, FX_640, C9_512, C9_1024"
            exit 1
            ;;
    esac
}

get_batch_size() {
    case $1 in
        FX_160|C9_128|C9_256)
            echo "16"
            ;;
        FX_640|C9_512)
            echo "8"
            ;;
        C9_1024)
            echo "4"
            ;;
        *)
            echo "16"
            ;;
    esac
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

# Get chunk and batch sizes for this dataset
CHUNK_SIZE=$(get_chunk_size "$DATASET_NAME")
BATCH_SIZE=$(get_batch_size "$DATASET_NAME")

echo "Custom CTC Data Preparation"
echo "==========================="
echo "Dataset: $DATASET_NAME"
echo "Input POD5: $INPUT_POD5"
echo "Reference: $REFERENCE"
echo "Output directory: $OUTPUT_DIR"
echo "Chunk size: $CHUNK_SIZE samples"
echo "Batch size: $BATCH_SIZE"
echo "Min accuracy: 0.8"
echo ""

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Create temporary directory for this dataset
TEMP_DIR="$OUTPUT_DIR/${DATASET_NAME}_temp"
mkdir -p "$TEMP_DIR"

# Copy POD5 file to temp directory
cp "$INPUT_POD5" "$TEMP_DIR/training_data.pod5"

echo "Generating CTC data with custom chunk size..."

# Run bonito basecaller to generate CTC data
bonito basecaller dna_r10.4.1_e8.2_400bps_hac@v4.3.0 \
    "$TEMP_DIR" \
    --device cuda \
    --reference "$REFERENCE" \
    --min-accuracy 0.8 \
    --save-ctc \
    --chunksize "$CHUNK_SIZE" \
    --batchsize "$BATCH_SIZE" \
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
chunks = np.load('$OUTPUT_DIR/chunks.npy')
refs = np.load('$OUTPUT_DIR/references.npy')
ref_lengths = np.load('$OUTPUT_DIR/reference_lengths.npy')
print(f'CTC dataset for $DATASET_NAME:')
print(f'  Chunks: {chunks.shape}')
print(f'  References: {refs.shape}')
print(f'  Reference lengths: {ref_lengths.shape}')
print(f'  Chunk size: {chunks.shape[1]} (expected: $CHUNK_SIZE)')
"
else
    echo "Error: No CTC files generated"
    echo "Check log file: $OUTPUT_DIR/${DATASET_NAME}_ctc_log.txt"
    exit 1
fi

# Clean up
rm -rf "$TEMP_DIR"

echo ""
echo "Custom CTC preparation complete for $DATASET_NAME"
echo "Files saved in: $OUTPUT_DIR"