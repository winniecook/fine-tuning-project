#!/bin/bash

# Train model with standard 10K chunks
# 5 epochs, learning rate 5e-4, batch size based on dataset

set -e

usage() {
    echo "Usage: $0 --ctc_dir <dir> --output_model <dir> --dataset_name <name>"
    echo ""
    echo "Required:"
    echo "  --ctc_dir      Directory with CTC training data (chunks.npy, etc.)"
    echo "  --output_model Output directory for trained model"
    echo "  --dataset_name Dataset name (e.g. C9_128, FX_640)"
    exit 1
}

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --ctc_dir)
            CTC_DIR="$2"
            shift 2
            ;;
        --output_model)
            OUTPUT_MODEL="$2"
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

if [[ -z "$CTC_DIR" || -z "$OUTPUT_MODEL" || -z "$DATASET_NAME" ]]; then
    echo "Error: Missing required arguments"
    usage
fi

# Set batch size based on dataset (standard 10K chunks use batch 32)
BATCH_SIZE=32

echo "Training Standard Model"
echo "======================"
echo "Dataset: $DATASET_NAME"
echo "CTC directory: $CTC_DIR"
echo "Output model: $OUTPUT_MODEL"
echo "Batch size: $BATCH_SIZE"
echo ""

# Create output directory
mkdir -p "$OUTPUT_MODEL"

# Train model
bonito train -f \
    --device cuda \
    --epochs 5 \
    --lr 5e-4 \
    --batch "$BATCH_SIZE" \
    --pretrained dna_r10.4.1_e8.2_400bps_hac@v4.3.0 \
    --directory "$CTC_DIR" \
    "$OUTPUT_MODEL"

echo ""
echo "Training complete for $DATASET_NAME"
echo "Model saved in: $OUTPUT_MODEL"