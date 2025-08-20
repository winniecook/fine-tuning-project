#!/bin/bash

# Evaluate custom chunk fine-tuned models
# Each model tested with its corresponding custom chunk size

set -e

usage() {
    echo "Usage: $0 --models_dir <dir> --test_datasets_dir <dir> --output_dir <dir>"
    echo ""
    echo "Evaluates custom chunk fine-tuned models"
    echo ""
    echo "Required:"
    echo "  --models_dir         Directory containing custom fine-tuned models"
    echo "  --test_datasets_dir  Directory containing test POD5 files"
    echo "  --output_dir         Output directory for results"
    echo ""
    echo "Expected models: FX_160_custom, FX_640_custom, C9_128_custom, etc."
    echo "Chunk sizes: 35K (FX_160,C9_128,C9_256), 55K (FX_640,C9_512), 75K (C9_1024)"
    exit 1
}

get_repeat_unit() {
    case $1 in
        FX_*)
            echo "CGG"
            ;;
        C9_*)
            echo "GGGGCC"
            ;;
        *)
            echo "Error: Unknown dataset type $1"
            exit 1
            ;;
    esac
}

get_flank_file() {
    case $1 in
        FX_*)
            echo "pJAZZ_FX_invariant_parts.fa"
            ;;
        C9_*)
            echo "pJAZZ_C9_invariant_parts.fa"
            ;;
        *)
            echo "Error: Unknown dataset type $1"
            exit 1
            ;;
    esac
}

get_custom_chunk_size() {
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
            echo "35000"  # default
            ;;
    esac
}

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --models_dir)
            MODELS_DIR="$2"
            shift 2
            ;;
        --test_datasets_dir)
            TEST_DIR="$2"
            shift 2
            ;;
        --output_dir)
            OUTPUT_DIR="$2"
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

if [[ -z "$MODELS_DIR" || -z "$TEST_DIR" || -z "$OUTPUT_DIR" ]]; then
    echo "Error: Missing required arguments"
    usage
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "Custom Training Evaluation"
echo "========================="
echo "Models directory: $MODELS_DIR"
echo "Test datasets: $TEST_DIR"
echo "Output directory: $OUTPUT_DIR"
echo ""

# Find all custom models
custom_models=($(find "$MODELS_DIR" -name "*_custom" -type d -exec basename {} \; | sort))

echo "Found custom models: ${custom_models[*]}"
echo ""

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Evaluate each custom model on its corresponding dataset
for model in "${custom_models[@]}"; do
    echo "Processing model: $model"
    
    # Extract dataset name (remove _custom suffix)
    dataset_name="${model%_custom}"
    
    # Get parameters
    repeat_unit=$(get_repeat_unit "$dataset_name")
    flank_file=$(get_flank_file "$dataset_name")
    chunk_size=$(get_custom_chunk_size "$dataset_name")
    
    # Set up paths
    model_path="$MODELS_DIR/$model"
    test_pod5="$TEST_DIR/${dataset_name}_test.pod5"
    result_dir="$OUTPUT_DIR/custom_${dataset_name}"
    fastq_file="$result_dir/basecalls.fastq"
    
    echo "  Dataset: $dataset_name, Chunk size: $chunk_size"
    
    # Skip if model or test dataset doesn't exist
    if [[ ! -d "$model_path" ]]; then
        echo "  Skipping - model not found: $model_path"
        continue
    fi
    
    if [[ ! -f "$test_pod5" ]]; then
        echo "  Skipping - test file not found: $test_pod5"
        continue
    fi
    
    # Create result directory
    mkdir -p "$result_dir"
    
    echo "  Basecalling with custom fine-tuned model..."
    
    # Basecall with custom fine-tuned model (dataset-specific chunk size)
    bonito basecaller \
        "$model_path" \
        "$test_pod5" \
        --device cuda \
        --chunksize "$chunk_size" \
        > "$fastq_file"
    
    echo "  Running evaluation pipeline..."
    
    # Run Nextflow evaluation
    cd "$SCRIPT_DIR"
    nextflow run get_per-read_repeat_metrics.nf \
        --reads_fq "$fastq_file" \
        --flank_fasta "$flank_file" \
        --publish_dir "$result_dir/" \
        --threads 4 \
        --repeat_expansion_unit "$repeat_unit" \
        --tolerance 5
    
    echo "  Complete: $result_dir"
done

echo ""
echo "Custom training evaluation complete!"
echo ""
echo "Results in: $OUTPUT_DIR"
echo "Directories: custom_FX_160, custom_FX_640, custom_C9_128, etc."