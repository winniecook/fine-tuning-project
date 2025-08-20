#!/bin/bash

# Evaluate baseline ONT model on all test datasets
# Uses standard 10K chunk size for all evaluations

set -e

usage() {
    echo "Usage: $0 --test_datasets_dir <dir> --output_dir <dir>"
    echo ""
    echo "Evaluates baseline ONT model on all test datasets"
    echo ""
    echo "Required:"
    echo "  --test_datasets_dir  Directory containing test POD5 files"
    echo "  --output_dir         Output directory for results"
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

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
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

if [[ -z "$TEST_DIR" || -z "$OUTPUT_DIR" ]]; then
    echo "Error: Missing required arguments"
    usage
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "Baseline Model Evaluation"
echo "========================"
echo "Test datasets: $TEST_DIR"
echo "Output directory: $OUTPUT_DIR"
echo "Chunk size: 10,000 (standard)"
echo ""

# Find all test datasets
test_datasets=($(find "$TEST_DIR" -name "*_test.pod5" -exec basename {} .pod5 \; | sed 's/_test$//' | sort))

echo "Found test datasets: ${test_datasets[*]}"
echo ""

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Evaluate baseline on each dataset
for dataset in "${test_datasets[@]}"; do
    echo "Processing dataset: $dataset"
    
    # Get parameters
    repeat_unit=$(get_repeat_unit "$dataset")
    flank_file=$(get_flank_file "$dataset")
    
    # Set up paths
    test_pod5="$TEST_DIR/${dataset}_test.pod5"
    result_dir="$OUTPUT_DIR/baseline_${dataset}"
    fastq_file="$result_dir/basecalls.fastq"
    
    # Skip if test dataset doesn't exist
    if [[ ! -f "$test_pod5" ]]; then
        echo "  Skipping - test file not found: $test_pod5"
        continue
    fi
    
    # Create result directory
    mkdir -p "$result_dir"
    
    echo "  Basecalling with baseline model..."
    
    # Basecall with baseline ONT model (10K chunks)
    bonito basecaller \
        dna_r10.4.1_e8.2_400bps_hac@v4.3.0 \
        "$test_pod5" \
        --device cuda \
        --chunksize 10000 \
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
echo "Baseline evaluation complete!"
echo ""
echo "Results in: $OUTPUT_DIR"
echo "Directories: baseline_FX_160, baseline_FX_640, baseline_C9_128, etc."