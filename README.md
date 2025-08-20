# Enhancing Nanopore Basecalling Accuracy for Pathogenic Repeat Expansions Through Fine-Tuning

**Author:** Winnie Cook  
**Institution:** King's College London

## Project Overview

This repository contains the codebase for investigating Oxford Nanopore basecalling performance on repeat expansion sequences. The project explores custom chunking strategies and fine-tuning approaches to address challenges in basecalling repetitive genomic regions.

The work investigates the sequence context hypothesis - that standard 10,000-sample chunking may limit complete sequence representation during training for longer repeat expansions. Through systematic analysis of training data integrity and custom chunk size implementation, we explored whether providing complete sequence context during training could improve repeat expansion basecalling.

Results showed improvements in repeat length accuracy (5-36 percentage points) and sequence purity (99.41-100.00%) on the datasets tested on, whilst maintaining accuracy on standard, genome-wide datasets (HG002) within 1.0% of baseline, suggesting that sequence context during training may influence repeat enumeration accuracy.

## 1. Standard Genome Processing (HG002)

### Filter BAM alignments
**Script:** `01_data_preparation/HG002_processing/filter_enhanced.py`  
**Methods section:** 2.6 - Standard Genomic Benchmarking Dataset (HG002)  
**Requirements:** Python 3.8+, pysam v0.22.0  
- **Input:** Aligned BAM file
- **Output:** Filtered BAM, read IDs list, stats
```bash
python filter_enhanced.py input.bam -o filtered_output/ -p hg002
```

### Compare POD5 with filtered reads
**Script:** `01_data_preparation/HG002_processing/compare_pod5_filtered.py`  
**Methods section:** 2.6 - Standard Genomic Benchmarking Dataset (HG002)  
**Requirements:** Python 3.8+  
- **Input:** Filtered read IDs, POD5 read IDs
- **Output:** Verification report
```bash
python compare_pod5_filtered.py filtered_read_ids.txt pod5_read_ids.txt -o output/
```

### HG002 processing guide
**File:** `01_data_preparation/HG002_processing/hg002_commands.md`  
**Methods section:** 2.6 - Standard Genomic Benchmarking Dataset (HG002)

## 2. Data Preparation

### Filter POD5 files
**Pipeline:** `pipelines/dual_flanking_filtering/filter_pod5_by_dual_match_reads.nf`  
**Methods section:** 2.7.4 - Initial Basecalling and Filtering (dual-flank strategy)  
**Requirements:** Nextflow v21.04+, BLAST+ v2.10+, SeqKit v2.0+  
- **Input:** Raw POD5 files, flank sequences
- **Output:** Filtered POD5 files
```bash
nextflow run filter_pod5_by_dual_match_reads.nf \
    --pod5_dir raw_pod5/ \
    --flank_fasta pJAZZ_FX_invariant_parts.fa \
    --publish_dir filtered/
```

### Split POD5 data
**Script:** `01_data_preparation/01_pod5_splitting.py`  
**Methods section:** 2.13.2 - Data Splitting Protocol  
**Requirements:** Python 3.8+, pod5 library v3.0.0  
- **Input:** Filtered POD5 file
- **Output:** 90% train / 10% test split
```bash
python 01_pod5_splitting.py \
    --input filtered_data.pod5 \
    --output_dir split_data/ \
    --dataset_name FX_640
```

## 3. Training Data Preparation

### Standard chunks (10K samples)
**Script:** `01_data_preparation/02_prepare_ctc_standard.sh`  
**Methods section:** 2.9.1 - Training Data Preparation  
**Requirements:** Bonito v0.7.3, CUDA GPU  
- **Input:** Training POD5, reference FASTA
- **Output:** CTC training arrays
```bash
bash 02_prepare_ctc_standard.sh \
    --input_pod5 train_data.pod5 \
    --reference reference.fa \
    --output_dir ctc_standard/
```

### Custom chunks (35K/55K/75K samples)
**Script:** `01_data_preparation/02_prepare_ctc_custom.sh`  
**Methods section:** 2.9.1 - Training Data Preparation, 2.12.2 - Chunk Size Calculation  
**Requirements:** Bonito v0.7.3, CUDA GPU  
- **Input:** Training POD5, reference FASTA
- **Output:** CTC training arrays
```bash
bash 02_prepare_ctc_custom.sh \
    --input_pod5 train_data.pod5 \
    --reference reference.fa \
    --output_dir ctc_custom/
```

## 4. Model Training

### Standard training
**Script:** `02_training/train_standard_model.sh`  
**Methods section:** 2.13.1 - Fine-tuning Methodology  
**Requirements:** Bonito v0.7.3, CUDA GPU  
- **Input:** Standard CTC data
- **Output:** Fine-tuned model
```bash
bash train_standard_model.sh \
    --ctc_dir ctc_standard/ \
    --output_model models/standard_FX_640/
```

### Custom training
**Script:** `02_training/training_custom/train_custom_model.sh`  
**Methods section:** 2.13.1 - Fine-tuning Methodology, 2.12.1 - Selection Strategy  
**Requirements:** Bonito v0.7.3, CUDA GPU  
- **Input:** Custom CTC data
- **Output:** Fine-tuned model
```bash
bash training_custom/train_custom_model.sh \
    --ctc_dir ctc_custom/ \
    --output_model models/custom_FX_640/
```

## 5. Evaluation

### Baseline evaluation
**Script:** `03_evaluation/baseline_evaluation.sh`  
**Methods section:** 4.1.1 - Bonito HAC Demonstrates Strong Genome-Wide Performance  
**Requirements:** Bonito v0.7.3, Nextflow v21.04+, BLAST+ v2.10+  
- **Input:** Test POD5 files
- **Output:** Baseline performance metrics
```bash
bash baseline_evaluation.sh \
    --test_datasets_dir test_data/ \
    --output_dir baseline_results/
```

### Standard model evaluation
**Script:** `03_evaluation/evaluation_standard/standard_evaluations.sh`  
**Methods section:** 4.2 - Initial Fine-Tuning Shows Improvements Limited to Short Expansions  
**Requirements:** Bonito v0.7.3, Nextflow v21.04+, BLAST+ v2.10+  
- **Input:** Standard models, test POD5
- **Output:** Standard model metrics
```bash
bash evaluation_standard/standard_evaluations.sh \
    --models_dir models/standard/ \
    --test_datasets_dir test_data/ \
    --output_dir standard_results/
```

### Custom model evaluation
**Script:** `03_evaluation/evaluation_custom/custom_evaluation.sh`  
**Methods section:** 4.4.1 - Custom Chunk Size Training Provides the Sequence Context Necessary  
**Requirements:** Bonito v0.7.3, Nextflow v21.04+, BLAST+ v2.10+  
- **Input:** Custom models, test POD5
- **Output:** Custom model metrics
```bash
bash evaluation_custom/custom_evaluation.sh \
    --models_dir models/custom/ \
    --test_datasets_dir test_data/ \
    --output_dir custom_results/
```

### Repeat metrics pipeline
**Pipeline:** `pipelines/repeat_metrics_pipeline/get_per-read_repeat_metrics.nf`  
**Methods section:** 2.8 - Evaluation Metrics and Pipeline for Performance on Repeat Expansion Dataset  
**Requirements:** Nextflow v21.04+, BLAST+ v2.10+, SeqKit v2.0+  
- **Input:** FASTQ basecalls, flank sequences
- **Output:** Per-read repeat metrics
```bash
nextflow run get_per-read_repeat_metrics.nf \
    --reads_fq basecalls.fastq \
    --flank_fasta pJAZZ_C9_invariant_parts.fa \
    --repeat_expansion_unit GGGGCC
```

## 6. Analysis

### Reference length analysis
**Script:** `04_analysis/reference_length_analysis.py`  
**Methods section:** 2.10.2 - Reference Length Coverage Analysis  
**Requirements:** Python 3.8+, NumPy  
- **Input:** CTC training directories
- **Output:** Reference length statistics
```bash
python reference_length_analysis.py \
    --base_path ctc_training_data/ \
    --datasets C9_128_ctc C9_256_ctc FX_160_ctc
```

### Read fragmentation analysis
**Script:** `04_analysis/read_fragmentation/read_fragmentation.py`  
**Methods section:** 2.10.3 - Read Fragmentation and Segmentation Analysis  
**Requirements:** Python 3.8+, pandas, NumPy  
- **Input:** CTC training data
- **Output:** Fragmentation statistics
```bash
python read_fragmentation/read_fragmentation.py \
    --ctc_dir ctc_training_data/ \
    --output fragmentation_stats.csv
```

### Missing chunk analysis
**Script:** `04_analysis/missing_chunk_analysis/missing_chunk_analysis.py`  
**Methods section:** 2.10.4 - Missing Chunk Detection and Data Loss Quantification  
**Requirements:** Python 3.8+, pandas, NumPy  
- **Input:** CTC training directories
- **Output:** Data completeness statistics
```bash
python missing_chunk_analysis/missing_chunk_analysis.py \
    --base_path ctc_training_data/ \
    --datasets C9_128_ctc C9_256_ctc FX_160_ctc
```

## 7. Visualisation

### Performance comparison plots
**Script:** `05_visualisation/10k_vs_custom.R`  
**Methods section:** Results section - Figures 4.1, 4.2  
**Requirements:** R 4.0+, ggplot2, dplyr  
- **Input:** Analysis CSV files
- **Output:** Performance comparison plots
```bash
Rscript 10k_vs_custom.R
```

### Missing chunk visualisation
**Script:** `05_visualisation/missing_chunk_vis.R`  
**Methods section:** Results section - Figure 4.2  
**Requirements:** R 4.0+, ggplot2  
- **Input:** Missing chunk analysis CSV
- **Output:** Data loss plots
```bash
Rscript missing_chunk_vis.R
```

### Chunk distribution plots
**Script:** `05_visualisation/plot_chunk_distribution.R`  
**Methods section:** Results section - Figure 4.1  
**Requirements:** R 4.0+, ggplot2  
- **Input:** Fragmentation analysis CSV
- **Output:** Chunk distribution plots
```bash
Rscript plot_chunk_distribution.R
```

## Workflow Summary

1. **Prepare data:** Filter POD5 → Split train/test → Generate CTC data
2. **Train models:** Standard (10K) vs Custom (35K/55K/75K) chunks
3. **Evaluate:** Compare baseline vs standard vs custom approaches
4. **Analyse:** Quantify sequence coverage and data completeness improvements
5. **Visualise:** Generate plots comparing all approaches

## Reference Files

- `pJAZZ_C9_invariant_parts.fa` - C9orf72 flanking sequences
- `pJAZZ_FX_invariant_parts.fa` - Fragile X flanking sequences
