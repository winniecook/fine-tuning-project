# HG002 Processing Commands

Simple command reference for HG002 benchmarking dataset preparation (Section 2.6).

## 1. Download HG002 Data

```bash
# Download HG002 reference genome
wget https://github.com/marbl/HG002/releases/download/v1.1/hg002v1.1.fasta.gz
gunzip hg002v1.1.fasta.gz

# Download R10.4.1 sequencing data from Garvan Institute
wget s3://gtgseq/ont-r10-5khz-dna/NA24385/raw/PGXXXX230339_reads.blow5
wget s3://gtgseq/ont-r10-5khz-dna/NA24385/raw/PGXXXX230339_reads.blow5.idx
```

## 2. Convert BLOW5 to POD5

```bash
# Convert using blue-crab with 16 threads
blue-crab s2p PGXXXX230339_reads.blow5 -d pod5_data -p 16
```

## 3. Basecall with Dorado SUP

```bash
# Initial basecalling with highest quality model
dorado basecaller dna_r10.4.1_e8.2_400bps_sup@v5.0.0 \
    pod5_data > PGXXXX230339_reads_dorado-sup-v5.fastq
```

## 4. Align to HG002 Reference

```bash
# Align with minimap2 using nanopore preset
minimap2 -ax map-ont -t 50 \
    hg002v1.1.fasta \
    PGXXXX230339_reads_dorado-sup-v5.fastq > aligned_reads.sam
```

## 5. Convert and Sort BAM

```bash
# Convert SAM to BAM and sort
samtools sort -@ 150 aligned_reads.sam > aligned_reads.bam

# Generate alignment statistics  
samtools flagstat aligned_reads.bam > alignment_stats.txt
```

## 6. Quality Filtering

```bash
# Remove secondary, supplementary, and duplicate alignments
python filter_enhanced.py aligned_reads.bam -o filtered_output/ -p hg002
```

## 7. Verify Data Integrity

```bash
# Cross-reference filtered reads with original POD5 data
python compare_pod5_filtered.py \
    filtered_output/hg002_*_filtered_read_ids.txt \
    pod5_read_ids.txt \
    -o verification_output/
```

## 8. Final Validation

```bash
# Confirm 145,425 high-quality reads with 100% primary alignment
samtools flagstat filtered_output/hg002_*_filtered.bam
```

## Expected Results

- **Initial reads:** 39,701,986 (99.10% mapping rate)
- **After filtering:** 418,905 reads  
- **Final dataset:** 145,425 reads (100% primary alignment)
- **Match rate with POD5:** 99.42%

---

**Note:** These commands reflect the exact methodology from the thesis. Adjust thread counts (`-t`, `-@`) based on your system.