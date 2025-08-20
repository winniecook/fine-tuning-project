#!/usr/bin/env python3
"""
Enhanced filtering script for nanopore sequencing BAM files.
Removes secondary alignments, supplementary alignments, and duplicate reads
to retain only high-quality, uniquely aligned reads.

This script is designed for processing long-read sequencing data where
duplicate removal and alignment filtering are crucial for downstream analysis.

Author: Winnie Cook
Version: 1.1
Date: 2025
"""

import sys
import os
import argparse
import pysam
from collections import defaultdict

def create_output_paths(input_path, output_dir, prefix):
    """
    Generate output file paths based on input parameters.
    
    Args:
        input_path (str): Path to input BAM file
        output_dir (str): Directory for output files
        prefix (str): Custom prefix for output files (optional)
    
    Returns:
        tuple: (output_bam_path, read_ids_file_path, log_file_path)
    """
    # Extract base name from input file
    base_name = os.path.splitext(os.path.basename(input_path))[0]
    
    # Apply custom prefix if provided
    if prefix:
        base_name = f"{prefix}_{base_name}"
    
    # Create output directory if it doesn't exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Generate full paths
    output_bam_path = os.path.join(output_dir, f"{base_name}_filtered.bam")
    read_ids_file_path = os.path.join(output_dir, f"{base_name}_filtered_read_ids.txt")
    log_file_path = os.path.join(output_dir, f"{base_name}_filtering_log.txt")
    
    return output_bam_path, read_ids_file_path, log_file_path

def write_log_file(log_path, stats_dict, input_path, output_path):
    """
    Write detailed filtering statistics to a log file.
    
    Args:
        log_path (str): Path to output log file
        stats_dict (dict): Dictionary containing filtering statistics
        input_path (str): Path to input BAM file
        output_path (str): Path to output BAM file
    """
    try:
        with open(log_path, 'w') as log_file:
            log_file.write("BAM Filtering Log\n")
            log_file.write("="*50 + "\n")
            log_file.write(f"Input file: {input_path}\n")
            log_file.write(f"Output file: {output_path}\n")
            log_file.write(f"Filter date: {stats_dict['timestamp']}\n\n")
            
            log_file.write("Filtering Parameters:\n")
            log_file.write("- Remove secondary alignments (SAM flag 256)\n")
            log_file.write("- Remove supplementary alignments (SAM flag 2048)\n")
            log_file.write("- Remove duplicate reads based on query name\n\n")
            
            log_file.write("Results Summary:\n")
            log_file.write(f"Total reads processed: {stats_dict['total_reads']:,}\n")
            log_file.write(f"Secondary alignments removed: {stats_dict['secondary_removed']:,}\n")
            log_file.write(f"Supplementary alignments removed: {stats_dict['supplementary_removed']:,}\n")
            log_file.write(f"Duplicate reads removed: {stats_dict['duplicate_removed']:,}\n")
            log_file.write(f"Reads retained: {stats_dict['filtered_reads']:,}\n")
            log_file.write(f"Retention rate: {stats_dict['retention_rate']:.2f}%\n")
            
    except Exception as e:
        print(f"Warning: Could not write log file: {e}")

def filter_bam(input_bam_path, output_bam_path, read_ids_file_path, log_file_path, verbose=True):
    """
    Filter BAM file to remove secondary, supplementary, and duplicate alignments.
    
    This function processes a BAM file and removes reads that meet any of the
    following criteria:
    - Secondary alignments (reads with multiple mapping locations)
    - Supplementary alignments (chimeric reads split across locations)
    - Duplicate reads (multiple reads with identical query names)
    
    Args:
        input_bam_path (str): Path to input BAM file
        output_bam_path (str): Path for filtered BAM output
        read_ids_file_path (str): Path for read IDs text file
        log_file_path (str): Path for filtering log file
        verbose (bool): Whether to print progress updates
    
    Returns:
        dict: Dictionary containing filtering statistics
    """
    
    if verbose:
        print(f"Processing: {input_bam_path}")
        print(f"Output BAM: {output_bam_path}")
        print(f"Read IDs file: {read_ids_file_path}")
        print(f"Log file: {log_file_path}")
    
    # Open input BAM file
    try:
        input_bam = pysam.AlignmentFile(input_bam_path, "rb")
    except Exception as e:
        raise IOError(f"Error opening input BAM file: {e}")
    
    # Open output BAM file
    try:
        output_bam = pysam.AlignmentFile(output_bam_path, "wb", template=input_bam)
    except Exception as e:
        input_bam.close()
        raise IOError(f"Error creating output BAM file: {e}")
    
    # Initialise statistics counters
    total_reads = 0
    filtered_reads = 0
    secondary_removed = 0
    supplementary_removed = 0
    duplicate_removed = 0
    
    # Store filtered read IDs for output
    filtered_read_ids = []
    
    # Track read names to identify duplicates
    seen_read_names = set()
    
    if verbose:
        print("Filtering reads...")
    
    # Process each read in the input BAM file
    for read in input_bam:
        total_reads += 1
        
        # Skip secondary alignments (SAM flag 256)
        # These represent alternative alignments for reads that map to multiple locations
        if read.is_secondary:
            secondary_removed += 1
            continue
        
        # Skip supplementary alignments (SAM flag 2048)
        # These represent additional alignments for chimeric reads
        if read.is_supplementary:
            supplementary_removed += 1
            continue
        
        # Skip duplicate reads based on query name
        # In nanopore data, true duplicates are rare, but this catches any present
        if read.query_name in seen_read_names:
            duplicate_removed += 1
            continue
        
        # Read passes all filtering criteria - retain it
        seen_read_names.add(read.query_name)
        output_bam.write(read)
        filtered_read_ids.append(read.query_name)
        filtered_reads += 1
        
        # Print progress updates for large files
        if verbose and total_reads % 100000 == 0:
            print(f"Processed {total_reads:,} reads, retained {filtered_reads:,}")
    
    # Close BAM files
    input_bam.close()
    output_bam.close()
    
    # Write filtered read IDs to text file
    # This allows tracking of which specific reads passed filtering
    try:
        with open(read_ids_file_path, 'w') as f:
            for read_id in filtered_read_ids:
                f.write(f"{read_id}\n")
    except Exception as e:
        raise IOError(f"Error writing read IDs file: {e}")
    
    # Create BAM index for the filtered file
    try:
        if verbose:
            print("Indexing output BAM file...")
        pysam.index(output_bam_path)
    except Exception as e:
        print(f"Warning: Could not index output BAM file: {e}")
    
    # Calculate statistics
    retention_rate = (filtered_reads/total_reads)*100 if total_reads > 0 else 0
    
    # Prepare statistics dictionary
    import datetime
    stats = {
        'total_reads': total_reads,
        'filtered_reads': filtered_reads,
        'secondary_removed': secondary_removed,
        'supplementary_removed': supplementary_removed,
        'duplicate_removed': duplicate_removed,
        'retention_rate': retention_rate,
        'timestamp': datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    }
    
    # Write log file
    write_log_file(log_file_path, stats, input_bam_path, output_bam_path)
    
    # Print summary statistics to console
    if verbose:
        print("\n" + "="*50)
        print("FILTERING SUMMARY")
        print("="*50)
        print(f"Total reads processed: {total_reads:,}")
        print(f"Secondary alignments removed: {secondary_removed:,}")
        print(f"Supplementary alignments removed: {supplementary_removed:,}")
        print(f"Duplicate reads removed: {duplicate_removed:,}")
        print(f"Reads retained: {filtered_reads:,}")
        print(f"Retention rate: {retention_rate:.2f}%")
        print("="*50)
    
    return stats

def main():
    """
    Main function to handle command line arguments and execute BAM filtering.
    
    Parses command line arguments and calls the filtering function with
    appropriate parameters.
    """
    
    # Set up argument parser
    parser = argparse.ArgumentParser(
        description='Filter BAM files to remove secondary, supplementary, and duplicate alignments',
        epilog='''
Example usage:
  python3 filter-enhanced.py input.bam
  python3 filter-enhanced.py input.bam -o /path/to/output/
  python3 filter-enhanced.py input.bam -o /path/to/output/ -p sample01
  python3 filter-enhanced.py input.bam -o /path/to/output/ -p sample01 --quiet
        ''',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    # Required argument
    parser.add_argument('input_bam', 
                       help='Input BAM file path')
    
    # Optional arguments
    parser.add_argument('-o', '--output-dir', 
                       default='.', 
                       help='Output directory (default: current directory)')
    
    parser.add_argument('-p', '--prefix', 
                       default='', 
                       help='Prefix for output files (optional)')
    
    parser.add_argument('--quiet', 
                       action='store_true', 
                       help='Suppress progress output')
    
    parser.add_argument('--version', 
                       action='version', 
                       version='%(prog)s 1.1')
    
    # Parse arguments
    args = parser.parse_args()
    
    # Validate input file
    if not os.path.exists(args.input_bam):
        print(f"Error: Input file '{args.input_bam}' does not exist.")
        sys.exit(1)
    
    # Check file extension (warning only)
    if not args.input_bam.lower().endswith('.bam'):
        print(f"Warning: Input file '{args.input_bam}' does not have .bam extension.")
    
    # Generate output paths
    try:
        output_bam, read_ids_file, log_file = create_output_paths(
            args.input_bam, args.output_dir, args.prefix
        )
    except Exception as e:
        print(f"Error creating output paths: {e}")
        sys.exit(1)
    
    # Execute filtering
    try:
        stats = filter_bam(
            args.input_bam, 
            output_bam, 
            read_ids_file, 
            log_file, 
            verbose=not args.quiet
        )
        
        if not args.quiet:
            print(f"\nFiltering completed successfully!")
            print(f"Output BAM: {output_bam}")
            print(f"Read IDs file: {read_ids_file}")
            print(f"Log file: {log_file}")
            print(f"Filtered reads: {stats['filtered_reads']:,}")
        
    except Exception as e:
        print(f"Error during filtering: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()