#!/usr/bin/env python3
"""
Compare filtered read IDs with Pod5 file read IDs to verify data integrity.
This script cross-references filtered BAM read IDs with the original Pod5 data
to ensure that filtered reads are available in the raw sequencing data.

Author: Winnie Cook 
Version: 1.0
Date: 2025
"""

import sys
import os
import argparse

def read_ids_from_file(filepath):
    """
    Read read IDs from a text file, one per line.
    
    Args:
        filepath (str): Path to file containing read IDs
        
    Returns:
        set: Set of unique read IDs
    """
    try:
        with open(filepath, 'r') as f:
            ids = set(line.strip() for line in f if line.strip())
        return ids
    except Exception as e:
        raise IOError(f"Error reading file {filepath}: {e}")

def write_report(output_path, filtered_count, pod5_count, found_count, missing_count, missing_ids):
    """
    Write comparison report to file.
    
    Args:
        output_path (str): Path for output report
        filtered_count (int): Number of filtered read IDs
        pod5_count (int): Number of Pod5 read IDs
        found_count (int): Number of matched read IDs
        missing_count (int): Number of missing read IDs
        missing_ids (set): Set of missing read IDs
    """
    found_percent = (found_count * 100.0 / filtered_count) if filtered_count > 0 else 0
    missing_percent = (missing_count * 100.0 / filtered_count) if filtered_count > 0 else 0
    
    try:
        with open(output_path, 'w') as f:
            f.write("Pod5 Read ID Verification Report\n")
            f.write("="*40 + "\n\n")
            f.write(f"Total filtered read IDs: {filtered_count:,}\n")
            f.write(f"Total Pod5 read IDs: {pod5_count:,}\n")
            f.write(f"Filtered IDs found in Pod5: {found_count:,} ({found_percent:.2f}%)\n")
            f.write(f"Filtered IDs missing from Pod5: {missing_count:,} ({missing_percent:.2f}%)\n\n")
            
            if missing_ids:
                f.write("Sample of missing read IDs (first 10):\n")
                f.write("-" * 40 + "\n")
                for read_id in list(missing_ids)[:10]:
                    f.write(f"{read_id}\n")
                    
        print(f"Report written to: {output_path}")
        
    except Exception as e:
        raise IOError(f"Error writing report file: {e}")

def write_matched_ids(output_path, matched_ids):
    """
    Write matched read IDs to file.
    
    Args:
        output_path (str): Path for matched IDs output
        matched_ids (set): Set of matched read IDs
    """
    try:
        with open(output_path, 'w') as f:
            for read_id in sorted(matched_ids):
                f.write(f"{read_id}\n")
                
        print(f"Matched read IDs written to: {output_path}")
        
    except Exception as e:
        raise IOError(f"Error writing matched IDs file: {e}")

def compare_read_ids(filtered_ids_file, pod5_ids_file, output_report, matched_ids_file, verbose=True):
    """
    Compare filtered read IDs with Pod5 read IDs and generate reports.
    
    Args:
        filtered_ids_file (str): Path to filtered read IDs file
        pod5_ids_file (str): Path to Pod5 read IDs file
        output_report (str): Path for output report
        matched_ids_file (str): Path for matched IDs output
        verbose (bool): Whether to print progress messages
        
    Returns:
        dict: Dictionary containing comparison statistics
    """
    
    # Read filtered read IDs
    if verbose:
        print("Reading filtered read IDs...")
    filtered_ids = read_ids_from_file(filtered_ids_file)
    if verbose:
        print(f"Found {len(filtered_ids):,} unique read IDs in filtered set")
    
    # Read Pod5 read IDs
    if verbose:
        print("Reading Pod5 read IDs...")
    pod5_ids = read_ids_from_file(pod5_ids_file)
    if verbose:
        print(f"Found {len(pod5_ids):,} unique read IDs in Pod5 file")
    
    # Find intersection and differences
    if verbose:
        print("Comparing read ID sets...")
    found_ids = filtered_ids.intersection(pod5_ids)
    missing_ids = filtered_ids - pod5_ids
    
    # Calculate percentages
    found_percent = (len(found_ids) * 100.0 / len(filtered_ids)) if filtered_ids else 0
    missing_percent = (len(missing_ids) * 100.0 / len(filtered_ids)) if filtered_ids else 0
    
    if verbose:
        print(f"Found {len(found_ids):,} filtered IDs in the Pod5 file ({found_percent:.2f}%)")
        print(f"Missing {len(missing_ids):,} filtered IDs from the Pod5 file ({missing_percent:.2f}%)")
    
    # Write output files
    write_report(output_report, len(filtered_ids), len(pod5_ids), 
                len(found_ids), len(missing_ids), missing_ids)
    write_matched_ids(matched_ids_file, found_ids)
    
    # Return statistics
    return {
        'filtered_count': len(filtered_ids),
        'pod5_count': len(pod5_ids),
        'found_count': len(found_ids),
        'missing_count': len(missing_ids),
        'found_percent': found_percent,
        'missing_percent': missing_percent
    }

def main():
    """
    Main function to handle command line arguments and execute comparison.
    """
    
    parser = argparse.ArgumentParser(
        description='Compare filtered read IDs with Pod5 file read IDs for data integrity verification',
        epilog='''
Example usage:
  python3 compare_pod5_filtered.py filtered_ids.txt pod5_ids.txt
  python3 compare_pod5_filtered.py filtered_ids.txt pod5_ids.txt -o /output/dir/ -p sample01
        ''',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    # Required arguments
    parser.add_argument('filtered_ids', 
                       help='Path to filtered read IDs file')
    parser.add_argument('pod5_ids', 
                       help='Path to Pod5 read IDs file')
    
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
                       version='%(prog)s 1.0')
    
    args = parser.parse_args()
    
    # Validate input files
    for filepath in [args.filtered_ids, args.pod5_ids]:
        if not os.path.exists(filepath):
            print(f"Error: Input file '{filepath}' does not exist.")
            sys.exit(1)
    
    # Generate output paths
    if not os.path.exists(args.output_dir):
        try:
            os.makedirs(args.output_dir)
        except Exception as e:
            print(f"Error creating output directory: {e}")
            sys.exit(1)
    
    # Create output filenames
    prefix = f"{args.prefix}_" if args.prefix else ""
    output_report = os.path.join(args.output_dir, f"{prefix}pod5_verification_report.txt")
    matched_ids_file = os.path.join(args.output_dir, f"{prefix}matched_pod5_read_ids.txt")
    
    # Execute comparison
    try:
        stats = compare_read_ids(
            args.filtered_ids,
            args.pod5_ids,
            output_report,
            matched_ids_file,
            verbose=not args.quiet
        )
        
        if not args.quiet:
            print(f"\nComparison completed successfully!")
            print(f"Match rate: {stats['found_percent']:.2f}% ({stats['found_count']:,}/{stats['filtered_count']:,})")
        
    except Exception as e:
        print(f"Error during comparison: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()