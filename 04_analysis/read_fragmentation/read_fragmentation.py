#!/usr/bin/env python3
"""
Read fragmentation analysis script for CTC training data.
Analyses how reads are split into chunks and quantifies fragmentation patterns.

Based on Section 2.10.3 - Read Fragmentation and Segmentation Analysis
"""

import numpy as np
import pandas as pd
import os
import argparse
import glob

def load_dataset_data(dataset_path):
    """Load TSV file or basecalls summary for a dataset"""
    # Try different possible file names
    possible_files = [
        "ctc_generation_log_summary.tsv",
        "basecalls_summary.tsv", 
        "summary.tsv",
        "ctc_log.tsv"
    ]
    
    for filename in possible_files:
        tsv_path = os.path.join(dataset_path, filename)
        if os.path.exists(tsv_path):
            try:
                return pd.read_csv(tsv_path, sep='\t')
            except Exception as e:
                print(f"  Error reading {filename}: {e}")
                continue
    
    return None

def analyze_read_fragmentation(df):
    """
    Analyze how reads are fragmented into chunks.
    Expects read_id format: 'read_id:chunk_id:total_chunks'
    """
    
    # Check if read_id contains chunk information
    if ':' not in str(df['read_id'].iloc[0]):
        print("  No chunk information found in read_id format")
        return None, None
    
    # Split read_id into components
    read_parts = df['read_id'].str.extract(r'(?P<read_id_base>[^:]+):(?P<chunk_id>\d+):(?P<total_chunks>\d+)')
    read_parts[['chunk_id', 'total_chunks']] = read_parts[['chunk_id', 'total_chunks']].astype(int)
    
    # Add parsed info back to dataframe
    df_with_parts = df.copy()
    df_with_parts = df_with_parts.join(read_parts)
    
    # Analyze fragmentation per read
    read_stats = df_with_parts.groupby('read_id_base').agg(
        observed_chunks=('chunk_id', 'count'),
        reported_total_chunks=('total_chunks', 'max'),
        duration_max=('duration', 'max') if 'duration' in df.columns else ('chunk_id', 'count')
    ).reset_index()
    
    # Calculate chunk distribution
    chunk_counts = read_stats['reported_total_chunks'].value_counts().sort_index()
    total_reads = len(read_stats)
    
    # Create fragmentation summary
    fragmentation_stats = {}
    for chunks, count in chunk_counts.items():
        percentage = (count / total_reads) * 100
        fragmentation_stats[chunks] = {
            'count': count, 
            'percentage': percentage,
            'is_single_chunk': chunks == 1
        }
    
    # Calculate key metrics
    single_chunk_reads = chunk_counts.get(1, 0)
    multi_chunk_reads = total_reads - single_chunk_reads
    
    single_chunk_percentage = (single_chunk_reads / total_reads) * 100
    multi_chunk_percentage = (multi_chunk_reads / total_reads) * 100
    
    # Find missing chunks (data loss)
    read_stats['missing_chunks'] = read_stats['reported_total_chunks'] - read_stats['observed_chunks']
    missing_chunk_reads = read_stats[read_stats['missing_chunks'] > 0]
    
    summary = {
        'total_reads': total_reads,
        'single_chunk_reads': single_chunk_reads,
        'multi_chunk_reads': multi_chunk_reads,
        'single_chunk_percentage': single_chunk_percentage,
        'multi_chunk_percentage': multi_chunk_percentage,
        'average_chunks_per_read': read_stats['reported_total_chunks'].mean(),
        'max_chunks_per_read': read_stats['reported_total_chunks'].max(),
        'reads_with_missing_chunks': len(missing_chunk_reads),
        'missing_chunk_percentage': (len(missing_chunk_reads) / total_reads) * 100
    }
    
    return fragmentation_stats, summary

def analyze_ctc_directory(ctc_dir, chunk_type="unknown"):
    """
    Analyze all datasets in a CTC directory for read fragmentation.
    """
    
    # Find dataset directories
    dataset_dirs = []
    for item in os.listdir(ctc_dir):
        item_path = os.path.join(ctc_dir, item)
        if os.path.isdir(item_path) and ('ctc' in item.lower() or chunk_type in item):
            dataset_dirs.append(item_path)
    
    dataset_dirs.sort()
    
    if not dataset_dirs:
        print(f"No dataset directories found in {ctc_dir}")
        return {}
    
    print(f"READ FRAGMENTATION ANALYSIS - {chunk_type.upper()} CHUNKS")
    print("=" * 60)
    
    all_results = {}
    
    for dataset_dir in dataset_dirs:
        dataset_name = os.path.basename(dataset_dir)
        
        print(f"\n{dataset_name}:")
        print("-" * 40)
        
        # Load data
        df = load_dataset_data(dataset_dir)
        
        if df is None:
            print("  No TSV file found")
            continue
        
        print(f"  Total entries: {len(df)}")
        
        # Analyze fragmentation
        try:
            fragmentation_stats, summary = analyze_read_fragmentation(df)
            
            if fragmentation_stats is None:
                continue
            
            # Display results
            print(f"  Total unique reads: {summary['total_reads']}")
            print(f"  Average chunks per read: {summary['average_chunks_per_read']:.1f}")
            print(f"  Max chunks per read: {summary['max_chunks_per_read']}")
            
            print(f"  Fragmentation breakdown:")
            print(f"    Single chunk reads: {summary['single_chunk_percentage']:5.1f}% ({summary['single_chunk_reads']} reads)")
            print(f"    Multi-chunk reads:  {summary['multi_chunk_percentage']:5.1f}% ({summary['multi_chunk_reads']} reads)")
            
            if summary['reads_with_missing_chunks'] > 0:
                print(f"    Missing chunks:     {summary['missing_chunk_percentage']:5.1f}% ({summary['reads_with_missing_chunks']} reads)")
            
            # Show detailed chunk distribution
            print(f"  Detailed chunk distribution:")
            for chunks in sorted(fragmentation_stats.keys()):
                count = fragmentation_stats[chunks]['count']
                percentage = fragmentation_stats[chunks]['percentage']
                print(f"    {percentage:5.1f}% of reads have {chunks} chunks ({count} reads)")
            
            # Store results
            dataset_short = dataset_name.replace('_ctc_10k', '').replace('_ctc_35k', '').replace('_ctc_55k', '').replace('_ctc_75k', '').replace('_ctc', '')
            all_results[dataset_short] = {
                'summary': summary,
                'fragmentation_stats': fragmentation_stats
            }
            
        except Exception as e:
            print(f"  Error analyzing fragmentation: {e}")
    
    return all_results

def compare_fragmentation(standard_dir, custom_dir):
    """
    Compare fragmentation between standard and custom chunk training data.
    """
    
    print("\nFRAGMENTATION COMPARISON: STANDARD vs CUSTOM")
    print("=" * 70)
    
    # Analyze both directories
    print("STANDARD 10K CHUNKS:")
    print("-" * 30)
    standard_results = analyze_ctc_directory(standard_dir, "10k")
    
    print("\nCUSTOM CHUNKS:")
    print("-" * 30)
    custom_results = analyze_ctc_directory(custom_dir, "custom")
    
    # Create comparison table
    if standard_results and custom_results:
        print("\nFRAGMENTATION COMPARISON SUMMARY:")
        print("-" * 70)
        print(f"{'Dataset':<10} {'Standard':<15} {'Custom':<15} {'Improvement':<15}")
        print(f"{'':10} {'Single%':<7} {'Multi%':<7} {'Single%':<7} {'Multi%':<7} {'Single pp':<15}")
        print("-" * 70)
        
        # Match datasets between standard and custom
        common_datasets = set(standard_results.keys()) & set(custom_results.keys())
        
        for dataset in sorted(common_datasets):
            std_single = standard_results[dataset]['summary']['single_chunk_percentage']
            std_multi = standard_results[dataset]['summary']['multi_chunk_percentage']
            cust_single = custom_results[dataset]['summary']['single_chunk_percentage']
            cust_multi = custom_results[dataset]['summary']['multi_chunk_percentage']
            
            improvement = cust_single - std_single
            
            print(f"{dataset:<10} {std_single:6.1f}  {std_multi:6.1f}  {cust_single:6.1f}  {cust_multi:6.1f}  {improvement:+6.1f}pp")

def main():
    parser = argparse.ArgumentParser(
        description='Analyze read fragmentation in CTC training data',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Analyze standard 10K chunk fragmentation
  python fragmentation_analysis.py --ctc_dir /path/to/ctc_10k_data/ --chunk_type 10k
  
  # Analyze custom chunk fragmentation
  python fragmentation_analysis.py --ctc_dir /path/to/ctc_custom_data/ --chunk_type custom
  
  # Compare standard vs custom fragmentation
  python fragmentation_analysis.py --standard_dir /path/to/ctc_10k/ --custom_dir /path/to/ctc_custom/
        """
    )
    
    parser.add_argument('--ctc_dir', 
                        help='Path to CTC training data directory')
    parser.add_argument('--chunk_type', 
                        choices=['10k', '35k', '55k', '75k', 'custom', 'standard'],
                        default='10k',
                        help='Type of chunks to analyze (default: 10k)')
    parser.add_argument('--standard_dir',
                        help='Path to standard (10K) CTC training data for comparison')
    parser.add_argument('--custom_dir',
                        help='Path to custom chunk CTC training data for comparison')
    parser.add_argument('--output',
                        help='Output CSV file for results (optional)')
    
    args = parser.parse_args()
    
    # Validate arguments
    if args.standard_dir and args.custom_dir:
        # Comparison mode
        if not os.path.exists(args.standard_dir):
            print(f"Error: Standard directory not found: {args.standard_dir}")
            return 1
        if not os.path.exists(args.custom_dir):
            print(f"Error: Custom directory not found: {args.custom_dir}")
            return 1
        
        compare_fragmentation(args.standard_dir, args.custom_dir)
        
    elif args.ctc_dir:
        # Single analysis mode
        if not os.path.exists(args.ctc_dir):
            print(f"Error: CTC directory not found: {args.ctc_dir}")
            return 1
        
        results = analyze_ctc_directory(args.ctc_dir, args.chunk_type)
        
        # Create summary table
        if results:
            print(f"\nSUMMARY TABLE:")
            print("-" * 80)
            print(f"{'Dataset':<10} {'Total':<8} {'Single%':<8} {'Multi%':<8} {'Avg Chunks':<10} {'Max Chunks':<10}")
            print("-" * 80)
            
            for dataset, data in results.items():
                summary = data['summary']
                print(f"{dataset:<10} {summary['total_reads']:<8} {summary['single_chunk_percentage']:6.1f}%  {summary['multi_chunk_percentage']:6.1f}%  {summary['average_chunks_per_read']:8.1f}  {summary['max_chunks_per_read']:8}")
        
        # Save results if output file specified
        if args.output and results:
            summary_data = []
            for dataset, data in results.items():
                summary = data['summary']
                summary['dataset'] = dataset
                summary_data.append(summary)
            
            df_output = pd.DataFrame(summary_data)
            df_output.to_csv(args.output, index=False)
            print(f"\nResults saved to: {args.output}")
            
    else:
        print("Error: Must specify either --ctc_dir or both --standard_dir and --custom_dir")
        parser.print_help()
        return 1
    
    return 0

if __name__ == "__main__":
    exit(main())