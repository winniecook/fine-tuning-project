#!/usr/bin/env python3
"""
Simple reference length analysis for CTC training data.
Analyzes mean reference lengths to validate sequence context hypothesis.
"""

import numpy as np
import os
import argparse

def analyze_reference_lengths(base_path, datasets):
    """
    Analyze mean reference lengths from CTC training data.
    """
    
    print("MEAN REFERENCE LENGTHS")
    print("=" * 40)
    
    results = []
    
    for dataset in datasets:
        dataset_path = os.path.join(base_path, dataset)
        ref_lengths_path = os.path.join(dataset_path, "reference_lengths.npy")
        
        if os.path.exists(ref_lengths_path):
            ref_lengths = np.load(ref_lengths_path)
            mean_length = np.mean(ref_lengths)
            dataset_short = dataset.replace('_ctc_10k', '').replace('_ctc', '')
            print(f"{dataset_short:<8}: {mean_length:.1f}")
            
            results.append({
                'dataset': dataset_short,
                'mean_length': mean_length,
                'total_chunks': len(ref_lengths)
            })
        else:
            print(f"{dataset}: File not found")
    
    return results

def main():
    parser = argparse.ArgumentParser(description='Analyze reference lengths from CTC training data')
    
    parser.add_argument('--base_path', required=True,
                        help='Base path to CTC training data')
    parser.add_argument('--datasets', nargs='+', required=True,
                        help='List of dataset directories to analyze')
    
    args = parser.parse_args()
    
    if not os.path.exists(args.base_path):
        print(f"Error: Base path not found: {args.base_path}")
        return 1
    
    results = analyze_reference_lengths(args.base_path, args.datasets)
    
    if results:
        print(f"\nSummary:")
        print("-" * 30)
        for result in results:
            print(f"{result['dataset']}: {result['mean_length']:.1f} bp (from {result['total_chunks']} chunks)")
    
    return 0

if __name__ == "__main__":
    exit(main())