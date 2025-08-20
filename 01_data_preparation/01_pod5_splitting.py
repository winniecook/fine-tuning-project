#!/usr/bin/env python3
"""
Split POD5 files into 90% train / 10% test sets.
Creates train and test POD5 files from input data.
"""

import argparse
import os
import random
import pod5


def main():
    parser = argparse.ArgumentParser(description='Split POD5 into 90% train / 10% test')
    
    parser.add_argument('--input', required=True, help='Input POD5 file')
    parser.add_argument('--output_dir', required=True, help='Output directory')
    parser.add_argument('--dataset_name', required=True, help='Dataset name (e.g. C9_128)')
    parser.add_argument('--seed', type=int, default=42, help='Random seed')
    
    args = parser.parse_args()
    
    print(f"Splitting {args.dataset_name}...")
    
    # Extract all read IDs
    read_ids = []
    with pod5.Reader(args.input) as reader:
        for read in reader.reads():
            read_ids.append(str(read.read_id))
    
    print(f"Found {len(read_ids)} total reads")
    
    # Shuffle and split 90/10
    random.seed(args.seed)
    random.shuffle(read_ids)
    
    split_point = int(len(read_ids) * 0.9)
    train_ids = set(read_ids[:split_point])
    test_ids = set(read_ids[split_point:])
    
    print(f"Train: {len(train_ids)} reads, Test: {len(test_ids)} reads")
    
    # Set up output paths
    os.makedirs(args.output_dir, exist_ok=True)
    train_file = os.path.join(args.output_dir, f"{args.dataset_name}_train.pod5")
    test_file = os.path.join(args.output_dir, f"{args.dataset_name}_test.pod5")
    
    # Create train and test files
    train_count = 0
    test_count = 0
    
    with pod5.Reader(args.input) as reader:
        with pod5.Writer(train_file) as train_writer, pod5.Writer(test_file) as test_writer:
            for read in reader.reads():
                read_id = str(read.read_id)
                if read_id in train_ids:
                    train_writer.add_read(read)
                    train_count += 1
                elif read_id in test_ids:
                    test_writer.add_read(read)
                    test_count += 1
    
    print(f"Created {train_file} with {train_count} reads")
    print(f"Created {test_file} with {test_count} reads")


if __name__ == '__main__':
    main()