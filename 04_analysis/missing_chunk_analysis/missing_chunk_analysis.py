import pandas as pd
import os

def analyse_chunks(dataset_path, dataset_name):
    """Analyse missing chunks for a single dataset"""
    
    # Look for the TSV file
    tsv_file = os.path.join(dataset_path, "ctc_generation_log_summary.tsv")
    
    if not os.path.exists(tsv_file):
        print(f"{dataset_name}: No TSV file found")
        return None
    
    # Load the data
    df = pd.read_csv(tsv_file, sep='\t')
    print(f"{dataset_name}: Loaded {len(df)} entries")
    
    # Check if we have chunk info in read_id
    if ':' not in str(df['read_id'].iloc[0]):
        print(f"{dataset_name}: No chunk info found")
        return None
    
    # Parse the read_id format: 'read_id:chunk_id:total_chunks'
    parts = df['read_id'].str.split(':', expand=True)
    df['base_read_id'] = parts[0]
    df['chunk_id'] = parts[1].astype(int)
    df['total_chunks'] = parts[2].astype(int)
    
    # Group by base read ID
    read_summary = df.groupby('base_read_id').agg({
        'chunk_id': 'count',        # how many chunks we actually have
        'total_chunks': 'max'       # how many chunks there should be
    }).rename(columns={'chunk_id': 'observed_chunks'})
    
    # Calculate missing chunks
    read_summary['missing_chunks'] = read_summary['total_chunks'] - read_summary['observed_chunks']
    
    # Get statistics
    total_reads = len(read_summary)
    complete_reads = len(read_summary[read_summary['missing_chunks'] == 0])
    missing_reads = total_reads - complete_reads
    
    # Show chunk distribution
    chunk_counts = read_summary['total_chunks'].value_counts().sort_index()
    print(f"  Chunk distribution:")
    for chunks, count in chunk_counts.items():
        pct = count / total_reads * 100
        print(f"    {chunks} chunks: {count} reads ({pct:.1f}%)")
    
    # Show missing chunk info
    if missing_reads > 0:
        missing_summary = read_summary[read_summary['missing_chunks'] > 0]
        missing_counts = missing_summary['missing_chunks'].value_counts().sort_index()
        
        print(f"  Missing chunks:")
        for missing, count in missing_counts.items():
            pct = count / total_reads * 100
            print(f"    {missing} missing: {count} reads ({pct:.1f}%)")
    
    return {
        'name': dataset_name,
        'total_reads': total_reads,
        'complete_reads': complete_reads,
        'missing_reads': missing_reads,
        'completion_rate': complete_reads / total_reads * 100
    }

def main():
    # Define your dataset paths here
    datasets = {
        'C9_128': '/path/to/C9_128_ctc',
        'C9_256': '/path/to/C9_256_ctc', 
        'C9_512': '/path/to/C9_512_ctc',
        'C9_1024': '/path/to/C9_1024_ctc',
        'FX_160': '/path/to/FX_160_ctc',
        'FX_640': '/path/to/FX_640_ctc'
    }
    
    print("Chunk Analysis Results")
    print("=" * 50)
    
    results = []
    
    # Analyse each dataset
    for name, path in datasets.items():
        print(f"\n{name}:")
        print("-" * 30)
        
        result = analyse_chunks(path, name)
        if result:
            results.append(result)
    
    # Summary table
    if results:
        print(f"\n\nSummary")
        print("=" * 50)
        print(f"{'Dataset':<10} {'Total':<8} {'Complete':<10} {'Missing':<8} {'% Complete':<12}")
        print("-" * 50)
        
        for r in results:
            print(f"{r['name']:<10} {r['total_reads']:<8} {r['complete_reads']:<10} {r['missing_reads']:<8} {r['completion_rate']:<12.1f}")
        
        # Best and worst
        best = max(results, key=lambda x: x['completion_rate'])
        worst = min(results, key=lambda x: x['completion_rate'])
        
        print(f"\nBest: {best['name']} ({best['completion_rate']:.1f}%)")
        print(f"Worst: {worst['name']} ({worst['completion_rate']:.1f}%)")

if __name__ == "__main__":
    main()