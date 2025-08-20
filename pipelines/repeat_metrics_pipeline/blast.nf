#!/usr/bin/env nextflow
nextflow.enable.dsl=2   

params.threads = 4  // Set default to 4 threads if not provided
def required_params = [
    flank_fasta: "the path to the flank fasta file",
    reads_fq: "the path to the reads fastq file",
    publish_dir: "the path to the output directory",
    threads: "the number of threads to use",
]

// Validate required parameters
required_params.each { param, description ->
    if (params."$param" == null) {
        println "Please enter ${description} using the '--${param}' parameter"
        System.exit(0)
    }
}

// Creates a nucleotide BLAST database named 'flank' from the provided flank FASTA file.
process MAKEBLASTDB_FLANK {
    tag "makeblastdb_flank"

    input:
    path flank_fasta

    output:
    path "flank.*"

    script:
    """
    makeblastdb -in ${flank_fasta} -dbtype nucl -out flank
    """
}

// Converts the input reads from FASTQ format to FASTA format using seqtk.
process FASTQ_TO_FASTA {
    tag "fastq_to_fasta"

    input:
    path reads_fq

    output:
    path "reads.fasta"

    publishDir params.publish_dir, mode: 'copy'

    script:
    """
    seqkit fq2fa ${reads_fq} -o reads.fasta
    """
}

// Split the reads.fasta into N chunks
process SPLIT_FASTA {
    tag "split_fasta"

    input:
    path reads_fasta
    val n_chunks

    output:
    path "reads.part*.fasta"

    script:
    """
    seqkit split -p ${n_chunks} ${reads_fasta} -O "."
    """
}

// Run BLASTN on each chunk
process BLASTN_CHUNK {
    tag "blastn_chunk"

    input:
    path chunk_fasta
    path blastdb_flank

    output:
    path "*.blast.txt"

    script:
    """
    outname=\$(basename ${chunk_fasta} .fasta).blast.txt
    blastn -query ${chunk_fasta} -db flank -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" -num_threads 1 > \$outname
    """
}
// Concatenate all BLAST results
process CONCAT_RESULTS {
    tag "concat_results"

    input:
    path blast_results_files

    output:
    path "blast_results_all.txt"

    publishDir params.publish_dir, mode: 'copy'

    script:
    """
    cat ${blast_results_files} > blast_results_all.txt
    """
}

workflow blast {
    reads_fasta_channel = FASTQ_TO_FASTA(params.reads_fq)
    blastdb_flank_channel = MAKEBLASTDB_FLANK(params.flank_fasta)

    // Split the fasta file into chunks
    split_chunks = SPLIT_FASTA(reads_fasta_channel, params.threads).flatten()

    // Run BLASTN on each chunk in parallel
    blast_results = BLASTN_CHUNK(split_chunks, blastdb_flank_channel)

    // Concatenate all results
    concat_results = CONCAT_RESULTS(blast_results.collect())

    // Emit reads_fasta and concat_results
    emit: 
    reads_fasta = reads_fasta_channel
    blast_results = concat_results
}