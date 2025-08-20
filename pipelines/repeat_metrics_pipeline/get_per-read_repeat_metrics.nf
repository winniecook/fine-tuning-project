#!/usr/bin/env nextflow
nextflow.enable.dsl=2  

include {blast} from './blast.nf'

params.tolerance = 5 // Set default tolerance to 5 if not provided
def required_params = [
    repeat_expansion_unit: "the repeat expansion unit kmer",
]

// Validate required parameters
required_params.each { param, description ->
    if (params."$param" == null) {
        println "Please enter ${description} using the '--${param}' parameter"
        System.exit(0)
    }
}

process GET_METRICS_FOR_DUAL_MATCH_READS {

    tag "get_metrics_for_dual_match_reads"

    input:
    path reads_fasta
    path blast_output

    output:
    path "per-read_repeat_metrics.txt"
    path "repeat_metric_summary.txt"  

    publishDir params.publish_dir, mode: 'copy'

    script:
    """
    get_metrics_for_dual_match_reads.py -r ${reads_fasta} -b ${blast_output} -o per-read_repeat_metrics.txt -s repeat_metric_summary.txt -k ${params.repeat_expansion_unit} -t ${params.tolerance}
    """
}

workflow {
    blast()  // Run everything from blast.nf
    GET_METRICS_FOR_DUAL_MATCH_READS(blast.out.reads_fasta, blast.out.blast_results)
}