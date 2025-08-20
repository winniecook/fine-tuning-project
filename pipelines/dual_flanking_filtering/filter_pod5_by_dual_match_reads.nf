#!/usr/bin/env nextflow
nextflow.enable.dsl=2  

include {blast} from './blast.nf'

params.tolerance = 10 // Set default tolerance to 10 if not provided
def required_params = [
    pod5_dir: "the path to the pod5 directory to filter",
]

// Validate required parameters
required_params.each { param, description ->
    if (params."$param" == null) {
        println "Please enter ${description} using the '--${param}' parameter"
        System.exit(0)
    }
}

process GET_DUAL_MATCH_READ_IDS {

    tag "get_dual_match_read_ids"

    input:
    path reads_fasta
    path blast_output

    output:
    path "dual_match_read_ids.txt"

    publishDir params.publish_dir, mode: 'copy'

    script:
    """
    get_read-id-file_for_dual_match_reads.py -r ${reads_fasta} -b ${blast_output} -o dual_match_read_ids.txt -t ${params.tolerance}
    """
}

process FILTER_POD5_BY_DUAL_MATCH_READS {

    tag "filter_pod5_by_dual_match_reads"

    input:
    path dual_match_read_ids

    output:
    path "filtered_pod5.pod5"

    publishDir params.publish_dir, mode: 'copy'

    script:
    """
    pod5 filter -t ${params.threads} -i ${dual_match_read_ids} -M -D ${params.pod5_dir} -o filtered_pod5.pod5
    """
}

workflow {
    blast()  // Run everything from blast.nf
    read_id_ch = GET_DUAL_MATCH_READ_IDS(blast.out.reads_fasta, blast.out.blast_results)
    FILTER_POD5_BY_DUAL_MATCH_READS(read_id_ch)
}