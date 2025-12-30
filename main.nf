#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.threads = 8


process RUN_SHOVILL {
    tag "$sample_id"
    publishDir "${params.outdir}", mode: 'copy'


    input:
    tuple val(sample_id), path(reads_1), path(reads_2)

    output:
    tuple val(sample_id), path("${sample_id}_contigs.fa")

    script:
    """
    shovill \
        --R1 ${reads_1} \
        --R2 ${reads_2} \
	    --assembler ${params.assembler} \
        --outdir . \
        --cpus ${task.cpus} \
        --force
    cp contigs.fa ${sample_id}_contigs.fa
    """
}

workflow {
    if (params.samplesheet_input != 'NO_FILE') {	    
	    ch_fastq = Channel.fromPath(params.samplesheet_input).splitCsv(header: true).map{ it -> [it['ID'], it['R1'], it['R2']] }
	
    } else {

	    ch_fastq = Channel.fromFilePairs( params.fastq_search_path, flat: true ).map{ it -> [it[0].split('_')[0], it[1], it[2]] }.unique{ it -> it[0]}

    }
    //Ch_rd = Channel.fromFilePairs(params.reads, flat: true).map{ it -> [it[0].split('_')[0], it[1], it[2]] }.unique{ it -> it[0] }
    RUN_SHOVILL(ch_fastq)
}
