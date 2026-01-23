#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.threads = 8


process RUN_SHOVILL {
    tag "$sample_id"
    publishDir "${params.outdir}/${sample_id}", mode: 'copy'


    input:
    tuple val(sample_id), path(reads_1), path(reads_2)

    output:
    tuple val(sample_id), path("${sample_id}_contigs.fa"), emit: contigs
    tuple val(sample_id), path("${sample_id}_shovill_provenance.yml"), emit: provenance
    tuple val(sample_id), path("${sample_id}_shovill.log"), emit: log

    script:
    """
    printf -- "- process_name: shovill\\n"                                                 >> ${sample_id}_shovill_provenance.yml
    printf -- "  tools:\\n"                                                                  >> ${sample_id}_shovill_provenance.yml
    printf -- "    - tool_name: shovill\\n"                                                >> ${sample_id}_shovill_provenance.yml
    printf -- "      tool_version: \$(shovill --version | cut -d ' ' -f 2 | tr -d 'v')\\n" >> ${sample_id}_shovill_provenance.yml

    shovill \
        --R1 ${reads_1} \
        --R2 ${reads_2} \
	    --assembler ${params.assembler} \
        --outdir . \
        --cpus ${task.cpus} \
        --force > ${sample_id}_shovill.log 2>&1
    cp contigs.fa ${sample_id}_contigs.fa
    """
}

process pipeline_provenance {

    tag { pipeline_name + " / " + pipeline_version }

    executor 'local'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple val(session_id), val(run_name), val(pipeline_name), val(pipeline_version), val(timestamp_analysis_start)

    output:
    file("pipeline_provenance.yml")

    script:
    """
    printf -- "- pipeline_name: ${pipeline_name}\\n"                       >> pipeline_provenance.yml
    printf -- "  pipeline_version: ${pipeline_version}\\n"                 >> pipeline_provenance.yml
    printf -- "  nextflow_session_id: ${session_id}\\n"                    >> pipeline_provenance.yml
    printf -- "  nextflow_run_name: ${run_name}\\n"                        >> pipeline_provenance.yml
    printf -- "  timestamp_analysis_start: ${timestamp_analysis_start}\\n" >> pipeline_provenance.yml
    """
}

workflow {
    ch_workflow_metadata = Channel.value([
	workflow.sessionId,
	workflow.runName,
	workflow.manifest.name,
	workflow.manifest.version,
	workflow.start,
    ])
    
    ch_pipeline_provenance = pipeline_provenance(ch_workflow_metadata)

    if (params.samplesheet_input != 'NO_FILE') {	    
	    ch_fastq = Channel.fromPath(params.samplesheet_input).splitCsv(header: true).map{ it -> [it['ID'], it['R1'], it['R2']] }
	
    } else {

	    ch_fastq = Channel.fromFilePairs( params.fastq_search_path, flat: true ).map{ it -> [it[0].split('_')[0], it[1], it[2]] }.unique{ it -> it[0]}

    }
    //Ch_rd = Channel.fromFilePairs(params.reads, flat: true).map{ it -> [it[0].split('_')[0], it[1], it[2]] }.unique{ it -> it[0] }
    RUN_SHOVILL(ch_fastq)

}
