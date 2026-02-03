process fastp {

    tag { sample_id }

    //publishDir "${params.outdir}/${sample_id}", pattern: "${sample_id}*.trim.fastq.gz", mode:'copy'
    publishDir "${params.outdir}/${sample_id}", pattern: "${sample_id}_fastp.*", mode:'copy'

    input:
    tuple val(sample_id), path(reads_1), path(reads_2)

    output:
    tuple val(sample_id), path("${sample_id}_trim_R1.fastq.gz"), path("${sample_id}_trim_R2.fastq.gz"), emit: trimmed_reads
    tuple val(sample_id), path("${sample_id}_fastp.json"), emit: json
    tuple val(sample_id), path("${sample_id}_fastp_provenance.yml"), emit: provenance
    tuple val(sample_id), path("${sample_id}_fastp.csv")            , emit: metrics
    tuple val(sample_id), path("${sample_id}_fastp.json")           , emit: report_json
    tuple val(sample_id), path("${sample_id}_fastp.html")           , emit: report_html

    script:
    """      
    printf -- "- process_name: fastp\\n" > ${sample_id}_fastp_provenance.yml
    printf -- "  tool_name: fastp\\n  tool_version: \$(fastp --version 2>&1 | cut -d ' ' -f 2)\\n" >> ${sample_id}_fastp_provenance.yml

    fastp \
      -t ${task.cpus} \
      -i ${reads_1} \
      -I ${reads_2} \
      -o ${sample_id}_trim_R1.fastq.gz \
      -O ${sample_id}_trim_R2.fastq.gz \
      --detect_adapter_for_pe \
      --trim_poly_g \
      --overrepresentation_analysis \
      --report_title "fastp report: ${sample_id}" \
      --json ${sample_id}_fastp.json \
      --html ${sample_id}_fastp.html
    
    fastp_json_to_csv.py -s ${sample_id} ${sample_id}_fastp.json > ${sample_id}_fastp.csv

    """
}

process cutadapter {

    tag { sample_id }

    //publishDir "${params.outdir}/${sample_id}", pattern: "${sample_id}*.out.fastq.gz", mode:'copy'
    publishDir "${params.outdir}/${sample_id}", pattern: "${sample_id}.cutadapt.log", mode:'copy'

    input:
    tuple val(sample_id), path(reads_1), path(reads_2)

    output:
    tuple val(sample_id), path("${sample_id}_out_R1.fastq.gz"), path("${sample_id}_out_R2.fastq.gz"), emit: out_reads
    path("${sample_id}.cutadapt.log"), emit: log
    tuple val(sample_id), path("${sample_id}_cutadapt_provenance.yml"), emit: provenance

    script:
    """
    printf -- "- process_name: cutadapt\\n" > ${sample_id}_cutadapt_provenance.yml
    printf -- "  tool_name: cutadapt\\n  tool_version: \$(cutadapt --version 2>&1 | cut -d ' ' -f 2)\\n" >> ${sample_id}_cutadapt_provenance.yml

    cutadapt \
      -j ${task.cpus} \
      -g GGGGGGGGGGGGGGGGGGGG -G GGGGGGGGGGGGGGGGGGGG \
      -n 10 \
      --nextseq-trim=20 -m 1 \
      -o ${sample_id}_out_R1.fastq.gz \
      -p ${sample_id}_out_R2.fastq.gz \
      ${reads_1}\
      ${reads_2}\
      > ${sample_id}.cutadapt.log
      
    """
}

