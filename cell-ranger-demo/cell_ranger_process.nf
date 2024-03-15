#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.patient_id = "4261"
params.current_path = "/cache/jwzhang/phs002748.v1.p1/RNA_TIL"
params.transcriptome_ref = "/cache/jwzhang/phs002748.v1.p1/refdata-gex-GRCh38-2020-A"

params.samples = [
    [params.patient_id, params.current_path, params.transcriptome_ref]
]

process runCellRangerCount {
    executor "local"
    cpus 3
    container "jiaweizhang1997/cellranger:8.0.0"
    publishDir "${params.current_path}/output", mode: 'move'

    input:
    tuple val(patient_id), path(current_path), path(transcriptome_ref)

    output:
    path("${current_path}/output")

    shell:
    """
    mkdir -p ${current_path}/output
    /opt/cellranger/cellranger-8.0.0/cellranger count \\
        --id "${patient_id}_GEX" \\
        --sample ${patient_id}_FRTU_TIL_GEX \\
        --fastqs ${current_path}/${patient_id} \\
        --transcriptome ${transcriptome_ref} \\
        --disable-ui \\
        --localcores ${task.cpus} \\
        --nosecondary \\
        --create-bam false
    """
}

workflow {
    fastq_pairs_ch = Channel.from(params.samples)
        .map { sample_info -> tuple(sample_info[0], sample_info[1], sample_info[2]) }

    cellranger_results = runCellRangerCount(fastq_pairs_ch)
}