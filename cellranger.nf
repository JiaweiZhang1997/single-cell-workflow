#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.patient_id = "4261"
params.current_path = "/cache/jwzhang/phs002748.v1.p1/RNA_TIL"
params.transcriptome_ref = "/cache/jwzhang/phs002748.v1.p1/refdata-gex-GRCh38-2020-A"
params.tcr_reference = "/cache/jwzhang/phs002748.v1.p1/refdata-cellranger-vdj-GRCh38-alts-ensembl-5.0.0"

params.samples = [
    [params.patient_id, params.current_path, params.transcriptome_ref, params.tcr_reference]
]

process runCellRangerCount {
    executor "local"
    cpus 4
    container "jiaweizhang1997/cellranger:8.0.0"
    publishDir "${params.current_path}/output", mode: 'copy'

    input:
    tuple val(patient_id), path(current_path), path(transcriptome_ref)

    output:
    file("${patient_id}_GEX/outs/filtered_feature_bc_matrix/*.gz")

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

process runCellRangerVDJ {
    executor "local"
    cpus 4
    container "jiaweizhang1997/cellranger:8.0.0"
    publishDir "${params.current_path}/output", mode: 'copy'

    input:
    tuple val(patient_id), path(current_path), path(tcr_reference)

    output:
    file("${patient_id}_TCR/outs/filtered_contig_annotations.csv")

    shell:
    """
    /opt/cellranger/cellranger-8.0.0/cellranger vdj \\
        --id ${patient_id}_TCR \\
        --sample ${patient_id}_FRTU_TIL_TCR \\
        --fastqs ${current_path}/${patient_id} \\
        --reference ${tcr_reference} \\
        --disable-ui \\
        --localcores ${task.cpus}
    """
}

workflow {
    fastq_Count_ch = Channel.from(params.samples)
        .map { sample_info -> tuple(sample_info[0], sample_info[1], sample_info[2])}

    fastq_VDJ_ch = Channel.from(params.samples)
        .map { sample_info -> tuple(sample_info[0], sample_info[1], sample_info[3])}

    cellranger_Count = runCellRangerCount(fastq_Count_ch)
    cellranger_VDJ = runCellRangerVDJ(fastq_VDJ_ch)

}
