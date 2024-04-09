#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.patient_id = "4261"
params.current_path = "/cache/jwzhang/phs002748.v1.p1/RNA_TIL"
params.transcriptome_ref = "/cache/jwzhang/phs002748.v1.p1/refdata-gex-GRCh38-2020-A"
params.tcr_reference = "/cache/jwzhang/phs002748.v1.p1/refdata-cellranger-vdj-GRCh38-alts-ensembl-5.0.0"
params.r_process = "/cache/jwzhang/phs002748.v1.p1/RNA_TIL/nextflow/process.R"
params.r_preprocess = "/cache/jwzhang/phs002748.v1.p1/RNA_TIL/nextflow/preprocess_tcr.R"
params.mt_thresh = "10"
params.max_features = "2500"
params.tcr_react = "/cache/jwzhang/phs002748.v1.p1/RNA_TIL/output/react_TCR_list.txt"
// params.cellranger = "false"
// params.process_TCR = "false"
// params.end = "false"


params.samples = [
    [params.patient_id, params.current_path, params.transcriptome_ref, params.tcr_reference, params.mt_thresh, params.max_features, params.tcr_react]
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

process preprocess_TCR {
    executor "local"
    cpus 4
    container "jiaweizhang1997/single-cell-r:4.3.3"
    publishDir "${params.current_path}/output", mode: 'copy'

    input:  
    tuple val(patient_id), path(current_path)

    output:
    file("${current_path}/output/${patient_id}_TCR/outs/filtered_contig_annotations.csv.format.csv")

    script:
    """
    Rscript ${params.r_preprocess} ${current_path}/output/${patient_id}_TCR/outs/filtered_contig_annotations.csv
    """
}

process r_process {
    executor "local"
    cpus 4
    container "jiaweizhang1997/single-cell-r:4.3.3"
    publishDir "${params.current_path}/output", mode: 'copy'

    input:
    tuple val(patient_id), path(current_path), val(mt_thresh), val(max_features), path(tcr_react)

    // output:
    // file("*.pdf")
    // when:
    // // params.cellranger == "true" && params.process_TCR == "true"
    // params.process_TCR == "true"

    script:
    """
    Rscript ${params.r_process} \\
    ${current_path}/output/${patient_id}_GEX/outs \\
    ${patient_id} \\
    ${mt_thresh} \\
    ${max_features} \\
    ${params.tcr_react} \\
    ${current_path}/output/${patient_id}_TCR/outs/filtered_contig_annotations.csv.format.csv
    """
}


// preprocess_TCR.after(runCellRangerVDJ)
// r_process.after(runCellRangerCount, preprocess_TCR)

workflow {
    fastq_Count_ch = Channel.from(params.samples)
        .map { sample_info -> tuple(sample_info[0], sample_info[1], sample_info[2])}

    fastq_VDJ_ch = Channel.from(params.samples)
        .map { sample_info -> tuple(sample_info[0], sample_info[1], sample_info[3])}

    cellranger_Count = runCellRangerCount(fastq_Count_ch)
    cellranger_VDJ = runCellRangerVDJ(fastq_VDJ_ch)
    // preprocess_ch = Channel.from(params.samples)
    //     .map { sample_info -> tuple(sample_info[0], sample_info[1])}
    // result_TCR = preprocess_TCR(preprocess_ch)
    // r_ch = Channel.from(params.samples)
    //     .map { sample_info -> tuple(sample_info[0], sample_info[1] ,sample_info[4], sample_info[5], sample_info[6])}
    // r_process(r_ch)
}


// workflow flow1{
//     fastq_Count_ch = Channel.from(params.samples)
//         .map { sample_info -> tuple(sample_info[0], sample_info[1], sample_info[2])}

//     fastq_VDJ_ch = Channel.from(params.samples)
//         .map { sample_info -> tuple(sample_info[0], sample_info[1], sample_info[3])}

//     cellranger_Count = runCellRangerCount(fastq_Count_ch)
//     cellranger_VDJ = runCellRangerVDJ(fastq_VDJ_ch)
//     params.cellranger = "true"
//     }

// workflow flow2{
//     preprocess_ch = Channel.from(params.samples)
//         .map { sample_info -> tuple(sample_info[0], sample_info[1])}
//     result_TCR = preprocess_TCR(preprocess_ch)
//     params.process_TCR = "true"
//     }

// workflow flow3{
//     r_ch = Channel.from(params.samples)
//         .map { sample_info -> tuple(sample_info[0], sample_info[1] ,sample_info[4], sample_info[5], sample_info[6])}
//     r_process(r_ch)
//     params.end = "true"

//     }