#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.patient_id = "4261"
params.current_path = "/cache/jwzhang/phs002748.v1.p1/RNA_TIL"
params.r_process = "/cache/jwzhang/phs002748.v1.p1/RNA_TIL/nextflow/process.R"
params.r_preprocess = "/cache/jwzhang/phs002748.v1.p1/RNA_TIL/nextflow/preprocess_tcr.R"
params.mt_thresh = "10"
params.max_features = "2500"
params.tcr_react = "/cache/jwzhang/phs002748.v1.p1/RNA_TIL/output/react_TCR_list.txt"



params.samples = [
    [params.patient_id, params.current_path, params.mt_thresh, params.max_features, params.tcr_react]
]

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

    output:
    file("${current_path}/output/${patient_id}_GEX/outs/*.rds")

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


workflow {
    preprocess_ch = Channel.from(params.samples)
        .map { sample_info -> tuple(sample_info[0], sample_info[1])}
    result_TCR = preprocess_TCR(preprocess_ch)
    r_ch = Channel.from(params.samples)
        .map { sample_info -> tuple(sample_info[0], sample_info[1] ,sample_info[2], sample_info[3], sample_info[4])}
    r_process(r_ch)
}

