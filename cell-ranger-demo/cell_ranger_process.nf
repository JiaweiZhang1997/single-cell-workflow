// nextflow_seurat_workflow.nf

params {
    patient_id = "4261"
    current_path = "/share/home/xxwang/data/phs002748.v1.p1/RNA_TIL"

    transcriptome_ref = "/share/home/xxwang/software/cellranger-6.1.2/refdata-gex-GRCh38-2020-A"
    gex_output_prefix = "${params.patient_id}_GEX"
    
    tcr_reference = "/share/home/xxwang/software/cellranger-6.1.2/refdata-cellranger-vdj-GRCh38-alts-ensembl-5.0.0"
    tcr_output_prefix = "${params.patient_id}_TCR"
}

process runCellRangerCount {
    executor "local"
    cpus 3
    container "cellranger"

    input:
    path(fastqs) from file("${fastqs_dir}/*.fastq.gz")

    output:
    dir("${gex_output_prefix}") into cellranger_output

    shell:
    """
    cellranger count \\
        --id "${params.patient_id}_GEX" \\
        --sample ${params.patient_id}_FRTU_TIL_GEX \\
        --fastqs ${params.current_path}/${params.patient_id} \\
        --transcriptome ${params.transcriptome_ref} \\
        --disable-ui \\
        --localcores ${task.cpus} \\
        --nosecondary
    """
}

process runCellRangerVDJ {
    executor "local"
    cpus 3
    container "cellranger"

    input:
    path(fastqs) from file("${tcr_fastqs_dir}/*.fastq.gz")

    output:
    dir("${params.patient_id}_TCR") into cellranger_vdj_output

    shell:
    """
    cellranger vdj \\
        --id ${params.patient_id}_TCR \\
        --sample ${params.patient_id}_FRTU_TIL_TCR \\
        --fastqs ${params.current_path}/${params.patient_id} \\
        --reference ${tcr_reference} \\
        --disable-ui \\
        --localcores ${task.cpus}
    """
}


workflow {
    cellranger_results = Channel.fromPath("${fastqs_dir}/*.fastq.gz").map { fastqs ->
        runCellRangerCount(fastqs)
    }

    cellranger_vdj_results = Channel.fromPath("${tcr_fastqs_dir}/*.fastq.gz").map { fastqs ->
        runCellRangerVDJ(fastqs)
    }
}


