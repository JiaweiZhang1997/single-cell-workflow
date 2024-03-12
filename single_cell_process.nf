// nextflow_seurat_workflow.nf

params {
    patient_id = "4261"
    current_path = "/share/home/xxwang/data/phs002748.v1.p1/RNA_TIL"

    transcriptome_ref = "/share/home/xxwang/software/cellranger-6.1.2/refdata-gex-GRCh38-2020-A"
    gex_output_prefix = "${params.patient_id}_GEX"
    
    tcr_reference = "/share/home/xxwang/software/cellranger-6.1.2/refdata-cellranger-vdj-GRCh38-alts-ensembl-5.0.0"
    tcr_output_prefix = "${params.patient_id}_TCR"

    tcr_react = "./react_TCR_list.txt"

    r_seurat = "./process_seurat.R"
    r_QCandanalyze = "./qc_and_analyze.R"
    r_addTCR = "./addTCRMetadata.R"
}

process runCellRangerCount {
    executor "local"
    cpus 3
    container "cellranger:6.1.2"

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
    container "cellranger:6.1.2"

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

process processSeurat {
    executor "local"
    container "rocker/r-ver:4.1.2"

    input:
    path filtered_matrix from cellranger_output.flatten().map { dir -> "${dir}/filtered_feature_bc_matrix" }

    output:
    file "seurat_obj.rds" into seurat_objs

    script:
    """
    Rscript "${params.r_seurat}" \\
        --filtered-matrix "${filtered_matrix}" \\
        --sample-name "${params.patient_id}"
    """
}

process qcAndAnalyze {
    executor "local"
    container "rocker/r-ver:4.1.2" 
    input:
    file seurat_obj from seurat_objs

    output:
    file "${sample_name}_QC_befor.pdf"
    file "${sample_name}_QC_after.pdf"
    file "qc_filtered_seurat_obj.rds"

    script:
    """
    Rscript "${params.r_QCandanalyze}" \\
        --seurat-obj "${seurat_obj}" \\
        --sample-name "${params.patient_id}" \\
        --mt-thresh 10 \\
        --max-features 2500
    """
}

process addTCRMetadata {
    executor "local"
    container "rocker/r-ver:4.1.2" 

    input:
    file seurat_obj from qc_results.out.seurat_obj
    path tcr_annotations from "${params.patient_id}_TCR/outs/filtered_contig_annotations.csv.format.csv"
    file react_list from "./react_TCR_list.txt"

    output:
    file "${params.patient_id}_QC_filtered_seurat_obj_with_TCR.rds"

    script:
    """
    Rscript "${params.r_addTCR}" \\
        --tcr_react "${params.tcr_react}"
        --tcr_dir "./${params.patient_id}_TCR/outs/filtered_contig_annotations.csv.format.csv"
    """

}

workflow {
    cellranger_results = Channel.fromPath("${fastqs_dir}/*.fastq.gz").map { fastqs ->
        runCellRangerCount(fastqs)
    }

    seurat_objs = cellranger_results.flatMap { result ->
        processSeurat(result)
    }
    
    // 可选：收集并查看Seurat对象的输出路径
    seurat_objs.view()

    qc_results = seurat_objs.map { seurat_obj_path ->
        sample_name = params.patient_id
        qcAndAnalyze(seurat_obj: seurat_obj_path)
    }

    qc_results.view()

    cellranger_vdj_results = Channel.fromPath("${tcr_fastqs_dir}/*.fastq.gz").map { fastqs ->
        runCellRangerVDJ(fastqs)
    }

    // 合并TCR元数据到已进行QC和分析的Seurat对象
    final_seurat_objs = qc_results.flatMap { result ->
        addTCRMetadata(
            seurat_obj: result.out.seurat_obj,
            tcr_annotations: "${params.patient_id}_TCR/outs/filtered_contig_annotations.csv.format.csv",
            react_list: "./react_TCR_list.txt"
        )
    }
}


