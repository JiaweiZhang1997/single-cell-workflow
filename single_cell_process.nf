// nextflow_seurat_workflow.nf

params {
    patient_id = '4261'
    fastqs_dir = "./${params.patient_id}_GEX/outs/fastq_path"
    transcriptome_ref = '/share/home/xxwang/software/cellranger-6.1.2/refdata-gex-GRCh38-2020-A'
    output_prefix = "${params.patient_id}_GEX"
    r_script = './process_seurat.R' // 假设存在一个名为process_seurat.R的R脚本文件
}

process runCellRangerCount {
    executor 'local'
    cpus 3
    container 'cellranger:6.1.2'

    input:
    path(fastqs) from file("${fastqs_dir}/*.fastq.gz")

    output:
    dir("${output_prefix}") into cellranger_output

    shell:
    """
    cellranger count \\
        --id ${output_prefix} \\
        --sample ${params.patient_id}_FRTU_TIL_GEX \\
        --fastqs ${fastqs} \\
        --transcriptome ${transcriptome_ref} \\
        --disable-ui \\
        --localcores ${task.cpus}
    """
}

process processSeurat {
    executor 'local'
    container 'rocker/r-ver:4.1.2' // 假设使用包含R和必要包的Docker镜像

    input:
    path filtered_matrix from cellranger_output.flatten().map { dir -> "${dir}/filtered_feature_bc_matrix" }

    output:
    file "seurat_obj.rds" into seurat_objs

    script:
    """
    Rscript '${params.r_script}' \\
        --filtered-matrix '${filtered_matrix}' \\
        --sample-name '${params.patient_id}'
    """
}

process qcAndAnalyze {
    executor 'local'
    container 'rocker/r-ver:4.1.2' // 假设使用包含R和必要包的Docker镜像

    input:
    file seurat_obj from seurat_objs

    output:
    file "${sample_name}_QC_befor.pdf"
    file "${sample_name}_QC_after.pdf"
    file "qc_filtered_seurat_obj.rds"

    script:
    """
    Rscript '${workflow.projectDir}/qc_and_analyze.R' \\
        --seurat-obj '${seurat_obj}' \\
        --sample-name '${params.patient_id}' \\
        --mt-thresh 10 \\
        --max-features 2500
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

    // QC并分析处理
    qc_results = seurat_objs.map { seurat_obj_path ->
        sample_name = params.patient_id
        qcAndAnalyze(seurat_obj: seurat_obj_path)
    }

    // 可选：收集并查看QC结果文件路径
    qc_results.view()
}


