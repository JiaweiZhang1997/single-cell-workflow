// nextflow_seurat_workflow.nf

params {
    patient_id = '4261'
    fastqs_dir = "./${params.patient_id}_GEX/outs/fastq_path"
    transcriptome_ref = '/share/home/xxwang/software/cellranger-6.1.2/refdata-gex-GRCh38-2020-A'
    output_prefix = "${params.patient_id}_GEX"
    r_script = './process_seurat.R'
    tcr_fastqs_dir = "./${params.patient_id}_TCR"
    tcr_reference = '/share/home/xxwang/software/cellranger-6.1.2/refdata-cellranger-vdj-GRCh38-alts-ensembl-5.0.0'
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
    container 'rocker/r-ver:4.1.2'

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
    container 'rocker/r-ver:4.1.2' 
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

process runCellRangerVDJ {
    executor 'local'
    cpus 3
    container 'cellranger:6.1.2'

    input:
    path(fastqs) from file("${tcr_fastqs_dir}/*.fastq.gz")

    output:
    dir("${params.patient_id}_TCR") into cellranger_vdj_output

    shell:
    """
    cellranger vdj \\
        --id ${params.patient_id}_TCR \\
        --sample ${params.patient_id}_FRTU_TIL_TCR \\
        --fastqs ${fastqs} \\
        --reference ${tcr_reference} \\
        --disable-ui \\
        --localcores ${task.cpus}
    """
}

process addTCRMetadata {
    executor 'local'
    container 'rocker/r-ver:4.1.2' 

    input:
    file seurat_obj from qc_results.out.seurat_obj
    path tcr_annotations from "${params.patient_id}_TCR/outs/filtered_contig_annotations.csv.format.csv"
    file react_list from './react_TCR_list.txt'

    output:
    file "${params.patient_id}_QC_filtered_seurat_obj_with_TCR.rds"

    script:
    """
    Rscript -e '
        library(dplyr)

        tcr_4261 <- read.csv(\"${tcr_annotations}\") %>%
            select(cdr3_aa1, cdr3_aa2, CTaa, barcode_raw)

        tcr_react_list <- read.table(\"${react_list}\", sep = \"\\t\", header = TRUE) %>%
            select(CD4.CD8, CDR3A_B, Archival.Prospective, Tumor.ID)

        metadata_add_TCR <- function(seurat_obj, tcr_df, tcr_react_df) {
            patient_id <- seurat_obj@meta.data$orig.ident[1]
            temp_mtx <- tcr_df %>%
                remove_rownames() %>%
                column_to_rownames(var = "barcode_raw")
            temp_combined <- merge(seurat_obj@meta.data, temp_mtx, by = 0, sort = FALSE, all.x = TRUE)
            temp_combined <- merge(temp_combined, subset(tcr_react_df, Tumor.ID == patient_id), by.x = "CTaa", by.y = "CDR3A_B", sort = FALSE, all.x = TRUE) %>%
                remove_rownames() %>%
                column_to_rownames(var = "Row.names")
            raw_barcode <- row.names(seurat_obj@meta.data)
            temp_combined <- temp_combined[raw_barcode, ]
            seurat_obj@meta.data <- temp_combined
            return(seurat_obj)
        }

        sc_4261_qc <- readRDS(\"${seurat_obj}\")
        sc_4261_qc_TCR <- metadata_add_TCR(sc_4261_qc, tcr_4261, tcr_react_list)
        saveRDS(sc_4261_qc_TCR, file=\"${params.patient_id}_QC_filtered_seurat_obj_with_TCR.rds\")
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
            react_list: './react_TCR_list.txt'
        )
    }
}


