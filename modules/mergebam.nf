#!/usr/bin/env nextflow
nextflow.enable.dsl=2

//---------------------------------------------------------------------
// Mergebam Pipeline: Merges multiple BAM files per sample and extracts regions of interest
//---------------------------------------------------------------------

// Load sample information from file mapping Sample_ID to Flow_cell_ID
def sample_info = [:]
file(params.bam_sample_id_file).splitEachLine(~/\t| /) { fields ->
    if (fields.size() == 2) {
        if (fields[1] != '') { // Check if Flow Cell ID is not empty
            sample_info[fields[0]] = fields[1] // Map Sample_ID to Flow_cell_ID
        } else {
            println "Skipping ${fields[0]}: Flow Cell ID is empty"
        }
    }
}

//---------------------------------------------------------------------
// Process Definitions
//---------------------------------------------------------------------

process merge_bam_files {
    label 'samtools_merge'
    publishDir "${params.merge_bam_dir}", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), path(bam_files) // Expecting a tuple of sample_id and a list of BAM file paths

    output:
    tuple val(sample_id), path("${sample_id}.merged.bam"), path("${sample_id}.merged.bam.bai"), emit: mergebamout

    script:
    """
    ulimit -n 50000
    echo "Merging BAM files for sample: ${sample_id}"
    printf '%s\\n' ${bam_files.join(' ')} > bam_list.txt
    echo "Total BAM files: \$(wc -l < bam_list.txt)"
    samtools merge -@ ${params.threads} -f -b bam_list.txt ${sample_id}.merged.bam
    samtools index -@ ${params.threads} ${sample_id}.merged.bam
    """
}

process extract_roi {
    label 'roi_extraction'
    publishDir "${params.roi_bam_dir}", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai), path(roi_bed)

    output:
    tuple val(sample_id), path("${sample_id}.roi.bam"), path("${sample_id}.roi.bam.bai"), emit: roi_bam

    script:
    """
    set -euo pipefail
    samtools view -@ ${params.threads} -b -L ${roi_bed} ${bam} \\
    | samtools sort -@ ${params.threads} -o ${sample_id}.roi.bam
    samtools index -@ ${params.threads} ${sample_id}.roi.bam
    """
}

//---------------------------------------------------------------------
// Workflow Definition
//---------------------------------------------------------------------

workflow mergebam {
    main:
        // Store start time
        start_time = new Date()

        // Print sample information when workflow is executed
        sample_info.each { sample_id, flow_cell_id ->
            println "Sample_id ${sample_id} with Flow Cell ID: ${flow_cell_id}"
        }

        // Process all samples listed in sample_info.txt
        bam_files = Channel.fromPath("${params.input_dir}/**/bam_pass/*.bam")
            .map { bam ->
                def sample_id = bam.getParent().getParent().getParent().getBaseName()
                def flow_cell_id = sample_info[sample_id]

                if (flow_cell_id) {
                    return tuple(sample_id, bam.toAbsolutePath())
                } else {
                    return null
                }
            }
            .filter { it != null }

        // Create unique_bam_files channel
        unique_bam_files = bam_files
            .groupTuple()
            .map { sample_id, bam_files ->
                def unique_bams = [:]
                bam_files.each { bam ->
                    def base_name = bam.getBaseName()
                    unique_bams[base_name] = bam
                }
                return tuple(sample_id, unique_bams.values().toList())
            }

        // Run merge_bam_files process
        merged_results = merge_bam_files(unique_bam_files)

        // Run extract_roi process
        grouped_files_out = merge_bam_files.out.mergebamout
            .map { sample_id, bam, bai ->
                tuple(sample_id, bam, bai, file(params.roi_bed))
            }

        extract_roi_results = extract_roi(grouped_files_out)

    emit:
        // Emit the merged BAM files and a completion flag
        merged_bams = merge_bam_files.out.mergebamout
        roi_bams = extract_roi.out.roi_bam
        complete = Channel.of(true)
}

