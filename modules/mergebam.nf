#!/usr/bin/env nextflow
nextflow.enable.dsl=2
// Define parameters
// params.input_dir = "/home/chbope/extension/testdata/single_bam_folder" // Base directory containing Sample_ID folders
// params.merge_bam_dir = "/home/chbope/extension/testdata/single_bam_folder/merge_bam" // Directory to store merged BAM files
// params.occ_bam_dir = "/home/chbope/extension/testdata/single_bam_folder/occ_bam"
// params.threads = 12 // Number of threads for samtools operations
// params.roi_bed = file("/home/chbope/extension/data/reference/OCC.protein_coding.bed")
params.bam_sample_id_file = file("/home/chbope/extension/testdata/sample_ids_bam.txt")

///
//params.sample_flowcell_id = file("/data/pipeline/data/200_GBMs_samples10.txt")

// Load sample information from a file
//def sample_info = [:]
//file("/data/pipeline/results/trash/sample_info.txt").splitEachLine('\t') { fields ->
//    if (fields.size() == 2) {
//        sample_info[fields[0]] = fields[1] // Map Sample_ID to Flow_cell_ID
//    }
//}

def sample_info = [:]
//file(params.bam_sample_id_file).splitEachLine('\t') { fields ->
file(params.bam_sample_id_file).splitEachLine('\t') { fields ->
    if (fields.size() == 2) {
        if (fields[1] != '') { // Check if Flow Cell ID is not empty
            sample_info[fields[0]] = fields[1] // Map Sample_ID to Flow_cell_ID
            println "Sample_id ${fields[0]} with Flow Cell ID: ${fields[1]}" 
        } else {
            println "Skipping ${fields[0]}: Flow Cell ID is empty"
        }
    }
}

process merge_bam_files {
    cpus 12
    memory '12 GB'
    label 'samtools_merge'
    publishDir "${params.merge_bam_dir}", mode: 'copy'

    input:
    tuple val(sample_id), path(bam_files) // Expecting a tuple of sample_id and a list of BAM file paths

    output:
    tuple val(sample_id), path("${sample_id}.merge.bam"), path("${sample_id}.merge.bam.bai"), emit: mergebamout

    script:
    """
    echo "Merging BAM files for sample: ${sample_id}"
    for bam in ${bam_files}; do
        echo "Processing file: \${bam}" # Print each BAM file being merged
    done
    samtools merge -@ ${params.threads} -f ${sample_id}.merge.bam ${bam_files.join(' ')}
    samtools index -@ ${params.threads} ${sample_id}.merge.bam
    """
}

process extract_roi {
    cpus 12
    memory '12GB'
    label 'roi_extraction'
    publishDir "${params.occ_bam_dir}", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai), path(roi_bed)

    output:
    path("${sample_id}.occ.bam")
    path("${sample_id}.occ.bam.bai")

    script:
    """
    intersectBed -a ${bam} -b ${roi_bed} > ${sample_id}.occ.bam
    samtools index -@ ${params.threads} ${sample_id}.occ.bam
    """
}

workflow mergebam {
    main:
        // Store start time
        start_time = new Date()

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
                tuple(sample_id, bam, bai, params.roi_bed)
            }
        
        extract_roi_results = extract_roi(grouped_files_out)

    emit:
        // Emit the merged BAM files and a completion flag
        merged_bams = merge_bam_files.out.mergebamout
        complete = Channel.of(true)
}

