#!/usr/bin/env nextflow
nextflow.enable.dsl=2


// Include pipeline modules conditionally
if (params.run_mode_mergebam || params.run_mode_order) {
    include { mergebam } from './modules/mergebam.nf'
}
if (params.run_mode_epi2me || params.run_mode_order || params.run_mode_epiannotation) {
    include { epi2me }   from './modules/epi2me.nf'
}
if (params.run_mode_annotation || params.run_mode_order || params.run_mode_epiannotation) {
    include { annotation } from './modules/annotation.nf'
}

workflow {
    if (params.run_mode_order) {
        // Set run_mode_annotation to 'all' when using run_mode_order to ensure all analyses run
        params.run_mode_annotation = 'all'
        
        log.info """
        Running pipelines sequentially in strict order:
        1. Mergebam Pipeline
        2. Epi2me Pipeline
        3. Annotation Pipeline
        """

        // Step 1: Run mergebam
        log.info "=== Starting Mergebam Pipeline ==="
        def mergebam_results = mergebam()

        // Step 2: Run epi2me (add reference genome files to channel)
        log.info "=== Starting Epi2me Pipeline ==="
        def epi2me_input = mergebam_results.merged_bams
            .map { sample_id, bam, bai ->
                tuple(
                    sample_id,
                    bam,
                    bai,
                    file(params.reference_genome),
                    file(params.reference_genome_bai)
                )
            }
        // epi2me runs extract_roi internally; roi_bam is emitted back
        def epi2me_results = epi2me(epi2me_input)

        // Step 3: Prepare analysis input using epi2me outputs
        def annotation_input = epi2me_results.results
            .join(epi2me_results.roi_bam, by: 0)
            .map { joined ->
                println "DEBUG: joined: ${joined}"
                def sample_id = joined[0]
                def modkit = joined[1]
                def segs_bed = joined[2]
                def bins_bed = joined[3]
                def segs_vcf = joined[4]
                def rds_file = joined[5]
                def sv = joined[6]
                def roi_bam = joined[7]
                def roi_bai = joined[8]

                def ref = file(params.reference_genome)
                def ref_bai = file(params.reference_genome_bai)
                def tr_bed = file(params.tr_bed_file)

                tuple(
                    sample_id,
                    roi_bam,
                    roi_bai,
                    ref,
                    ref_bai,
                    tr_bed,
                    modkit,
                    segs_bed,
                    bins_bed,
                    segs_vcf,
                    rds_file,
                    sv
                )
            }

        // Step 4: Run analysis
        log.info "=== Starting Annotation Pipeline ==="
        def annotation_results = annotation(annotation_input)

    } else if (params.run_mode_epiannotation) {
        // Combined mode: run epi2me then analysis (assumes merged BAM files already exist)
        // Set run_mode_annotation to 'all' to ensure all analyses run
        params.run_mode_annotation = 'all'

        log.info """
        Running combined Epi2me + Annotation Pipeline sequentially:
        1. Epi2me Pipeline (using existing merged BAM files)
        2. Annotation Pipeline (using epi2me outputs)
        """

        // Step 1: Create input channel for epi2me from existing merged BAM files
        // Read sample IDs from the epi2me sample_ids file to filter which samples to process
        def merged_bam_channel = Channel
            .from(file(params.epi2me_sample_id_file).readLines())
            .map { line ->
                def fields = line.trim().split(/\s+/) as List
                def sample_id = fields[0].trim()

                // Try exact match first
                def bam = file("${params.merge_bam_folder}/${sample_id}.merged.bam")
                def bai = file("${params.merge_bam_folder}/${sample_id}.merged.bam.bai")

                // If exact match doesn't exist, try wildcard pattern
                if (!bam.exists()) {
                    def bam_files = file("${params.merge_bam_folder}/${sample_id}.*.bam")
                    bam = bam_files instanceof List ? bam_files[0] : bam_files
                }
                if (!bai.exists()) {
                    def bai_files = file("${params.merge_bam_folder}/${sample_id}.*.bam.bai")
                    bai = bai_files instanceof List ? bai_files[0] : bai_files
                }

                if (!bam || !bai || !bam.exists() || !bai.exists()) {
                    error "BAM file or index file not found for sample ID: ${sample_id}"
                }

                tuple(
                    sample_id,
                    file(bam),
                    file(bai),
                    file(params.reference_genome),
                    file(params.reference_genome_bai)
                )
            }

        // Step 2: Run epi2me with merged BAM files (runs extract_roi internally)
        log.info "=== Starting Epi2me Pipeline ==="
        def epi2me_results = epi2me(merged_bam_channel)

        // Step 3: Combine epi2me results with OCC BAM files (generated by extract_roi inside epi2me)
        def annotation_input = epi2me_results.results
            .join(epi2me_results.roi_bam, by: 0)
            .map { joined ->
                def sample_id = joined[0]
                def modkit = joined[1]
                def segs_bed = joined[2]
                def bins_bed = joined[3]
                def segs_vcf = joined[4]
                def rds_file = joined[5]
                def sv = joined[6]
                def roi_bam = joined[7]
                def roi_bai = joined[8]

                def ref = file(params.reference_genome)
                def ref_bai = file(params.reference_genome_bai)
                def tr_bed = file(params.tr_bed_file)

                tuple(
                    sample_id,
                    roi_bam,
                    roi_bai,
                    ref,
                    ref_bai,
                    tr_bed,
                    modkit,
                    segs_bed,
                    bins_bed,
                    segs_vcf,
                    rds_file,
                    sv
                )
            }

        // Step 3: Run analysis
        log.info "=== Starting Annotation Pipeline ==="
        def annotation_results = annotation(annotation_input)

    } else {
        if (params.run_mode_mergebam) mergebam()
        if (params.run_mode_epi2me) epi2me(Channel.empty())
        if (params.run_mode_annotation) annotation(Channel.empty())
    }
}

// Single workflow completion handler
workflow.onComplete {
    def msg = """
        Pipeline execution summary
        ---------------------------
        Completed at : ${workflow.complete}
        Duration    : ${workflow.duration}
        Success     : ${workflow.success}
        workDir     : ${workflow.workDir}
        exit status : ${workflow.exitStatus}
        """
    if (workflow.success) {
        log.info msg
        // Note: Don't print "RMD report generated successfully" here as it may be misleading
        // when processes are cached. The run script will validate actual results.
    } else {
        log.error msg
    }
}

workflow.onError {
    log.error """
        Pipeline execution failed
        ---------------------------
        Error message: ${workflow.errorMessage}
        """
}
