#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Include pipeline modules
include { mergebam } from './modules/mergebam.nf'
include { epi2me }   from './modules/epi2me.nf'
include { analysis } from './modules/analysis.nf'

workflow {
    if (params.run_order_mode) {
        log.info """
        Running pipelines sequentially in strict order:
        1. Mergebam Pipeline
        2. Epi2me Pipeline
        3. Analysis Pipeline
        """

        // Step 1: Run mergebam
        log.info "=== Starting Mergebam Pipeline ==="
        def mergebam_results = mergebam()

        // Step 2: Run epi2me
        log.info "=== Starting Epi2me Pipeline ==="
        def epi2me_results = epi2me(mergebam_results.merged_bams)

        // Step 3: Prepare analysis input
        def analysis_input = epi2me_results.results.map { result ->
            def (sample_id, modkit, segs_bed, bins_bed, segs_vcf, sv) = result
            
            // Verify required files exist
            def bam = file("${params.occ_bam_folder}/${sample_id}.occ.bam")
            def bai = file("${params.occ_bam_folder}/${sample_id}.occ.bam.bai")
            def ref = file(params.reference_genome)
            def ref_bai = file(params.reference_genome_bai)
            def tr_bed = file(params.tr_bed_file)
            def bedmethyl = file("${params.bedmethyl_folder}/${sample_id}.wf_mods.bedmethyl.gz")

            assert bam.exists() : "BAM file not found: ${bam}"
            assert bai.exists() : "BAI file not found: ${bai}"
            assert bedmethyl.exists() : "Bedmethyl file not found: ${bedmethyl}"

            tuple(
                sample_id, 
                bam, 
                bai, 
                ref, 
                ref_bai, 
                tr_bed, 
                modkit, 
                segs_bed, 
                bins_bed, 
                segs_vcf, 
                sv
            )
        }

        // Step 4: Run analysis
        log.info "=== Starting Analysis Pipeline ==="
        analysis(analysis_input)
    } else {
        if (params.run_mode_mergebam) mergebam()
        if (params.run_mode_epi2me) epi2me(Channel.empty())
        if (params.run_mode_analysis) analysis(Channel.empty())
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
        if (params.run_mode_analysis == 'rmd') {
            log.info "RMD report generated successfully"
        }
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