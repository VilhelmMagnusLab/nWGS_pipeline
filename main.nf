#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Include pipeline modules
include { mergebam } from './modules/mergebam.nf'
include { epi2me }   from './modules/epi2me.nf'
include { analysis } from './modules/analysis.nf'

workflow {
    if (params.run_order_mode) {
        println """
        Running pipelines sequentially in strict order:
        1. Mergebam Pipeline
        2. Epi2me Pipeline
        3. Analysis Pipeline
        """
        
        // Run mergebam first
        println "=== Starting Mergebam Pipeline ==="
        def mergebam_results = mergebam()

        // Run epi2me after mergebam completes
        println "=== Starting Epi2me Pipeline ==="
        def epi2me_results = epi2me(mergebam_results.merged_bams)

        // Wait for epi2me processes to complete
        epi2me_results.results.subscribe { results ->
            def (sample_id, bam, bai, ref, ref_bai, tr_bed, modkit, segs_bed, bins_bed, segs_vcf, sv) = results
            println "Epi2me processes completed for ${sample_id}"
        }

        // Create channel for analysis
        def analysis_input = epi2me_results.results
            .map { results ->
                def (sample_id, bam, bai, ref, ref_bai, tr_bed, modkit, segs_bed, bins_bed, segs_vcf, sv) = results
                
                // Get paths from config
                def cnv_files = [
                    segs_vcf: file("${params.epicnv}/qdna_seq/${sample_id}_segs.vcf"),
                    segs_bed: file("${params.epicnv}/qdna_seq/${sample_id}_segs.bed"),
                    bins_bed: file("${params.epicnv}/qdna_seq/${sample_id}_bins.bed"),
                    calls_bed: file("${params.epicnv}/qdna_seq/${sample_id}_calls.bed")
                ]
                
                def sv_file = file("${params.episv}/${sample_id}_wf_sv.vcf.gz")
                def modkit_file = file("${params.epimodkit}/${sample_id}_wf_mods.bedmethyl.gz")
                
                // Wait for files to be ready
                while (!sv_file.exists() || !modkit_file.exists() || !cnv_files.every { k, f -> f.exists() }) {
                    sleep(300000) // Wait 5 minutes
                    println "Waiting for output files to be ready for ${sample_id}..."
                }

                // Verify file sizes
                def verifyFileSize = { file ->
                    if (file.exists() && file.size() == 0) {
                        error "File ${file} exists but is empty"
                    }
                }
                
                verifyFileSize(sv_file)
                verifyFileSize(modkit_file)
                cnv_files.each { k, f -> verifyFileSize(f) }
                
                println "All files verified for ${sample_id}, proceeding with analysis"
                results
            }

        // Run analysis after verification
        println "=== Starting Analysis Pipeline ==="
        analysis(analysis_input)
        
    } else {
        // Run workflows independently based on mode parameters
        if (params.run_mode_mergebam) mergebam()
        if (params.run_mode_epi2me) epi2me()
        if (params.run_mode_analysis) analysis()
    }
}
