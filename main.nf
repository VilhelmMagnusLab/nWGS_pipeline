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
        
        // Run mergebam first and capture its output
        log.info "=== Starting Mergebam Pipeline ==="
        def mergebam_results = mergebam()
        
        // Run epi2me after mergebam completes
        log.info "=== Starting Epi2me Pipeline ==="
        // Wait for mergebam to complete before starting epi2me
        mergebam_results.complete.view { "Mergebam completed" }
        def epi2me_results = mergebam_results.merged_bams
            | epi2me
        
        // Run analysis after epi2me completes
        log.info "=== Starting Analysis Pipeline ==="
        // Wait for epi2me to complete before starting analysis
        epi2me_results.complete.view { "Epi2me completed" }
        epi2me_results.results
            | analysis
    }
    else {
        if (params.run_mode_mergebam) {
            mergebam()
        }
        if (params.run_mode_epi2me) {
            epi2me(Channel.empty())
        }
        if (params.run_mode_analysis) {
            analysis(Channel.empty())
        }
    }
}
