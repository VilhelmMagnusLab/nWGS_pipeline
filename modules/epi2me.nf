#!/usr/bin/env nextflow
nextflow.enable.dsl=2

//---------------------------------------------------------------------
// Epi2me Pipeline: Modified base calling, structural variant calling, and copy number variation analysis
//---------------------------------------------------------------------

//---------------------------------------------------------------------
// Process Definitions
//---------------------------------------------------------------------

process run_epi2me_modkit {
    label 'modkit'
    publishDir "${params.output_path}/routine_epi2me/${sample_id}", mode: "copy", overwrite: true

    input:
    tuple val(sample_id), path(bam), path(bai), path(reference_genome), path(reference_genome_bai)

    output:
    tuple val(sample_id), path("${sample_id}.wf_mods.bedmethyl.gz")

    script:
    """
    export PATH=/opt/custflow/epi2meuser/conda/bin:\$PATH
    
    # Check if modkit is available
    which modkit || echo "modkit not found in PATH"
    
    # Check if input files exist
    if [ ! -f "${bam}" ]; then
        echo "ERROR: BAM file not found: ${bam}"
        exit 1
    fi
    
    if [ ! -f "${reference_genome}" ]; then
        echo "ERROR: Reference genome not found: ${reference_genome}"
        exit 1
    fi
    
    modkit pileup \
      ${bam} \
      ${sample_id}.wf_mods.bedmethyl \
      --ref ${reference_genome} \
      --interval-size ${params.interval_size} \
      --threads ${task.cpus} \
      --log-filepath ${sample_id}_modkit.log \
      --cpg \
      --combine-strands

    gzip ${sample_id}.wf_mods.bedmethyl

    
    # Check if output was created
    ls -la ${sample_id}.wf_mods.bedmethyl.gz
    """
}

process run_epi2me_sv {
    label 'pipeline1'
    publishDir "${params.output_path}/routine_epi2me/${sample_id}", mode: "copy", overwrite: true

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    tuple val(sample_id), path("${sample_id}.sniffles.vcf.gz"), path("${sample_id}.sniffles.vcf.gz.tbi"), emit: svvcf

    script:
    """
    # Check if sniffles is available
    which sniffles || echo "sniffles not found in PATH"
    
    # Check if input file exists
    if [ ! -f "${bam}" ]; then
        echo "ERROR: BAM file not found: ${bam}"
        exit 1
    fi
    
    # Run Sniffles2
    sniffles --input ${bam} --vcf ${sample_id}.sniffles.vcf.gz

    # Check if output was created
    ls -la ${sample_id}.sniffles.vcf.gz
    """
}

process run_epi2me_cnv {
    label 'epi2me'
    publishDir "${params.output_path}/routine_epi2me/${sample_id}", mode: "copy", overwrite: true

    input:
    tuple val(sample_id), path(bam), path(bai), path(reference_genome), path(reference_genome_bai)

    output:
    tuple val(sample_id), path("${sample_id}_segs.bed"), path("${sample_id}_bins.bed"), path("${sample_id}_segs.vcf"), path("${sample_id}_copyNumbersCalled.rds"), path("${sample_id}_calls.bed"), path("${sample_id}_calls.vcf")
    
    script:
    """
    # Create a .Renviron file to forcefully override R library paths
    # This prevents R from using host libraries even if they're accessible
    cat > .Renviron <<'RENVEOF'
R_LIBS_USER=
R_LIBS_SITE=
R_USER_CACHE_DIR=
R_PROFILE_USER=
R_ENVIRON_USER=
R_HOME=
RENVEOF

    # Clear all R environment variables to prevent host library conflicts
    # Use unset to completely remove the variables
    unset R_HOME R_LIBS R_LIBS_USER R_LIBS_SITE R_USER_CACHE_DIR R_PROFILE_USER R_ENVIRON_USER

    echo "R environment cleared for container isolation"

    # Run QDNAseq R script with explicit --vanilla flag to ignore R environment
    Rscript --no-site-file --no-init-file --no-environ $baseDir/bin/run_qdnaseq_rds.r \
        --bam ${bam} \
        --binsize ${params.binsize} \
        --out_prefix ${sample_id}

    # Check if output files were created
    ls -la ${sample_id}_*
    """
}

// Extract regions of interest from merged BAM — must complete before run_clair3 and run_clairs_to
process extract_roi {
    label 'roi_extraction'
    publishDir "${params.roi_bam_dir}", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), path(bam), path(bai), path(roi_bed)

    output:
    tuple val(sample_id), path("${sample_id}.roi.bam"), path("${sample_id}.roi.bam.bai"), emit: roi_bam

    script:
    """
    set -euo pipefail
    samtools view -@ ${task.cpus} -b -L ${roi_bed} ${bam} \\
    | samtools sort -@ ${task.cpus} -o ${sample_id}.roi.bam
    samtools index -@ ${task.cpus} ${sample_id}.roi.bam
    """
}

// SNV calling using Clair3 for OCC (regions of interest) regions
process run_clair3 {
    label 'clair3'
    publishDir "${params.output_path}/routine_epi2me/${sample_id}", mode: "copy", overwrite: true

    input:
    tuple val(sample_id), path(roi_bam), path(roi_bam_bai), path(reference_genome), path(reference_genome_bai), path(refGene), path(hg38_refGeneMrna), path(clinvar), path(clinvarindex), path(hg38_cosmic100), path(hg38_cosmic100index)

    output:
    tuple val(sample_id), path('output_clair3/'), emit: clair3_output_dir
    tuple val(sample_id), path('output_clair3/pileup.vcf.gz'), emit: pileup_vcf
    tuple val(sample_id), path('output_clair3/merge_output.vcf.gz'), emit: merge_vcf

    script:
    """
    # Activate conda environment for Clair3
    source /opt/conda/etc/profile.d/conda.sh
    conda activate clair3_v2 2>/dev/null || conda activate clair3 2>/dev/null || true

    /opt/bin/run_clair3.sh \
        --bam_fn=$roi_bam \
        --ref_fn=$reference_genome  \
        --threads=8 \
        --var_pct_full=1 \
        --ref_pct_full=1 \
        --var_pct_phasing=1 \
        --platform="ont" \
        --no_phasing_for_fa \
        --model_path=${params.clair3_model == "hac" ? params.clair3_model_path_hac : params.clair3_model_path} \
        --output=output_clair3

    # remove tmp folder
    rm -rf output_clair3/tmp*
    """
}

// SNV calling using ClairS-TO for somatic variants
process run_clairs_to {
    label 'clairsto'
    publishDir "${params.output_path}/routine_epi2me/${sample_id}", mode: "copy", overwrite: true

    input:
    tuple val(sample_id), path(roi_bam), path(roi_bam_bai), path(reference_genome), path(reference_genome_bai), path(refGene), path(hg38_refGeneMrna), path(clinvar), path(clinvarindex), path(hg38_cosmic100), path(hg38_cosmic100index), path(roi_protein_coding_bed)

    output:
    tuple val(sample_id), path('clairsto_output/'), emit: clairsto_output_dir
    tuple val(sample_id), path('clairsto_output/snv.vcf.gz'), emit: snv_vcf
    tuple val(sample_id), path('clairsto_output/indel.vcf.gz'), emit: indel_vcf

    script:
    """
    # Set micromamba environment variables
    export MAMBA_ROOT_PREFIX=/opt/micromamba
    export MAMBA_EXE=/opt/micromamba/bin/micromamba

    # Set MKL variables before activation
    export MKL_INTERFACE_LAYER=LP64
    export MKL_THREADING_LAYER=SEQUENTIAL

    # Activate conda environment
    source /opt/micromamba/etc/profile.d/micromamba.sh
    micromamba activate clairs-to

    # Run ClairS-TO (may produce empty VCFs if no variants found - this is normal)
    /opt/bin/run_clairs_to \
        --tumor_bam_fn=${roi_bam} \
        --ref_fn=${reference_genome} \
        --threads=${task.cpus} \
        --platform="ont_r10_dorado_4khz" \
        --output_dir=clairsto_output \
        --bed_fn=${roi_protein_coding_bed} \
        --conda_prefix /opt/micromamba/envs/clairs-to || echo "ClairS-TO completed with warnings (possibly no variants found)"

    # Ensure output files exist (create empty gzipped VCFs if missing)
    if [ ! -f clairsto_output/snv.vcf.gz ]; then
        echo "##fileformat=VCFv4.2" | bgzip > clairsto_output/snv.vcf.gz
    fi
    if [ ! -f clairsto_output/indel.vcf.gz ]; then
        echo "##fileformat=VCFv4.2" | bgzip > clairsto_output/indel.vcf.gz
    fi

    # Create index files for the VCF files (required for downstream annotation)
    tabix -p vcf clairsto_output/snv.vcf.gz || true
    tabix -p vcf clairsto_output/indel.vcf.gz || true

    # remove tmp folder
    rm -rf clairsto_output/tmp*
    """
}

// Quality assessment and statistics using Cramino
process cramino_report {
    label 'epic'
    publishDir "${params.output_path}/routine_epi2me/${sample_id}/cramino", mode: "copy", overwrite: true

    input:
    tuple val(sample_id), path(merge_bam), path(merge_bam_bai), path(reference_genome), path(reference_genome_bai)

    output:
    tuple val(sample_id), path("${sample_id}_cramino_statistics.txt"), emit:craminostatout

    script:
    """
    #!/bin/bash
    set -e

    # Check if we're in a container and use appropriate conda setup
    if [ -f "/opt/conda/etc/profile.d/conda.sh" ]; then
        # Container environment
        source /opt/conda/etc/profile.d/conda.sh
        conda activate annotatecnv_env
    else
        # Local environment
        source activate annotatecnv_env
    fi
    cramino $merge_bam --reference $reference_genome > ${sample_id}_cramino_statistics.txt
    """
}

//---------------------------------------------------------------------
// Workflow Definition
//---------------------------------------------------------------------

workflow epi2me {
    take:
        merged_data

    main:
        params.run_mode = params.run_mode_epi2me ?: 'all'
        println "epi2me run mode: ${params.run_mode}"

        // Validate required parameters
        if (!params.reference_genome) {
            error "ERROR: Reference genome path not specified (params.reference_genome)"
        }
        if (!params.reference_genome_bai) {
            error "ERROR: Reference genome index path not specified (params.reference_genome_bai)"
        }

        // Create file objects for reference files
        reference_genome = file(params.reference_genome)
        reference_genome_bai = file(params.reference_genome_bai)
        episv = file(params.episv)
        epimodkit = file(params.epimodkit)
        epicnv = file(params.epicnv)    

        // Validate reference files exist
        if (!reference_genome.exists()) {
            error "ERROR: Reference genome file not found: ${params.reference_genome}"
        }
        if (!reference_genome_bai.exists()) {
            error "ERROR: Reference genome index file not found: ${params.reference_genome_bai}"
        }

        // Validate run_mode parameter
        if (!['modkit', 'cnv', 'sv', 'snv', 'stat', 'all'].contains(params.run_mode)) {
            error "ERROR: Invalid run_mode '${params.run_mode}' for epi2me. Valid modes: modkit, cnv, sv, snv, stat, all."
        }

        // Create input channel based on run mode
        // Use merged_data input for both run_mode_order and run_mode_epiannotation
        input_channel = (params.run_mode_order || params.run_mode_epiannotation) ?
            merged_data.map { sid, bam, bai, ref, ref_bai ->
                tuple(
                    sid,
                    bam,
                    bai,
                    ref,
                    ref_bai
                )
            } : Channel
            .from(file(params.epi2me_sample_id_file).readLines())
            .map { line ->
                def fields = line.trim().split(/\s+/) as List
                def sample_id = fields[0].trim()
                // Try exact match first, then wildcard pattern
                def bam = file("${params.merge_bam_folder}/${sample_id}.merged.bam")
                def bai = file("${params.merge_bam_folder}/${sample_id}.merged.bam.bai")

                // If exact match doesn't exist, try wildcard pattern
                if (!bam.exists()) {
                    bam = file("${params.merge_bam_folder}/${sample_id}.*.bam")
                    bam = bam.find()
                }
                if (!bai.exists()) {
                    bai = file("${params.merge_bam_folder}/${sample_id}.*.bam.bai")
                    bai = bai.find()
                }

                if (!bam || !bai || !bam.exists() || !bai.exists()) {
                    error "BAM file or index file not found for sample ID: ${sample_id}. Tried both exact match (${sample_id}.merged.bam) and wildcard pattern (${sample_id}.*.bam)"
                }

                return tuple(
                    sample_id,
                    bam,
                    bai,
                    reference_genome,
                    reference_genome_bai
                )
            }

        // Run processes based on mode
        modkit_ch = Channel.empty()
        cnv_ch = Channel.empty()
        sv_ch = Channel.empty()
        clair3_ch = Channel.empty()
        clairsto_ch = Channel.empty()
        cramino_ch = Channel.empty()
        roi_bam_ch = Channel.empty()

        // extract_roi runs only for snv and all modes — output feeds into run_clair3 and run_clairs_to
        if (params.run_mode in ['snv', 'all']) {
            roi_bam_ch = extract_roi(
                input_channel.map { sid, bam, bai, ref, ref_bai ->
                    tuple(sid, bam, bai, file(params.roi_bed))
                }
            ).roi_bam
        }

        roi_input_channel = roi_bam_ch
            .map { sid, roi_bam, roi_bai ->
                tuple(sid, roi_bam, roi_bai, file(params.reference_genome), file(params.reference_genome_bai))
            }

        if (params.run_mode in ['modkit', 'all']) {
            println "Running modkit..."
            modkit_ch = run_epi2me_modkit(input_channel)
        }

        if (params.run_mode in ['cnv', 'all']) {
            println "Running cnv..."
            cnv_ch = run_epi2me_cnv(input_channel)
        }

        if (params.run_mode in ['sv', 'all']) {
            println "Running sv..."
            sv_ch = input_channel.map { sid, bam, bai, ref, ref_bai ->
                tuple(
                    sid,
                    bam,
                    bai
                )
            } | run_epi2me_sv
        }

        // SNV calling with Clair3 and ClairS-TO
        // Uses OCC/ROI BAM files and annotation databases
        if (params.run_mode in ['snv', 'all']) {
            println "Running SNV calling (Clair3 and ClairS-TO)..."

            // Load annotation files as channels
            def refgene_ch = Channel.value(file(params.refgene))
            def hg38_refgenemrna_ch = Channel.value(file(params.hg38_refgenemrna))
            def clinvar_ch = Channel.value(file(params.clinvar))
            def clinvarindex_ch = Channel.value(file(params.clinvarindex))
            def hg38_cosmic100_ch = Channel.value(file(params.hg38_cosmic100))
            def hg38_cosmic100index_ch = Channel.value(file(params.hg38_cosmic100index))
            def roi_protein_coding_bed_ch = Channel.value(file(params.roi_protein_coding_bed))

            // Prepare input for Clair3 (OCC BAM + annotation files)
            def clair3_input = roi_input_channel
                .combine(refgene_ch)
                .combine(hg38_refgenemrna_ch)
                .combine(clinvar_ch)
                .combine(clinvarindex_ch)
                .combine(hg38_cosmic100_ch)
                .combine(hg38_cosmic100index_ch)

            // Prepare input for ClairS-TO (OCC BAM + annotation files + OCC BED)
            def clairsto_input = roi_input_channel
                .combine(refgene_ch)
                .combine(hg38_refgenemrna_ch)
                .combine(clinvar_ch)
                .combine(clinvarindex_ch)
                .combine(hg38_cosmic100_ch)
                .combine(hg38_cosmic100index_ch)
                .combine(roi_protein_coding_bed_ch)

            // Run variant calling processes
            def clair3_result = run_clair3(clair3_input)
            clair3_ch = clair3_result.clair3_output_dir  // Use one of the outputs for dependency tracking

            def clairsto_result = run_clairs_to(clairsto_input)
            clairsto_ch = clairsto_result.clairsto_output_dir  // Use one of the outputs for dependency tracking
        }

        // Cramino statistics (runs for 'stat' mode or 'all' mode)
        if (params.run_mode in ['stat', 'all']) {
            println "Running Cramino statistics..."

            // Cramino uses merged BAM files
            def cramino_input = input_channel
                .view { "Cramino input: $it" }

            def cramino_result = cramino_report(cramino_input)
            cramino_ch = cramino_result.craminostatout  // Use the output for dependency tracking
        }

        // Create default channels for empty processes
        default_modkit = input_channel.map { sid, bam, bai, ref, ref_bai ->
            tuple(sid, file("${params.epimodkit}/${sid}.wf_mods.bedmethyl.gz"))
        }

        default_cnv = input_channel.map { sid, bam, bai, ref, ref_bai ->
            tuple(
                sid,
                file("${params.epicnv}/${sid}_segs.bed"),
                file("${params.epicnv}/${sid}_bins.bed"),
                file("${params.epicnv}/${sid}_segs.vcf")
            )
        }

        default_sv = input_channel.map { sid, bam, bai, ref, ref_bai ->
            tuple(sid, file("${params.episv}/${sid}.sniffles.vcf.gz"))
        }

        default_clair3 = input_channel.map { sid, bam, bai, ref, ref_bai ->
            tuple(sid, file("${params.output_path}/routine_epi2me/${sid}/output_clair3"))
        }

        default_clairsto = input_channel.map { sid, bam, bai, ref, ref_bai ->
            tuple(sid, file("${params.output_path}/routine_epi2me/${sid}/clairsto_output"))
        }

        default_cramino = input_channel.map { sid, bam, bai, ref, ref_bai ->
            tuple(sid, file("${params.output_path}/routine_epi2me/${sid}/cramino/${sid}_cramino_statistics.txt"))
        }

        // Mix actual results with defaults (but only for standalone mode, not for run_mode_order/epiannotation)
        // In run_mode_order and run_mode_epiannotation, we want to wait for actual process outputs
        modkit_results = (params.run_mode_order || params.run_mode_epiannotation) ?
            modkit_ch : modkit_ch.mix(default_modkit)
        cnv_results = (params.run_mode_order || params.run_mode_epiannotation) ?
            cnv_ch : cnv_ch.mix(default_cnv)
        sv_results = (params.run_mode_order || params.run_mode_epiannotation) ?
            sv_ch : sv_ch.mix(default_sv)
        clair3_results = (params.run_mode_order || params.run_mode_epiannotation) ?
            clair3_ch : clair3_ch.mix(default_clair3)
        clairsto_results = (params.run_mode_order || params.run_mode_epiannotation) ?
            clairsto_ch : clairsto_ch.mix(default_clairsto)
        cramino_results = (params.run_mode_order || params.run_mode_epiannotation) ?
            cramino_ch : cramino_ch.mix(default_cramino)

        // Combine all results
        results_ch = input_channel
            .join(modkit_results, by: 0)
            .join(cnv_results, by: 0)
            .join(sv_results, by: 0)

        // In epiannotation/order mode, wait for SNV and cramino to complete
        if (params.run_mode_epiannotation || params.run_mode_order) {
            // Collect all SNV and cramino outputs first as a barrier signal
            // NOTE: must use .combine() not .cross() — .cross() consumes the barrier item once,
            // so only the first sample passes; all remaining samples are silently dropped.
            // .combine() creates a cartesian product: every sample pairs with the single barrier.
            def snv_cramino_barrier = clair3_ch
                .mix(clairsto_ch)
                .mix(cramino_ch)
                .collect()
                .map { true }  // single "all done" signal

            results_ch = results_ch
                .combine(snv_cramino_barrier)          // appends barrier flag to each tuple
                .map { args -> args[0..(args.size()-2)] }  // drop the barrier flag
        }

        results_ch = results_ch
            .map { args ->
                def sample_id = args[0]
                def bam = args[1]
                def bai = args[2]
                def ref = args[3]
                def ref_bai = args[4]
                def modkit = args[5]

                // CNV outputs (6 files): segs_bed, bins_bed, segs_vcf, copyNumbersCalled.rds, calls.bed, calls.vcf
                def segs_bed = args[6]
                def bins_bed = args[7]
                def segs_vcf = args[8]
                def rds_file = args[9]  // copyNumbersCalled.rds - needed for ACE
                def calls_bed = args[10]
                def calls_vcf = args[11]

                // SV outputs: sv.vcf.gz and sv.vcf.gz.tbi
                def sv = args[12]
                def sv_index = args[13]

                // SNV outputs (Clair3, ClairS-TO) - args[18], args[19]
                // Cramino output - args[20]
                // These are included for dependency tracking but not returned

                log.info "Processing completed for sample: ${sample_id}"
                tuple(
                    sample_id,
                    modkit,
                    segs_bed,
                    bins_bed,
                    segs_vcf,
                    rds_file,
                    sv
                )
            }

    emit:
        results = results_ch
            .map { results ->
                log.info "Epi2me processing completed for sample: ${results[0]}"
                results
            }
        roi_bam = roi_bam_ch
}