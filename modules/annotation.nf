#!/usr/bin/env nextflow
nextflow.enable.dsl=2

//---------------------------------------------------------------------
// Annotation Pipeline: Comprehensive annotation of glioblastoma samples
// Includes MGMT methylation, structural variant annotation, CNV analysis,
// SNV annotation, and report generation
//---------------------------------------------------------------------

def start_time = new Date()

//---------------------------------------------------------------------
// Helper function definitions
//---------------------------------------------------------------------

def validateParameters() {
    params.run_mode = params.run_mode_annotation ?: 'rmd'
    println "Annotation run mode: ${params.run_mode}"
    if (!['mgmt', 'svannasv', 'cnv', 'roi', 'tertp', 'rmd', 'all'].contains(params.run_mode)) {
        error "ERROR: Invalid run_mode '${params.run_mode}' for annotation. Valid modes: mgmt, svannasv, cnv, roi, tertp, rmd, all"
    }
}

def loadSampleThresholds() {
    return file(params.analyse_sample_id_file).readLines()
        .collectEntries { line ->
            def tokens = line.trim().split(/\s+/)  // Handle any whitespace
            if (tokens.size() == 2) {
                // User provided tumor content: [sample_id, tumor_content]
                def threshold = tokens[1].toFloat()
                if (threshold < 0 || threshold > 1) {
                    throw new IllegalArgumentException("Tumor content must be between 0-1: ${line}")
                }
                [(tokens[0]) : threshold]
            } else if (tokens.size() == 1) {
                // User needs ACE calculation: [sample_id]
                [(tokens[0]) : null]  // null indicates need ACE calculation
            } else {
                throw new IllegalArgumentException("Invalid line format in sample_id_file: ${line}")
            }
        }
}

//---------------------------------------------------------------------
// Process Definitions
//---------------------------------------------------------------------

// Extract EPIC methylation sites and MGMT CpG islands from bedmethyl data
process extract_epic {
    label 'epic'
    tag "${sample_id}"
    publishDir "${params.output_path}/routine_annotation/${sample_id}/methylation/", mode: "copy", overwrite: true
    publishDir "${params.path}/routine_results/${sample_id}", mode: "copy", overwrite: true, pattern: "*_mnpflex_input.bed"

    input:
    tuple val(sample_id), file(bedmethyl), file(epicsites), file(mgmt_cpg_island_hg38)

    output:
    tuple val(sample_id), path("${sample_id}_EpicSelect_header.bed"), emit: epicselectnanodxinput
    path("${sample_id}_MGMT.bed")
    tuple val(sample_id), path("${sample_id}_MGMT_header.bed"), emit: MGMTheaderout
    tuple val(sample_id), path("${sample_id}_wf_mods.bedmethyl_intersect.bed"), emit: sturgeonbedinput
    tuple val(sample_id), path("${sample_id}_EpicSelect_m.bed")
    tuple val(sample_id), path("${sample_id}_mnpflex_input.bed"), emit: mnpflex_bed

    script:
    """
    which intersectBed
    intersectBed -a $bedmethyl -b $epicsites -wb | \
    awk -v OFS="\\t" '\$1=\$1' | awk -F'\\t' 'BEGIN{ OFS="\\t" }{print \$1,\$2,\$3,\$4,\$5,\$11,\$23}' > ${sample_id}_EpicSelect.bed

    intersectBed -a $bedmethyl -b $epicsites -wb | awk -v OFS="\\t" '\$1=\$1' | awk -F'\\t' 'BEGIN{ OFS="\\t" } {print \$1, \$2, \$3, \$4, \$5, \$6, \$7, \$8, \$9, \$10, \$11, \$12, \$13, \$14, \$15, \$16, \$17, \$18}' > ${sample_id}_wf_mods.bedmethyl_intersect.bed

    awk 'BEGIN {print "Chromosome\\tStart\\tEnd\\tmodBase\\tCoverage\\tMethylation_frequency\\tIllumina_ID"} 1' ${sample_id}_EpicSelect.bed > ${sample_id}_EpicSelect_header.bed

    grep -w 'm'  ${sample_id}_EpicSelect_header.bed >  ${sample_id}_EpicSelect_m.bed

    awk 'BEGIN{
    OFS="\\t";
    print "\\"chr\\"","\\"start\\"","\\"end\\"","\\"coverage\\"","\\"methylation_percentage\\"","\\"IlmnID\\""
    }
    {
    if (\$1 ~ /^chr/) {
        printf "\\"%s\\"\\t%s\t%s\\t%s\\t%s\\t\\"%s\\"\\n", \$1, \$2, \$3, \$5, \$6, \$7
    }
    }   ' ${sample_id}_EpicSelect_m.bed > ${sample_id}_mnpflex_input.bed


    intersectBed -a $bedmethyl -b $mgmt_cpg_island_hg38 | \
    awk -v OFS="\\t" '\$1=\$1' | awk -F'\\t' 'BEGIN{ OFS="\\t" }{print \$1,\$2,\$3,\$4,\$5,\$11,\$12,\$13,\$14,\$15,\$16}'  > ${sample_id}_MGMT.bed

    awk 'BEGIN {print "Chrom\\tStart\\tEnd\\tmodBase\\tDepth\\tMethylation\\tNmod\\tNcanon\\tNother\\tNdelete\\tNfail"} 1' ${sample_id}_MGMT.bed > ${sample_id}_MGMT_header.bed
    """
}


// Prepare nanoDX input for methylation data processing and CpG site intersection
process nanodx {
    label 'epic'
    publishDir "${params.output_path}/routine_annotation/${sample_id}/classifier/nanodx", mode: "copy", overwrite: true

    input:
    tuple val(sample_id), path(nanodx_bed), path(hg19_450model)

    output:
    tuple val(sample_id), path("${sample_id}_nanodx_bedmethyl.bed"), emit: nanodx450out

    script:
    """
    nanodx450intersectdataframe.py $hg19_450model $nanodx_bed ${sample_id}_output_cpg.bed ${sample_id}_nanodx_bedmethyl.bed ${sample_id}_nanodx_bedmethylfilter.bed
    """
}

//Sturgeon classifier
process sturgeon {
    label 'epic'
    publishDir "${params.output_path}/routine_annotation/${sample_id}/classifier/sturgeon", mode: "copy", overwrite: true
    publishDir "${params.path}/routine_results/${sample_id}", mode: "copy", overwrite: true, pattern: "*_bedmethyl_sturgeon_general.pdf"

    input:
    tuple val(sample_id), path(sturgeon_bed), path(sturgeon_model)

    output:
    path("${sample_id}_bedmethyl_sturgeon.bed")
    tuple val(sample_id), path("${sample_id}_bedmethyl_sturgeon_general.pdf"), optional: true, emit: sturgeon_pdf
    path("${sample_id}_bedmethyl_sturgeon_general.csv"), optional: true


    """
    /sturgeon/venv/bin/sturgeon inputtobed -i $sturgeon_bed  -o ${sample_id}_bedmethyl_sturgeon.bed  -s modkit_pileup  --reference-genome hg38

    /sturgeon/venv/bin/sturgeon predict -i ${sample_id}_bedmethyl_sturgeon.bed   -o  ${sample_id}_bedmethyl_sturgeon --model-files $sturgeon_model  --plot-results

    # Copy the PDF and CSV files to the work directory for publishDir to find them
    if [ -f "${sample_id}_bedmethyl_sturgeon/${sample_id}_bedmethyl_sturgeon_general.pdf" ]; then
        cp "${sample_id}_bedmethyl_sturgeon/${sample_id}_bedmethyl_sturgeon_general.pdf" .
    fi
    if [ -f "${sample_id}_bedmethyl_sturgeon/${sample_id}_bedmethyl_sturgeon_general.csv" ]; then
        cp "${sample_id}_bedmethyl_sturgeon/${sample_id}_bedmethyl_sturgeon_general.csv" .
    fi
    # Note: ${sample_id}_bedmethyl_sturgeon.bed is already in the work directory (created at line 134)
    rm -rf ${sample_id}_bedmethyl_sturgeon
    """
    
}


// Neural network classification for tumor type prediction using NanoDx CrossNN classifier
process crossNN {
    label 'nanodx'
    publishDir "${params.output_path}/routine_annotation/${sample_id}/classifier/nanodx", mode: "copy", overwrite: true
    
    input:
    tuple val(sample_id), path(bed_file), path(model_file), path(snakefile), path(nn_model)
    
    output:
    tuple val(sample_id), path("${sample_id}_nanodx_classifier.txt")
    tuple val(sample_id), path("${sample_id}_nanodx_classifier.tsv"), emit: rmdnanodx
    
    script:
    """
    #!/bin/bash
    export TMPDIR="\${PWD}/tmp/"
    export XDG_CACHE_HOME="\${TMPDIR}.cache"
    mkdir -p "\$TMPDIR/.cache"
    
    # Use container's conda environment
    source /opt/conda/etc/profile.d/conda.sh
    conda activate base
    conda activate nanodx_env2feb

    # Remove symlinked Snakefile to avoid overwriting the reference file
    rm -f Snakefile

    # Create Snakefile with correct paths and reduced memory
    cat << EOF > Snakefile
rule all:
    input:
        "${sample_id}_nanodx_classifier.txt",
        "${sample_id}_nanodx_classifier.tsv"

rule NN_classifier:
    input:
        bed = "${bed_file}",
        model = "${model_file}"
    output:
        txt = "${sample_id}_nanodx_classifier.txt",
        votes = "${sample_id}_nanodx_classifier.tsv"
    threads: 2
    resources:
        mem_mb = 2048
    script: "${params.nanodx_workflow_dir}/scripts/classify_NN_bedMethyl.py"
EOF

    # Run snakemake with error handling
    if ! snakemake --cores ${task.cpus} --verbose NN_classifier; then
        echo "Snakemake failed, creating empty output files"
        touch "${sample_id}_nanodx_classifier.txt"
        touch "${sample_id}_nanodx_classifier.tsv"
        echo "NN classifier failed - pickle compatibility issue" > "${sample_id}_nanodx_classifier.txt"
    fi
    """
}

// Dimensionality reduction plot generation using t-SNE/UMAP
process tsne_plot {
    label 'tsne'
    stageInMode 'copy'
    publishDir "${params.output_path}/routine_annotation/${sample_id}/classifier/nanodx", mode: "copy", overwrite: true
    publishDir "${params.path}/routine_results/${sample_id}", mode: "copy", overwrite: true, pattern: "*_tsne_plot.html"

    input:
    tuple val(sample_id), path(epic_bed)
    path(color_map)
    path(training_set)

    output:
    tuple val(sample_id), path("${sample_id}_tsne_plot.pdf"), emit: tsne_out
    tuple val(sample_id), path("${sample_id}_tsne_plot.html"), emit: tsne_html

    script:
    """
    #!/bin/bash
    set -e

    # Activate conda environment to access Rscript
    source /opt/conda/etc/profile.d/conda.sh
    conda activate tsneenv

    # Verify Rscript is available
    if ! command -v Rscript >/dev/null 2>&1; then
        echo "ERROR: Rscript not found after activating tsneenv"
        echo "Available conda environments:"
        conda env list
        echo "PATH: \$PATH"
        exit 1
    fi

    echo "Using Rscript from: \$(which Rscript)"

    # Run the memory-optimized t-SNE script (uses 30k probes instead of 100k to reduce RAM usage)
    Rscript ${params.diana_dir}/bin/crossnn_tsne_fixed.R \\
        --color-map ${color_map} \\
        --bed ${epic_bed} \\
        --trainingset ${training_set} \\
        --method umap \\
        --umap-n-neighbours 10 \\
        --umap-min-dist 0.5 \\
        --umap-pca-dim 100 \\
        --pdf ${sample_id}_tsne_plot.pdf \\
        --html ${sample_id}_tsne_plot.html

    echo "t-SNE plot generated for ${sample_id}"
    """
}

// Pan-cancer CrossNN classifier — always runs alongside Capper et al.
process crossNN_pancan {
    label 'nanodx'
    publishDir "${params.output_path}/routine_annotation/${sample_id}/classifier/nanodx", mode: "copy", overwrite: true

    input:
    tuple val(sample_id), path(bed_file), path(model_file), path(snakefile), path(nn_model)

    output:
    tuple val(sample_id), path("${sample_id}_nanodx_classifier_pancan.txt")
    tuple val(sample_id), path("${sample_id}_nanodx_classifier_pancan.tsv"), emit: rmdnanodx_pancan

    script:
    """
    #!/bin/bash
    export TMPDIR="\${PWD}/tmp/"
    export XDG_CACHE_HOME="\${TMPDIR}.cache"
    mkdir -p "\$TMPDIR/.cache"

    source /opt/conda/etc/profile.d/conda.sh
    conda activate base
    conda activate nanodx_env2feb

    rm -f Snakefile

    cat << EOF > Snakefile
rule all:
    input:
        "${sample_id}_nanodx_classifier_pancan.txt",
        "${sample_id}_nanodx_classifier_pancan.tsv"

rule NN_classifier:
    input:
        bed = "${bed_file}",
        model = "${model_file}"
    output:
        txt = "${sample_id}_nanodx_classifier_pancan.txt",
        votes = "${sample_id}_nanodx_classifier_pancan.tsv"
    threads: 2
    resources:
        mem_mb = 2048
    script: "${params.nanodx_workflow_dir}/scripts/classify_NN_bedMethyl.py"
EOF

    if ! snakemake --cores ${task.cpus} --verbose NN_classifier; then
        echo "Snakemake failed, creating empty output files"
        touch "${sample_id}_nanodx_classifier_pancan.txt"
        touch "${sample_id}_nanodx_classifier_pancan.tsv"
        echo "NN pancan classifier failed - pickle compatibility issue" > "${sample_id}_nanodx_classifier_pancan.txt"
    fi
    """
}

// Pan-cancer t-SNE/UMAP plot — always runs alongside Capper et al.
process tsne_plot_pancan {
    label 'tsne'
    stageInMode 'copy'
    publishDir "${params.output_path}/routine_annotation/${sample_id}/classifier/nanodx", mode: "copy", overwrite: true
    publishDir "${params.path}/routine_results/${sample_id}", mode: "copy", overwrite: true, pattern: "*_tsne_plot_pancan.html"

    input:
    tuple val(sample_id), path(epic_bed)
    path(color_map)
    path(training_set)

    output:
    tuple val(sample_id), path("${sample_id}_tsne_plot_pancan.pdf"), emit: tsne_pancan_out
    tuple val(sample_id), path("${sample_id}_tsne_plot_pancan.html"), emit: tsne_pancan_html

    script:
    """
    #!/bin/bash

    source /opt/conda/etc/profile.d/conda.sh
    conda activate tsneenv

    if ! command -v Rscript >/dev/null 2>&1; then
        echo "WARNING: Rscript not found — creating placeholder outputs"
        touch "${sample_id}_tsne_plot_pancan.pdf"
        touch "${sample_id}_tsne_plot_pancan.html"
        exit 0
    fi

    echo "Using Rscript from: \$(which Rscript)"

    if ! Rscript ${params.diana_dir}/bin/crossnn_tsne_fixed.R \\
        --color-map ${color_map} \\
        --bed ${epic_bed} \\
        --trainingset ${training_set} \\
        --method umap \\
        --umap-n-neighbours 10 \\
        --umap-min-dist 0.5 \\
        --umap-pca-dim 100 \\
        --pdf ${sample_id}_tsne_plot_pancan.pdf \\
        --html ${sample_id}_tsne_plot_pancan.html; then
        echo "WARNING: Pan-cancer t-SNE failed (training set may be corrupt) — creating placeholder outputs"
        touch "${sample_id}_tsne_plot_pancan.pdf"
        touch "${sample_id}_tsne_plot_pancan.html"
    else
        echo "Pan-cancer t-SNE plot generated for ${sample_id}"
    fi
    """
}

// MGMT promoter methylation analysis and quantification
process mgmt_promoter {
    label 'epic'
    publishDir "${params.output_path}/routine_annotation/${sample_id}/methylation/", mode: "copy", overwrite: true

    input:
    tuple val(sample_id), path(mgmt_bed)

    output:
    tuple val(sample_id), path("${sample_id}_MGMT_results.csv"), emit: mgmtresultsout
    
    script:
  """
    #!/bin/bash
    set -e  # Exit on error
    
    if [ ! -f "${mgmt_bed}" ]; then
        echo "Error: Input file ${mgmt_bed} not found"
        exit 1
    fi
    
    MGMT_Prospective2.R ${mgmt_bed} "${sample_id}_MGMT_results.csv"
    
    if [ ! -f "${sample_id}_MGMT_results.csv" ]; then
        echo "Error: Output file not generated"
        exit 1
    fi
    """
}

// Structural variant annotation using Svanna for region of interest genes (occ genes)
process svannasv {

  label 'svannasv'
   publishDir "${params.output_path}/routine_annotation/${sample_id}/structure_variant/svannasv/", mode: "copy", overwrite: true
   publishDir "${params.path}/routine_results/${sample_id}", mode: "copy", overwrite: true, pattern: "*_occ_svanna_annotation.html"

   input:
   tuple val(sample_id), path(wf_sv), path(wf_sv_tbi),path(roi_protein_coding_bed)

   output:
   
   tuple val(sample_id), path("${sample_id}_occ_svanna_annotation.html"), emit:rmdsvannahtml
   tuple val(sample_id), path("${sample_id}_occ_svanna_annotation.vcf.gz"), emit: occsvannaannotationannotationvcf
   tuple val(sample_id), path("${sample_id}_sniffles2_under30mb.vcf.gz")

   script:
   """
   # Debug: List input files
   echo "Input files:"
   echo "SV file: $wf_sv"
   echo "SV index: $wf_sv_tbi"
   echo "OCC protein coding bed: $roi_protein_coding_bed"
   ls -la

   bcftools view -O z -o ${sample_id}_sniffles2_under30mb.vcf.gz   -i '(INFO/SVTYPE="DUP" || INFO/SVTYPE="INV" || INFO/SVTYPE="INS") && INFO/SVLEN < 30000000 || (INFO/SVTYPE="DEL" && INFO/SVLEN > -30000000)'   $wf_sv
   bcftools index -t ${sample_id}_sniffles2_under30mb.vcf.gz


   java -jar ${params.svanna_bin_dir}/svanna-cli-1.0.4.jar prioritize  \
   -d ${params.ref_dir}/svanna-data  \
   --vcf  $wf_sv \
   --phenotype-term HP:0100836 \
   --output-format html,vcf \
   --prefix ${sample_id}_occ_svanna_annotation

  # cp "${sample_id}_occ_svanna_annotation.html" "${params.output_path}/report/${sample_id}_svanna.html"

   """
}

// Fusion event analysis and filtering from structural variants
process svannasv_fusion_events {
    label 'svannasv'
    publishDir "${params.output_path}/routine_annotation/${sample_id}/structure_variant/svannasv/", mode: "copy", overwrite: true

    input:
    tuple val(sample_id), path(occ_svannavcf), path(genecode_bed), path(occ_fusions_genes)

    output:
    tuple val(sample_id), path("${sample_id}_filter_fusion_event.tsv"), emit: filterfusioneventout
    tuple val(sample_id), path("${sample_id}_filter_fusion_event_detailed.tsv")
    tuple val(sample_id), path("${sample_id}_filter_fusion_event_ensembl_only.tsv"), optional: true
    tuple val(sample_id), path("${sample_id}_filter_fusion_event_detailed_ensembl_only.tsv"), optional: true

    script:

    """
    # Create enhanced GFF3 with both exons and introns
    gunzip -c $occ_svannavcf > ${sample_id}_occ_svanna_annotation.vcf
   
    create_gff3_with_introns.py --gff3 $genecode_bed --out ${sample_id}_genecode_with_introns.gff3

    breaking_point_bed_translocation_exon.py --vcf ${sample_id}_occ_svanna_annotation.vcf --out ${sample_id}_breaking_bedpoints.bed

    awk 'BEGIN{OFS="\t"} {if (\$1 !~ /^chr/) \$1 = "chr"\$1; print}' ${sample_id}_breaking_bedpoints.bed > ${sample_id}_breaking_bedpoints_sort.bed

    # Intersect with enhanced GFF3 that includes introns
    intersectBed -a ${sample_id}_breaking_bedpoints_sort.bed  -b ${sample_id}_genecode_with_introns.gff3  -wb  > ${sample_id}_breaking_bedpoints_genecode.bed

    #remove duplicate bed points and add exon/intron annotations

    remove_duplicate_report_exon.py --in  ${sample_id}_breaking_bedpoints_genecode.bed  \
            --formatted ${sample_id}_breaking_bedpoints_genecode_format.bed \
             --out ${sample_id}_breaking_bedpoints_genecode_clean.bed    \
             --paired ${sample_id}_breaking_bedpoints_genecode_clean_paired.bed  \
             --gene-list $occ_fusions_genes \
             --filtered ${sample_id}_filter_fusion_event_detailed_temp.tsv

    # Add intergenic annotations for breakpoints that don't overlap any features
    annotate_intergenic_breakpoints.py --original-bed ${sample_id}_breaking_bedpoints_sort.bed \
             --annotated ${sample_id}_filter_fusion_event_detailed_temp.tsv \
             --out ${sample_id}_filter_fusion_event_detailed.tsv

    #summarize exon/intron/intergenic features into compact format

    summarize_fusion_features.py --in ${sample_id}_filter_fusion_event_detailed.tsv \
             --out ${sample_id}_filter_fusion_event_temp.tsv

    # Filter to keep only complete fusion pairs (IDs with both start and end breakpoints)
    # Then keep only one representative fusion per unique gene pair
    awk 'NR==1 {header=\$0; next}
         {id=\$4; breaking=\$6; gene=\$7;
          data[id,breaking]=\$0; ids[id]++;
          genes[id,breaking]=gene;
          if(breaking=="start") has_start[id]=1;
          if(breaking=="end") has_end[id]=1}
         END {print header;
              for(id in ids) {
                if(has_start[id] && has_end[id]) {
                  gene_start=genes[id,"start"];
                  gene_end=genes[id,"end"];
                  gene_pair=(gene_start < gene_end) ? gene_start"-"gene_end : gene_end"-"gene_start;
                  if(!seen_pair[gene_pair]++) {
                    if((id,"start") in data) print data[id,"start"];
                    if((id,"end") in data) print data[id,"end"]
                  }
                }
              }
         }' ${sample_id}_filter_fusion_event_temp.tsv > ${sample_id}_filter_fusion_event.tsv

    # Filter fusion events to keep only those with official Ensembl gene IDs (ENSG...)
    echo "Filtering for fusions with complete Ensembl gene IDs..."
    filter_fusion_complete_ensembl.py \
        --input ${sample_id}_filter_fusion_event_detailed.tsv \
        --output ${sample_id}_filter_fusion_event_detailed_ensembl_only.tsv \
        --gff3 $genecode_bed \
        --stats ${sample_id}_ensembl_filter_stats.txt

    # Annotate fusion breakpoints with exon coordinates and coding phase
    if [ -s ${sample_id}_filter_fusion_event_detailed_ensembl_only.tsv ]; then
        echo "Annotating breakpoints with exon coordinates and phase information..."
        annotate_fusion_exon_phase.py \
            --input ${sample_id}_filter_fusion_event_detailed_ensembl_only.tsv \
            --output ${sample_id}_filter_fusion_event_detailed_ensembl_only_exon_phase.tsv \
            --gff3 $genecode_bed \
            --stats ${sample_id}_exon_phase_stats.txt

        # Use the exon/phase annotated version as the detailed output
        mv ${sample_id}_filter_fusion_event_detailed_ensembl_only_exon_phase.tsv ${sample_id}_filter_fusion_event_detailed_ensembl_only.tsv

        # Create summarized version of Ensembl-filtered fusions
        summarize_fusion_features.py --in ${sample_id}_filter_fusion_event_detailed_ensembl_only.tsv \
             --out ${sample_id}_filter_fusion_event_ensembl_only_temp.tsv

        # Filter to keep only complete fusion pairs (IDs with both start and end breakpoints)
        # Then keep only one representative fusion per unique gene pair
        awk 'NR==1 {header=\$0; next}
             {id=\$4; breaking=\$6; gene=\$7;
              data[id,breaking]=\$0; ids[id]++;
              genes[id,breaking]=gene;
              if(breaking=="start") has_start[id]=1;
              if(breaking=="end") has_end[id]=1}
             END {print header;
                  for(id in ids) {
                    if(has_start[id] && has_end[id]) {
                      gene_start=genes[id,"start"];
                      gene_end=genes[id,"end"];
                      gene_pair=(gene_start < gene_end) ? gene_start"-"gene_end : gene_end"-"gene_start;
                      if(!seen_pair[gene_pair]++) {
                        if((id,"start") in data) print data[id,"start"];
                        if((id,"end") in data) print data[id,"end"]
                      }
                    }
                  }
             }' ${sample_id}_filter_fusion_event_ensembl_only_temp.tsv > ${sample_id}_filter_fusion_event_ensembl_only.tsv
    fi

    """

}

// Circos plot generation for structural variant visualization
process circosplot {
   label 'circos'
   publishDir "${params.output_path}/routine_annotation/${sample_id}/structure_variant/svannasv/", mode: "copy", overwrite: true
   
   input:
   tuple val(sample_id), path(svanna_output), path(vcf2circos_json)

   output:
   tuple val(sample_id), path("${sample_id}_vcf2circo.html"), optional: true, emit: circosout

   script:
   """
   # Check if file is empty (excluding header)
   # Handle both gzipped and uncompressed VCF files
   if [[ ${svanna_output} == *.gz ]]; then
      # For gzipped files, use zcat
      VARIANT_COUNT=\$(zcat ${svanna_output} | grep -v '^#' | wc -l)
   else
      # For uncompressed files, use grep directly
      VARIANT_COUNT=\$(grep -v '^#' ${svanna_output} | wc -l)
   fi

   if [ \$VARIANT_COUNT -eq 0 ]; then
      echo "Warning: ${svanna_output} has no variants (0 non-header lines). Skipping vcf2circos plot generation."
      touch ${sample_id}_vcf2circo.html
      exit 0
   else
      echo "Found \$VARIANT_COUNT variants in ${svanna_output}. Generating circos plot..."
      # Try to generate circos plot, but create empty file if it fails
      set +e  # Don't exit on error
      vcf2circos -i $svanna_output -o ${sample_id}_vcf2circo.html -p $vcf2circos_json -a hg38
      CIRCOS_EXIT=\$?
      set -e  # Re-enable exit on error

      if [ \$CIRCOS_EXIT -ne 0 ]; then
          echo "Warning: vcf2circos failed with exit code \$CIRCOS_EXIT. Creating empty placeholder file."
          touch ${sample_id}_vcf2circo.html
      else
          echo "Circos plot generated successfully."
      fi
   fi
   """
}

// Copy number variant annotation and analysis with ACE
process annotatecnv {
   label 'annotatecnv'
    publishDir "${params.output_path}/routine_annotation/${sample_id}/cnv/", mode: "copy", overwrite: true

   input:
    tuple val(sample_id),
          path(vcf_file),
          path(roi_protein_coding_bed),
          path(calls_bed),
          path(seg_bed),
          val(threshold),  // Now explicitly receiving threshold
          path(cnv_genes_tuned)  // CNV genes annotation file

   output:
   tuple val(sample_id), path("${sample_id}_calls_fixed.vcf"), emit: callsfixedout
   tuple val(sample_id), path("${sample_id}_annotatedcnv.csv"), emit:annotatedcnvcsvout
   tuple val(sample_id), path("${sample_id}_annotatedcnv_filter.csv"), emit:annotatedcnvfiltercsvout
   tuple val(sample_id), path("${sample_id}_annotatedcnv_filter_header.csv"), emit:rmdannotatedcnvfilter
   tuple val(sample_id), path("${sample_id}_tumor_copy.txt"), path("${sample_id}_bins_filter.bed"), emit:tumorcopyandbinsfilterout
   tuple val(sample_id), path("${sample_id}_cnv_plot_full.pdf"), path("${sample_id}_tumor_copy_number.txt"), path("${sample_id}_annotatedcnv_filter_header.csv"), path("${sample_id}_cnv_chr9.pdf"), path("${sample_id}_cnv_chr7.pdf"), emit: rmdcnvtumornumber
   // TIFF outputs from v2 script (for high-resolution images)
   ///tuple val(sample_id), path("${sample_id}_cnv_plot_full.tiff"), emit: cnvplotfulltiff
   ///tuple val(sample_id), path("${sample_id}_cnv_chr9.tiff"), emit: cnvchr9tiff
   ///tuple val(sample_id), path("${sample_id}_cnv_chr7.tiff"), emit: cnvchr7tiff

   script:
   """
    #!/bin/bash
    ##source /opt/conda/etc/profile.d/conda.sh
    ##conda activate annotatecnv_env

    # Check if we're in a container and use appropriate conda setup
    if [ -f "/opt/conda/etc/profile.d/conda.sh" ]; then
        # Container environment
   source /opt/conda/etc/profile.d/conda.sh
   conda activate annotatecnv_env
    else
        # Local environment
        source activate annotatecnv_env
    fi

    # Debug output
    echo "Processing sample: ${sample_id}"
    echo "Using threshold: ${threshold}"
    
    # Process VCF file
    awk 'OFS="\\t" {if (NR > 13) \$1="chr"\$1; print}' $vcf_file > ${sample_id}_calls_fixed.vcf
    
    # Intersect and process
    intersectBed -a ${sample_id}_calls_fixed.vcf -b $roi_protein_coding_bed -wa -wb | \
        cut -f1,2,5,8,20 | awk '/protein_coding/' | \
        awk -v OFS=";" '\$1=\$1' | \
        awk 'BEGIN { FS=";"; OFS="\\t"} {\$1=\$1; print}' | \
   cut -f1,2,3,5,6,8,9,13 > ${sample_id}_annotatedcnv.csv

    # Generate plots and reports
   #cnv_html.R $calls_bed ${sample_id}_annotatedcnv.csv ${sample_id}_CNV_plot.pdf ${sample_id}_CNV_plot.html $sample_id

    CNV_function_new_update.R $calls_bed $cnv_genes_tuned $seg_bed \
        ${sample_id}_cnv_plot_full.pdf ${sample_id}_cnv_chr9.pdf ${sample_id}_cnv_chr7.pdf $sample_id

    # Process annotation files
    awk 'BEGIN { OFS="," } { gsub(/<[^>]+>/, substr(\$3, 2, length(\$3) - 2), \$3); print \$0 }' \
        ${sample_id}_annotatedcnv.csv > ${sample_id}_annotatedcnv_filter.csv

    awk 'BEGIN {print "Chrom,Start,Type,End,SVLEN,Score,LOG2CNT,Gene"} 1' \
        ${sample_id}_annotatedcnv_filter.csv > ${sample_id}_annotatedcnv_filter_header.csv

    # Run CNV mapping with threshold
    cnv_mapping_occfusion_update.py $seg_bed $roi_protein_coding_bed \
        ${sample_id}_tumor_copy.txt ${sample_id}_bins_filter.bed ${threshold}
    
    cnv_mapping_occfusion_update_nofilter.py $seg_bed \
        ${sample_id}_tumor_copy_number.txt ${threshold}
    """
}

// Annotation of Clair3 variant calling results
process clair3_annotate {
    label 'clair3'
    publishDir "${params.output_path}/routine_annotation/${sample_id}/snv_annotation", mode: "copy", overwrite: true

    input:
    tuple val(sample_id), path(clair3_output_dir), path(pileup_vcf), path(merge_vcf)

    output:
    tuple val(sample_id), path("${sample_id}_occ_pileup_snvs_avinput")
    tuple val(sample_id), path("${sample_id}_occ_pileup_annotateandfilter.csv"), emit:occpileupannotateandfilterout
    tuple val(sample_id), path("${sample_id}_occ_merge_snv_avinpt")
    tuple val(sample_id), path("${sample_id}_occ_merge.hg38_multianno.txt")
    tuple val(sample_id), path("${sample_id}_merge_annotateandfilter.csv"), emit:mergeannotateandfilterout
    tuple val(sample_id), path("${sample_id}_occ_pileup_annotateandfilter.csv"), path("${sample_id}_merge_annotateandfilter.csv"), emit:clair3output

    script:
    """
    convert2annovar.pl ${pileup_vcf} \
        --format vcf4 \
        --withfreq \
        --filter pass \
        --fraction 0.1 \
        --includeinfo \
        --outfile ${sample_id}_occ_pileup_snvs_avinput

    table_annovar.pl  ${sample_id}_occ_pileup_snvs_avinput \
        -outfile occ_pileup \
        -buildver hg38 -protocol refGene,clinvar_20240611,cosmic100coding2024\
        -operation g,f,f \
        ${params.humandb_dir} \
        -otherinfo

    awk '/exonic/ && /nonsynonymous/ && !/Benign/ && !/Likely_benign/|| /upstream/ || /Func.refGene/ || /splicing/ && !/Benign/ && !/Likely_benign/ || /Pathogenic/' occ_pileup.hg38_multianno.txt \
        | awk '/exonic/ || /TERT/ || /Func.refGene/ || /Pathogenic/'  \
        | awk '!/dist=166/' \
        | cut -f1-16,26,28,29 > ${sample_id}_occ_pileup_annotateandfilter.csv

    convert2annovar.pl ${merge_vcf} \
        --format vcf4 \
        --withfreq \
        --filter pass \
        --fraction 0.1 \
        --includeinfo \
        --outfile ${sample_id}_occ_merge_snv_avinpt

    table_annovar.pl ${sample_id}_occ_merge_snv_avinpt \
        -outfile occ_merge \
        -buildver hg38 -protocol refGene,clinvar_20240611,cosmic100coding2024\
        -operation g,f,f \
        ${params.humandb_dir} \
        -otherinfo

    awk '/exonic/ && /nonsynonymous/ && !/Benign/ && !/Likely_benign/|| /upstream/ || /Func.refGene/ || /splicing/ && !/Benign/ && !/Likely_benign/ || /frameshift/ && !/Benign/ && !/Likely_benign/ || /stopgain/ && !/Benign/ && !/Likely_benign/ || /Pathogenic/' \
        occ_merge.hg38_multianno.txt \
        | awk '/exonic/ || /TERT/ || /Func.refGene/ || /Pathogenic/'  \
        | awk '!/dist=166/' \
        | cut -f1-16,26,28,29 \
        > ${sample_id}_merge_annotateandfilter.csv

    cp occ_merge.hg38_multianno.txt ${sample_id}_occ_merge.hg38_multianno.txt
    """
}



// Somatic variant calling using ClairS-TO for tumor-only samples
// Annotation of ClairS-TO variant calling results
process clairs_to_annotate {
    label 'clairsto'
    publishDir "${params.output_path}/routine_annotation/${sample_id}/snv_annotation", mode: "copy", overwrite: true

    input:
    tuple val(sample_id), path('clairsto_output'), path(snv_vcf), path(indel_vcf)

    output:
    //tuple val(sample_id), path('clairsto_output/')
    tuple val(sample_id), path("${sample_id}_clairS_To_snv_avinput")
    tuple val(sample_id), path("${sample_id}_ClairS_TO_snv.hg38_multianno.txt")
    tuple val(sample_id), path("${sample_id}_annotateandfilter_clairsto.csv"), emit:annotateandfilter_clairstoout
    tuple val(sample_id), path("${sample_id}_merge_snv_indel_claisto.vcf.gz"), emit: clairsto_merged_vcf

    script:
    """
    # VCF index files (.tbi) should already exist from the epi2me:run_clairs_to process
    # bcftools concat requires these index files to be present
    # The directory is staged as 'clairsto_output' to access both VCF files and their .tbi index files
    bcftools concat -a -d all clairsto_output/snv.vcf.gz clairsto_output/indel.vcf.gz -Oz -o ${sample_id}_merge_snv_indel_claisto.vcf.gz

    convert2annovar.pl ${sample_id}_merge_snv_indel_claisto.vcf.gz \
        --format vcf4 \
        --filter pass \
        --includeinfo \
        --outfile  ${sample_id}_clairS_To_snv_avinput

    table_annovar.pl ${sample_id}_clairS_To_snv_avinput \
        -outfile ClairS_TO_snv \
        -buildver hg38 -protocol refGene,clinvar_20240611,cosmic100coding2024\
        -operation g,f,f \
        ${params.humandb_dir} \
        -otherinfo

    awk '/exonic/ && /nonsynonymous/ && !/Benign/ && !/Likely_benign/|| /upstream/ || /Func.refGene/ || /splicing/ && !/Benign/ && !/Likely_benign/ || /frameshift/ && !/Benign/ && !/Likely_benign/ || /stopgain/ && !/Benign/ && !/Likely_benign/ || /Pathogenic/' \
        ClairS_TO_snv.hg38_multianno.txt \
        | awk '/exonic/ || /TERT/ || /Func.refGene/ || /Pathogenic/'  \
        | awk '!/dist=166/' \
        | cut -f1-16,25,26  > ${sample_id}_annotateandfilter_clairsto.csv

    cp ClairS_TO_snv.hg38_multianno.txt ${sample_id}_ClairS_TO_snv.hg38_multianno.txt
    """
}

// Merge and filter annotations from Clair3 and ClairS-TO results
process merge_annotation {
    debug true
    label 'merge_annotation'
    publishDir "${params.output_path}/routine_annotation/${sample_id}/merge_annot_clair3andclairsto/", mode: "copy", overwrite: true

    input:
    tuple val(sample_id), path(merge_file), path(pileup_file), path(clairsto_file), path(occ_genes)
    
    output:
    //tuple val(sample_id), path("${sample_id}_merge_annotation_filter_snvs_allcall_filter.csv")
    tuple val(sample_id), path("${sample_id}_merge_annotation_filter_snvs_allcall.csv"), emit: occmergeout
    

    script:
    """
    #!/bin/bash
    set -e  # Exit on error

    echo "Input files:"
    echo "Sample ID: ${sample_id}"
    echo "Merge file: ${merge_file}"
    echo "Pileup file: ${pileup_file}"
    echo "ClairSTo file: ${clairsto_file}"
    echo "OCC genes file: ${occ_genes}"

    merge_annotations_prospective.R \\
        "${merge_file}" \\
        "${pileup_file}" \\
        "${clairsto_file}" \\
        "${sample_id}_merge_annotation_filter_snvs_allcall.csv" \\
        "${occ_genes}" \\
        ${sample_id}_merge_annotation_filter_snvs_allcall_filter.csv

    """
}


// TERT promoter variant visualization using IGV tools
process igv_tools {
    label 'epic'
    publishDir "${params.output_path}/routine_annotation/${sample_id}/coverage", mode: "copy", overwrite: true

    input:
    tuple val(sample_id), path(occ_bam), path(occ_bam_bai), path(tertp_variants), path(ncbirefseq), path(reference_genome), path(reference_genome_bai)

    output:
    tuple val(sample_id), file("${sample_id}_tertp_id1.html"), emit: tertp_out_igv

    script:
    """
    # Activate conda environment in the container
    source /opt/conda/etc/profile.d/conda.sh
    conda activate base
    
    export CURL_CA_BUNDLE=/etc/ssl/certs/ca-certificates.crt
    create_report $tertp_variants --fasta $reference_genome --flanking 1000 --tracks $tertp_variants $occ_bam $ncbirefseq --output ${sample_id}_tertp_id1.html
    ##cp "${sample_id}_tertp_id1.html" "${params.output_path}/report/${sample_id}_tertp_id1.html"
    """
}

// Quality assessment and statistics using Cramino
// Genomic region coverage plotting for EGFR, IDH1, and TERTp
process plot_genomic_regions {
    publishDir "${params.output_path}/routine_annotation/${sample_id}/coverage", mode: 'copy'
    label 'gviz'
    
    input:
    tuple val(sample_id), 
          path(gviz_data),
          path(bam_file),
          path(bam_index),
          path(cytoband_file)

    output:
    tuple val(sample_id), 
          path("${sample_id}_egfr_coverage.pdf"),
          path("${sample_id}_idh1_coverage.pdf"),
          path("${sample_id}_tertp_coverage.pdf"),
          path("${sample_id}_idh2_coverage.pdf"),
          emit: plot_genomic_regions_out

    script:
    """
    #!/bin/bash
    ls -l
    echo "Rscript path: \$(which Rscript)"
    echo "GViz file: ${gviz_data}"
    eval \"\$(micromamba shell hook --shell bash)\"
    micromamba activate gviz_env
    plot_genomic_regions.R \
        "${gviz_data}" \
        "${sample_id}" \
        "${bam_file}" \
        "${sample_id}_egfr_coverage.pdf" \
        "${sample_id}_idh1_coverage.pdf" \
        "${sample_id}_idh2_coverage.pdf" \
        "${sample_id}_tertp_coverage.pdf" \
        "${cytoband_file}"
    """
}

// Comprehensive PDF report generation using R Markdown
process markdown_report {
    publishDir "${params.result_path}/${sample_id}", mode: "copy", overwrite: true, pattern: "*.pdf"

    input:
    tuple val(sample_id),
          path(craminoreport),
          val(sample_id_file),
          path(dictionaire),
          path(logo),
          path(cnv_plot),
          path(tumor_number),
          path(annotatecnv),
          path(cnv_chr9),
          path(cnv_chr7),
          path(mgmt_results),
          path(merge_results),
          path(fusion_events),
          path(svannahtml),
          path(tertphtml),
          path(egfr_coverage),
          path(idh1_coverage),
          path(idh2_coverage),
          path(tertp_coverage),
          path(tsne_plot_file),
          path(nanodx_classifier),
          path(snv_target_genes),
          path(protein_coding_bed),
          path(rmd_template),
          path(warning_img),
          path(nanodx_classifier_pancan),
          path(dictionaire_pancan),
          path(tsne_plot_pancan_file)

    output:
    file("${sample_id}_markdown_pipeline_report.pdf")

    script:
    """
    # Output PDF path
    output_file="${sample_id}_markdown_pipeline_report.pdf"

    # Handle different run modes
    if [ "${params.run_mode_order}" = "true" ]; then
        # For run_mode_order: Create sample_file.txt with sample_id and threshold_value
        THRESHOLD_VALUE=""
        if [ -f "${params.output_path}/routine_annotation/${sample_id}/cnv/ace/${sample_id}_ace_results/threshold_value.txt" ]; then
            THRESHOLD_VALUE=\$(cat "${params.output_path}/routine_annotation/${sample_id}/cnv/ace/${sample_id}_ace_results/threshold_value.txt")
            # Create sample_file.txt
            echo -e "${sample_id}\\t\${THRESHOLD_VALUE}" > sample_file.txt
            echo "Created sample_file.txt for run_mode_order with: ${sample_id} \${THRESHOLD_VALUE}"
            SAMPLE_FILE="\${PWD}/sample_file.txt"
        else
            echo "ERROR: ACE threshold file not found for ${sample_id} at ${params.output_path}/routine_annotation/${sample_id}/cnv/ace/${sample_id}_ace_results/threshold_value.txt"
            exit 1
        fi
    else
        # For run_mode_annotation: Check if user provided tumor content first (takes priority)
        if [ ! -f "${sample_id_file}" ]; then
            echo "ERROR: Sample ID file not found: ${sample_id_file}"
            exit 1
        fi

        # Count columns in the sample_id_file
        NUM_COLS=\$(awk '{print NF; exit}' "${sample_id_file}")
        echo "Detected \$NUM_COLS columns in sample_id_file: ${sample_id_file}"

        if [ "\$NUM_COLS" -ge 2 ]; then
            # User provided 2 columns (sample_id + tumor content) - use user-provided value
            # Extract the line for this sample_id
            grep "^${sample_id}[[:space:]]" "${sample_id_file}" > sample_file.txt || {
                echo "WARNING: Sample ${sample_id} not found in ${sample_id_file}, creating file with sample_id only"
                echo "${sample_id}" > sample_file.txt
            }
            SAMPLE_FILE="\${PWD}/sample_file.txt"
            echo "Created local sample_file.txt with user-provided tumor content:"
            cat "\${SAMPLE_FILE}"
        elif [ -f "${params.output_path}/routine_annotation/${sample_id}/cnv/ace/${sample_id}_ace_results/threshold_value.txt" ]; then
            # User provided only 1 column, but ACE results available - use ACE-calculated value
            THRESHOLD_VALUE=\$(cat "${params.output_path}/routine_annotation/${sample_id}/cnv/ace/${sample_id}_ace_results/threshold_value.txt")
            echo -e "${sample_id}\\t\${THRESHOLD_VALUE}" > sample_file.txt
            echo "Created sample_file.txt with ACE-calculated tumor content: ${sample_id} \${THRESHOLD_VALUE}"
            SAMPLE_FILE="\${PWD}/sample_file.txt"
        else
            # Only 1 column and no ACE results - copy as-is
            cp "${sample_id_file}" sample_file.txt
            SAMPLE_FILE="\${PWD}/sample_file.txt"
            echo "Created local sample_file.txt (single column, no tumor content)"
            cat "\${SAMPLE_FILE}"
        fi
    fi

    # Now call the Rscript with the updated Rmd file using absolute paths
    # Activate conda environment and find Rscript
    source /opt/conda/etc/profile.d/conda.sh
    conda activate markdown_env
    
    # Try to find Rscript
    RSCRIPT_PATH=""
    if command -v Rscript >/dev/null 2>&1; then
        RSCRIPT_PATH="Rscript"
    elif [ -f "/opt/conda/envs/markdown_env/bin/Rscript" ]; then
        RSCRIPT_PATH="/opt/conda/envs/markdown_env/bin/Rscript"
    elif [ -f "/opt/conda/bin/Rscript" ]; then
        RSCRIPT_PATH="/opt/conda/bin/Rscript"
    elif [ -f "/usr/local/bin/Rscript" ]; then
        RSCRIPT_PATH="/usr/local/bin/Rscript"
    elif [ -f "/usr/bin/Rscript" ]; then
        RSCRIPT_PATH="/usr/bin/Rscript"
    else
        echo "ERROR: Rscript not found. Available R-related files:"
        find /opt/conda -name "*R*" 2>/dev/null | head -10
        find /usr -name "*Rscript*" 2>/dev/null | head -10
        exit 1
    fi
    
    echo "Using Rscript at: \$RSCRIPT_PATH"

    \$RSCRIPT_PATH -e "rmarkdown::render('\${PWD}/${rmd_template}', output_file=commandArgs(trailingOnly=TRUE)[24])" \
      "${sample_id}" \
      "\${PWD}/${craminoreport}" \
      "\${SAMPLE_FILE}" \
      "\${PWD}/${nanodx_classifier}" \
      "\${PWD}/${dictionaire}" \
      "\${PWD}/${logo}" \
      "\${PWD}/${cnv_plot}" \
      "\${PWD}/${tumor_number}" \
      "\${PWD}/${annotatecnv}" \
      "\${PWD}/${cnv_chr9}" \
      "\${PWD}/${cnv_chr7}" \
      "\${PWD}/${mgmt_results}" \
      "\${PWD}/${merge_results}" \
      "\${PWD}/${fusion_events}" \
      "\${PWD}/${tertphtml}" \
      "\${PWD}/${svannahtml}" \
      "\${PWD}/${egfr_coverage}" \
      "\${PWD}/${idh1_coverage}" \
      "\${PWD}/${idh2_coverage}" \
      "\${PWD}/${tertp_coverage}" \
      "\${PWD}/${tsne_plot_file}" \
      "\${PWD}/${snv_target_genes}" \
      "\${PWD}/${protein_coding_bed}" \
      "\${PWD}/\${output_file}" \
      "${workflow.manifest.version}" \
      "${params.snv_depth_threshold}" \
      "${params.snv_gq_threshold}" \
      "\${PWD}/${warning_img}" \
      "\${PWD}/${nanodx_classifier_pancan}" \
      "\${PWD}/${dictionaire_pancan}" \
      "\${PWD}/${tsne_plot_pancan_file}"

    # Flatten PDF with Ghostscript for maximum printer compatibility (soft — skips if gs absent)
    if command -v gs >/dev/null 2>&1; then
        echo "Ghostscript found — flattening PDF to PDF 1.4 for printer compatibility..."
        gs -q -dBATCH -dNOPAUSE -dSAFER \
           -sDEVICE=pdfwrite \
           -dCompatibilityLevel=1.4 \
           -dEmbedAllFonts=true \
           -dSubsetFonts=true \
           -dPrinted=false \
           -sOutputFile="\${output_file}_flat.pdf" "\${output_file}" \
        && mv "\${output_file}_flat.pdf" "\${output_file}" \
        && echo "PDF flattened successfully." \
        || echo "WARNING: Ghostscript flattening failed — keeping original PDF."
    else
        echo "Ghostscript not found on host — skipping PDF flattening (PDF 1.4 pdflatex flag still active)."
    fi

    # Clean up intermediate files and folders created by R Markdown (belt and suspenders)
    echo "Cleaning up intermediate R Markdown files..."
    rm -rf *_files *.html *.log *.tex *_cache 2>/dev/null || true

    echo "RMD report generated successfully"


    # This removes any intermediate files that might have been published before pattern filter
    echo "Final cleanup of report directory..."
    find ${params.output_path}/report -type d -name "*_files" -exec rm -rf {} + 2>/dev/null || true
    find ${params.output_path}/report -type f -name "*.tex" -delete 2>/dev/null || true
    find ${params.output_path}/report -type f -name "*.html" -delete 2>/dev/null || true
    find ${params.output_path}/report -type f -name "*.log" -delete 2>/dev/null || true
    echo "Cleanup complete - only PDF files remain"
    """
}

// Copy result files to routine_results folder for easy access
process copy_results_to_summary {
    publishDir "${params.result_path}/${sample_id}", mode: "copy", overwrite: true, pattern: "${sample_id}_*"
    publishDir "${params.result_path}/${sample_id}", mode: "copy", overwrite: true, pattern: "roi.protein_coding_*.bed"

    input:
    tuple val(sample_id), path(mnpflex_bed)
    tuple val(sample_id), path(sturgeon_pdf)
    tuple val(sample_id), path(tsne_html)
    tuple val(sample_id), path(svanna_html)
    tuple val(sample_id), path(tsne_pancan_html)
    tuple val(sample_id), path(tsne_pancan_pdf)
    path(roi_bed)

    output:
    path("${sample_id}_mnpflex_input.bed"), optional: true
    path("${sample_id}_bedmethyl_sturgeon_general.pdf"), optional: true
    path("${sample_id}_tsne_plot.html"), optional: true
    path("${sample_id}_occ_svanna_annotation.html"), optional: true
    path("${sample_id}_tsne_plot_pancan.html"), optional: true
    path("${sample_id}_tsne_plot_pancan.pdf"), optional: true
    path("roi.protein_coding_*.bed"), optional: true

    script:
    """
    # Copy files dereferencing symlinks
    # Use temp names first, then rename to avoid "same file" error
    if [ -f "${mnpflex_bed}" ]; then
        cp -L ${mnpflex_bed} temp_mnpflex.bed
        mv temp_mnpflex.bed ${sample_id}_mnpflex_input.bed
        echo "Created ${sample_id}_mnpflex_input.bed"
    fi

    if [ -f "${sturgeon_pdf}" ]; then
        cp -L ${sturgeon_pdf} temp_sturgeon.pdf
        mv temp_sturgeon.pdf ${sample_id}_bedmethyl_sturgeon_general.pdf
        echo "Created ${sample_id}_bedmethyl_sturgeon_general.pdf"
    fi

    if [ -f "${tsne_html}" ]; then
        cp -L ${tsne_html} temp_tsne.html
        mv temp_tsne.html ${sample_id}_tsne_plot.html
        echo "Created ${sample_id}_tsne_plot.html"
    fi

    if [ -f "${svanna_html}" ]; then
        cp -L ${svanna_html} temp_svanna.html
        mv temp_svanna.html ${sample_id}_occ_svanna_annotation.html
        echo "Created ${sample_id}_occ_svanna_annotation.html"
    fi

    if [ -f "${tsne_pancan_html}" ]; then
        cp -L ${tsne_pancan_html} temp_tsne_pancan.html
        mv temp_tsne_pancan.html ${sample_id}_tsne_plot_pancan.html
        echo "Created ${sample_id}_tsne_plot_pancan.html"
    fi

    if [ -f "${tsne_pancan_pdf}" ]; then
        cp -L ${tsne_pancan_pdf} temp_tsne_pancan.pdf
        mv temp_tsne_pancan.pdf ${sample_id}_tsne_plot_pancan.pdf
        echo "Created ${sample_id}_tsne_plot_pancan.pdf"
    fi

    # Copy ROI BED file with timestamp for run traceability
    TIMESTAMP=\$(date +%d_%m_%Y_%H%M)
    if [ -f "${roi_bed}" ]; then
        cp -L ${roi_bed} "roi.protein_coding_\${TIMESTAMP}.bed"
        echo "Created roi.protein_coding_\${TIMESTAMP}.bed"
    fi

    echo "Files ready for publishing:"
    ls -lh ${sample_id}_* roi.protein_coding_*.bed 2>/dev/null || echo "No output files created"
    """
}

// ACE tumor content calculation and copy number analysis
process ace_tmc {
    label 'ace_tmc'
    publishDir "${params.output_path}/routine_annotation/${sample_id}/cnv/ace/", mode: "copy", overwrite: true
    
    input:
    tuple val(sample_id), path(rds_file)
    
    output:
        tuple val(sample_id), path("${sample_id}_ace_results"), emit: aceresults
    tuple val(sample_id), env(threshold_value), emit: threshold_value
    
    script:
    """
    #!/bin/bash
    set -e

    # Check if we're in a container and use appropriate conda setup
    if [ -f "/opt/conda/etc/profile.d/conda.sh" ]; then
        source /opt/conda/etc/profile.d/conda.sh
        conda activate ace_env
    else
    source activate ace_env
    fi
    
    # Debug info
    echo "Processing sample: ${sample_id}"
    echo "Input RDS file: ${rds_file}"
    ls -l ${rds_file}

    # Create output directory
    mkdir -p ${sample_id}_ace_results
    
    # Run ACE TMC analysis
    ace_tmc.R "${rds_file}" "${sample_id}_ace_results" "${sample_id}"
    
    # Read and export the threshold value
    threshold_value=\$(cat "${sample_id}_ace_results/threshold_value.txt")
    echo "Threshold value for ${sample_id}: \$threshold_value"
    """
}

//---------------------------------------------------------------------
// Static Reference File Channels (Module Level - Created Once)
//---------------------------------------------------------------------
// Create value channels for static reference files at module level
// to prevent cache invalidation. These are created ONCE when the module
// is loaded, not every time the workflow runs for a new sample.

def roi_protein_coding_bed_ch = Channel.value(file(params.roi_protein_coding_bed))
def nanodx_450model_ch = Channel.value(file(params.nanodx_450model))
def snakefile_nanodx_ch = Channel.value(file(params.snakefile_nanodx))
def nn_model_ch = Channel.value(file(params.nn_model))
def nanodxcolormap_ch = Channel.value(file(params.nanodxcolormap))
def nanodxh5_ch = Channel.value(file(params.nanodxh5))
def nanodx_pancan_model_ch = Channel.value(file(params.nanodx_pancan_model))
def nanodxcolormap_pancan_ch = Channel.value(file(params.nanodxcolormap_pancan))
def nanodxh5_pancan_ch = Channel.value(file(params.nanodxh5_pancan))
def hg19_450model_ch = Channel.value(file(params.hg19_450model))
def vcf2circos_json_ch = Channel.value(file(params.vcf2circos_json))
def genecode_bed_ch = Channel.value(file(params.genecode_bed))
def occ_fusion_genes_list_ch = Channel.value(file(params.occ_fusion_genes_list))
def refgene_ch = Channel.value(file(params.refgene))
def hg38_refgenemrna_ch = Channel.value(file(params.hg38_refgenemrna))
def clinvar_ch = Channel.value(file(params.clinvar))
def clinvarindex_ch = Channel.value(file(params.clinvarindex))
def hg38_cosmic100_ch = Channel.value(file(params.hg38_cosmic100))
def hg38_cosmic100index_ch = Channel.value(file(params.hg38_cosmic100index))
def tertp_variants_ch = Channel.value(file(params.tertp_variants))
def ncbirefseq_ch = Channel.value(file(params.ncbirefseq))
def gviz_data_ch = Channel.value(file(params.gviz_data))
def cytoband_file_ch = Channel.value(file(params.cytoband_file))
def epicsites_ch = Channel.value(file(params.epicsites))
def mgmt_cpg_island_hg38_ch = Channel.value(file(params.mgmt_cpg_island_hg38))
def _sturgeon_file = file(params.sturgeon_model)
def sturgeon_available = _sturgeon_file.exists()
def sturgeon_model_ch = sturgeon_available \
    ? Channel.value(_sturgeon_file) \
    : Channel.empty()
if (!sturgeon_available) {
    log.warn "Sturgeon model not found (${params.sturgeon_model}) — Sturgeon classification will be skipped. Download general.zip from Zenodo and place it in ${params.ref_dir}/"
}

//---------------------------------------------------------------------
// Workflow definition
//---------------------------------------------------------------------

workflow annotation {
    take:
        input_data

    main:
        validateParameters()

        // Define fusion events channel conditionally to avoid undefined output errors
        def fusion_events_channel = Channel.empty()


        // Initialize channels as empty by default
        def annotatecnv_results = null
        def annotatecnv_out = Channel.empty()
        def merge_annotation_out = Channel.empty()
        def crossNN_out = Channel.empty()
        def mgmt_promoter_out = Channel.empty()
        def svannasv_out = Channel.empty()
        def igv_tools_out = Channel.empty()
        def cramino_report_out = Channel.empty()

        // NOTE: Static reference file channels are now defined at module level (before workflow block)
        // to ensure they're created once and reused across all samples for optimal caching

        // Initialize sample_thresholds based on run mode
        def sample_thresholds
        if (params.run_mode_epiannotation) {
            // For epiannotation mode, load samples from epi2me_sample_id_file
            sample_thresholds = file(params.epi2me_sample_id_file).readLines()
                .collectEntries { line ->
                    def sample_id = line.trim()
                    [(sample_id): null]  // null indicates ACE calculation needed
                }
            println "Sample thresholds (from epi2me_sample_id_file): ${sample_thresholds}"
        } else if (params.run_mode_order) {
            // Load sample IDs from bam_sample_id_file and initialize with null thresholds
            def sample_ids = file(params.bam_sample_id_file).readLines().collect { line ->
                // Extract only the sample ID part (first column), removing flow cell ID
                line.trim().split(/\s+/)[0]
            }
            // Initialize sample_thresholds with null values for ACE calculation
            sample_thresholds = sample_ids.collectEntries { [(it): null] }
            println "Sample thresholds (run_mode_order): ${sample_thresholds}"
        } else {
            sample_thresholds = loadSampleThresholds()
            println "Sample thresholds: ${sample_thresholds}"
        }

        // Create segsfromepi2me channel based on mode
        boosts_segsfromepi2me_channel = (params.run_mode_order || params.run_mode_epiannotation) ?
            input_data.combine(roi_protein_coding_bed_ch).map { args ->
                def sample_id = args[0]
                def bam = args[1]
                def bai = args[2]
                def ref = args[3]
                def ref_bai = args[4]
                def tr_bed = args[5]
                def modkit = args[6]
                def segs_bed = args[7]
                def bins_bed = args[8]
                def segs_vcf = args[9]
                def rds_file = args[10]
                def sv = args[11]
                def occ_bed = args[12]

                tuple(
                    sample_id,
                    segs_vcf,
                    occ_bed,
                    bins_bed,
                    segs_bed,
                    sample_thresholds[sample_id]
                )
            } :
            Channel.fromList(sample_thresholds.keySet().collect())
                .combine(roi_protein_coding_bed_ch)
                .map { sample_id, occ_bed ->
                    tuple(
                        sample_id,
                        file("${params.segsfromepi2me_folder}/${sample_id}/${sample_id}_segs.vcf"),
                        occ_bed,
                        file("${params.segsfromepi2me_folder}/${sample_id}/${sample_id}_bins.bed"),
                        file("${params.segsfromepi2me_folder}/${sample_id}/${sample_id}_segs.bed"),
                        sample_thresholds[sample_id]
                    )
                }

        boosts_svanna_channel = (params.run_mode_order || params.run_mode_epiannotation) ?
            input_data.combine(roi_protein_coding_bed_ch).map { args ->
                def sample_id = args[0]
                def bam = args[1]
                def bai = args[2]
                def ref = args[3]
                def ref_bai = args[4]
                def tr_bed = args[5]
                def modkit = args[6]
                def segs_bed = args[7]
                def bins_bed = args[8]
                def segs_vcf = args[9]
                def rds_file = args[10]
                def sv = args[11]
                def occ_bed = args[12]

                //log.info "Creating Svanna input for sample: ${sample_id} (order mode)"
                //log.info "SV file path: ${sv}"

                // Create index file path from the published SV file location
                def sv_file = file("${params.output_path}/routine_epi2me/${sample_id}/${sample_id}.sniffles.vcf.gz")
                def sv_index = file("${params.output_path}/routine_epi2me/${sample_id}/${sample_id}.sniffles.vcf.gz.tbi")

                tuple(
                    sample_id,
                    sv_file,               // SV VCF file from epi2me
                    sv_index,              // Index file
                    occ_bed
                )
            } :
            Channel.fromList(sample_thresholds.keySet().collect())
                .combine(roi_protein_coding_bed_ch)
                .map { sample_id, occ_bed ->
                    //log.info "Creating Svanna input for sample: ${sample_id} (standalone mode)"
                    def sv_file = file("${params.sv_folder}/${sample_id}/${sample_id}.sniffles.vcf.gz")

                    //if (!sv_file.exists()) {
                    //    error "SV file not found: ${sv_file}"
                    //}

                    tuple(
                        sample_id,
                        sv_file,
                        file("${sv_file}.tbi"),
                        occ_bed
                    )
                }

        boosts_clair3_channel = (params.run_mode_order || params.run_mode_epiannotation) ?
            input_data.combine(refgene_ch).combine(hg38_refgenemrna_ch).combine(clinvar_ch)
                .combine(clinvarindex_ch).combine(hg38_cosmic100_ch).combine(hg38_cosmic100index_ch)
                .map { args ->
                def sample_id = args[0]
                def bam = args[1]
                def bai = args[2]
                def ref = args[3]
                def ref_bai = args[4]
                def tr_bed = args[5]
                def modkit = args[6]
                def segs_bed = args[7]
                def bins_bed = args[8]
                def segs_vcf = args[9]
                def rds_file = args[10]
                def sv = args[11]
                def refgene = args[12]
                def refgenemrna = args[13]
                def clinvar = args[14]
                def clinvarindex = args[15]
                def cosmic100 = args[16]
                def cosmic100index = args[17]

                tuple(
                    sample_id,
                    bam,
                    bai,
                    ref,
                    ref_bai,
                    refgene,
                    refgenemrna,
                    clinvar,
                    clinvarindex,
                    cosmic100,
                    cosmic100index
                )
            } :
            Channel.fromList(sample_thresholds.keySet().collect())
                .combine(refgene_ch).combine(hg38_refgenemrna_ch).combine(clinvar_ch)
                .combine(clinvarindex_ch).combine(hg38_cosmic100_ch).combine(hg38_cosmic100index_ch)
                .map { sample_id, refgene, refgenemrna, clinvar, clinvarindex, cosmic100, cosmic100index ->
                    tuple(
                        sample_id,
                        file("${params.roi_bam_folder}/${sample_id}.roi.bam"),
                        file("${params.roi_bam_folder}/${sample_id}.roi.bam.bai"),
                        file(params.reference_genome),
                        file(params.reference_genome_bai),
                        refgene,
                        refgenemrna,
                        clinvar,
                        clinvarindex,
                        cosmic100,
                        cosmic100index
                    )
                }

        boosts_clairSTo_channel = (params.run_mode_order || params.run_mode_epiannotation) ?
            input_data.combine(refgene_ch).combine(hg38_refgenemrna_ch).combine(clinvar_ch)
                .combine(clinvarindex_ch).combine(hg38_cosmic100_ch).combine(hg38_cosmic100index_ch)
                .combine(roi_protein_coding_bed_ch).map { args ->
                def sample_id = args[0]
                def bam = args[1]
                def bai = args[2]
                def ref = args[3]
                def ref_bai = args[4]
                def tr_bed = args[5]
                def modkit = args[6]
                def segs_bed = args[7]
                def bins_bed = args[8]
                def segs_vcf = args[9]
                def rds_file = args[10]
                def sv = args[11]
                def refgene = args[12]
                def refgenemrna = args[13]
                def clinvar = args[14]
                def clinvarindex = args[15]
                def cosmic100 = args[16]
                def cosmic100index = args[17]
                def occ_bed = args[18]

                tuple(
                    sample_id,
                    bam,
                    bai,
                    ref,
                    ref_bai,
                    refgene,
                    refgenemrna,
                    clinvar,
                    clinvarindex,
                    cosmic100,
                    cosmic100index,
                    occ_bed
                )
            } :
            Channel.fromList(sample_thresholds.keySet().collect())
                .combine(refgene_ch).combine(hg38_refgenemrna_ch).combine(clinvar_ch)
                .combine(clinvarindex_ch).combine(hg38_cosmic100_ch).combine(hg38_cosmic100index_ch)
                .combine(roi_protein_coding_bed_ch)
                .map { sample_id, refgene, refgenemrna, clinvar, clinvarindex, cosmic100, cosmic100index, occ_bed ->
                    tuple(
                        sample_id,
                        file("${params.roi_bam_folder}/${sample_id}.roi.bam"),
                        file("${params.roi_bam_folder}/${sample_id}.roi.bam.bai"),
                        file(params.reference_genome),
                        file(params.reference_genome_bai),
                        refgene,
                        refgenemrna,
                        clinvar,
                        clinvarindex,
                        cosmic100,
                        cosmic100index,
                        occ_bed

                    )
                }

        boosts_igv_channel = (params.run_mode_order || params.run_mode_epiannotation) ?
            input_data.combine(tertp_variants_ch).combine(ncbirefseq_ch).map { args ->
                def sample_id = args[0]
                def bam = args[1]
                def bai = args[2]
                def ref = args[3]
                def ref_bai = args[4]
                def tr_bed = args[5]
                def modkit = args[6]
                def segs_bed = args[7]
                def bins_bed = args[8]
                def segs_vcf = args[9]
                def rds_file = args[10]
                def sv = args[11]
                def tertp_var = args[12]
                def ncbi = args[13]

                tuple(sample_id, bam, bai, tertp_var, ncbi, ref, ref_bai)
            } :
            Channel.fromList(sample_thresholds.keySet().collect())
                .combine(tertp_variants_ch).combine(ncbirefseq_ch)
                .map { sample_id, tertp_var, ncbi ->
                    tuple(sample_id, file("${params.roi_bam_folder}/${sample_id}.roi.bam"),
                          file("${params.roi_bam_folder}/${sample_id}.roi.bam.bai"),
                          tertp_var,
                          ncbi,
                          file(params.reference_genome),
                          file("${params.reference_genome}.fai"))
                }

        boosts_plot_genomic_regions_channel = (params.run_mode_order || params.run_mode_epiannotation) ?
            input_data.combine(gviz_data_ch).combine(cytoband_file_ch).map { args ->
                def sample_id = args[0]
                def bam = args[1]
                def bai = args[2]
                def ref = args[3]
                def ref_bai = args[4]
                def tr_bed = args[5]
                def modkit = args[6]
                def segs_bed = args[7]
                def bins_bed = args[8]
                def segs_vcf = args[9]
                def rds_file = args[10]
                def sv = args[11]
                def gviz = args[12]
                def cytoband = args[13]

                tuple(
                    sample_id,
                    gviz,
                    bam,
                    bai,
                    cytoband
                )
            } :
            Channel.fromList(sample_thresholds.keySet().collect())
                .combine(gviz_data_ch).combine(cytoband_file_ch)
                .map { sample_id, gviz, cytoband ->
                    tuple(
                        sample_id,
                        gviz,
                        file("${params.roi_bam_folder}/${sample_id}.roi.bam"),
                        file("${params.roi_bam_folder}/${sample_id}.roi.bam.bai"),
                        cytoband
                    )
                }

        boosts_cramino = (params.run_mode_order || params.run_mode_epiannotation) ?
            input_data.map { args ->
                def sample_id = args[0]
                def occ_bam = args[1]    
                def occ_bai = args[2]
                def ref = args[3]
                def ref_bai = args[4]

                // IMPORTANT: Cramino should analyze the FULL merged BAM, not the ROI BAM
                // Read merged BAM from published directory
                def bam_file = file("${params.merge_bam_folder}/${sample_id}.merged.bam")
                def bai_file = file("${params.merge_bam_folder}/${sample_id}.merged.bam.bai")

                // If exact match doesn't exist, try .merge.bam specifically (legacy)
                if (!bam_file.exists()) {
                    bam_file = file("${params.merge_bam_folder}/${sample_id}.merged.bam")
                    bai_file = file("${params.merge_bam_folder}/${sample_id}.merged.bam.bai")
                }

                // If still not found, try wildcard but exclude .roi.bam files
                if (!bam_file.exists()) {
                    def pattern_files = file("${params.merge_bam_folder}/${sample_id}.*.bam")
                    def potential_bams = pattern_files instanceof List ? pattern_files : [pattern_files]
                    // Filter out .roi.bam files
                    def filtered_bams = potential_bams.findAll { !it.name.contains('.roi.') }
                    if (filtered_bams.size() > 0) {
                        bam_file = filtered_bams[0]
                        bai_file = file("${bam_file}.bai")
                    }
                }

                // Only return valid tuples
                if (bam_file.exists() && bai_file.exists()) {
                    tuple(
                        sample_id,
                        bam_file,    // <- Now using merged BAM
                        bai_file,
                        ref,
                        ref_bai
                    )
                } else {
                    println "WARNING: Cramino - Merged BAM file not found for ${sample_id}. Tried: ${sample_id}.merged.bam, ${sample_id}.merge.bam, ${sample_id}.*.bam"
                    null
                }
            }
            .filter { it != null } :
            Channel.fromList(sample_thresholds.keySet().collect())
                .map { sample_id -> 
                    // Try exact match first, then wildcard pattern (avoid .roi.bam files)
                    def bam_file = file("${params.merge_bam_folder}/${sample_id}.merged.bam")
                    def bai_file = file("${params.merge_bam_folder}/${sample_id}.merged.bam.bai")

                    // If exact match doesn't exist, try .merge.bam specifically (legacy)
                    if (!bam_file.exists()) {
                        bam_file = file("${params.merge_bam_folder}/${sample_id}.merge.bam")
                        bai_file = file("${params.merge_bam_folder}/${sample_id}.merge.bam.bai")
                    }

                    // If still not found, try wildcard but exclude .roi.bam files
                    if (!bam_file.exists()) {
                        def pattern_files = file("${params.merge_bam_folder}/${sample_id}.*.bam")
                        def potential_bams = pattern_files instanceof List ? pattern_files : [pattern_files]
                        // Filter out .roi.bam files
                        def filtered_bams = potential_bams.findAll { !it.name.contains('.roi.') }
                        if (filtered_bams.size() > 0) {
                            bam_file = filtered_bams[0]
                            bai_file = file("${bam_file}.bai")
                        }
                    }

                    // Only process samples that have corresponding BAM files
                    if (bam_file.exists() && bai_file.exists()) {
                        tuple(
                            sample_id,
                            bam_file,
                            bai_file,
                            file(params.reference_genome, checkIfExists: true),
                            file("${params.reference_genome}.fai", checkIfExists: true)
                        )
                    } else {
                        println "WARNING: Skipping sample ${sample_id} - BAM file not found. Tried: ${sample_id}.merged.bam, ${sample_id}.merge.bam, ${sample_id}.*.bam (excluding .roi.bam)"
                        null
                    }
                }
                .filter { it != null }

        // MGMT analysis
        if (params.run_mode in ['mgmt', 'rmd', 'all'] || params.run_mode_order) {
            println "Running MGMT Analysis..."
            
            // Create MGMT channel that works for combined modes (order and epiannotation)
            mgmt_ch = (params.run_mode_order || params.run_mode_epiannotation) ?
                input_data.combine(epicsites_ch).combine(mgmt_cpg_island_hg38_ch).map { args ->
                    def sample_id = args[0]
                    def bam = args[1]
                    def bai = args[2]
                    def ref = args[3]
                    def ref_bai = args[4]
                    def tr_bed = args[5]
                    def modkit = args[6]
                    def segs_bed = args[7]
                    def bins_bed = args[8]
                    def segs_vcf = args[9]
                def rds_file = args[10]
                    def sv = args[11]
                    def epicsites = args[12]
                    def mgmt_cpg = args[13]

                    tuple(
                        sample_id,
                        modkit,
                        epicsites,
                        mgmt_cpg
                    )
                } :
                Channel.fromList(sample_thresholds.keySet().collect())
                    .combine(epicsites_ch).combine(mgmt_cpg_island_hg38_ch)
                    .map { sample_id, epicsites, mgmt_cpg ->
                        tuple(
                            sample_id,
                            file("${params.bedmethyl_folder}/${sample_id}/${sample_id}.wf_mods.bedmethyl.gz"),
                            epicsites,
                            mgmt_cpg
                        )
                    }
            
            // Run MGMT related processes
            extract_epic(mgmt_ch)
            
            // Create channels for downstream processes
            MGMT_output = extract_epic.out.MGMTheaderout
            
            MGMT_sturgeon = extract_epic.out.sturgeonbedinput
                 .combine(sturgeon_model_ch)
                 .map { sample_id, sturgeoninput, sturgeon_model ->
                     tuple(sample_id, sturgeoninput, sturgeon_model)
                 }

            mgmt_nanodx = extract_epic.out.epicselectnanodxinput
                .combine(hg19_450model_ch)
                .map { sample_id, epicselectnanodxinput, hg19_model ->
                    tuple(sample_id, epicselectnanodxinput, hg19_model)
                }

            // Run the processes
            if (sturgeon_available) { sturgeon(MGMT_sturgeon) }
            mgmt_promoter(MGMT_output)
            nanodx(mgmt_nanodx)

            nanodx_out = nanodx.out.nanodx450out
                .combine(nanodx_450model_ch).combine(snakefile_nanodx_ch).combine(nn_model_ch)
                .map { sample_id, nanodx450out, nanodx_model, snakefile, nn_model ->
                    tuple(
                        sample_id,
                        nanodx450out,
                        nanodx_model,
                        snakefile,
                        nn_model
                    )
                }

            crossNN(nanodx_out)

            // Pan-cancer classifier (always runs alongside Capper)
            nanodx_pancan_out = nanodx.out.nanodx450out
                .combine(nanodx_pancan_model_ch).combine(snakefile_nanodx_ch).combine(nn_model_ch)
                .map { sample_id, nanodx450out, nanodx_model, snakefile, nn_model ->
                    tuple(sample_id, nanodx450out, nanodx_model, snakefile, nn_model)
                }
            crossNN_pancan(nanodx_pancan_out)

            // Generate t-SNE plots for both classifiers
            tsne_plot(
                extract_epic.out.epicselectnanodxinput,
                nanodxcolormap_ch,
                nanodxh5_ch
            )
            tsne_plot_pancan(
                extract_epic.out.epicselectnanodxinput,
                nanodxcolormap_pancan_ch,
                nanodxh5_pancan_ch
            )
        }

        // Svanna analysis
        if (params.run_mode in ['svannasv', 'rmd', 'all'] || params.run_mode_order) {
            println "Running Svanna Analysis..."
            
            svannasv(boosts_svanna_channel)
            svannasv_out = svannasv.out.occsvannaannotationannotationvcf
                .combine(vcf2circos_json_ch)
                .map{ sample_id, svannavcfoutput, vcf2circos ->
                [sample_id, svannavcfoutput, vcf2circos]
            }
            circosplot(svannasv_out)
            circosplot_out=circosplot.out.circosout
            svannaoutfusion_events= svannasv.out.occsvannaannotationannotationvcf
                .combine(genecode_bed_ch).combine(occ_fusion_genes_list_ch)
                .map{ sample_id, occsvannaannotationannotationvcf, genecode, fusion_genes ->
                [sample_id, occsvannaannotationannotationvcf, genecode, fusion_genes]
            }
            svannasv_fusion_events(svannaoutfusion_events)
            
            // Assign fusion events channel when SV analysis is run
            fusion_events_channel = svannasv_fusion_events.out.filterfusioneventout

            // Copy result files to routine_results for easy access
            if (params.run_mode_order || params.run_mode_epiannotation || params.run_mode_annotation) {
                def sturgeon_pdf_ch = sturgeon_available
                    ? sturgeon.out.sturgeon_pdf
                    : extract_epic.out.mnpflex_bed.map { sid, f -> tuple(sid, file("NO_STURGEON_PDF")) }
                copy_results_to_summary(
                    extract_epic.out.mnpflex_bed,
                    sturgeon_pdf_ch,
                    tsne_plot.out.tsne_html,
                    svannasv.out.rmdsvannahtml,
                    tsne_plot_pancan.out.tsne_pancan_html,
                    tsne_plot_pancan.out.tsne_pancan_out,
                    file(params.roi_protein_coding_bed)
                )
            }
        }

        // OCC analysis
        // NOTE: Variant calling (Clair3/ClairS-TO) is now handled by epi2me workflow
        // This section reads pre-existing VCF files from epi2me output and runs annotation
        if (params.run_mode in ['roi', 'rmd', 'all'] || params.run_mode_order) {
            println "Running ROI Analysis (annotation of pre-existing VCFs from epi2me)..."

            // Step 1: Create channels for VCF files
            // In epiannotation/order mode, derive from input_data to ensure dependency on epi2me completion
            // In standalone mode, create from sample list
            def clair3_annot_input = (params.run_mode_epiannotation || params.run_mode_order) ?
                input_data.map { args ->
                    def sample_id = args[0]
                    def clair3_output_dir = file("${params.output_path}/routine_epi2me/${sample_id}/output_clair3")
                    def pileup_vcf = file("${params.output_path}/routine_epi2me/${sample_id}/output_clair3/pileup.vcf.gz")
                    def merge_vcf = file("${params.output_path}/routine_epi2me/${sample_id}/output_clair3/merge_output.vcf.gz")
                    tuple(sample_id, clair3_output_dir, pileup_vcf, merge_vcf)
                }.view { "Clair3 annotation input: $it" } :
                Channel.fromList(sample_thresholds.keySet().collect())
                    .map { sample_id ->
                        def clair3_output_dir = file("${params.output_path}/routine_epi2me/${sample_id}/output_clair3")
                        def pileup_vcf = file("${params.output_path}/routine_epi2me/${sample_id}/output_clair3/pileup.vcf.gz")
                        def merge_vcf = file("${params.output_path}/routine_epi2me/${sample_id}/output_clair3/merge_output.vcf.gz")

                        if (!clair3_output_dir.exists() || !pileup_vcf.exists() || !merge_vcf.exists()) {
                            error "ERROR: Clair3 VCF files not found for ${sample_id}. Please run epi2me workflow with SNV mode first."
                        }

                        tuple(sample_id, clair3_output_dir, pileup_vcf, merge_vcf)
                    }
                    .view { "Clair3 annotation input: $it" }

            def clairsto_annot_input = (params.run_mode_epiannotation || params.run_mode_order) ?
                input_data.map { args ->
                    def sample_id = args[0]
                    def clairsto_output_dir = file("${params.output_path}/routine_epi2me/${sample_id}/clairsto_output")
                    def snv_vcf = file("${params.output_path}/routine_epi2me/${sample_id}/clairsto_output/snv.vcf.gz")
                    def indel_vcf = file("${params.output_path}/routine_epi2me/${sample_id}/clairsto_output/indel.vcf.gz")
                    tuple(sample_id, clairsto_output_dir, snv_vcf, indel_vcf)
                } :
                Channel.fromList(sample_thresholds.keySet().collect())
                    .map { sample_id ->
                        def clairsto_output_dir = file("${params.output_path}/routine_epi2me/${sample_id}/clairsto_output")
                        def snv_vcf = file("${params.output_path}/routine_epi2me/${sample_id}/clairsto_output/snv.vcf.gz")
                        def indel_vcf = file("${params.output_path}/routine_epi2me/${sample_id}/clairsto_output/indel.vcf.gz")

                        if (!clairsto_output_dir.exists() || !snv_vcf.exists() || !indel_vcf.exists()) {
                            error "ERROR: ClairS-TO VCF files not found for ${sample_id}. Please run epi2me workflow with SNV mode first."
                        }

                        tuple(sample_id, clairsto_output_dir, snv_vcf, indel_vcf)
                    }
                .view { "ClairSTo annotation input: $it" }

            // Step 2: Run annotation processes
            clair3_annotate(clair3_annot_input)
            clairs_to_annotate(clairsto_annot_input)

            // Step 3: Create properly structured channels for combination
            def clair3_results = clair3_annotate.out.clair3output
                .map { args ->
                    def sample_id = args[0]
                    def pileup_file = args[1]
                    def merge_file = args[2]
                    tuple(sample_id, pileup_file, merge_file)
                }
                .view { "Clair3 annotated mapped: $it" }

            def clairsto_results = clairs_to_annotate.out.annotateandfilter_clairstoout
                .view { "ClairSTo annotated output: $it" }

            // Step 4: Combine results and create input for merge_annotation
            combine_file = clair3_results
                .combine(clairsto_results, by: 0)
                .combine(occ_fusion_genes_list_ch)
                .map { sample_id, pileup_file, merge_file, clairsto_file, occ_genes ->
                    println "Creating merge input for sample: $sample_id"
                    tuple(
                        sample_id,
                        merge_file,
                        pileup_file,
                        clairsto_file,
                        occ_genes
                    )
                }
                .view { "Merge annotation input: $it" }

            // Run merge annotation
            merge_annotation(combine_file)
         }

        // tertp analysis
        if (params.run_mode in ['tertp', 'rmd', 'all'] || params.run_mode_order) {
            println "Running tertp Analysis..."
            igv_tools(boosts_igv_channel)
        //    igv_tools.out.tertp_out_igv.view { "tertp output: $it" }
            plot_genomic_regions(boosts_plot_genomic_regions_channel)
        }

        // CNV analysis (now handles both CNV and RMD modes)
        if (params.run_mode in ['cnv', 'rmd', 'all'] || params.run_mode_order) {
            println "Running CNV Analysis..."
            
            // Handle sample_thresholds for different run modes
            def samples_needing_ace = []
            def samples_with_provided_threshold = [:]

            if (params.run_mode_order || params.run_mode_epiannotation) {
                // In run_mode_order or run_mode_epiannotation, we always compute ACE to get thresholds
                println "Using ${params.run_mode_order ? 'run_mode_order' : 'run_mode_epiannotation'} - computing ACE for all samples to get thresholds"

                // For run_mode_order, load from bam_sample_id_file
                // For run_mode_epiannotation, load from epi2me_sample_id_file
                if (params.run_mode_order) {
                    def sample_ids = file(params.bam_sample_id_file).readLines().collect { line ->
                        // Extract only the sample ID part (first column), removing flow cell ID
                        line.trim().split(/\s+/)[0]
                    }
                    samples_needing_ace = sample_ids.toSet()
                } else {
                    // For run_mode_epiannotation, load from epi2me_sample_id_file
                    def sample_ids = file(params.epi2me_sample_id_file).readLines().collect { line ->
                        // Extract only the sample ID part (first column)
                        line.trim().split(/\t/)[0]
                    }
                    samples_needing_ace = sample_ids.toSet()
                }
                println "Samples for ACE calculation: ${samples_needing_ace}"
            } else {
                // Separate samples that need ACE from those that don't
                samples_needing_ace = sample_thresholds.findAll { k, v -> v == null }.keySet()
                samples_with_provided_threshold = sample_thresholds.findAll { k, v -> v != null }
            }

            println "Samples needing ACE calculation: ${samples_needing_ace}"
            println "Samples with provided thresholds: ${samples_with_provided_threshold}"

            // Run ACE only for samples that need calculation (have null threshold)
            if (samples_needing_ace.size() > 0) {
        // For run_mode_epiannotation or run_mode_order, use RDS from input_data channel
        // For standalone mode, scan filesystem
        if (params.run_mode_order || params.run_mode_epiannotation) {
            ace_input = input_data
                .map { args ->
                    def sample_id = args[0]
                    def rds_file = args[10]

                    if (samples_needing_ace.contains(sample_id)) {
                        println "Found matching RDS file for sample ${sample_id} from epi2me output"
                        tuple(sample_id, rds_file)
                    } else {
                        null
                    }
                }
                .filter { it != null }
        } else {
            ace_input = Channel
                .fromList(samples_needing_ace.collect { sid ->
                    def rds_file = file("${params.cnv_rds}/${sid}/${sid}_copyNumbersCalled.rds")
                    if (rds_file.exists()) {
                        println "Found matching RDS file for sample ${sid}: ${rds_file}"
                        tuple(sid, rds_file)
                    } else {
                        println "WARNING: RDS file not found for sample ${sid}: ${rds_file}"
                        null
                    }
                }.findAll { it != null })
        }

                // Run ACE analysis
            ace_tmc(ace_input)
                
            ace_thresholds = ace_tmc.out.threshold_value
                .map { args -> 
                    def sample_id = args[0]
                    def threshold = args[1]
                    println "Calculated threshold for ${sample_id}: ${threshold}"
                    tuple(sample_id, threshold.toFloat())
                }
            } else {
                ace_thresholds = Channel.empty()
            }

            // Create final threshold mapping
            def final_thresholds = [:]
            samples_with_provided_threshold.each { sample_id, threshold ->
                final_thresholds[sample_id] = threshold
                println "Using provided threshold for ${sample_id}: ${threshold}"
            }

            println "Final threshold mapping: ${final_thresholds}"

            // Create channel for annotatecnv based on run mode
            if (params.run_mode_order || params.run_mode_epiannotation) {
                // Use epi2me output paths when in run_mode_order or run_mode_epiannotation
                annotatecnv_input = input_data
                    .map { args -> 
                        def sample_id = args[0]
                        def bam = args[1]
                        def bai = args[2]
                        def ref = args[3]
                        def ref_bai = args[4]
                        def tr_bed = args[5]
                        def modkit = args[6]
                        def segs_bed = args[7]
                        def bins_bed = args[8]
                        def segs_vcf = args[9]
                def rds_file = args[10]
                        def sv = args[11]
                        println "Processing epi2me results for sample: ${sample_id} from epi2me output"

                        // Use published epi2me output directory
                        def segs_vcf_path = "${params.output_path}/routine_epi2me/${sample_id}/${sample_id}_segs.vcf"
                        def bins_bed_path = "${params.output_path}/routine_epi2me/${sample_id}/${sample_id}_bins.bed"
                        def segs_bed_path = "${params.output_path}/routine_epi2me/${sample_id}/${sample_id}_segs.bed"

                        println "Checking file existence:"
                        println "  segs_vcf_path: ${segs_vcf_path}"
                        println "  bins_bed_path: ${bins_bed_path}"
                        println "  segs_bed_path: ${segs_bed_path}"

                        tuple(
                            sample_id,
                            file(segs_vcf_path),
                            file(params.roi_protein_coding_bed),
                            file(bins_bed_path),
                            file(segs_bed_path)
                        )
                    }
            } else {
                // Use configured input folders when in standalone mode
                annotatecnv_input = Channel.fromList(sample_thresholds.keySet().collect())
                    .map { sample_id ->
                        println "Processing sample: ${sample_id} from configured folders"
                        tuple(
                            sample_id,
                            file("${params.segsfromepi2me_folder}/${sample_id}/${sample_id}_segs.vcf"),
                            file(params.roi_protein_coding_bed),
                            file("${params.segsfromepi2me_folder}/${sample_id}/${sample_id}_bins.bed"),
                            file("${params.segsfromepi2me_folder}/${sample_id}/${sample_id}_segs.bed")
                        )
                    }
            }

            // Combine with thresholds and prepare final input
            // For run_mode_order and run_mode_epiannotation, we use calculated thresholds from ACE
            def annotatecnv_with_provided = (params.run_mode_order || params.run_mode_epiannotation) ?
                annotatecnv_input.combine(ace_thresholds, by: 0).map { args ->
                    def sample_id = args[0]
                    def segs_vcf = args[1]
                    def roi_protein_coding_bed = args[2]
                    def bins_bed = args[3]
                    def segs_bed = args[4]
                    def threshold = args[5]

                    println "Preparing annotatecnv input for ${sample_id} with ACE threshold: ${threshold}"
                    tuple(
                        sample_id,
                        segs_vcf,
                        roi_protein_coding_bed,
                        bins_bed,
                        segs_bed,
                        threshold.toString(),
                        file(params.cnv_genes_tuned)
                    )
                } :
                annotatecnv_input
                    .filter { args ->
                        def sample_id = args[0]
                        def segs_vcf = args[1]
                        def roi_protein_coding_bed = args[2]
                        def bins_bed = args[3]
                        def segs_bed = args[4]

                        // Check if this sample has a provided threshold
                        def has_threshold = final_thresholds.containsKey(sample_id)
                        println "Sample ${sample_id} has provided threshold: ${has_threshold}"
                        has_threshold
                    }
                    .map { args ->
                        def sample_id = args[0]
                        def segs_vcf = args[1]
                        def roi_protein_coding_bed = args[2]
                        def bins_bed = args[3]
                        def segs_bed = args[4]

                        // Add the threshold and cnv_genes_tuned to the tuple
                        tuple(
                            sample_id,
                            segs_vcf,
                            roi_protein_coding_bed,
                            bins_bed,
                            segs_bed,
                            final_thresholds[sample_id],
                            file(params.cnv_genes_tuned)
                        )
                    }

            // For samples with calculated thresholds, combine with ace_thresholds
            def annotatecnv_with_calculated = (params.run_mode_order || params.run_mode_epiannotation) ?
                Channel.empty() :  // Skip this in run_mode_order/run_mode_epiannotation since we handle it above
                annotatecnv_input
                    .filter { args ->
                        def sample_id = args[0]
                        def segs_vcf = args[1]
                        def roi_protein_coding_bed = args[2]
                        def bins_bed = args[3]
                        def segs_bed = args[4]
                        !final_thresholds.containsKey(sample_id)
                    }
                    .combine(ace_thresholds, by: 0)
                    .map { args ->
                        def sample_id = args[0]
                        def segs_vcf = args[1]
                        def roi_protein_coding_bed = args[2]
                        def bins_bed = args[3]
                        def segs_bed = args[4]
                        def threshold = args[5]
                        println "Preparing annotatecnv input for ${sample_id} with calculated tumor content ${threshold}"
                        tuple(
                            sample_id,              // sample_id
                            segs_vcf,               // segs_vcf
                            roi_protein_coding_bed, // roi_protein_coding_bed
                            bins_bed,               // bins_bed
                            segs_bed,               // segs_bed
                            threshold.toString(),   // threshold value as string
                            file(params.cnv_genes_tuned)  // CNV genes annotation file
                        )
                    }

            // Combine both channels
            annotatecnv_input = annotatecnv_with_provided.mix(annotatecnv_with_calculated)
                .view { "Annotatecnv input: $it" }

            // Run annotatecnv
            annotatecnv(annotatecnv_input)
            annotatecnv_results = annotatecnv.out

            // Run plot_genomic_regions for coverage analysis
            //plot_genomic_regions(boosts_plot_genomic_regions_channel)
        }

        // Statistics mode - cramino is now run by epi2me workflow
        // Analysis workflow reads the pre-generated cramino outputs
        // RMD report generation
        if (params.run_mode in ['rmd', 'all'] || params.run_mode_order || params.run_mode_epiannotation) {
            println "Running RMD Report Generation..."

            // annotatecnv_results is set by CNV section which always runs before RMD
            // No fallback needed - if CNV section ran, annotatecnv_results is defined

            // Reuse MGMT outputs from earlier analysis if available
            // If MGMT analysis was not run, we need to run it here
            if (!(params.run_mode in ['mgmt', 'rmd', 'all'])) {
                println "MGMT analysis not run earlier, running now for RMD report..."
                
                // Create channel for MGMT analysis
        mgmt_ch = (params.run_mode_order || params.run_mode_epiannotation) ? 
                input_data.map { args -> 
                    def sample_id = args[0]
                    def bam = args[1]
                    def bai = args[2]
                    def ref = args[3]
                    def ref_bai = args[4]
                    def tr_bed = args[5]
                    def modkit = args[6]
                    def segs_bed = args[7]
                    def bins_bed = args[8]
                    def segs_vcf = args[9]
                def rds_file = args[10]
                    def sv = args[11]
                    
                    tuple(
                        sample_id,
                        modkit,
                        epicsites,
                        mgmt_cpg
                    )
                } :
                Channel.fromList(sample_thresholds.keySet().collect())
                    .combine(epicsites_ch).combine(mgmt_cpg_island_hg38_ch)
                    .map { sample_id, epicsites, mgmt_cpg ->
                        tuple(
                            sample_id,
                            file("${params.bedmethyl_folder}/${sample_id}/${sample_id}.wf_mods.bedmethyl.gz"),
                            epicsites,
                            mgmt_cpg
                        )
                    }

            // Run MGMT related processes
            extract_epic(mgmt_ch)

            // Create channels for downstream processes
            MGMT_output = extract_epic.out.MGMTheaderout

            MGMT_sturgeon = extract_epic.out.sturgeonbedinput
                 .combine(sturgeon_model_ch)
                 .map { sample_id, sturgeoninput, sturgeon_model ->
                     tuple(sample_id, sturgeoninput, sturgeon_model)
                 }

            mgmt_nanodx = extract_epic.out.epicselectnanodxinput
                .combine(hg19_450model_ch)
                .map { sample_id, epicselectnanodxinput, hg19_model ->
                    tuple(sample_id, epicselectnanodxinput, hg19_model)
                }

            // Run the processes
            if (sturgeon_available) { sturgeon(MGMT_sturgeon) }
            mgmt_promoter(MGMT_output)
            nanodx(mgmt_nanodx)

            nanodx_out = nanodx.out.nanodx450out
                .combine(nanodx_450model_ch).combine(snakefile_nanodx_ch).combine(nn_model_ch)
                .map { sample_id, nanodx450out, nanodx_model, snakefile, nn_model ->
                    tuple(
                        sample_id,
                        nanodx450out,
                        nanodx_model,
                        snakefile,
                        nn_model
                    )
                }

            crossNN(nanodx_out)
            rmd_nanodx_out = crossNN.out.rmdnanodx

            // Pan-cancer classifier (always runs alongside Capper)
            nanodx_pancan_out = nanodx.out.nanodx450out
                .combine(nanodx_pancan_model_ch).combine(snakefile_nanodx_ch).combine(nn_model_ch)
                .map { sample_id, nanodx450out, nanodx_model, snakefile, nn_model ->
                    tuple(sample_id, nanodx450out, nanodx_model, snakefile, nn_model)
                }
            crossNN_pancan(nanodx_pancan_out)

            // Generate t-SNE plots for both classifiers
            tsne_plot(
                extract_epic.out.epicselectnanodxinput,
                nanodxcolormap_ch,
                nanodxh5_ch
            )
            tsne_plot_pancan(
                extract_epic.out.epicselectnanodxinput,
                nanodxcolormap_pancan_ch,
                nanodxh5_pancan_ch
            )
            } else {
                println "Reusing MGMT outputs from earlier analysis"
                // The outputs are already available from the MGMT section
            }

            // Svanna analysis - reuse outputs if already run
            if (!(params.run_mode in ['svannasv', 'rmd', 'all'])) {
                println "Svanna analysis not run earlier, running now for RMD report..."
        svannasv(boosts_svanna_channel)
        rmd_svanna_html = svannasv.out.rmdsvannahtml
            svannasv_out = svannasv.out.occsvannaannotationannotationvcf
                .combine(vcf2circos_json_ch)
                .map { sample_id, occsvannaannotationannotationvcf, vcf2circos ->
                    [sample_id, occsvannaannotationannotationvcf, vcf2circos]
                }
        circosplot(svannasv_out)

        // Also run fusion events analysis for RMD mode
        svannaoutfusion_events = svannasv.out.occsvannaannotationannotationvcf
                .combine(genecode_bed_ch).combine(occ_fusion_genes_list_ch)
                .map { sample_id, occsvannaannotationannotationvcf, genecode, fusion_genes ->
                [sample_id, occsvannaannotationannotationvcf, genecode, fusion_genes]
            }
        svannasv_fusion_events(svannaoutfusion_events)
        
        // Assign fusion events channel
        fusion_events_channel = svannasv_fusion_events.out.filterfusioneventout
            } else {
                println "Reusing Svanna outputs from earlier analysis"
                // Ensure fusion_events_channel is defined when reusing outputs
                fusion_events_channel = svannasv_fusion_events.out.filterfusioneventout
            }

            // ROI analysis - reuse outputs if already run
            if (!(params.run_mode in ['roi', 'rmd', 'all'])) {
                println "ROI analysis not run earlier, running annotation for RMD report..."
                // NOTE: Variant calling should have been done by epi2me workflow
                // Read pre-existing VCF files from epi2me output directory

                // Step 1: Read pre-existing VCF files
                def clair3_annot_input = Channel.fromList(sample_thresholds.keySet().collect())
                    .map { sample_id ->
                        def clair3_output_dir = file("${params.output_path}/routine_epi2me/${sample_id}/output_clair3")
                        def pileup_vcf = file("${params.output_path}/routine_epi2me/${sample_id}/output_clair3/pileup.vcf.gz")
                        def merge_vcf = file("${params.output_path}/routine_epi2me/${sample_id}/output_clair3/merge_output.vcf.gz")

                        // Skip existence check in epiannotation and run_mode_order (files created by epi2me workflow)
                        if (!params.run_mode_epiannotation && !params.run_mode_order && (!clair3_output_dir.exists() || !pileup_vcf.exists() || !merge_vcf.exists())) {
                            error "ERROR: Clair3 VCF files not found for ${sample_id}. Please run epi2me workflow with SNV mode first."
                        }

                        tuple(sample_id, clair3_output_dir, pileup_vcf, merge_vcf)
                    }

                def clairsto_annot_input = Channel.fromList(sample_thresholds.keySet().collect())
                    .map { sample_id ->
                        def clairsto_output_dir = file("${params.output_path}/routine_epi2me/${sample_id}/clairsto_output")
                        def snv_vcf = file("${params.output_path}/routine_epi2me/${sample_id}/clairsto_output/snv.vcf.gz")
                        def indel_vcf = file("${params.output_path}/routine_epi2me/${sample_id}/clairsto_output/indel.vcf.gz")

                        // Skip existence check in epiannotation and run_mode_order (files created by epi2me workflow)
                        if (!params.run_mode_epiannotation && !params.run_mode_order && (!clairsto_output_dir.exists() || !snv_vcf.exists() || !indel_vcf.exists())) {
                            error "ERROR: ClairS-TO VCF files not found for ${sample_id}. Please run epi2me workflow with SNV mode first."
                        }

                        tuple(sample_id, clairsto_output_dir, snv_vcf, indel_vcf)
                    }

                // Step 2: Run annotation
                clair3_annotate(clair3_annot_input)
                clairs_to_annotate(clairsto_annot_input)

                clair3_out = clair3_annotate.out.clair3output
                clairs_to_out = clairs_to_annotate.out.annotateandfilter_clairstoout
            combine_file = clair3_out.combine(clairs_to_out, by: 0)
                .combine(occ_fusion_genes_list_ch)
                .map { sample_id, pileup_file, merge_file, clairsto_file, occ_genes ->
                    tuple(sample_id, merge_file, pileup_file, clairsto_file, occ_genes)
    }
        merge_annotation(combine_file)
            } else {
                println "Reusing ROI outputs from earlier analysis"
            }

            // Other tools - reuse outputs if already run or run for RMD
            // For 'rmd' mode, igv_tools and plot_genomic_regions were already run in tertp section
            // For 'tertp' mode, they were also already run
            // For other modes (rmd, cnv, stat, etc), we need to run them now
            if (!(params.run_mode in ['tertp', 'rmd', 'all'])) {
                println "Running tertp analysis for RMD report..."
                igv_tools(boosts_igv_channel)
                plot_genomic_regions(boosts_plot_genomic_regions_channel)
            } else {
                println "Reusing tertp outputs from earlier analysis"
            }

            // Read cramino outputs from epi2me workflow
            println "Reading Cramino statistics from epi2me output..."
            def cramino_output_ch = Channel.fromList(sample_thresholds.keySet().collect())
                .map { sample_id ->
                    def cramino_file = file("${params.output_path}/routine_epi2me/${sample_id}/cramino/${sample_id}_cramino_statistics.txt")
                    // Skip existence check in epiannotation and run_mode_order (files created by epi2me workflow)
                    if (!params.run_mode_epiannotation && !params.run_mode_order && !cramino_file.exists()) {
                        error "ERROR: Cramino statistics file not found for ${sample_id}. Please run epi2me workflow with 'stat' or 'rmd' mode first."
                    }
                    tuple(sample_id, cramino_file)
                }

            // Combine all results for markdown report
            // Use the stored annotatecnv results from earlier CNV analysis
            if (!annotatecnv_results) {
                error "CNV analysis results not found. Make sure CNV analysis runs before RMD generation."
            }
            
            // Check if all required channels are available
            if (!merge_annotation.out.occmergeout) {
                error "Merge annotation results not found. Make sure OCC analysis runs before RMD generation."
            }
            
            if (!crossNN.out.rmdnanodx) {
                error "NN classifier results not found. Make sure MGMT analysis runs before RMD generation."
            }
            
            if (!mgmt_promoter.out.mgmtresultsout) {
                error "MGMT promoter results not found. Make sure MGMT analysis runs before RMD generation."
            }
            
            if (!svannasv.out.rmdsvannahtml) {
                error "Svanna HTML results not found. Make sure Svanna analysis runs before RMD generation."
            }
            
            if (!fusion_events_channel) {
                error "Fusion events results not found. Make sure Svanna analysis runs before RMD generation."
            }
            
            if (!igv_tools.out.tertp_out_igv) {
                error "tertp HTML results not found. Make sure tertp analysis runs before RMD generation."
            }

            if (!cramino_output_ch) {
                error "Cramino statistics not found. Make sure Cramino analysis runs before RMD generation."
            }

            if (!plot_genomic_regions.out.plot_genomic_regions_out) {
                error "Genomic regions plot results not found. Make sure tertp analysis runs before RMD generation."
            }

            mergecnv_out = annotatecnv_results.rmdcnvtumornumber
            .combine(merge_annotation.out.occmergeout, by:0)
            .combine(crossNN.out.rmdnanodx, by: 0)
            .combine(mgmt_promoter.out.mgmtresultsout, by:0)
            .combine(svannasv.out.rmdsvannahtml, by:0)
            .combine(fusion_events_channel, by:0)
            .combine(igv_tools.out.tertp_out_igv, by:0)
            .combine(cramino_output_ch, by:0)
            .combine(plot_genomic_regions.out.plot_genomic_regions_out, by:0)
            .combine(tsne_plot.out.tsne_out, by:0)
            .combine(crossNN_pancan.out.rmdnanodx_pancan, by:0)
            .combine(tsne_plot_pancan.out.tsne_pancan_out, by:0)
            // Create final map for markdown report
        mergecnv_out_map = mergecnv_out.map { args ->
                def sample_id = args[0]
                def cnv_plot = args[1]
                def tumor_copy_number = args[2]
                def annotatedcnv_filter_header = args[3]
                def cnv_chr9 = args[4]
                def cnv_chr7 = args[5]
                def merge_annotation_filter_snvs_allcall = args[6]
                def nanodx_classifier = args[7]
                def mgmt_results = args[8]
                def svannahtml = args[9]
                def fusion_events = args[10]
                def tertphtml = args[11]
                def craminoreport = args[12]
                def egfr_coverage = args[13]
                def idh1_coverage = args[14]
                def tertp_coverage = args[15]
                def idh2_coverage = args[16]
                def tsne_plot_file = args[17]
                def nanodx_classifier_pancan = args[18]
                def tsne_plot_pancan_file = args[19]

                // Use correct sample ID file based on run mode
                def sample_id_file = params.run_mode_order ? "placeholder" : params.analyse_sample_id_file

                [
                    sample_id,
                    craminoreport,
                    sample_id_file,
                    params.nanodx_dictinaire,
                    params.mardown_logo,
                    cnv_plot,
                    tumor_copy_number,
                    annotatedcnv_filter_header,
                    cnv_chr9,
                    cnv_chr7,
                    mgmt_results,
                    merge_annotation_filter_snvs_allcall,
                    fusion_events,
                    svannahtml,
                    tertphtml,
                    egfr_coverage,
                    idh1_coverage,
                    idh2_coverage,
                    tertp_coverage,
                    tsne_plot_file,
                    nanodx_classifier,
                    file("${params.ref_dir}/snv_target_genes.txt"),
                    file(params.roi_protein_coding_bed),
                    file("${params.diana_dir}/bin/nextflow_markdown_pipeline_update_final.Rmd"),
                    file("${params.diana_dir}/docs/warning.png"),
                    nanodx_classifier_pancan,
                    file(params.nanodx_pancan_dictinaire),
                    tsne_plot_pancan_file
                ]
            }.view()

            // Generate markdown report
            markdown_report(mergecnv_out_map)
        }

    emit:
        markdown_out = (params.run_mode_annotation == 'rmd' || params.run_mode_order) ? 
            markdown_report.out : 
            Channel.empty()
}

// Helper function to extract sample_id from BAM filename
//def extractSampleId(file) {
//    def filename = file.getName()
//    return filename.split('\\.')[0]  // Get everything before the first dot
//}
//}
