#!/usr/bin/env nextflow

nextflow.enable.dsl=2
def start_time = new Date()
//params.mgmt_promoter_r_script = "mnt/scripts/MGMT_Prospective2.R"
// Define the base path as a parameter
// params.path = "/home/chbope/extension"
// params.input_path = "${params.path}"
// params.output_path = "${params.path}/results"
// params.annotate_dir = "/home/chbope/extension/data/annotations"
// params.config_dir = "/home/chbope/Documents/nanopore/packages/knotannotsv/knotAnnotSV"
//   params.ref_dir ="/home/chbope/extension/data/reference"
// params.model_path="/home/chbope/extension/Data_for_Bope"
// params.clair3_dir="/home/chbope/Documents/nanopore/Data_for_Bope/results/sample_id1/callvariantclair3/clair3_output"
// params.humandb_dir="/home/chbope/extension/data/annovar/humandb"
// params.clairSTo_dir="/home/chbope/Documents/nanopore/Data_for_Bope/results/sample_id1/callvariantclairsto/clairsto_output"
// params.svanna_dir="/home/chbope/extension/data/svanna-cli-1.0.4/"
// params.bin_dir="/home/chbope/Documents/nanopore/nextflow/bin/"
// params.epi2me_dir="/home/chbope/Documents/nanopore/epi2me/wf-human-variation-master/"
// params.out_dir_epi2me ="/home/chbope/extension/out_dir_epi2me"

// // input reference files for all samples

// params.reference_genome = file("${params.ref_dir}/GCF_000001405.39_GRCh38.p13_genomic_chr_only_plus_mt.fa")
// params.reference_genome_bai = file("${params.ref_dir}/GCF_000001405.39_GRCh38.p13_genomic_chr_only_plus_mt.fa.fai")

// params.occ_fusions=file("${params.ref_dir}/OCC.fusions.bed")
// params.epicsites=file("${params.ref_dir}/EPIC_sites_NEW.bed")
// params.mgmt_cpg_island_hg38=file("${params.ref_dir}/MGMT_CpG_Island.hg38.bed")

// params.occ_snv_screening=file("${params.ref_dir}/OCC.SNV.screening.bed")
// params.tertp_variants=file("${params.ref_dir}/TERTp_variants.bed")
// params.ncbirefseq=file("${params.ref_dir}/ncbiRefSeq.txt.gz")
// params.refgene=file("${params.humandb_dir}/hg38_refGene.txt")
// params.hg38_refgenemrna=file("${params.humandb_dir}/hg38_refGeneMrna.fa")
// params.clinvar=file("${params.humandb_dir}/hg38_clinvar_20240611.txt")
// params.clinvarindex=file("${params.humandb_dir}/hg38_clinvar_20240611.txt.idx")
// params.hg38_cosmic100=file("${params.humandb_dir}/hg38_cosmic100coding2024.txt")
// params.hg38_cosmic100index=file("${params.humandb_dir}/hg38_cosmic100coding2024.txt.idx")
// params.knotannotsv_conf=file("${params.config_dir}/config_AnnotSV.yaml")
// params.occ_fusion_genes_list=file("${params.ref_dir}/occ_fusions_genes.txt")
// //params.bins_bed=file("/home/chbope/Documents/nanopore/Data_for_Bope/segs/T21-058_bins.bed")
// params.convertvariant_config = file("${params.ref_dir}/annotsv2.json")
// params.occ_genes = file("${params.ref_dir}/OCCgenes.rds")
// params.vcf2circos_json = file("${params.ref_dir}/options.json")
// params.sturgeon_model = file("${params.ref_dir}/general.zip")
// params.hg19_450model = file("${params.ref_dir}/nanoDx/static/hg19_450model.bed")
// params.nanodx_450model = file("${params.ref_dir}/nanoDx/static/Capper_et_al_NN.pkl")
//   params.nanodx_dictinaire = file("${params.ref_dir}/nanoDx/static/Capper_et_al_dictionary.txt")
// params.nn_model = file("${params.ref_dir}/nanoDx/workflow/envs/NN_model.yaml")
// params.snakefile_nanodx = file("${params.ref_dir}/nanoDx/workflow/Snakefile")
// params.sub_snakefile = file("${params.ref_dir}/nanoDx/workflow/methylation.snakefile")
// params.tr_bed_file = file("${params.ref_dir}/human_GRCh38_trf.bed")
// params.epi2me_config_file = file("/home/chbope/Documents/nanopore/epi2me/wf-human-variation-master/nextflow.config")
// params.epi2me_base_config_file = file("/home/chbope/Documents/nanopore/epi2me/wf-human-variation-master/base.config")
// params.sample_id_file = file("${params.path}/Data_for_Bope/sample_ids1.txt")
// params.mardown_logo = file("/home/chbope/Downloads/log_update.pdf")
// //params.ma = file("/home/chbope/extension/Data_for_Bope/sample_ids.txt")
// //input individual samples files

// //input file folder
// //merge_bam_folder='/home/chbope/extension/Data_for_Bope/merge_bams/'
// //vcf_folder='/home/chbope/extension/Data_for_Bope/vcfs/'
// params.bams_folder='/home/chbope/extension/Data_for_Bope/bams/'
// params.bedmethyl_folder='/home/chbope/extension/Data_for_Bope/bedmethyl/'
// params.segsfromepi2me_folder='/home/chbope/extension/Data_for_Bope/segs/'
// params.sv_folder='/home/chbope/extension/Data_for_Bope/sv/'
// params.nanodx_workflow_dir='/home/chbope/Documents/nanopore/packages/nanoDx/workflow'



// Define the available run modes
//def mode_descriptions = [
//    'methylation': 'Run methylation processes.',
//    'annotsv': 'Run structure variant annotation process (AnnotSV and Svanna).',
//    'cnv' : 'Run a copy number variation process.', 
//    'occ' : 'Run OCC processes (Clair3 and Clais_to).', 
//    'terp' : 'Run igv_tool process.',
//    'rmd'   : 'Run the rmd report',
//    'stat' : 'Generate Sequencing stats',
//    'all': 'Run all the pipeline processes.'
//    ]

// Define the run_mode parameter (default is "both")
//params.run_mode = params.run_mode ?: 'all'


// Show available modes and descriptions before starting the pipeline
//println "Available run modes with descriptions:"
//mode_descriptions.each { mode, description ->
//    println "- ${mode}: ${description}"
//}
//println "Selected run mode: ${params.run_mode}"

//if (!mode_descriptions.keySet().contains(params.run_mode)) {
//    println "ERROR: Invalid run mode '${params.run_mode}'. Please select from ${mode_descriptions.keySet().toList().join(', ')}."
//    exit 1  // Exit with error status
//}


//##############################3
//## Extract overlapping EPIC sites for methylation based classification 
//#############################

//# wf_mods_bedmethyl generated by modkit in Epi2me-labs on the P24
//# EPIC_sites_NEW.bed contains chromosomal coordinates and names of EPIC 850 probes
//# Uses Bedtools intersect, awk and sed


//epi2me pipeline to run methyl, cnv and sv



process extract_epic {
    cpus 4
    memory '2 GB'
    label 'epic'
    tag "${sample_id}"

    publishDir "${params.output_path}/methylation/", mode: "copy", overwrite: true

    input:
    tuple val(sample_id), file(bedmethyl), file(epicsites), file(mgmt_cpg_island_hg38)

    output:
   // path("${sample_id}_EpicSelect.bed")
    tuple val(sample_id), path("${sample_id}_EpicSelect_header.bed"), emit: epicselectnanodxinput
    path("${sample_id}_MGMT.bed")
    tuple val(sample_id), path("${sample_id}_MGMT_header.bed"), emit: MGMTheaderout
    tuple val(sample_id), path("${sample_id}_wf_mods.bedmethyl_intersect.bed"), emit: sturgeonbedinput

    script:
    """
    # Step 1: Run intersectBed and save the output to EpicSelect.bed
    which intersectBed 
    intersectBed -a $bedmethyl -b $epicsites -wb | \
    awk -v OFS="\\t" '\$1=\$1' | awk -F'\\t' 'BEGIN{ OFS="\\t" }{print \$1,\$2,\$3,\$4,\$5,\$11,\$23}' > ${sample_id}_EpicSelect.bed

    intersectBed -a $bedmethyl  -b $epicsites   -wb | awk -v OFS="\\t" '\$1=\$1' | awk -F'\\t' 'BEGIN{ OFS="\\t" } {print \$1, \$2, \$3, \$4, \$5, \$6, \$7, \$8, \$9, \$10, \$11, \$12, \$13, \$14, \$15, \$16, \$17, \$18}' > ${sample_id}_wf_mods.bedmethyl_intersect.bed

    # Step 2: Add header to EpicSelect.bed
    awk 'BEGIN {print "Chromosome\\tStart\\tEnd\\tmodBase\\tCoverage\\tMethylation_frequency\\tIllumina_ID"} 1' ${sample_id}_EpicSelect.bed > ${sample_id}_EpicSelect_header.bed

    # Step 3: Run intersectBed with MGMTp_hg38
    intersectBed -a $bedmethyl  -b $mgmt_cpg_island_hg38 | \
    awk -v OFS="\\t" '\$1=\$1' | awk -F'\\t' 'BEGIN{ OFS="\\t" }{print \$1,\$2,\$3,\$4,\$5,\$11,\$12,\$13,\$14,\$15,\$16}'  > ${sample_id}_MGMT.bed

    # Step 4: Add header to MGMT.bed
    awk 'BEGIN {print "Chrom\\tStart\\tEnd\\tmodBase\\tDepth\\tMethylation\\tNmod\\tNcanon\\tNother\\tNdelete\\tNfail"} 1' ${sample_id}_MGMT.bed > ${sample_id}_MGMT_header.bed
    """
}

//Sturgeon classifier
process sturgeon {
    cpus 4
    memory '2 GB'
    label 'epic'
    publishDir "${params.output_path}/classifier/sturgeon", mode: "copy", overwrite: true

    input:
    tuple val(sample_id), path(sturgeon_bed), path(sturgeon_model)

    output:
    tuple path("${sample_id}_bedmethyl_sturgeon.bed"), path("${sample_id}_bedmethyl_sturgeon")


    """
    /sturgeon/venv/bin/sturgeon inputtobed -i $sturgeon_bed  -o ${sample_id}_bedmethyl_sturgeon.bed  -s modkit_pileup  --reference-genome hg38
   
    /sturgeon/venv/bin/sturgeon predict -i ${sample_id}_bedmethyl_sturgeon.bed   -o  ${sample_id}_bedmethyl_sturgeon --model-files $sturgeon_model  --plot-results

    """
}

//NanoDX classifier

process nanodx {
    cpus 4
    memory '16 GB'
    label 'epic'
    publishDir "${params.output_path}/classifier/nanodx", mode: "copy", overwrite: true

    input:
    tuple val(sample_id), path(nanodx_bed), path(hg19_450model)

    output:
    tuple val(sample_id), path("${sample_id}_nanodx_bedmethyl.bed"), emit: nanodx450out

    script:
    """
    nanodx450intersectdataframe.py $hg19_450model $nanodx_bed ${sample_id}_output_cpg.bed ${sample_id}_nanodx_bedmethyl.bed ${sample_id}_nanodx_bedmethylfilter.bed
    """
}

//nanodx classifier for 45
process run_nn_classifier {
    label 'snakemake'
    publishDir "${params.output_path}/classifier/nanodx", mode: 'copy'
    
    cpus 4
    memory '12 GB'

    input:
    tuple val(sample_id), path(bed_file), path(model_file), path(nn_model), path(snakefile_nanodx)

    output:
    tuple val(sample_id), path("${sample_id}_nanodx_classifier.txt")
    tuple val(sample_id), path("${sample_id}_nanodx_classifier.tsv"), emit:rmdnanodx

    containerOptions = "-B /home/chbope/extension/trash/tmp:/home/chbope/extension/trash/tmp/"

    script:
    """
    # Set TMPDIR to the provided path
    export TMPDIR="/home/chbope/extension/trash/tmp/"
    mkdir -p \$TMPDIR

    # Verify that the directory is writable
    echo "Testing write permissions in \$TMPDIR"
    touch \$TMPDIR/test_file
    if [ -f \$TMPDIR/test_file ]; then
        echo "Write permissions are working."
    else
        echo "Write permissions are NOT working."
    fi

    # Source the Conda environment
    source /opt/conda/etc/profile.d/conda.sh
    conda activate nanodx_environ
    conda info --envs
    conda env list 
    echo \$CONDA_PREFIX
    which conda
    
    # Verify environment setup
    which python
    python -c "import pandas; print(pandas.__version__)"
    
    # Create a temporary Snakefile with the dynamic inputs
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
    threads: 4
    resources: 
        mem_mb = 16384
    script: "${params.nanodx_workflow_dir}/scripts/classify_NN_bedMethyl.py"
EOF

    # Run snakemake within the container
    snakemake --use-conda \
        --conda-prefix \$TMPDIR/.snakemake/conda \
        --cores ${task.cpus} \
        --verbose \
        NN_classifier
    """
}



// Extract MGMT promoter Ccd pG-island #extracting MGMT prmoter sites
process mgmt_promoter {
    cpus 4
    memory '2 GB'
    label 'epic'
    publishDir "${params.output_path}/methylation/", mode: "copy", overwrite: true

    input:
    tuple val(sample_id), path(mgmt_bed)

    output:
    tuple val(sample_id), path("${sample_id}_MGMT_results.csv"), emit:mgmtresultsout
    
    script:
    MGMT_results = "${sample_id}_MGMT_results.csv"

  //# calculate methylation in MGMT promoter

  """

   MGMT_Prospective2.R ${mgmt_bed} ${sample_id}_MGMT_results.csv

  """
}

//####################################
//## Annotate structural variants
//####################

//# _wf_sv.vcf.gz greated by Sniffles2 via Epi2melabs on the P24
//# Fusion gene candidate regions in OCC.fusions.bed
//# Variants annotated with AnnotSV

process svannasv {
    cpus 4
    memory '2 GB'
    label 'svannasv'
    publishDir "${params.output_path}/structure_variant/svannasv/", mode: "copy", overwrite: true
   
    input:
    tuple val(sample_id), 
          path(sv_file), 
          path(sv_file_tbi), 
          path(occ_fusions)

    output:
    tuple val(sample_id), path("${sample_id}_occ_svanna_annotation.html"), emit:rmdsvannahtml
    tuple val(sample_id), path("${sample_id}_occ_svanna_annotation.vcf.gz"), emit: occsvannaannotationannotationvcf
    
    script:
    """
    # Handle gzipped files
    if [[ "${sv_file}" == *.gz && "${sv_file_tbi}" != "NO_INDEX_NEEDED" ]]; then
        gunzip -c ${sv_file} > tmp.vcf
        SV_INPUT="tmp.vcf"
    else
        SV_INPUT="${sv_file}"
    fi

    # Create intersection file
    intersectBed -a \$SV_INPUT -b ${occ_fusions} -header > ${sample_id}_OCC_SVs.vcf

    # Check if intersection file is empty (excluding header)
    if [ \$(grep -v '^#' ${sample_id}_OCC_SVs.vcf | wc -l) -eq 0 ]; then
        # If empty, use original wf_sv file
        INPUT_FILE=\$SV_INPUT
    else
        # If not empty, use intersection file
        INPUT_FILE=${sample_id}_OCC_SVs.vcf
    fi
 
    # Run svanna prioritization (redirect stderr to avoid treating warnings as errors)
    /usr/lib/jvm/java-17-openjdk-amd64/bin/java -jar ${params.bin_dir}/svanna-cli-1.0.3.jar prioritize \\
        -d ${params.svanna_dir}/svanna-data \\
        --vcf \$INPUT_FILE \\
        --phenotype-term HP:0100836 \\
        --output-format html,vcf \\
        --prefix ${sample_id}_occ_svanna_annotation

    # Copy report to output directory (continue even if file doesn't exist)
    cp "${sample_id}_occ_svanna_annotation.html" "${params.output_path}/report/${sample_id}_occ_svanna_annotation.html"

    """
}



process annotesv {
    cpus 4
    memory '2 GB'
    label 'annotesv'
    publishDir "${params.output_path}/structure_variant/annotsv/", mode: "copy", overwrite: true
   
    input:
    tuple val(sample_id), 
          path(sv_file), 
          path(sv_file_tbi), 
          path(occ_fusions), 
          path(occ_fusion_genes_list)

    output:
    //tuple val(sample_id), path("${sample_id}_annotated_variants.tsv"), emit: annotatedvariantsout
    path("${sample_id}_OCC_SVs.vcf")
    tuple val(sample_id), path("annotated_variants.tsv"), emit: annotatedvariantsout
    tuple val(sample_id), path("${sample_id}_annotSV_fusion_extraction.csv"), emit: annotsvfusion
    path("${sample_id}_annotated_variants.tsv")
    
    script:
    """
    # Create intersection file
    if [[ "${sv_file}" == *.gz && "${sv_file_tbi}" != "NO_INDEX_NEEDED" ]]; then
        gunzip -c ${sv_file} > tmp.vcf
        SV_INPUT="tmp.vcf"
    else
        SV_INPUT="${sv_file}"
    fi

    intersectBed -a \$SV_INPUT -b $occ_fusions -header > ${sample_id}_OCC_SVs.vcf

    # Check if intersection file is empty (excluding header)
    if [ \$(grep -v '^#' ${sample_id}_OCC_SVs.vcf | wc -l) -eq 0 ]; then
        # If empty, use original wf_sv file
        INPUT_FILE=\$SV_INPUT
    else
  
        INPUT_FILE=${sample_id}_OCC_SVs.vcf
    fi

    # Run AnnotSV with selected input file
    AnnotSV -SVinputFile \$INPUT_FILE \
        -annotationsDir ${params.annotate_dir} \
        -hpo "HPO:0100836" \
        -vcf 1 \
        -genomeBuild GRCh38 \
        -outputFile ./annotated_variants.tsv
            
    annotsv_fusion_filter.py ./annotated_variants.tsv $occ_fusion_genes_list ${sample_id}_annotSV_fusion_extraction.csv

   
    cp ./annotated_variants.tsv ${sample_id}_annotated_variants.tsv

    # Verify the copy operation
    if [ ! -f "${sample_id}_annotated_variants.tsv" ]; then
        echo "Error: Copy failed. Debugging information:"
        ls -la
        pwd
        exit 1
    fi
    """
}


//visualize AnnotSV result using knotAnnotSV
process knotannotsv {
    cpus 4
    memory '2 GB'
    label 'knotannotsv'
    publishDir "${params.output_path}/structure_variant/annotsv/", mode: "copy", overwrite: true
   
    input:
    tuple val(sample_id), path(annotsv_output), path(knotannov_conf)

    output:
    tuple val(sample_id), path("${sample_id}_annotated_variants.html"), emit: rmdannotsvhtml
    path("${sample_id}_annotated_variants.xlsm")

    script:
    """
    # Run knotAnnotSV for HTML output
    knotAnnotSV.pl --annotSVfile ${annotsv_output} --configFile ${knotannov_conf} --outDir 

    # Run knotAnnotSV for Excel output
    knotAnnotSV2XL.pl --annotSVfile ${annotsv_output} --configFile ${knotannov_conf} --outDir 

    # Create empty files if they don't exist
    #touch annotated_variants.html
    
    
    #touch annotated_variants.xlsm

    # Copy and rename output files (with error handling)
    cp annotated_variants.html ${sample_id}_annotated_variants.html 
    #|| touch ${sample_id}_annotated_variants.html
    cp annotated_variants.xlsm ${sample_id}_annotated_variants.xlsm 
    #|| touch ${sample_id}_annotated_variants.xlsm
    """
}

//visualized Structure variant from AnnotSv using vfc2circos packages


process circosplot {
    cpus 4
    memory '2 GB'
    label 'circos'
    publishDir "${params.output_path}/structure_variant/svannasv/", mode: "copy", overwrite: true
   
    input:
    tuple val(sample_id), path(annotsv_output), path(vcf2circos_json)

    output:
    tuple val(sample_id), path("${sample_id}_vcf2circo.html"), optional: true, emit: circosout

    script:
    """
    # Check if file is empty (excluding header)
    if [ \$(grep -v '^#' ${annotsv_output} | wc -l) -eq 0 ]; then
        echo "Warning: ${annotsv_output} is empty. Skipping vcf2circos plot generation."
        touch ${sample_id}_vcf2circo.html
        exit 0
    else
        vcf2circos -i $annotsv_output -o ${sample_id}_vcf2circo.html -p $vcf2circos_json -a hg38
    fi
    """
}

//# Create an annotated CNV plot from the outputs of Epi2me qDNAseq
// # Uses _segs.vcf which is created by qDNAseq via Epi2melabs on P24
// # Custom Rscript to add annotation labels to CNV plot where gene amplifications or deletions are detected

// Process definition


process annotatecnv {
    cpus 4
    memory '2 GB'
    label 'annotatecnv'
    publishDir "${params.output_path}/copy_number_variation/", mode: "copy", overwrite: true

    input:
    tuple val(sample_id), path(vcf_file), path(occ_fusions), path(calls_bed), path(seg_bed), val(threshold)

    output:
    tuple val(sample_id), path("${sample_id}_calls_fixed.vcf"), emit: callsfixedout
   path("${sample_id}_annotatedcnv.csv")
   path("${sample_id}_annotatedcnv_filter.csv")
   path("${sample_id}_CNV_plot.pdf")
   path("${sample_id}_annotatedcnv_filter_header.csv"), emit:rmdannotatedcnvfilter
   path("${sample_id}_CNV_plot.html")
   tuple path("${sample_id}_tumor_copy.txt"), path("${sample_id}_bins_filter.bed")
   tuple path("${sample_id}_CNV_plot.pdf"), path("${sample_id}_annotatedcnv_filter.csv"), emit:cnvpdfandcsvout
   tuple val(sample_id), path("${sample_id}_cnv_plot_full.pdf"), path("${sample_id}_tumor_copy_number.txt"), path("${sample_id}_annotatedcnv_filter_header.csv"), path("${sample_id}_cnv_chr9.pdf"), path("${sample_id}_cnv_chr7.pdf"), emit: rmdcnvtumornumber
//   source /opt/conda/etc/profile.d/conda.sh
   script:
   """
   source /opt/conda/etc/profile.d/conda.sh
   conda init
   conda activate annotatecnv_env
   awk 'OFS="\\t" {if (NR > 13) \$1="chr"\$1; print}' $vcf_file  > ${sample_id}_calls_fixed.vcf

   intersectBed -a ${sample_id}_calls_fixed.vcf  -b $occ_fusions -wa -wb | \
   cut -f1,2,5,8,20 | awk '/protein_coding/'| awk -v OFS=";" '\$1=\$1'| awk 'BEGIN { FS=";"; OFS="\t"} {\$1=\$1; print}' | \
   cut -f1,2,3,5,6,8,9,13 > ${sample_id}_annotatedcnv.csv
   
   cnv_html.R $calls_bed ${sample_id}_annotatedcnv.csv ${sample_id}_CNV_plot.pdf ${sample_id}_CNV_plot.html $sample_id
   CNV_function_new_update.R $calls_bed ${sample_id}_annotatedcnv.csv $seg_bed ${sample_id}_cnv_plot_full.pdf ${sample_id}_cnv_chr9.pdf ${sample_id}_cnv_chr7.pdf $sample_id 

   awk 'BEGIN { OFS="," } {
    gsub(/<[^>]+>/, substr(\$3, 2, length(\$3) - 2), \$3);
    print \$0
   }' ${sample_id}_annotatedcnv.csv > ${sample_id}_annotatedcnv_filter.csv

   awk 'BEGIN {print "Chrom,Start,Type,End,SVLEN,Score,LOG2CNT,Gene"} 1' ${sample_id}_annotatedcnv_filter.csv > ${sample_id}_annotatedcnv_filter_header.csv

   cnv_mapping_occfusion_update.py  $calls_bed $occ_fusions ${sample_id}_tumor_copy.txt ${sample_id}_bins_filter.bed  $threshold

   cnv_mapping_occfusion_update_nofilter.py $calls_bed ${sample_id}_tumor_copy_number.txt   $threshold
   """
}

//#######################
//# Call small variants in ROI.bam file with clair3
//###############
//# ROI.bam file created on P24 by cross-referencing merged bam file with list of OCC genes
//# ROI.bam file transfered to office PC for variant calling and annotation
//# clair3 installed from https://github.com/HKU-BAL/Clair3

process clair3 {
    cpus 4
    memory '5 GB'
    label 'clair3'
    publishDir "${params.output_path}/OCC/$sample_id", mode: "copy", overwrite: true
   
    input:
    tuple val(sample_id), path(occ_bam), path(occ_bam_bai), path(reference_genome), path(reference_genome_bai),  path(refGene), path(hg38_refGeneMrna), path(clinvar), path(clinvarindex),path(hg38_cosmic100),path(hg38_cosmic100index)

    output:
    tuple val(sample_id), path('output_clair3/')
    tuple val(sample_id), path('occ_pileup_snvs_avinput')
    path("${sample_id}_occ_pileup_annotateandfilter.csv")
    path("${sample_id}_occ_pileup_annotateandfilter.csv"), emit:occpileupannotateandfilterout
    path('occ_merge_snv_avinpt')
    path('occ_merge.hg38_multianno.txt')
    path("${sample_id}_merge_annotateandfilter.csv"), emit:mergeannotateandfilterout
    tuple path("${sample_id}_occ_pileup_annotateandfilter.csv"), path("${sample_id}_merge_annotateandfilter.csv"), emit:clair3output 


    script:
   
   """ 

   run_clair3.sh \
    --bam_fn=$occ_bam \
    --ref_fn=$reference_genome  \
    --threads=8 \
    --var_pct_full=1 \
    --ref_pct_full=1 \
    --var_pct_phasing=1 \
    --platform="ont" \
    --no_phasing_for_fa \
    --model_path=${params.ref_dir}/r1041_e82_400bps_sup_v420 \
    --output=output_clair3
 
 convert2annovar.pl output_clair3/pileup.vcf.gz \
    --format vcf4 \
	--withfreq \
	--filter pass \
	--fraction 0.1 \
	--includeinfo \
	--outfile occ_pileup_snvs_avinput

   
   table_annovar.pl occ_pileup_snvs_avinput \
         -outfile occ_pileup \
         -buildver hg38 -protocol refGene,clinvar_20240611,cosmic100coding2024\
         -operation g,f,f \
         ${params.humandb_dir} \
         -otherinfo
      
    awk '/exonic/ && /nonsynonymous/ && !/Benign/ && !/Likely_benign/|| /upstream/ || /Func.refGene/ || /splicing/ && !/Benign/ && !/Likely_benign/' occ_pileup.hg38_multianno.txt \
| awk '/exonic/ || /TERT/ || /Func.refGene/'  \
| awk '!/dist=166/' \
| cut -f1-16,26,28,29 > ${sample_id}_occ_pileup_annotateandfilter.csv

convert2annovar.pl \
    output_clair3/merge_output.vcf.gz \
    --format vcf4 \
    --withfreq \
    --filter pass \
    --fraction 0.1 \
    --includeinfo \
    --outfile occ_merge_snv_avinpt

table_annovar.pl occ_merge_snv_avinpt \
    -outfile occ_merge \
    -buildver hg38 -protocol refGene,clinvar_20240611,cosmic100coding2024\
    -operation g,f,f \
    ${params.humandb_dir} \
    -otherinfo

    awk '/exonic/ && /nonsynonymous/ && !/Benign/ && !/Likely_benign/|| /upstream/ || /Func.refGene/ || /splicing/ && !/Benign/ && !/Likely_benign/ || /    frameshift/ && !/Benign/ && !/Likely_benign/' \
        occ_merge.hg38_multianno.txt \
    | awk '/exonic/ || /TERT/ || /Func.refGene/'  \
    | awk '!/dist=166/' \
    | cut -f1-16,26,28,29 \
    > ${sample_id}_merge_annotateandfilter.csv 

    """
   }


//#################################
//#### ClairS-TO
//#################################
//# ClairS-TO is a recent development to specifically call somatic variants in Tumor-only samples
//# It's run separate from Clair3
// # installed from https://github.com/HKU-BAL/ClairS-TO via micromamba

process clairs_to {
    cpus 4
    memory '2 GB'
    label 'clairsto'
    publishDir "${params.output_path}/OCC/$sample_id", mode: "copy", overwrite: true

    input:
    tuple val(sample_id), path(occ_bam), path(occ_bam_bai), path(reference_genome), path(reference_genome_bai),  path(refGene), path(hg38_refGeneMrna), path(clinvar), path(clinvarindex),path(hg38_cosmic100),path(hg38_cosmic100index), path(occ_snv_screening)
    
    output:
    tuple val(sample_id), path('clairsto_output/')
    tuple val(sample_id), path('clairS_To_snv_avinput')
    path('ClairS_TO_snv.hg38_multianno.txt')
    tuple val(sample_id), path("${sample_id}_annotateandfilter_clairsto.csv"), emit:annotateandfilter_clairstoout
    path("${sample_id}_merge_snv_indel_claisto.vcf.gz")


    script:

   """

   run_clairs_to \
    --tumor_bam_fn=$occ_bam \
    --ref_fn=$reference_genome  \
    --threads=8 \
    --platform="ont_r10_dorado_4khz" \
    --output_dir=clairsto_output \
    --bed_fn=$occ_snv_screening \


    bcftools merge --force-samples clairsto_output/snv.vcf.gz clairsto_output/indel.vcf.gz -o ${sample_id}_merge_snv_indel_claisto.vcf.gz

    convert2annovar.pl ${sample_id}_merge_snv_indel_claisto.vcf.gz \
   --format vcf4 \
   --filter pass \
   --includeinfo \
   --outfile  clairS_To_snv_avinput


  table_annovar.pl clairS_To_snv_avinput \
   -outfile ClairS_TO_snv \
   -buildver hg38 -protocol refGene,clinvar_20240611,cosmic100coding2024\
   -operation g,f,f \
    ${params.humandb_dir} \
   -otherinfo  

   awk '/exonic/ && /nonsynonymous/ && !/Benign/ || /upstream/ || /Func.refGene/' \
   ClairS_TO_snv.hg38_multianno.txt \
   | awk '/exonic/ || /TERT/ || /Func.refGene/'  \
  | awk '!/dist=166/' \
  | cut -f1-16,25,26  > ${sample_id}_annotateandfilter_clairsto.csv


    """
   }


// Simple R-script to merge annotation calls

process merge_annotation {
    cpus 4
    memory '2 GB'
    label 'merge_annotation'
    publishDir "${params.output_path}/merge_annot_clair3andclairsto/", mode: "copy", overwrite: true

    input:
    tuple val(sample_id), path(mergeannotateandfilterout), path(occpileupannotateandfilterout), path(annotateandfilter_clairstoout), path(occ_genes)
    
    output:
    tuple val(sample_id), path("${sample_id}_merge_annotation_filter_snvs_allcall.csv"), emit: occmergeout

    """

    merge_annotations_prospective.R $mergeannotateandfilterout $occpileupannotateandfilterout $annotateandfilter_clairstoout ${sample_id}_merge_annotation_filter_snvs_allcall.csv $occ_genes

    """
}
//# Uses igv_tools to extract 3 regions with recurring mutations.
process igv_tools {
    cpus 4
    memory '2 GB'
    label 'epic'


    publishDir "${params.output_path}/terp", mode: "copy", overwrite: true


    input:
    tuple val(sample_id),path(occ_bam), path(occ_bam_bai),  path(tertp_variants),path(ncbirefseq), path(reference_genome), path(reference_genome_bai)

    output:
    tuple val(sample_id), file("${sample_id}_tertp_id1.html"), emit: tertp_out_igv

    script:

    """
    export CURL_CA_BUNDLE=/etc/ssl/certs/ca-certificates.crt
    create_report $tertp_variants \
        --fasta $reference_genome\
        --flanking 1000 \
        --tracks $tertp_variants $occ_bam  $ncbirefseq \
        --output ${sample_id}_tertp_id1.html



    cp "${sample_id}_tertp_id1.html" "${params.output_path}/report/${sample_id}_tertp_id1.html"

    """
}

    process cramino_report {
        cpus 4
        memory '2 GB'
        label 'epic'
        publishDir "${params.output_path}/cramino", mode: "copy", overwrite: true


    input:
    tuple val(sample_id),path(merge_bam), path(merge_bam_bai), path(reference_genome), path(reference_genome_bai)

    output:
    tuple val(sample_id), file("${sample_id}_cramino_statistics.txt"), emit:craminostatout
    
    script:
    """
    source /opt/conda/etc/profile.d/conda.sh
    conda init
    conda activate annotatecnv_env

    cramino $merge_bam --reference $reference_genome > ${sample_id}_cramino_statistics.txt
   """

    }

    process markdown_report {
    cpus 4
    memory '2 GB'
    publishDir "${params.output_path}/report", mode: "copy", overwrite: true

    input:
    tuple val(sample_id), 
          val(craminoreport),
          val(nanodx), 
          val(dictionaire), 
          val(logo), 
          val(cnv_plot), 
          val(tumor_number), 
          val(annotatecnv),
          val(cnv_chr9),
          val(cnv_chr7), 
          val(mgmt_results),
          val(merge_results),
          val(annotSV_fusion), 
          path(terphtml),
          path(svannahtml), 
          path(annotsvhtml)

    output:
    file("${sample_id}_markdown_pipeline_report.pdf")


    script:
    """
    Rscript --version
    Rscript -e "rmarkdown::render('/home/chbope/Documents/nanopore/clone/nextflow/pipeline/bin/nextflow_markdown_pipeline3.Rmd',output_file=commandArgs(trailingOnly=TRUE)[17])" "${sample_id}" "${craminoreport}" "${nanodx}" "${dictionaire}" "${logo}" "${cnv_plot}" "${tumor_number}" "${annotatecnv}" "${cnv_chr9}" "${cnv_chr7}" "${mgmt_results}" "${merge_results}" "${annotSV_fusion}"  "${terphtml}" "${svannahtml}" "${annotsvhtml}" "\${PWD}/${sample_id}_markdown_pipeline_report.pdf"

    """
}


workflow analysis {
    take:
        epi2me_results  // Input channel from epi2me workflow
        
    main:
        params.run_mode = params.run_mode_analysis ?: 'all'
        println "Analysis run mode: ${params.run_mode}"

        // Validate run_mode parameter
        if (!['methylation', 'annotsv', 'cnv', 'occ', 'terp', 'mgmt', 'rmd', 'all'].contains(params.run_mode)) {
            error "ERROR: Invalid run_mode '${params.run_mode}' for analysis. Valid modes: methylation, annotsv, cnv, occ, terp, mgmt, rmd, all"
        }

        println "Running Analysis in mode: ${params.run_mode}"

        // Load sample thresholds only for analysis mode
        def sample_thresholds = params.run_order_mode ? [:] : 
            file(params.analyse_sample_id_file).readLines()
            .collectEntries { line ->
                def tokens = line.tokenize("\t")
                if (tokens.size() < 2) {
                    throw new IllegalArgumentException("Invalid line format in sample_id_file: ${line}")
                }
                def sample_id = tokens[0]
                def threshold = tokens[1]

                // Verify sample files exist before adding to thresholds
                def segs_bed = file("${params.segsfromepi2me_folder}/${sample_id}_segs.bed")
                def bins_bed = file("${params.segsfromepi2me_folder}/${sample_id}_bins.bed")
                def vcf_file = file("${params.segsfromepi2me_folder}/${sample_id}_segs.vcf")
                
                if (!segs_bed.exists() || !bins_bed.exists() || !vcf_file.exists()) {
                    println """
                    Warning: Skipping sample ${sample_id} - missing required files:
                    segs_bed: ${segs_bed} (exists: ${segs_bed.exists()})
                    bins_bed: ${bins_bed} (exists: ${bins_bed.exists()})
                    vcf_file: ${vcf_file} (exists: ${vcf_file.exists()})
                    """
                    return null
                }
                [(sample_id): threshold]
            }
            .findAll { it != null }  // Remove null entries
        
        if (!params.run_order_mode) {
            println "Sample thresholds: ${sample_thresholds}"
        }

        // Create bam_files channel based on run mode
        bam_files = params.run_order_mode ? 
            epi2me_results.map { sid, bam, bai, ref, ref_bai, tr_bed, modkit, segs_bed, bins_bed, segs_vcf, sv ->
                tuple(sid, bam, bai)
            } : Channel
            .fromList(sample_thresholds.keySet().collect())
            .map { sample_id -> 
                // Find all matching BAM files
                def bam_files = file("${params.merge_bam_folder}/${sample_id}*.bam")
                if (bam_files.isEmpty()) {
                    error "No BAM file found for sample ${sample_id} in ${params.bams_folder}"
                }
                if (bam_files.size() > 1) {
                    println "Warning: Multiple BAM files found for sample ${sample_id}, using: ${bam_files[0]}"
                }
                def bam = bam_files[0]
                def bai = file("${bam}.bai")
                
                if (!bai.exists()) {
                    throw new FileNotFoundException("Missing index file for ${bam}")
                }
                return tuple(sample_id, bam, bai)
            }

        // Create OCC BAM channel that can come from either mergebam output or direct input
        occ_bam_files = params.run_order_mode ? 
            Channel.empty() : Channel
            .fromList(sample_thresholds.keySet().collect())
            .map { sample_id ->
                // Try to find .occ.bam first, then fall back to regular .bam
                def bam = file("${params.occ_bam_folder}/${sample_id}.occ.bam")
                if (!bam.exists()) {
                    bam = file("${params.occ_bam_folder}/${sample_id}.bam")
                }
                if (!bam.exists()) {
                    log.warn("No BAM file found for sample ${sample_id} in ${params.occ_bam_folder}")
                    return null
                }
                def bai = file(bam.toString() + '.bai')
                if (!bai.exists()) {
                    log.warn("Missing index file for ${bam}")
                    return null
                }
                return tuple(sample_id, bam, bai)
            }
            .filter { it != null }

        // Create channels for all processes
        boosts_cramino = bam_files.map { sample_id, bam, bai -> 
            tuple(sample_id, bam, bai, file(params.reference_genome), file(params.reference_genome_bai))
        }

        boosts_clair3_channel = occ_bam_files.map { sample_id, bam, bai -> 
            tuple(sample_id, bam, bai, file(params.reference_genome), file(params.reference_genome_bai), 
                  file(params.refgene), file(params.hg38_refgenemrna), file(params.clinvar), 
                  file(params.clinvarindex), file(params.hg38_cosmic100), file(params.hg38_cosmic100index))
        }

        boosts_clairSTo_channel = occ_bam_files.map { sample_id, bam, bai -> 
            tuple(sample_id, bam, bai, file(params.reference_genome), file(params.reference_genome_bai),
                  file(params.refgene), file(params.hg38_refgenemrna), file(params.clinvar), 
                  file(params.clinvarindex), file(params.hg38_cosmic100), file(params.hg38_cosmic100index),
                  file(params.occ_snv_screening))
        }

        // Create SV channels based on run mode
        boosts_annotsv_channel = params.run_order_mode ?
            epi2me_results.map { sid, bam, bai, ref, ref_bai, tr_bed, modkit, segs_bed, bins_bed, segs_vcf, sv ->
                tuple(
                    sid,
                    sv,  // This is the _wf_sv.vcf.gz file from epi2me
                    file(sv.toString() + '.tbi'),
                    file(params.occ_fusions),
                    file(params.occ_fusion_genes_list)
                )
            } : Channel
            .fromPath("${params.sv_folder}/*.{vcf,vcf.gz}")
            .map { vcf -> 
                def sample_id = vcf.getName().replaceAll(/\.wf_sv\.(vcf|vcf\.gz)$/, '')
                def tbi = vcf.toString().endsWith('.gz') ? 
                    file(vcf.toString() + '.tbi') : 
                    file('NO_INDEX_NEEDED')
                tuple(
                    sample_id,
                    vcf,
                    tbi,
                    file(params.occ_fusions),
                    file(params.occ_fusion_genes_list)
                )
            }
            .filter { it[0] in sample_thresholds }

        boosts_svanna_channel = params.run_order_mode ?
            epi2me_results.map { sid, bam, bai, ref, ref_bai, tr_bed, modkit, segs_bed, bins_bed, segs_vcf, sv ->
                tuple(
                    sid,
                    sv,  // This is the _wf_sv.vcf.gz file from epi2me
                    file(sv.toString() + '.tbi'),
                    file(params.occ_fusions)
                )
            } : Channel
            .fromPath("${params.sv_folder}/*.{vcf,vcf.gz}")
            .map { vcf -> 
                def sample_id = vcf.getName().replaceAll(/\.wf_sv\.(vcf|vcf\.gz)$/, '')
                def tbi = vcf.toString().endsWith('.gz') ? 
                    file(vcf.toString() + '.tbi') : 
                    file('NO_INDEX_NEEDED')
                tuple(
                    sample_id,
                    vcf,
                    tbi,
                    file(params.occ_fusions)
                )
            }
            .filter { it[0] in sample_thresholds }

        boosts_epic_channel = bam_files.map { sample_id, bam, bai -> 
            tuple(sample_id, bam, bai, file(params.epicsites), file(params.mgmt_cpg_island_hg38))
        }

        boosts_segsfromepi2me_channel = params.run_order_mode ?
            epi2me_results.map { sid, bam, bai, ref, ref_bai, tr_bed, modkit, segs_bed, bins_bed, segs_vcf, sv ->
                tuple(sid, segs_vcf, segs_bed, bins_bed, file(params.occ_fusions), thresholdValue)
            } : Channel
            .fromList(sample_thresholds.keySet().collect())
            .map { sample_id -> 
                // Check if files exist with different naming patterns
                def segs_bed = file("${params.segsfromepi2me_folder}/${sample_id}_segs.bed")
                def bins_bed = file("${params.segsfromepi2me_folder}/${sample_id}_bins.bed")
                def vcf_file = file("${params.segsfromepi2me_folder}/${sample_id}_segs.vcf")
                
                // Verify files exist
                if (!segs_bed.exists() || !bins_bed.exists() || !vcf_file.exists()) {
                    error """
                    Missing required files for sample ${sample_id}:
                    segs_bed: ${segs_bed} (exists: ${segs_bed.exists()})
                    bins_bed: ${bins_bed} (exists: ${bins_bed.exists()})
                    vcf_file: ${vcf_file} (exists: ${vcf_file.exists()})
                    """
                }

                def threshold = sample_thresholds[sample_id]
                if (threshold == null) {
                    error "No threshold found for sample ${sample_id}"
                }

                tuple(
                    sample_id, 
                    vcf_file, 
                    file(params.occ_fusions),
                    bins_bed,
                    segs_bed, 
                    threshold.toFloat()
                )
            }

        // Create IGV channel
        boosts_igv_channel = occ_bam_files.map { sample_id, bam, bai -> 
            tuple(
                sample_id, 
                bam, 
                bai, 
                file(params.tertp_variants),
                file(params.ncbirefseq),
                file(params.reference_genome),
                file(params.reference_genome_bai)
            )
        }

        // Initialize channels
        terp_ch = Channel.empty()
        mgmt_ch = Channel.empty()
        markdown_ch = Channel.empty()

        if (params.run_mode in ['terp', 'all']) {
            println "Running TERP Analysis..."
            terp_ch = boosts_igv_channel | igv_tools
        }

        if (params.run_mode in ['mgmt', 'all']) {
            println "Running MGMT Analysis..."
            if (params.run_order_mode) {
                mgmt_ch = epi2me_results.map { sid, bam, bai, ref, ref_bai, tr_bed, modkit, segs_bed, bins_bed, segs_vcf, sv ->
                    tuple(sid, modkit, file(params.epicsites), file(params.mgmt_cpg_island_hg38))
                }
            } else {
                mgmt_ch = Channel
                    .fromPath("${params.bedmethyl_folder}/*.gz")
                    .map { gz ->
                        def sample_id = gz.getBaseName().split("\\.")[0]
                        tuple(sample_id, gz, file(params.epicsites), file(params.mgmt_cpg_island_hg38))
                    }
                    .filter { sample_id, gz, epicsites, mgmt_cpg ->
                        if (!sample_thresholds.containsKey(sample_id)) {
                            println "Warning: Sample ${sample_id} not found in sample_id_file, skipping..."
                            return false
                        }
                        return true
                    }
            }
            //mgmt_ch | mgmt_analysis
                extract_epic(mgmt_ch)
                MGMT_output = extract_epic.out.MGMTheaderout
                MGMT_sturgeon = extract_epic.out.sturgeonbedinput.map{sample_id, sturgeoninput -> [sample_id, sturgeoninput, params.sturgeon_model]}
                sturgeon(MGMT_sturgeon)
                mgmt_promoter(MGMT_output)
                mgmt_promoter_out=mgmt_promoter.out.mgmtresultsout //.view()
                mgmt_nanodx = extract_epic.out.epicselectnanodxinput.map{sample_id, epicselectnanodxinput -> [sample_id, epicselectnanodxinput, params.hg19_450model]}
                nanodx(mgmt_nanodx)
                nanodx_out = nanodx.out.nanodx450out.map{sample_id, nanodx450out -> [sample_id, nanodx450out, params.nanodx_450model, params.snakefile_nanodx, params.nn_model]}.view()
                run_nn_classifier(nanodx_out)
        }

        // Generate markdown report only if rmd mode
        if (params.run_mode in ['rmd', 'all']) {
            // Create MGMT channel for RMD mode
            mgmt_ch = Channel
                .fromPath("${params.bedmethyl_folder}/*.gz")
                .map { gz ->
                    def sample_id = gz.getBaseName().split("\\.")[0]
                    tuple(sample_id, gz, file(params.epicsites), file(params.mgmt_cpg_island_hg38))
                }
                .filter { it[0] in sample_thresholds }

            // Run MGMT analysis pipeline
            extract_epic(mgmt_ch)
            MGMT_output = extract_epic.out.MGMTheaderout
            
            // Run Sturgeon analysis
            MGMT_sturgeon = extract_epic.out.sturgeonbedinput
                .map { sample_id, sturgeoninput -> 
                    tuple(sample_id, sturgeoninput, params.sturgeon_model)
                }
            sturgeon(MGMT_sturgeon)
            
            // Run MGMT promoter analysis
            mgmt_promoter(MGMT_output)
            mgmt_promoter_out = mgmt_promoter.out.mgmtresultsout
            
            // Run NanoDx analysis
            mgmt_nanodx = extract_epic.out.epicselectnanodxinput
                .map { sample_id, epicselectnanodxinput -> 
                    tuple(sample_id, epicselectnanodxinput, params.hg19_450model)
                }
            nanodx(mgmt_nanodx)
            
            // Run neural network classifier
            nanodx_out = nanodx.out.nanodx450out
                .map { sample_id, nanodx450out -> 
                    tuple(sample_id, nanodx450out, params.nanodx_450model, params.snakefile_nanodx, params.nn_model)
                }
            run_nn_classifier(nanodx_out)
            rmd_nanodx_out = run_nn_classifier.out.rmdnanodx

            svannasv(boosts_svanna_channel)
            rmd_svanna_html = svannasv.out.rmdsvannahtml
            svannasv_out = svannasv.out.occsvannaannotationannotationvcf.map { sample_id, svannavcfoutput -> [ sample_id, svannavcfoutput, params.vcf2circos_json ]
    }
            circosplot(svannasv_out)
            circosplot_out = circosplot.out.circosout //.view()

            annotesv(boosts_annotsv_channel)
            annotsv_output = annotesv.out.annotatedvariantsout.map { sample_id, annotated_variants -> [ sample_id, annotated_variants, params.knotannotsv_conf ]
    }
            knotannotsv(annotsv_output)
            rmd_knotannotsv_html = knotannotsv.out.rmdannotsvhtml

            annotatecnv(boosts_segsfromepi2me_channel)
            rmd_annotatecnv_report = annotatecnv.out.rmdcnvtumornumber 

            clair3(boosts_clair3_channel)
            clair3_out = clair3.out.clair3output
            clairs_to(boosts_clairSTo_channel)
            clairs_to_out = clairs_to.out.annotateandfilter_clairstoout
            combine_file = clair3_out.combine(clairs_to_out).map { occ_pileup_annotateandfilter, merge_annotateandfilter, sample_id, annotateandfilter_clairsto -> [ sample_id, merge_annotateandfilter, occ_pileup_annotateandfilter, annotateandfilter_clairsto, params.occ_genes ]
        }
        merge_annotation(combine_file)
//    occmergefile_out = merge_annotation.out.occmergefile
        igv_tools(boosts_igv_channel)
        cramino_report(boosts_cramino)
        mergecnv_out = annotatecnv.out.rmdcnvtumornumber
                .combine(merge_annotation.out.occmergeout, by:0)
                .combine(run_nn_classifier.out.rmdnanodx, by: 0 )
                .combine(mgmt_promoter.out.mgmtresultsout, by:0)
                .combine(annotesv.out.annotsvfusion, by:0)
                .combine(svannasv.out.rmdsvannahtml, by:0)
                .combine(knotannotsv.out.rmdannotsvhtml, by:0)
                .combine(igv_tools.out.tertp_out_igv, by:0)
                .combine(cramino_report.out.craminostatout, by:0)
        mergecnv_out_map = mergecnv_out.map{def sample_id, cnv_plot, tumor_copy_number, annotatedcnv_filter_header, cnv_chr9, cnv_chr7, merge_annotation_filter_snvs_allcall, nanodx_classifier, mgmt_results, annotSV_fusion_extraction,svannahtml, annotsvhtml, terphtml, craminoreport -> [sample_id, craminoreport, nanodx_classifier,  params.nanodx_dictinaire, params.mardown_logo, cnv_plot, tumor_copy_number, annotatedcnv_filter_header, cnv_chr9, cnv_chr7, mgmt_results,
                        merge_annotation_filter_snvs_allcall, annotSV_fusion_extraction, terphtml, svannahtml, annotsvhtml]}.view()

        markdown_report(mergecnv_out_map).view()

    }

        // Add process execution blocks for each mode
        if (params.run_mode in ['cnv', 'all']) {
            println "Running CNV Analysis..."
            annotatecnv(boosts_segsfromepi2me_channel)
            cramino_report(boosts_cramino)
        }

        if (params.run_mode in ['annotsv', 'all']) {
            println "Running AnnotSV Analysis..."
            // Run SvannaSV analysis
            svannasv(boosts_svanna_channel)
            
            // Prepare output for circos plot
            svannasv_out = svannasv.out.occsvannaannotationannotationvcf
                .map { sample_id, svannavcfoutput -> 
                    tuple(sample_id, svannavcfoutput, file(params.vcf2circos_json))
                }
            
            // Generate circos plot
            circosplot(svannasv_out)
            
            // Run AnnotSV analysis
            annotesv(boosts_annotsv_channel)
            
            // Prepare output for knotAnnotSV
            annotsv_output = annotesv.out.annotatedvariantsout.view()
            .map { sample_id, annotated_variants -> 
                tuple(sample_id, annotated_variants, file(params.knotannotsv_conf))
                }
            
            // Run knotAnnotSV
            knotannotsv(annotsv_output)
        }

        if (params.run_mode in ['occ', 'all']) {
            println "Running OCC Analysis..."
            clair3(boosts_clair3_channel)
            clairs_to(boosts_clairSTo_channel)
            combine_file = clair3.out.clair3output
                .combine(clairs_to.out.annotateandfilter_clairstoout)
                .map { occ_pileup_annotateandfilter, merge_annotateandfilter, sample_id, annotateandfilter_clairsto -> 
                    tuple(sample_id, merge_annotateandfilter, occ_pileup_annotateandfilter, annotateandfilter_clairsto, file(params.occ_genes))
                }
            merge_annotation(combine_file)
        }

    emit:
        complete = Channel.of(true)
        results = terp_ch
}

// Keep only the necessary helper functions
def validateParameters() {
    params.run_mode = params.run_mode_analysis ?: 'all'
    println "Analysis run mode: ${params.run_mode}"

    if (!['methylation', 'annotsv', 'cnv', 'occ', 'terp', 'mgmt', 'rmd', 'all'].contains(params.run_mode)) {
        error "ERROR: Invalid run_mode '${params.run_mode}' for analysis. Valid modes: methylation, annotsv, cnv, occ, terp, mgmt, rmd, all"
    }
}

def loadSampleThresholds() {
    return file(params.analyse_sample_id_file).readLines()
        .collectEntries { line ->
            def tokens = line.tokenize("\t")
            if (tokens.size() < 2) {
                throw new IllegalArgumentException("Invalid line format in sample_id_file: ${line}")
            }
            [(tokens[0]) : tokens[1]]
        }
}

