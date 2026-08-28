#!/bin/bash

#==============================================================================
# Diana Pipeline Report Generator (Singularity Version)
#==============================================================================
#
# DESCRIPTION:
#   This script generates comprehensive PDF reports for nanopore whole genome 
#   sequencing (Diana) analysis results using R Markdown within a Singularity 
#   container. It processes multiple samples and creates detailed reports 
#   containing methylation analysis, structural variant annotation, copy number 
#   variation analysis, SNV calling results, and quality assessment metrics. 
#   This script must be used after the pipeline has finished running and the 
#   user can re-run specific process to generate the report.
#
# USAGE:
#   ./generate_report_singularity.sh [path_to_singularity_image]
#
# ARGUMENTS:
#   path_to_singularity_image: Optional path to the Singularity image file
#                             Defaults to: markdown_images_28feb2025_latest.sif
#
# REQUIREMENTS:
#   - Singularity container with R and required packages
#   - Sample IDs file with tumor content information
#   - All analysis results from the Diana pipeline
#   - R Markdown template file
#
# INPUT FILES:
#   - sample_ids.txt: Two-column file with sample ID and tumor content (decimal)
#     Location: /data/routine_diana/sample_ids.txt
#   - Various result files from routine_annotation/{sample_id}/ and routine_epi2me/{sample_id}/ directories
#   - The PATH for each result file is configured in the script.
#
# OUTPUT:
#   - PDF reports for each sample in routine_results/{sample_id}/ directory
#   - Format: {sample_id}_markdown_pipeline_report.pdf
#
# DEPENDENCIES:
#   - Singularity container with R and required packages
#   - Container must include: rmarkdown, data.table, kableExtra, and other R packages
#   - Default image: markdown_images_28feb2025_latest.sif
#
# AUTHOR: Diana Pipeline Development Team
# DATE: 2024
#==============================================================================

# Set the Singularity image path
SINGULARITY_IMAGE="${1:-markdown_images_28feb2025_latest.sif}"

# Check if Singularity image exists
if [ ! -f "$SINGULARITY_IMAGE" ]; then
    echo "Error: Singularity image '$SINGULARITY_IMAGE' not found!"
    echo "Please provide the correct path to the Singularity image as an argument:"
    echo "  ./generate_report_singularity.sh /path/to/markdown_images_28feb2025_latest.sif"
    exit 1
fi

echo "Using Singularity image: $SINGULARITY_IMAGE"

# Define base paths based on current pipeline structure
PIPELINE_DIR="path/diana"  # To be updated by the user
REFERENCE_PATH="${PIPELINE_DIR}/data/reference"
OUTPUT_PATH="output_path" #To be updated by the user

# Pipeline version string embedded in the report footer (kept in sync with nextflow.config manifest.version)
PIPELINE_VERSION=$(grep -m1 "version" "${PIPELINE_DIR}/nextflow.config" | sed -E "s/.*=\s*'([^']+)'.*/\1/")
PIPELINE_VERSION="${PIPELINE_VERSION:-unknown}"

# SNV filtering thresholds (kept in sync with conf/annotation.config defaults)
SNV_DEPTH_THRESHOLD="${SNV_DEPTH_THRESHOLD:-10}"
SNV_GQ_THRESHOLD="${SNV_GQ_THRESHOLD:-10}"

# Sample IDs file path (hardcoded location for routine processing)
samples_file="${OUTPUT_PATH}/sample_ids.txt"

# RMarkdown template file path
rmd_template="${PIPELINE_DIR}/bin/nextflow_markdown_pipeline_update_final.Rmd"

# Reference files shared across all samples
warning_img="${PIPELINE_DIR}/docs/warning.png"
protein_coding_bed="${REFERENCE_PATH}/roi.protein_coding.bed"
pancan_dictionaire="${REFERENCE_PATH}/nanoDx/static/pancan_devel_v5i_dictionary.txt"

# Check if required files exist
if [ ! -f "$samples_file" ]; then
    echo "Error: Sample IDs file not found at: $samples_file"
    exit 1
fi

if [ ! -f "$rmd_template" ]; then
    echo "Error: R Markdown template not found at: $rmd_template"
    exit 1
fi

# Loop through each sample ID in the samples file
while read -r sample_id tumor_content; do
    echo "Processing sample: $sample_id"

    # Build dynamic input paths for this sample based on routine_annotation/routine_epi2me structure
    ANNOTATION_PATH="${OUTPUT_PATH}/routine_annotation/${sample_id}"
    EPI2ME_PATH="${OUTPUT_PATH}/routine_epi2me/${sample_id}"

    craminoreport="${EPI2ME_PATH}/cramino/${sample_id}_cramino_statistics.txt"

    # Build a per-sample sample_ids_file, mirroring the priority logic in
    # modules/annotation.nf's markdown_report process: user-provided tumor
    # content (2-column sample_ids.txt) takes priority over the ACE-calculated
    # value, which is only used as a fallback when no user value is given.
    ace_threshold_file="${ANNOTATION_PATH}/cnv/ace/${sample_id}_ace_results/threshold_value.txt"
    sample_ids_file="$(mktemp "${TMPDIR:-/tmp}/sample_file_${sample_id}.XXXXXX")"

    NUM_COLS=$(awk 'NF{print NF; exit}' "${samples_file}")
    if [ "${NUM_COLS:-0}" -ge 2 ]; then
        if grep "^${sample_id}[[:space:]]" "${samples_file}" > "${sample_ids_file}"; then
            echo "Using user-provided tumor content for ${sample_id} from ${samples_file}"
        else
            echo "WARNING: Sample ${sample_id} not found in ${samples_file}, creating file with sample_id only"
            echo "${sample_id}" > "${sample_ids_file}"
        fi
    elif [ -f "${ace_threshold_file}" ]; then
        THRESHOLD_VALUE=$(cat "${ace_threshold_file}")
        printf '%s\t%s\n' "${sample_id}" "${THRESHOLD_VALUE}" > "${sample_ids_file}"
        echo "Using ACE-calculated tumor content for ${sample_id}: ${THRESHOLD_VALUE}"
    else
        cp "${samples_file}" "${sample_ids_file}"
        echo "WARNING: No tumor content available for ${sample_id} (single-column sample list, no ACE results at ${ace_threshold_file})"
    fi

    nanodx="${ANNOTATION_PATH}/classifier/nanodx/${sample_id}_nanodx_classifier.tsv"
    dictionaire="${REFERENCE_PATH}/nanoDx/static/Capper_et_al_dictionary.txt"
    logo="${REFERENCE_PATH}/log_update.pdf"
    cnv_plot="${ANNOTATION_PATH}/cnv/${sample_id}_cnv_plot_full.pdf"
    tumor_number="${ANNOTATION_PATH}/cnv/${sample_id}_tumor_copy_number.txt"
    annotatecnv="${ANNOTATION_PATH}/cnv/${sample_id}_annotatedcnv_filter_header.csv"
    cnv_chr9="${ANNOTATION_PATH}/cnv/${sample_id}_cnv_chr9.pdf"
    cnv_chr7="${ANNOTATION_PATH}/cnv/${sample_id}_cnv_chr7.pdf"
    mgmt_results="${ANNOTATION_PATH}/methylation/${sample_id}_MGMT_results.csv"
    merge_results="${ANNOTATION_PATH}/merge_annot_clair3andclairsto/${sample_id}_merge_annotation_filter_snvs_allcall.csv"
    fusion_events="${ANNOTATION_PATH}/structure_variant/svannasv/${sample_id}.sv_fusions_both_gbm.tsv"
    gbm_protein_fusion_events="${ANNOTATION_PATH}/structure_variant/svannasv/${sample_id}.sv_fusions_gbm_protein.tsv"
    tertphtml="${ANNOTATION_PATH}/coverage/${sample_id}_tertp_id1.html"
    svannahtml="${ANNOTATION_PATH}/structure_variant/svannasv/${sample_id}_roi_svanna_annotation.html"
    egfr_coverage="${ANNOTATION_PATH}/coverage/${sample_id}_egfr_coverage.pdf"
    idh1_coverage="${ANNOTATION_PATH}/coverage/${sample_id}_idh1_coverage.pdf"
    idh2_coverage="${ANNOTATION_PATH}/coverage/${sample_id}_idh2_coverage.pdf"
    tertp_coverage="${ANNOTATION_PATH}/coverage/${sample_id}_tertp_coverage.pdf"
    tsneplot="${ANNOTATION_PATH}/classifier/nanodx/${sample_id}_tsne_plot.pdf"
    snv_target_genes="${REFERENCE_PATH}/snv_target_genes.txt"
    nanodx_pancan="${ANNOTATION_PATH}/classifier/nanodx/${sample_id}_nanodx_classifier_pancan.tsv"
    tsneplot_pancan="${ANNOTATION_PATH}/classifier/nanodx/${sample_id}_tsne_plot_pancan.pdf"

    # Output PDF path - routine_results for final reports
    output_file="${OUTPUT_PATH}/routine_results/${sample_id}/${sample_id}_markdown_pipeline_report.pdf"

    # Create output directory if it doesn't exist
    mkdir -p "$(dirname "$output_file")"

    # Now call the Rscript using Singularity container
    singularity exec --bind /data:/data "$SINGULARITY_IMAGE" \
  	Rscript -e "rmarkdown::render('${rmd_template}', output_file=commandArgs(trailingOnly=TRUE)[24])" \
      "${sample_id}" \
      "${craminoreport}" \
      "${sample_ids_file}" \
      "${nanodx}" \
      "${dictionaire}" \
      "${logo}" \
      "${cnv_plot}" \
      "${tumor_number}" \
      "${annotatecnv}" \
      "${cnv_chr9}" \
      "${cnv_chr7}" \
      "${mgmt_results}" \
      "${merge_results}" \
      "${fusion_events}" \
      "${tertphtml}" \
      "${svannahtml}" \
      "${egfr_coverage}" \
      "${idh1_coverage}" \
      "${idh2_coverage}" \
      "${tertp_coverage}" \
      "${tsneplot}" \
      "${snv_target_genes}" \
      "${protein_coding_bed}" \
      "${output_file}" \
      "${PIPELINE_VERSION}" \
      "${SNV_DEPTH_THRESHOLD}" \
      "${SNV_GQ_THRESHOLD}" \
      "${warning_img}" \
      "${nanodx_pancan}" \
      "${pancan_dictionaire}" \
      "${tsneplot_pancan}" \
      "${gbm_protein_fusion_events}"

    # Clean up temporary R files and the per-sample sample_ids_file
    rm -rf /tmp/Rtmp*
    rm -f "${sample_ids_file}"

    # Clean up RMarkdown temporary files and folders in the output directory
    output_dir="$(dirname "$output_file")"
    output_basename="${sample_id}_markdown_pipeline_report"

    # Remove temporary files created by RMarkdown
    rm -rf "${output_dir}/${output_basename}_files"
    rm -f "${output_dir}/${output_basename}.tex"
    rm -f "${output_dir}/${output_basename}.log"
    rm -f "${output_dir}/${output_basename}.aux"

    echo "Finished sample: $sample_id (temporary files cleaned up)"

done < "${samples_file}"

echo "All samples processed successfully!"
