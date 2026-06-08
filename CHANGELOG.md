# My Pipeline: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### `Added`
- Added `--clair3_model sup|hac` parameter to select Clair3/ClairS-TO basecalling model at runtime
  - `sup` (Super Accuracy, default) uses `r1041_e82_400bps_sup_v520`
  - `hac` (High Accuracy) uses `r1041_e82_400bps_hac_v520`
  - Both paths defined in `conf/epi2me.config`; selection resolved in `modules/epi2me.nf` at runtime
  - Default set in `nextflow.config` as `clair3_model = "sup"`
- Updated Clair3 container to PyTorch-based `hkubal/clair3:latest` (conda env `clair3_v2`)
  - Supports newer PyTorch-format model files (`.pt`) required by v520+ models
  - `run_clair3` conda activation now falls back gracefully: `clair3_v2` → `clair3` → continues
- Added timestamped ROI BED file copy to `routine_results/{sample_id}/` for run traceability
  - Filename format: `roi.protein_coding_dd_mm_YYYY_HHmm.bed`
  - Allows tracking which ROI version was used for any given sample run
  - `copy_results_to_summary` now runs for all annotation modes (`run_mode_order`, `run_mode_epiannotation`, `run_mode_annotation`)
- Added Sturgeon classifier as optional — pipeline continues without `general.zip` if absent
  - Warning logged when model not found; all other analyses unaffected
  - README updated with Dropbox download link and note that file must remain as a zip archive
- Added both Capper et al. and pan-cancer (pancan) classifiers running by default in every pipeline run
  - New `crossNN_pancan` process outputs `{sample_id}_nanodx_classifier_pancan.tsv/txt`
  - New `tsne_plot_pancan` process outputs `{sample_id}_tsne_plot_pancan.pdf/html`
  - Pancan files copied to `routine_results/{sample_id}/` alongside Capper files
  - Executive summary shows two classification tables: Pancan first, then Capper et al.
  - Dedicated container mappings added for `crossNN_pancan` and `tsne_plot_pancan` in `conf/annotation.config`
  - Graceful fallback: if pancan tsne fails (e.g. corrupted training set), pipeline continues and report still generates
- Added `setup_vep_cache()` to `setup_pipeline.sh` — automatically downloads Ensembl VEP cache (v115, GRCh38) from Ensembl FTP into `data/humandb/homo_sapiens_refseq/` during installation
  - Download is **non-fatal**: if Ensembl FTP is unreachable the installation continues normally with ANNOVAR as the default annotator
  - Prints manual download instructions when the download fails so users can install the cache later
  - Skips silently if `homo_sapiens_refseq/` already exists
- Added preprint citation to README (Bope CD et al. 2026, https://doi.org/10.64898/2026.03.25.714119)
- Added `extract_roi` process to `modules/epi2me.nf` — runs for `snv` and `all` modes, directly feeding `run_clair3` and `run_clairs_to` without requiring a separate mergebam step
- Added missing `occ_bam_dir` and `roi_bed` parameters to `conf/epi2me.config`
- Added `extract_roi` container and cpus configuration to `conf/epi2me.config` (Singularity and Docker profiles)
- Added quality caution warnings to executive summary in PDF report
  - Low Tumor Content warning shown when estimated tumor content is ≤15%
  - Low Coverage warning shown when mean sequencing coverage is <30x
  - Both warnings are independent and display simultaneously when both thresholds are missed
  - Warnings use `warning.png` icon staged via Nextflow path and rendered with LaTeX `\includegraphics`
  - Fallback plain-text warning rendered when image is unavailable
- Added PDF printer compatibility improvements
  - `\pdfminorversion=4` in LaTeX header forces PDF 1.4 output (universally printable)
  - Ghostscript post-processing step in `annotation:markdown_report` flattens embedded graphics and fonts
  - Ghostscript step is soft — skips silently if `gs` is not installed on the host, no extra install required

### `Changed`
- Bundled ANNOVAR Perl scripts into `bin/` with permission from the ANNOVAR authors
  - Scripts included: `annotate_variation.pl`, `coding_change.pl`, `convert2annovar.pl`, `table_annovar.pl`, `index_annovar.pl`, `prepare_annovar_user.pl`
  - Created `THIRD_PARTY_LICENSES.md` with ANNOVAR licence notice and required citation (Wang et al. 2010, *NAR*)
  - README ANNOVAR section updated: manual copy instructions removed, licence notice and link to `THIRD_PARTY_LICENSES.md` added
- Changed default SNV annotator from `"vep"` to `"annovar"` in `nextflow.config`
  - VEP remains available via `--snv_annotator vep` or by editing `nextflow.config`
  - README annotator table and code examples updated to reflect new default
- Updated Zenodo record from `19232427` to `20596458` (DOI: `10.5281/zenodo.20596458`)
  - Updated `setup_pipeline.sh` (header comment, `ZENODO_RECORD` variable, summary banner)
  - Updated all DOI links in `README.md`
  - Rebuilt Zenodo archives: `humandb.tar.gz`, `reference_core.tar.gz`, `diana_dummy.tar.gz`
  - `humandb.tar.gz` no longer bundles the VEP cache (`homo_sapiens_refseq/`) — now downloaded directly from Ensembl during setup (see below)
- README: COSMIC section updated — added login link (https://cancer.sanger.ac.uk/cosmic/login) and replaced hardcoded `v100` with generic `v{version}` so instructions stay valid for future releases
- Removed `docs/pancan_devel_v5i_dictionary.txt` from repository — updated dictionary now distributed via Zenodo reference archive; README note pointing to `docs/` removed
- Renamed `occ_fusions_genes.txt` → `roi_fusions_genes.txt` and all related identifiers
  - Config param: `occ_fusion_genes_list` → `roi_fusion_genes_list` in `annotation.config` and `example.config`
  - Process input variable: `occ_fusions_genes` → `roi_fusions_genes` in `modules/annotation.nf`
  - Channel variable: `occ_fusion_genes_list_ch` → `roi_fusion_genes_list_ch` in `modules/annotation.nf`
  - Physical file renamed on disk: `data/reference/occ_fusions_genes.txt` → `roi_fusions_genes.txt`
  - README directory tree and notes updated accordingly
- Removed pan-cancer dictionary from `docs/` — updated file now distributed via Zenodo reference archive
  - `docs/pancan_devel_v5i_dictionary.txt` deleted from repository
  - README note pointing to `docs/` removed
- Removed ANNOVAR Perl scripts from repository (license restricts redistribution)
  - `annotate_variation.pl`, `coding_change.pl`, `convert2annovar.pl`, `table_annovar.pl` removed from git and added to `.gitignore`
  - README updated with download instructions from https://annovar.openbioinformatics.org
- COSMIC database no longer distributed — bundled `hg38_cosmic100coding2024.txt` is now a placeholder with empty IDs
  - README updated with full preparation commands and note that institutional registration is required
- README: added Sturgeon, Clair3 model, ANNOVAR, COSMIC, and update script documentation sections
- Moved `extract_roi` from `modules/mergebam.nf` into `modules/epi2me.nf`; ROI extraction now happens inside the epi2me workflow and the result is emitted as `occ_bam`
- Updated `main.nf` to remove the `occ_bams` argument from `epi2me()` calls and join annotation input via `epi2me_results.occ_bam` instead
- `extract_roi` now uses `task.cpus` instead of `params.threads` (fixes undefined-parameter warning)
- `sample_ids_bam.txt` delimiter handling now accepts both tab and space in all parsers
  - `modules/mergebam.nf`: changed `splitEachLine('\t')` to `splitEachLine(~/\t| /)`
  - `modules/epi2me.nf` and `main.nf`: changed `tokenize("\t")` to `trim().split(/\s+/)`
- Fixed SNV executive summary CLNSIG filter to catch combined ClinVar pathogenic classifications
  - Changed from exact match (`CLNSIG == "Pathogenic"`) to case-insensitive substring match
  - Now catches values like `"Pathogenic/Likely_pathogenic/risk_factor"`, `"Likely_pathogenic"`, `"Likely_pathogenic|Affects"`, etc.
  - Only affects the executive summary table; full SNV table is unchanged
- Fixed CLNSIG simplification to handle compound conflicting classifications
  - Previously only the bare `"Conflicting_classifications_of_pathogenicity"` was shortened to `"conflicting"`
  - Compound values like `"Conflicting_classifications_of_pathogenicity|risk_factor"` were left long and incorrectly passed the pathogenic filter
  - Now uses `sub()` to replace the long prefix in all compound forms → `"conflicting|risk_factor"` etc.
  - These values no longer match the pathogenic filter and are correctly excluded from the executive summary
- Added `"conflicting: Conflicting classifications of pathogenicity (ClinVar)"` to the SNV table legend so users understand the abbreviation

### `Fixed`
- Fixed `smart_sample_monitor_v2.sh` loading 0 samples when `set -eo pipefail` is active
  - Root cause: `((line_count++))` returns exit code 1 when count is 0, killing the subshell silently
  - Fix: changed to pre-increment `((++line_count))` and `((++valid_count))`
  - Changed `load_samples` output from `echo "${samples[@]}"` to `printf '%s\n'` for safe newline-separated capture
  - Changed sample capture from word-split string to `mapfile -t samples < <(load_samples)`
- Fixed `smart_sample_monitor_v2.sh` only logging one sample as completed when pipeline processes multiple samples
  - Root cause: completion check skipped samples already marked `completed` in memory from a prior run
  - Fix: log and count all samples with a report file present, regardless of prior status
- Fixed annotation module only processing 1 of N samples in `run_mode_order` and `run_mode_epiannotation`
  - Root cause: SNV/cramino barrier used `.cross()` which consumes the barrier signal once, silently dropping all samples after the first
  - Fix: replaced `.cross()` with `.combine()` so all samples pair with the single barrier signal
- Fixed `annotation:extract_epic` failing with `Unable to open file *.wf_mods.bedmethyl.gz`
  - Root cause: `intersectBed` in the container does not support gzip input natively
  - Fix: added `zcat -f` to decompress the bedmethyl file before passing to all `intersectBed` calls


- Added automated Zenodo upload script (`upload_to_zenodo.sh`) for reference file distribution
  - Supports creating new deposits, new versions, or uploading to existing drafts
  - Implements bucket-based upload API with progress bars
  - Includes metadata management and publishing workflow
  - Fixed multiple issues through versions 5.0-5.5 (stdout/stderr, variable scoping, Content-Type header)
- Added packaging script (`package_for_zenodo.sh`) to prepare reference files for Zenodo
  - Creates properly formatted archives: Assembly.zip, general.zip, humandb.tar.gz, reference_core.tar.gz, etc.
  - Handles file verification and size reporting
  - Provides upload instructions for both automated and manual workflows
- Added unified setup script (`setup_pipeline.sh`) for streamlined installation
  - Single command setup for both Docker and Singularity/Apptainer
  - Automatically downloads all reference files from Zenodo (DOI: 10.5281/zenodo.18802824)
  - Handles container downloads and Nextflow installation
  - Eliminates need for manual reference file downloads
- Added configurable SNV filtering thresholds in annotation configuration
  - `snv_depth_threshold` parameter (default: 10) for minimum sequencing depth filtering
  - `snv_gq_threshold` parameter (default: 10) for minimum Genotype Quality filtering
  - Filters variants in Markdown PDF reports while preserving all variants in raw VCF files
  - GQ filtering supports multiple values from different callers (keeps variant if ANY value meets threshold)
- Added comprehensive software version documentation
  - `SOFTWARE_VERSIONS.md` - Human-readable documentation of all tools and versions
  - `versions.yml` - Machine-readable YAML format for automation
  - `CONTAINERS.md` - Quick reference guide for container usage and troubleshooting
- Added VCF indexing (.tbi files) to epi2me:run_clairs_to process
  - Automatically creates tabix index files for ClairS-TO VCF outputs
  - Ensures downstream annotation processes can access indexed VCF files
  - Handles empty VCF files gracefully when no variants are found
- Added automatic log file organization
  - All Nextflow logs now written to `logs/` directory
  - Log filenames include run mode and sample ID for easy identification
  - Format: `trace-{mode}_{sample_id}_{timestamp}.txt`
  - HTML reports: `execution_report-{mode}_{sample_id}_{timestamp}.html`

### `Changed`
- Updated Zenodo record reference from 17589248 (v4) to 18802824 (v1)
  - New DOI: 10.5281/zenodo.18802824
  - Concept DOI: 10.5281/zenodo.18802823
  - Updated in `setup_pipeline.sh` to download from new record
- Simplified reference file structure to avoid redundant file movements
  - Files are now extracted directly to their final locations
  - `general.zip` remains compressed (required for Sturgeon classifier)
  - Other archives (Assembly.zip, r1041_e82_400bps_sup_v420.zip, svanna-data.zip) are automatically extracted
- Updated README.md to reflect automated setup workflow
  - Documented unified `setup_pipeline.sh` usage for both Docker and Singularity
  - Removed external download requirement for gencode.v48.annotation.gff3 (now included in reference_core.tar.gz)
  - Updated Quick Start section with automated setup commands
  - Simplified Required Reference Data section to emphasize automation
- Removed deprecated documentation files
  - `ZENODO_FILES_GUIDE.md` (superseded by packaging script)
  - `ZENODO_UPLOAD_GUIDE.md` (superseded by upload script)
  - `smart_sample_monitor.sh` (unused utility)
- **BREAKING**: Renamed "analysis" to "annotation" throughout the entire pipeline
  - Module renamed: `modules/analysis.nf` → `modules/annotation.nf`
  - Config renamed: `conf/analysis.config` → `conf/annotation.config`
  - Run mode renamed: `--run_mode_analysis` → `--run_mode_annotation`
  - Run mode renamed: `--run_mode_epianalyse` → `--run_mode_epiannotation`
  - Directory renamed: `routine_analysis/` → `routine_annotation/`
  - Updated all documentation, scripts, and configuration files
- Updated SNV variant caller abbreviations in reports
  - Changed ClairS-TO abbreviation from "C" to "S_TO" for clarity
  - Updated legend: "P:Clair3 Pileup, M:Clair3 Merged and S_TO:ClairS-TO Somatic tumour-only"
- Renamed process `run_nn_classifier` to `crossNN` for clarity
  - Updated in `modules/annotation.nf` (7 occurrences)
  - Updated in `conf/annotation.config`
  - Updated in `conf/example.config`
- Updated smart_sample_monitor scripts for consistency
  - `smart_sample_monitor_v2.sh`: Changed variable names from `analysis_config` → `annotation_config`
  - `smart_sample_monitor.sh`: Changed variable names from `analysis_config` → `annotation_config`
- Updated setup scripts
  - `setup_docker.sh`: Updated pipeline mode examples and config file references
  - `setup_singularity.sh`: Updated default config, directory structure, and mode examples
  - `run_pipeline_conda.sh`: Updated variable names and config file references
- Enhanced README.md documentation
  - Added SNV Filtering Configuration section with examples
  - Updated all pipeline mode references to use "annotation" terminology
  - Added SNV mode to epi2me and annotation module options
  - Updated container and configuration file references
- Modified ClairS-TO VCF file access in annotation workflow
  - Changed from passing individual VCF files to using clairsto_output_dir
  - Ensures .tbi index files are accessible alongside VCF files
  - Fixes "could not retrieve index file" errors in run_mode_order

### `Fixed`
- Fixed "tabix: command not found" error in run_mode_order
  - Root cause: VCF files from epi2me lacked .tbi index files in run_mode_order
  - Solution: Added tabix indexing to epi2me:run_clairs_to process
  - Removed redundant tabix calls from annotation:clairs_to_annotate process
- Fixed bcftools concat index file access issues
  - Updated annotation:clairs_to_annotate to reference VCF files via clairsto_output_dir
  - Ensures both .vcf.gz and .vcf.gz.tbi files are accessible in same directory
- Fixed trace-*.txt files cluttering pipeline root directory
  - Implemented automatic log organization to `logs/` directory
  - Moved 197+ existing trace files to organized location

## [1.0.1] - 2024-11-14

### `Added`
- Added `--run_mode_epianalyse` pipeline execution mode to run Epi2me and Analysis modules sequentially when merged BAM files already exist
- Added comprehensive "Pipeline Run Modes" documentation table in README.md
- Added BED file format validation notes in README.md

### `Fixed`
- Fixed tumor content not appearing in PDF report when user provides 2-column `sample_ids.txt` file (sample_id + tumor_content)
  - Updated `modules/analysis.nf` (lines 862-895) to always create local `sample_file.txt` in work directory
  - Implemented priority logic: user-provided tumor content → ACE-calculated tumor content → no tumor content
  - Resolved issue where `sample_ids.txt` was passed as string value instead of staged file in `run_mode_analysis`
- Fixed circos plot always generating empty output files
  - Updated `modules/analysis.nf` (lines 437-468) to properly detect gzipped VCF files (.vcf.gz)
  - Added gzip detection: uses `zcat` for compressed files, `grep` for uncompressed files
  - Added graceful error handling when vcf2circos fails (creates empty placeholder instead of pipeline failure)
- Fixed BED file formatting errors in `OCC.protein_coding.bed`
  - Corrected line 204: Removed trailing tab character (was causing 11 fields error)
  - Corrected line 206: Converted spaces to tabs for proper BED format
  - Corrected line 208: Removed extra field (11th field)
  - Corrected line 209: Removed extra fields (11th and 12th fields)
  - All 208 lines now have proper 10-field BED format required by bedtools

### `Changed`
- Consolidated BED file usage: replaced `occ_snv_screening` parameter with `occ_protein_coding_bed` throughout pipeline
  - Updated `modules/analysis.nf`: Replaced all 4 occurrences of `occ_snv_screening` with `occ_protein_coding_bed`
  - Updated `conf/analysis.config`: Removed `occ_snv_screening` parameter definition
  - Updated `conf/example.config`: Removed `occ_snv_screening` parameter definition
  - Unified protein-coding region file usage across mergebam, SNV screening, and analysis modules
- Enhanced README.md documentation
  - Added pipeline execution mode flags to module headers
  - Added usage examples for `--run_mode_epianalyse` mode
  - Updated reference files section to document `OCC.protein_coding.bed` usage and requirements
  - Improved directory structure examples to reflect current file naming
  - Fixed typo: "Zenado" → "Zenodo"
  - Improved tool name capitalization consistency (Clair3, ClairS-TO)

### `Removed`
- Removed redundant `OCC.SNV.screening.bed` reference file parameter (consolidated into `OCC.protein_coding.bed`)

## [1.0dev] - 2025-01-16

### `Added`
- Fusion event exon/intron annotation pipeline with new Python scripts:
  - `bin/breaking_point_bed_translocation_exon.py` - Extract breakpoints from VCF to BED format
  - `bin/create_gff3_with_introns.py` - Calculate intron regions from exon boundaries in GFF3 files
  - `bin/remove_duplicate_report_exon.py` - Process intersectBed output to extract gene and feature information
  - `bin/summarize_fusion_features.py` - Consolidate fusion annotations with exon/intron/CDS/UTR/intergenic features
  - `bin/annotate_intergenic_breakpoints.py` - Add intergenic annotations for unannotated breakpoints
- Enhanced report generation (`bin/nextflow_markdown_pipeline_update_finalexecsummary.Rmd`) with:
  - Conditional EGFR section (only shown if EGFR is in CNV table)
  - Conditional SNV explanatory text (only shown if SNVs are detected)
  - Italic formatting for all gene names in tables (SNV, CNV, Fusion Events)
  - *TERT*p formatting for TERT promoter variants
  - LaTeX packages for preventing title/content page break separation (`needspace`, `etoolbox`, `titlesec`)
- Enhanced t-SNE plotting script (`bin/crossnn_tsne_fixedupdate.R`) with improved styling and bigger/bolder unknown cross symbols
- Singularity-based report generation script (`bin/generate_report_singularity.sh`)
- Crossnnumap container image support in setup scripts (`setup_singularity.sh`, `setup_docker.sh`)
- Dockerfile for crossnnumap container (`dockerfiles/Dockerfile_tsne`)
- Enhanced MGMT methylation table in reports with percentage values and improved headers
- Comprehensive `.gitignore` patterns to exclude large output directories and temporary files
- Added AA Change column to SNV summary table in executive summary report (`bin/nextflow_markdown_pipeline_update_finalexecsummary.Rmd`)
- Added placeholder regions for acrocentric chromosome p-arms (chr13, chr14, chr15, chr21, chr22) in CNV plots to display empty space instead of missing data

### `Fixed`
- Fixed t-SNE plot process failure in `--run_mode_analysis rmd` mode
- Corrected R environment handling in Nextflow configuration to prevent conda conflicts
- Fixed memory allocation for t-SNE plot process (increased to 75 GB)
- Corrected UMAP parameters in t-SNE plotting (`--umap-n-neighbours 10`, `--umap-min-dist 0.5`, `--umap-pca-dim 100`)
- Fixed data loading in enhanced t-SNE script to use correct probe IDs (`Illumina_ID`)
- Resolved plotting issues with unknown cross symbols in t-SNE visualizations

### `Changed`
- Updated `modules/analysis.nf`:
  - Enhanced `svannasv_fusion_events` process to add exon/intron annotation pipeline
  - Integrated new Python scripts for fusion feature extraction
  - Added GFF3 enhancement step to include calculated introns
  - Added fusion event filtering to report only complete events (with both start and end breakpoints)
- Updated `modules/epi2me.nf` with fusion annotation workflow
- Updated `modules/analysis.nf` to use enhanced t-SNE script with corrected parameters
- Modified `nextflow.config` to conditionally disable R environment clearing for RMD mode
- Updated MGMT table headers in reports:
  - "Mean Methylation Pyro" → "Mean Methylation, CpG 76–79 (%)"
  - "Mean Methylation Full" → "Mean Methylation, CpG 1–98 (%)"
  - "Classification by Pyro" → "Classification by CpG 76–79"
  - "Classification by Full" → "Classification by CpG 1–98"
- Converted MGMT methylation values to percentages in final reports
- Enhanced SNV summary table processing to extract and display protein changes (p.XXX notation) in AA Change column
- Improved CNV visualization (`bin/CNV_function_new_update.R`) to correctly represent acrocentric chromosomes with visual empty space for missing p-arm regions

### `Removed`
- Obsolete `bin/generate_report.sh` (replaced by Singularity version)
- Obsolete `bin/nextflow_markdown_pipeline_update_final9sep.Rmd` (replaced by updated version)

### `Dependencies`
- Added crossnnumap container image for enhanced t-SNE plotting
- Updated R package dependencies for improved visualization

### `Deprecated`

---

## Release Links

[Unreleased]: https://github.com/VilhelmMagnusLab/nWGS_pipeline/compare/v1.0.1...HEAD
[1.0.1]: https://github.com/VilhelmMagnusLab/nWGS_pipeline/compare/v1.0dev...v1.0.1
[1.0dev]: https://github.com/VilhelmMagnusLab/nWGS_pipeline/releases/tag/v1.0dev 