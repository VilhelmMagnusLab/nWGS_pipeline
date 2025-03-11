# nWGS_pipeline: Nanopore Whole Genome Sequencing Pipeline

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A522.10.1-23aa62.svg)](https://www.nextflow.io/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed.svg?labelColor=000000)](https://www.docker.com/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Release](https://img.shields.io/badge/release-v1.0.0-blue.svg)](https://github.com/yourusername/nWGS_pipeline/releases)

## Introduction

nWGS_pipeline is a comprehensive bioinformatics pipeline for analyzing glioblastoma samples using Oxford Nanopore sequencing data. The pipeline integrates multiple analyses including CNV detection, methylation profiling, structural variant calling, and MGMT promoter status determination.

## Pipeline Summary

The pipeline consists of three main modules that can be run independently or sequentially:

### 1. Mergebam Pipeline
- Merges multiple BAM files per sample
- Extracts regions of interest using OCC.protein_coding.bed
- Quality control of merged BAMs
- Outputs merged and indexed BAM files

### 2. Epi2me Pipeline
- CNV calling using QDNAseq
  - Bin size: 50kb
  - Output: segs.bed, bins.bed, segs.vcf
- Structural variant calling using Sniffles2
  - Mosaic mode enabled
  - Output: *_wf_sv.vcf.gz
- Modified base calling using Modkit
  - Output: *_wf_mods.bedmethyl.gz

### 3. Analysis Pipeline
- MGMT promoter methylation analysis
  - Uses EPIC array sites
  - Methylation level calculation
- Methylation-based classification
  - Sturgeon classifier
  - NanoDx neural network classifier
- Structural variant annotation
  - AnnotSV annotation
  - SvANNA pathogenicity prediction
  - Fusion gene detection
- CNV analysis
  - ACE threshold determination
  - Copy number annotation
  - Chromosome visualization
- Report generation
  - Interactive HTML reports
  - IGV snapshots
  - Circos plots

## Required Containers

```bash
# Epi2me Containers
wget https://[container-registry]/wf-human-variation_latest.sif
wget https://[container-registry]/snifflesv252_update_latest.sif
wget https://[container-registry]/modkit_latest.sif
wget https://[container-registry]/wf-cnv_latest.sif
wget https://[container-registry]/wf-human-variation-snp_latest.sif
wget https://[container-registry]/wf-human-variation-str_latest.sif
wget https://[container-registry]/snpeff_latest.sif

# Analysis Containers
wget https://[container-registry]/ace_images_6mars2025_latest.sif
wget https://[container-registry]/annotsv_3.3.4--py311hdfd78af_1.sif
wget https://[container-registry]/clair3_latest.sif
wget https://[container-registry]/clairs-to_latest.sif
wget https://[container-registry]/sturgeon_amd64_21jan_latest.sif
wget https://[container-registry]/igv_report_amd64_latest.sif
wget https://[container-registry]/vcf2circos_latest.sif
wget https://[container-registry]/nanodx_env_latest.sif
wget https://[container-registry]/markdown_images_28feb2025_latest.sif
wget https://[container-registry]/annotcnv_images_27feb1025_latest.sif
wget https://[container-registry]/mgmt_nanopipe_amd64_18feb2025_cramoni_latest.sif
```

## Required Reference Data

```bash
# Reference Genome
wget https://[reference-data]/GRCh38/GCF_000001405.39_GRCh38.p13_genomic_chr_only_plus_mt.fa

# Annotation Files
wget https://[reference-data]/OCC.fusions.bed
wget https://[reference-data]/EPIC_sites_NEW.bed
wget https://[reference-data]/MGMT_CpG_Island.hg38.bed
wget https://[reference-data]/OCC.SNV.screening.bed
wget https://[reference-data]/TERTp_variants.bed
wget https://[reference-data]/human_GRCh38_trf.bed

# Classifier Models
wget https://[reference-data]/general.zip  # Sturgeon model
wget https://[reference-data]/nanoDx/static/Capper_et_al_NN.pkl  # NanoDx model
```

## Quick Start

1. Install dependencies:
```bash
# Install Nextflow
curl -s https://get.nextflow.io | bash

# Install Singularity
# Follow instructions at: https://sylabs.io/guides/3.0/user-guide/installation.html
```

## Configuration

### 1. Singularity Configuration

The pipeline uses Singularity containers. You need to configure the bind paths in the following files:

```nextflow
// In nextflow.config
profiles {
    epi2me_auto_singularity {
        singularity {
            enabled = true
            autoMounts = true
            // Update these paths to match your system
            runOptions = "-B /path/to/conda/env/bin:/usr/local/bin -B /path/to/user/home/"
            cacheDir = "singularity-images"
        }
    }
}
```

### 2. Path Configuration

Update the following paths in the configuration files:

#### In conf/analysis.config:
```nextflow
params {
    // Base paths
    path = "/path/to/base/directory"
    input_path = "${params.path}"
    output_path = "${params.path}/results"

    // Reference and annotation directories
    annotate_dir = "/path/to/annotations"
    ref_dir = "/path/to/reference"
    humandb_dir = "/path/to/annovar/humandb"

    // Tool-specific directories
    svanna_dir = "/path/to/svanna-cli-1.0.4/"
    bin_dir = "/path/to/bin/"
    epi2me_dir = "/path/to/epi2me/wf-human-variation-master/"
}
```

#### In conf/epi2me.config:
```nextflow
params {
    // Base path
    path = "/path/to/base/directory"
    
    // Epi2me paths
    epi2me_dir = "/path/to/epi2me/wf-human-variation-master/"
    episv = "${path}/results/epi2me/episv"
    epimodkit = "${path}/results/epi2me/epimodkit"
    epicnv = "${path}/results/epi2me/epicnv"
    
    // Input/Output paths
    merge_bam_folder = "/path/to/merge_bam/folder"
    epi2me_sample_id_file = "${path}/sample_ids.txt"
}
```

### 3. Container Paths

Update container paths in both config files:

```nextflow
// In conf/analysis.config
process {
    withName:annotesv {
        container = '/path/to/containers/annotsv_3.3.4--py311hdfd78af_1.sif'
    }
    // ... other container paths
}

// In conf/epi2me.config
process {
    container = '/path/to/containers/wf-human-variation_latest.sif'
    // ... other container paths
}
```

### 4. Sample Configuration

Create a sample ID file (sample_ids.txt) with the following format:
```
sample_id1   flowcell_id1
sample_id2   flowcell_id2
```

The file paths should follow this structure:
```
/path/to/base/directory/
├── sample_id1/
│   └── bam_pass/
│       ├── file1.bam
│       └── file2.bam
└── sample_id2/
    └── bam_pass/
        ├── file1.bam
        └── file2.bam
```


### Important Notes:

1. All paths should be absolute paths
2. Ensure read/write permissions for all directories
3. Container bind paths must include all required directories
4. The singularity cache directory needs sufficient storage space
5. Memory and CPU requirements can be adjusted in nextflow.config

## Input Requirements

1. Sample ID file format (sample_ids.txt):
```
sample_id1   flowcell_id1
sample_id2   flowcell_id2
```

2. Directory structure:
```
input_dir/
├── sample_id1/
│   └── bam_pass/
│       ├── file1.bam
│       └── file2.bam
└── sample_id2/
    └── bam_pass/
        ├── file1.bam
        └── file2.bam
```

## Output Structure

```
results/
├── mergebam/
│   ├── merge_bam/
│   └── occ_bam/
├── epi2me/
│   ├── episv/
│   ├── epimodkit/
│   └── epicnv/
└── analysis/
    ├── cnv/
    ├── sv/
    ├── methylation/
    └── reports/
```

## 2. **Running Pipeline**

### Pipeline Modes
```bash
# Run complete workflow in order
nextflow run main.nf --run_order_mode

# Run individual modules
nextflow run main.nf --run_mode_mergebam
nextflow run main.nf --run_mode_epi2me all
nextflow run main.nf --run_mode_analysis rmd
```
#### Epi2me Pipeline Modes
```bash
# Run all Epi2me analyses
nextflow run main.nf --run_mode_epi2me all

# Run specific Epi2me analyses
nextflow run main.nf --run_mode_epi2me cnv    # Run only CNV analysis
nextflow run main.nf --run_mode_epi2me sv     # Run only structural variant calling
nextflow run main.nf --run_mode_epi2me modkit # Run only modified base calling
```

#### Analysis Pipeline Modes
```bash
# Run all analyses
nextflow run main.nf --run_mode_analysis all

# Run specific analyses
nextflow run main.nf --run_mode_analysis occ     # Run only OCC analysis
nextflow run main.nf --run_mode_analysis mgmt    # Run only MGMT analysis
nextflow run main.nf --run_mode_analysis annotsv # Run only AnnotSV analysis
nextflow run main.nf --run_mode_analysis cnv     # Run only CNV analysis
nextflow run main.nf --run_mode_analysis terp    # Run only TERP analysis
nextflow run main.nf --run_mode_analysis rmd     # Generate only the pdf report
```

## Acknowledgment and General Information

### Citing this Workflow
If you use this pipeline in your research, please cite:
```
[Citation details to be added]
```

### Acknowledgments
We would like to thank:
- The Nextflow community for their excellent framework
- Oxford Nanopore Technologies for their sequencing technology and tools
- All contributors and testers of this pipeline

### License
This project is licensed under the MIT License - see the LICENSE file for details.

### Funding
This work was supported by:
- [Funding details to be added]
- [Grant numbers to be added]

### Download
The latest version of the pipeline can be downloaded from:
```bash
git clone https://github.com/yourusername/nWGS_pipeline.git
```

## Questions and Feedback

For questions, bug reports, or feature requests, please contact:

**Maintainers:**
- Christian Bope (chbope@ous-hf.no / christianbope@gmail.com)
- Skabbi (skabbi@gmail.com)

You can also:
1. Open an issue on GitHub
2. Submit a pull request with improvements
3. Contact the maintainers directly via email

## Citations

If you use this pipeline, please cite:
- [Citations to be added]

# Mergebam Pipeline

Example of actual paths:
```
input_dir/
├── T10-01/
│   ├── 20231215_1340_3E_PAM69496_5c1d2ed7/bam_pass/
│   │   ├── PAM69496_pass_barcode01_*.bam
│   │   └── PAM69496_pass_barcode01_*.bam.bai
│   └── 20231216_1420_3E_PAM69496_7d4e9fc2/bam_pass/
│       ├── PAM69496_pass_barcode01_*.bam
│       └── PAM69496_pass_barcode01_*.bam.bai
└── T10-02/
    └── 20231217_1510_3E_PAM69497_8f3g1hj4/bam_pass/
        ├── PAM69497_pass_barcode02_*.bam
        └── PAM69497_pass_barcode02_*.bam.bai
```

## Configuration

In your `nextflow.config`, specify the input directory:

```nextflow
params {
    input_dir = "/path/to/nanopore/data"  // Directory containing sample folders
    // The pipeline will automatically find all BAM files in */bam_pass/ subdirectories
}
```

## Usage

Run the pipeline with:


The pipeline will:
1. Find all BAM files in the `*/bam_pass/` directories for each sample
2. Merge BAM files for each sample
3. Create index files for merged BAMs
4. Extract regions of interest

## Expected Output

```
results/
├── merged_bams/
│   ├── T10-01.merged.bam
│   ├── T10-01.merged.bam.bai
│   ├── T10-02.merged.bam
│   └── T10-02.merged.bam.bai
└── occ_bam/
    ├── T10-01_roi.bam
    ├── T10-01_roi.bam.bai
    ├── T10-02_roi.bam
    └── T10-02_roi.bam.bai
```
## Requirements

- Nextflow >= 21.04.0
- SAMtools >= 1.13
- Sufficient disk space for merged BAM files
- Memory requirements depend on BAM file sizes

## Input Directory Structure

The pipeline expects Nanopore sequencing data with BAM files in a specific structure:

```
input_dir/
├── sample_id1/
│   └── */bam_pass/
│       ├── *.bam
│       └── *.bam.bai
├── sample_id2/
│   └── */bam_pass/
│       ├── *.bam
│       └── *.bam.bai
└── sample_id3/
    └── */bam_pass/
        ├── *.bam
        └── *.bam.bai
```


