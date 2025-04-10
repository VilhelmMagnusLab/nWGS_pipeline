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
  - Output: *_wf_sv.vcf.gz
- Modified base calling using Modkit
  - Output: *_wf_mods.bedmethyl.gz

### 3. Analysis Pipeline
- MGMT promoter methylation analysis
  - Uses EPIC array sites
  - Methylation level calculation
- Methylation-based classification
  - NanoDx neural network classifier
- Structural variant annotation
  - AnnotSV annotation
  - Svanna pathogenicity prediction
  - Fusion gene detection
- CNV analysis
  - ACE to determine tumor content
  - Copy number annotation
  - Chromosome visualization
- Report generation
  - Interactive HTML reports
  - IGV snapshots
  - Circos plots

## Required Containers
### Epi2me Containers: 

```bash
# Epi2me Containers
# The following containers are required to run the Epi2me pipeline. The containers are downloaded from epi2me docker hub repository (https://hub.docker.com/r/epi2me/epi2me-workflows).
wget https://hub.docker.com/r/ontresearch/wf-common
wget https://hub.docker.com/r/ontresearch/wf-human-variation
wget https://hub.docker.com/r/ontresearch/wf-human-variation-snp
wget https://hub.docker.com/r/ontresearch/wf-human-variation-str
wget https://hub.docker.com/r/ontresearch/wf-human-variation-sv
wget https://hub.docker.com/r/ontresearch/modkit
wget https://hub.docker.com/r/ontresearch/wf-cnv

 ```

### Analysis Containers: 

.

```bash
# Analysis Containers
# The following containers are required to run the analysis pipeline. The containers are downloaded from the container registry (https://hub.docker.com/repositories/vilhelmmagnuslab/)
wget https://hub.docker.com/repositories/vilhelmmagnuslab/ace_1.24.0
wget https://hub.docker.com/repositories/vilhelmmagnuslab/annotsv_3.3.4--py311hdfd78af_1
wget https://hub.docker.com/repositories/vilhelmmagnuslab/clair3_amd64
wget https://hub.docker.com/r/hkubal/clairs-to
wget https://hub.docker.com/repositories/vilhelmmagnuslab/igv_report_amd64
wget https://hub.docker.com/repositories/vilhelmmagnuslab/vcf2circos
wget https://hub.docker.com/repositories/vilhelmmagnuslab/nanodx_env
wget https://hub.docker.com/repositories/vilhelmmagnuslab/markdown_images_28feb2025
wget https://hub.docker.com/repositories/vilhelmmagnuslab/annotcnv_images_27feb1025
wget https://hub.docker.com/repositories/vilhelmmagnuslab/mgmt_nanopipe_amd64_18feb2025_cramoni
wget https://hub.docker.com/repositories/vilhelmmagnuslab/nwgs_default_images_latest
```

## Required Reference Data

The required reference files are provided in the `refdata` folder:

```bash
refdata/
├── OCC.fusions.bed           # Fusion genes bed file
├── EPIC_sites_NEW.bed       # EPIC methylation sites
├── MGMT_CpG_Island.hg38.bed # MGMT CpG islands
├── OCC.SNV.screening.bed    # SNV screening regions
├── TERTp_variants.bed       # TERT promoter variants
├── hg38_refGene.txt         # RefGene annotation
├── hg38_refGeneMrna.fa      # RefGene mRNA sequences
├── hg38_clinvar_20240611.txt # ClinVar annotations
├── hg38_cosmic100coding2024.txt # Cosmic annotations
└── human_GRCh38_trf.bed    # Tandem repeat regions

reference_genome = "${params.ref_dir}/ref.fa" need to be provided or downloaded from the following link: https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/reference-genomes/human-reference-genomes/
```

These files are essential for:
- Methylation analysis (EPIC sites, MGMT)
- Structural variant analysis (Fusions)
- Copy number analysis
- SNV detection
- TERT promoter analysis

# Classifier Models

```bash
# The NanoDx model is downloaded from the following link: https://gitlab.com/pesk/nanoDx. The user need to download the model from the link and provide the path to the model in the analysis.config file (nanodx_workflow_dir). It is also required to make sure that the file Capper_et_al_NN.pkl is in the /$Path/nanoDx/static.

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
    bin_dir = "/path/to/bin/"
    epi2me_dir = "/path/to/epi2me/wf-human-variation-master/"
}
```

#### In conf/epi2me.config:
The epi2me folder need to download the epi2me workflow from the following link: https://github.com/epi2me-labs/wf-human-variation. In the epi2me.config file the user need to specify the path to the epi2me workflow.
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
        container = '/path/to/containers/annotsv*.sif'
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

#### Input structure for running directly from sample ID file:

1. Create a sample ID file with the following format: This file is used to run the analysis pipeline knowing the tumor content of the sample.
```
sample_id1   tumor content (float)
sample_id2   tumor content (float)
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


2. Sample ID file format used to merge the bam files: The sample ID and the flowcell ID is considered both to avoid mixing up the bam files from different flowcells. For this case the tumor content is calculated using ACE packages and inputted directly. The user can also analyse the tumur content generated by ACE and select the best fit and rerun the analyis with an updated tumor content has describes in step 1. 
```
sample_id1   flowcell_id1
sample_id2   flowcell_id2
```

#### Mergebam input structure:


```
input_dir/
├── V1001/
│   ├── 20231215_1340_3E_PAM69496_5c1d2ed7/bam_pass/
│   │   ├── PAM69496_pass_barcode01_*.bam
│   │   └── PAM69496_pass_barcode01_*.bam.bai
│   └── 20231216_1420_3E_PAM69496_7d4e9fc2/bam_pass/
│       ├── PAM69496_pass_barcode01_*.bam
│       └── PAM69496_pass_barcode01_*.bam.bai
└── V1002/
    └── 20231217_1510_3E_PAM69497_8f3g1hj4/bam_pass/
        ├── PAM69497_pass_barcode02_*.bam
        └── PAM69497_pass_barcode02_*.bam.bai
```

#### Expected from mergebam input:

```
results/
├── merged_bams/
│   ├── V1001.merged.bam
│   ├── V1001.merged.bam.bai
│   ├── V1002.merged.bam
│   └── V1002.merged.bam.bai
└── occ_bam/
    ├── V1001_roi.bam
    ├── V1001_roi.bam.bai
    ├── V1002_roi.bam
    └── V1002_roi.bam.bai
```

### Output results folder structure:

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
### Important Notes:

1. All paths should be absolute paths
2. Ensure read/write permissions for all directories
3. Container bind paths must include all required directories
4. The singularity cache directory needs sufficient storage space
5. Memory and CPU requirements can be adjusted in nextflow.config

## 2. **Running Pipeline**

### Pipeline Modes
```bash
# Run complete workflow in order
nextflow run main.nf --run_order_mode

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
nextflow run main.nf --run_mode_analysis tertp    # Run only TERTP analysis
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
git clone https://github.com/VilhelmMagnusLab/nWGS_pipeline.git
```

## Questions and Feedback

For questions, bug reports, or feature requests, please contact:

**Maintainers:**
- Christian Bope (chbope@ous-hf.no / christianbope@gmail.com)
- Skabbi (skahal@ous-hf.no / skabbi@gmail.com)

You can also:
1. Open an issue on GitHub
2. Submit a pull request with improvements
3. Contact the maintainers directly via email

## Citations

If you use this pipeline, please cite:
- [Citations to be added]



