# Software Versions - nWGS Pipeline

This document lists all software tools, packages, and their versions used in the nWGS (nanopore Whole Genome Sequencing) pipeline for CNS tumor analysis.

**Pipeline Version:** 1.0.1
**Last Updated:** 2026-02-26
**Authors:** Christian Domilongo Bope, Skarphéðinn Halldórsson, Richard Nagymihaly

---

## Pipeline Framework

| Tool | Version | Purpose |
|------|---------|---------|
| Nextflow | >=23.10.1 | Workflow management system |
| Singularity/Apptainer | Latest | Container runtime |

---

## Core Analysis Containers

### 1. SNV/Indel Calling

#### Clair3 (`clair3_amd64`)
**Purpose:** Long-read SNV and indel variant calling

| Component | Version |
|-----------|---------|
| Base Image | continuumio/miniconda3:latest |
| Python | 3.9.0 |
| samtools | 1.15.1 |
| whatshap | 1.7 |
| pypy3.6 | latest |
| tensorflow-cpu | 2.8.0 |
| tensorflow-addons | latest |
| pytables | latest |
| mpmath | 1.2.1 |
| cffi | 1.14.4 |
| parallel | 20191122 |
| pigz | latest |
| zstd | latest |

**Models:** Clair3 pre-trained models from http://www.bio8.cs.hku.hk/clair3/clair3_models/

#### ClairS-TO (`clairsto_amd64`)
**Purpose:** Somatic variant calling for tumor-only samples

| Component | Version |
|-----------|---------|
| Base Image | ubuntu:16.04 |
| micromamba | 1.5.1-2 |
| clair3 | latest |
| bcftools | latest |
| pytorch | latest |
| einops | latest |
| tqdm | latest |
| torchinfo | latest |
| scipy | latest (pip) |
| scikit-learn | latest (pip) |

**Models & Databases:**
- ClairS-TO models: http://www.bio8.cs.hku.hk/clairs-to/models/
- ClairS-TO databases: http://www.bio8.cs.hku.hk/clairs-to/databases/
- CNA reference files: http://www.bio8.cs.hku.hk/clairs-to/cna_data/

---

### 2. Structural Variant Calling

#### Sniffles v2 (`snifflesv252_update`)
**Purpose:** Structural variant detection from long-read alignments

| Component | Version |
|-----------|---------|
| Base Image | ubuntu:22.04 |
| Python | 3.10 |
| sniffles | latest (pip) |
| pysam | latest (pip) |
| truvari | latest (pip) |
| edlib | 1.3.9 |
| psutil | 5.9.4 |
| samtools | latest (apt) |
| bcftools | latest (apt) |
| tabix | latest (apt) |

---

### 3. Copy Number Variation Analysis

#### QDNAseq (`qdnaseq_amd64`)
**Purpose:** Copy number aberration detection from sequencing data

| Component | Version |
|-----------|---------|
| Base Image | rocker/r-ver:4.3.2 OR ubuntu:22.04 |
| R | 4.3.2 |
| Bioconductor | 3.18 |
| samtools | latest (apt) |

**R Packages:**
- argparser (CRAN)
- BiocManager (CRAN)
- QDNAseq (Bioconductor 3.18)
- QDNAseq.hg38 (Bioconductor 3.18 / GitHub: asntech/QDNAseq.hg38)
- QDNAseq.hg19 (Bioconductor 3.18)
- AnnotationHub (Bioconductor 3.18)
- AnnotationDbi (Bioconductor 3.18)
- remotes (CRAN)

#### ACE (`ace_1.24.0`)
**Purpose:** Absolute Copy number Estimation

| Component | Version |
|-----------|---------|
| Base Image | ubuntu:20.04 |
| R | 4.2 (conda) |
| Miniconda | Latest Linux x86_64 |

**Bioconductor Packages:**
- BiocManager
- genomicranges
- rtracklayer
- qdnaseq
- biobase
- QDNAseq
- QDNAseq.hg19
- ACE
- QDNAseq.hg38 (GitHub: asntech/QDNAseq.hg38@main)

**R Packages:**
- r-devtools
- r-remotes
- r-optparse
- r-data.table
- r-ggplot2
- r-reshape2
- r-gridextra

#### CNV Annotation (`annotcnv_images_27feb1025`)
**Purpose:** CNV annotation and visualization

| Component | Version |
|-----------|---------|
| Base Image | continuumio/miniconda3 |
| R | 4.2.3 |
| Python | 3.13.2 |
| cramino | 0.16.0 |
| bedtools | 2.31.1 |
| numpy | 2.2.3 |
| pandas | 2.2.3 |
| scipy | 1.15.2 |
| wkhtmltopdf | latest |

**Key R Packages (80+ total):**
- r-ggplot2: 3.5.1
- r-data.table: 1.15.4
- r-dplyr: 1.1.4
- r-plotly: 4.10.4
- r-rmarkdown: 2.27
- r-knitr: 1.47
- r-remotes
- r-ggpp: 0.5.8-1

---

### 4. DNA Methylation Analysis

#### Modkit (`modkit`)
**Purpose:** Modified base calling and methylation analysis

| Component | Version |
|-----------|---------|
| Container | modkit_latest.sif |
| Purpose | Extract 5mC methylation calls from BAM |

#### Sturgeon (`sturgeon_amd64_21jan`)
**Purpose:** DNA methylation-based CNS tumor classification

| Component | Version |
|-----------|---------|
| Base Image | python:3.9-slim |
| Python | 3.9 |
| sturgeon | latest (GitHub: marcpaga/sturgeon) |
| procps | latest |

---

### 5. Quality Control & Statistics

#### Cramino (`mgmt_nanopipe_amd64_18feb2025_cramoni`)
**Purpose:** CRAM/BAM quality metrics and statistics

| Component | Version |
|-----------|---------|
| Base Image | continuumio/miniconda3 |
| cramino | 0.16.0 |
| bedtools | 2.31.1 |
| r-base | 4.2.3 |
| Python | 3.13.2 |
| numpy | 2.2.3 |
| pandas | 2.2.3 |
| scipy | 1.15.2 |
| wkhtmltopdf | latest |

**R Packages (comprehensive list):**
- r-ggplot2: 3.5.1
- r-data.table: 1.15.4
- r-dplyr: 1.1.4
- r-plotly: 4.10.4
- r-rmarkdown: 2.27
- r-knitr: 1.47

---

### 6. Classification & Machine Learning

#### NanoDx (`nanodx_images_3feb25`)
**Purpose:** Neural network-based tumor classification

| Component | Version |
|-----------|---------|
| Base Image | continuumio/miniconda3 |
| Python | 3.10.2 |
| numpy | 1.22.4 |
| pandas | 1.5.1 |
| scikit-learn | 1.0.2 |
| pytorch | 1.12.1 |
| snakemake | latest |
| nextflow | latest |

#### t-SNE Visualization (`crossnnumap`)
**Purpose:** Dimensionality reduction and visualization

| Component | Version |
|-----------|---------|
| Base Image | ubuntu:24.04 |
| Miniconda | Latest Linux x86_64 |

**R Packages (via conda):**
- r-optparse
- r-data.table
- r-dplyr
- r-stringr
- bioconductor-rhdf5
- r-matrixstats
- r-ggplot2
- r-rtsne
- r-htmlwidgets
- r-r.utils
- r-ggtext
- r-uwot
- r-plotly

**Environment Variables:**
- HDF5_USE_FILE_LOCKING: FALSE
- HDF5_DISABLE_VERSION_CHECK: 1
- OMP_NUM_THREADS: 4
- R_MAX_VSIZE: 8GB

---

### 7. Visualization & Reporting

#### VCF2Circos (`vcf2circos`)
**Purpose:** Circos plot generation from VCF files

| Component | Version |
|-----------|---------|
| Base Image | python:3.9.16 |
| Python | 3.9.16 |
| vcf2circos | latest (GitHub: JbaptisteLam/vcf2circos@manuscript) |

**Configuration:** https://www.lbgi.fr/~lamouche/vcf2circos/config_vcf2circos_29032023.tar.gz

#### IGV Report (`igv_report_amd64`)
**Purpose:** Interactive Genome Viewer reports

| Component | Version |
|-----------|---------|
| Container | igv_report_amd64_latest.sif |
| Purpose | Generate IGV.js HTML reports |

#### Genomic Region Plots (`gviz_amd64ps`)
**Purpose:** Gene coverage visualization

| Component | Version |
|-----------|---------|
| Container | gviz_amd64ps_latest.sif |
| Purpose | Generate coverage plots for target genes |

#### Markdown Reports (`markdown_images_28feb2025`)
**Purpose:** R Markdown report generation

| Component | Version |
|-----------|---------|
| Base Image | continuumio/miniconda3 |
| System Tools | wkhtmltopdf, pandoc, texlive-full |

**R Packages:**
- r-base
- r-ggplot2
- r-scales
- r-data.table
- r-rmarkdown
- r-dplyr
- r-knitr
- r-kableextra
- r-webshot
- r-pagedown
- phantomjs
- pandoc

---

### 8. Default Tools Container (`nwgs_default_images`)

**Purpose:** General-purpose tools for various pipeline processes

| Component | Version |
|-----------|---------|
| Container | nwgs_default_images_latest.sif |
| Purpose | Standard bioinformatics utilities |

**Used by processes:**
- extract_epic
- sturgeon
- svannasv
- nanodx
- merge_annotation
- tertp

---

## Container Registry

All containers are available from:
- **Docker Hub:** `vilhelmmagnuslab/[container-name]:latest`
- **Local Singularity:** `${params.nWGS_dir}/containers/[container-name].sif`

---

## Reference Databases & Models

### Pre-trained Models
- **Clair3 Models:** http://www.bio8.cs.hku.hk/clair3/clair3_models/
- **ClairS-TO Models:** http://www.bio8.cs.hku.hk/clairs-to/models/
- **ClairS-TO Databases:** http://www.bio8.cs.hku.hk/clairs-to/databases/
- **CNA Reference Files:** http://www.bio8.cs.hku.hk/clairs-to/cna_data/

### Reference Genomes
- GRCh38 (hg38)
- GRCh37 (hg19)

### Annotation Resources
- QDNAseq reference bins (hg38, hg19)
- NanoDx classification models
- Capper et al. methylation dictionary
- Target gene lists for ROI analysis

---

## System Requirements

### Compute Resources (configurable)
- **CPUs:** Default 16 cores (max_cpus)
- **Memory:** Default 2 GB per process (max_memory)
- **Runtime:** Default 240 hours (max_time)

### Container System
- **Apptainer/Singularity:** Enabled by default
- **Docker:** Optional (via profile)
- **Conda:** Optional (via profile)

---

## Notes

1. **Version Notation:**
   - "latest" indicates the most recent stable version at build time
   - Specific versions are pinned where critical for reproducibility

2. **Container Updates:**
   - Containers are periodically rebuilt with updated dependencies
   - Date stamps in container names indicate last major update (e.g., `28feb2025`)

3. **Reproducibility:**
   - All containers are archived with specific tags
   - Dockerfile definitions available in `/dockerfiles/` directory
   - Conda environment files preserve exact dependency versions

4. **License:**
   - Pipeline: MIT License
   - Individual tools: See respective tool documentation

---

## Citation

If you use this pipeline, please cite:
- nWGS Pipeline: VilhelmMagnusLab/nWGS_pipeline
- Individual tools: See their respective publications

For questions or issues, visit: https://github.com/VilhelmMagnusLab/nWGS_pipeline
