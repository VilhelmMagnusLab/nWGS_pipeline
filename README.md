# DIANA: An integrated pipeline for analysis of long-read whole-genome sequencing data for molecular neuropathology

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.10.1-23aa62.svg)](https://www.nextflow.io/)
[![run with apptainer](https://img.shields.io/badge/run%20with-apptainer-1d355c.svg?labelColor=000000)](https://apptainer.org/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed.svg?labelColor=000000)](https://www.docker.com/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Release](https://img.shields.io/github/v/tag/VilhelmMagnusLab/DIANA?label=release&color=blue)](https://github.com/VilhelmMagnusLab/DIANA/releases)

## Overview

 Diagnostic Integrated Analytics for Neoplastic Alterations (DIANA) is a comprehensive bioinformatics pipeline for analyzing Neoplastic alteration. It integrates multiple analyses including CNV detection, methylation profiling, structural variant calling, and MGMT promoter status determination.

## Pipeline Schematic

DIANA pipeline follows a modular architecture with four main Nextflow modules that can be run independently or sequentially:

<div align="center">

![Diana Pipeline Schematic](nWGS.png)

</div>

*Pipeline workflow showing the flow from BAM files through Mergebam, Epi2me, and Annotation modules to final PDF reports.*

## Quick Start

### Prerequisites
- **Docker** (Desktop/Local) or **Singularity/Apptainer** (HPC)
- **Java 11–21** (auto-installed by setup script if missing; required by Nextflow)
- **Nextflow** (auto-installed by setup script)
- **Internet connection** for downloading reference files from Zenodo

### Automated Setup & Run

The pipeline now features a unified setup script that automatically downloads all reference files from Zenodo:

**For Docker (Desktop/Local):**
```bash
git clone https://github.com/VilhelmMagnusLab/Diana.git
cd Diana
./setup_pipeline.sh docker
./run_pipeline_docker.sh --run_mode_order --sample_id YOUR_SAMPLE_ID
```

**For Singularity/Apptainer (HPC):**
```bash
git clone https://github.com/VilhelmMagnusLab/Diana.git
cd Diana
./setup_pipeline.sh singularity
./run_pipeline_singularity.sh --run_mode_order --sample_id YOUR_SAMPLE_ID
```

**What the setup script does:**
- Checks for compatible Java (11–21) and installs it if missing
- Installs Nextflow and adds it to PATH via `.diana_env`
- Downloads all reference files from Zenodo (DOI: [10.5281/zenodo.20761496](https://doi.org/10.5281/zenodo.20761496))
- Extracts and organizes files into the correct directory structure
- Downloads and sets up Docker containers or Singularity images

**Note:** First-time setup downloads ~14 GB of reference data and may take 10-30 minutes depending on your internet connection.

### Test with Demo Data

A minimal test dataset (`diana_dummy`) is automatically downloaded and extracted into `data/diana_dummy/` by `setup_pipeline.sh`. It contains a single sample (`diana-001`) with a small BAM file and the required `final_summary` trigger file, mirroring the expected input structure.

After setup completes, run the pipeline on the demo data:

```bash
# Docker
bash smart_sample_monitor_v2.sh -d data/diana_dummy

# Singularity/Apptainer
bash smart_sample_monitor_v2.sh --singularity -d data/diana_dummy
```

The monitor will detect the `final_summary` file in `diana_dummy/diana-001/`, trigger the full pipeline, and write results to `~/routine_diana/routine_results/diana-001/`.

> The sample ID files created by `setup_pipeline.sh` already contain `diana-001` / `PBE00000` — no manual configuration needed.

## Updating the Pipeline

To get the latest version of DIANA, run the built-in update script from inside your pipeline directory:

```bash
cd /path/to/Diana
bash update.sh
```

**What it does:**
- Fetches and applies the latest changes from GitHub
- Backs up any modified config files (`conf/annotation.config`, `conf/epi2me.config`, `conf/example.config`) to `conf/backup_<timestamp>/` before updating
- Works for both cloned repositories and ZIP downloads
- After updating, re-apply your custom paths from the backup folder if needed

**For ZIP installs** (downloaded instead of cloned): `update.sh` automatically initialises a git repository and connects it to the upstream remote, so future updates are fast.

> **Note:** `update.sh` only updates pipeline code and configuration. Reference files and containers are not re-downloaded. If a new version requires updated reference files, the release notes will indicate this.

## Pipeline Modules

The pipeline consists of three main modules that can be run independently or sequentially:

### 1. **Mergebam Pipeline** (`--run_mode_mergebam`)
- Merges multiple BAM files per sample
- Extracts protein-coding regions of interest using `roi.protein_coding.bed`

### 2. **Epi2me Pipeline** (`--run_mode_epi2me`)
Four independent analysis types:

| Analysis | Tool | Purpose | Output |
|----------|------|---------|---------|
| **Modified Base Calling** | Modkit | DNA modifications (5mC, 5hmC) | `*_wf_mods.bedmethyl.gz` |
| **Structural Variants** | Sniffles2 | Structural variant detection | `*.sniffles.vcf.gz` |
| **Copy Number Variation** | QDNAseq | CNV detection | `*_segs.bed`, `*_bins.bed`, `*_segs.vcf` |
| **SNV Calling** | Clair3 (germline) + ClairS-TO (somatic) | Small variant calling on protein-coding ROI | `*.clair3.vcf.gz`, `*.clairsto.vcf.gz` |

#### PacBio HiFi Data: BAM Alignment Pre-processing

PacBio HiFi BAM files from the sequencer need to be **unaligned**. Before running the pipeline (specifically before modkit modified base calling), the BAM must be aligned to the reference genome. Use the following command:

```bash
samtools fastq -T MM,ML /path/to/input.hifi_reads.bam \
    | minimap2 -y -ax map-hifi -t 4 /path/to/GRCh38_reference.fa - \
    | samtools sort -@ 4 -o /path/to/output.hifi_reads.aligned.bam

# Then index the aligned BAM
samtools index /path/to/output.hifi_reads.aligned.bam
```

> **Important:** The `-T MM,ML` flag in `samtools fastq` is required to preserve the base modification tags (MM and ML) that encode 5mC methylation information. Without these tags, modkit will not be able to extract methylation calls.

### 3. **Annotation Pipeline** (`--run_mode_annotation`)
- **MGMT methylation analysis** using EPIC array sites
- **NanoDx neural network classification** with dual reference data:
  - **Capper CNS tumor reference data** - Optimized for brain tumors
  - **Pan-cancer reference data** - Broader tumor type coverage
- **Structural variant annotation** with Svanna
- **SNV annotation** of pre-called VCFs from the Epi2me module (ANNOVAR, ClinVar, COSMIC), filtered by configurable Depth and GQ thresholds for the report
- **CNV analysis** with ACE tumor content determination
- **Comprehensive reporting** (HTML, IGV snapshots, Circos plots, Markdown)

#### Methylation Classifier Selection

The pipeline supports two NanoDx methylation classifiers:

| Classifier | Recommended For | Description |
|------------|-----------------|-------------|
| **Capper et al.** | Brain tumors | Optimized for CNS tumor classification |
| **Pan-cancer v5i** | Broader tumor types | Extended classifier covering wider range of tumor types |

Both classifiers run by default. Example usage:
```bash
./run_pipeline_singularity.sh --run_mode_order --sample_id SAMPLE_001
```

## Pipeline Run Modes

The pipeline can be executed in different modes:

| Mode | Flag | Description | Use Case |
|------|------|-------------|----------|
| **Complete Pipeline** | `--run_mode_order` | Runs all three modules sequentially (Mergebam → Epi2me → Annotation) | Starting from raw BAM files |
| **Epi2me + Annotation** | `--run_mode_epiannotation` | Runs Epi2me and Annotation sequentially (assumes merged BAM files exist) | When BAM files are already merged |
| **Mergebam Only** | `--run_mode_mergebam` | Merges BAM files and extracts regions of interest | BAM preparation only |
| **Epi2me Only** | `--run_mode_epi2me [all\|modkit\|cnv\|sv\|snv]` | Runs specific Epi2me analyses | Methylation, CNV, SV, or SNV calling |
| **Annotation Only** | `--run_mode_annotation [all\|mgmt\|cnv\|svannasv\|terp\|snv\|rmd]` | Runs specific downstream analyses | Report generation or specific analyses |

## Container Systems

| Feature | Docker | Singularity/Apptainer |
|---------|--------|----------------------|
| **Best for** | Desktop/Local | HPC/Shared systems |
| **Setup Script** | `setup_docker.sh` | `setup_singularity.sh` |
| **Run Script** | `run_pipeline_docker.sh` | `run_pipeline_singularity.sh` |

All containers are automatically downloaded from [vilhelmmagnuslab Docker Hub](https://hub.docker.com/repositories/vilhelmmagnuslab).

> **Architecture:** All container images are built for **linux/amd64 (x86_64)**. ARM64 systems (e.g. Apple Silicon) are not currently validated — see [CONTAINERS.md](CONTAINERS.md) for details on emulation options.

## Usage Examples

### Complete Pipeline (Recommended)
```bash
# Docker - Full pipeline starting from raw BAM files
./run_pipeline_docker.sh --run_mode_order --sample_id T001

# Singularity/Apptainer - Full pipeline starting from raw BAM files
./run_pipeline_singularity.sh --run_mode_order --sample_id T001
```

### Epi2me + Annotation (When BAM files are already merged)
```bash
# Docker - Skip mergebam, run Epi2me and Annotation
./run_pipeline_docker.sh --run_mode_epiannotation --sample_id T001

# Singularity/Apptainer - Skip mergebam, run Epi2me and Annotation
./run_pipeline_singularity.sh --run_mode_epiannotation --sample_id T001
```

### Individual Modules

**Docker Commands:**
```bash
# Mergebam only
./run_pipeline_docker.sh --run_mode_mergebam

# Epi2me analyses
./run_pipeline_docker.sh --run_mode_epi2me all          # All Epi2me analyses
./run_pipeline_docker.sh --run_mode_epi2me stat         # QC statistics (cramino) only
./run_pipeline_docker.sh --run_mode_epi2me modkit       # Modified base calling only
./run_pipeline_docker.sh --run_mode_epi2me cnv          # CNV analysis only
./run_pipeline_docker.sh --run_mode_epi2me sv           # Structural variants only
./run_pipeline_docker.sh --run_mode_epi2me snv          # SNV calling (Clair3 + ClairS-TO) only

# Annotation modules
./run_pipeline_docker.sh --run_mode_annotation all        # All analyses
./run_pipeline_docker.sh --run_mode_annotation mgmt       # MGMT analysis only
./run_pipeline_docker.sh --run_mode_annotation cnv        # CNV analysis only
./run_pipeline_docker.sh --run_mode_annotation svannasv   # Svanna SV annotation only
./run_pipeline_docker.sh --run_mode_annotation terp       # TERTp promoter analysis only
./run_pipeline_docker.sh --run_mode_annotation snv        # SNV annotation (Clair3 + ClairS-TO) only
./run_pipeline_docker.sh --run_mode_annotation rmd        # Markdown report only
```

**Singularity/Apptainer Commands:**
```bash
# Mergebam only
./run_pipeline_singularity.sh --run_mode_mergebam

# Epi2me analyses
./run_pipeline_singularity.sh --run_mode_epi2me all          # All Epi2me analyses
./run_pipeline_singularity.sh --run_mode_epi2me stat         # QC statistics (cramino) only
./run_pipeline_singularity.sh --run_mode_epi2me modkit       # Modified base calling only
./run_pipeline_singularity.sh --run_mode_epi2me cnv          # CNV analysis only
./run_pipeline_singularity.sh --run_mode_epi2me sv           # Structural variants only
./run_pipeline_singularity.sh --run_mode_epi2me snv          # SNV calling (Clair3 + ClairS-TO) only

# Annotation modules
./run_pipeline_singularity.sh --run_mode_annotation all        # All analyses
./run_pipeline_singularity.sh --run_mode_annotation mgmt       # MGMT analysis only
./run_pipeline_singularity.sh --run_mode_annotation cnv        # CNV analysis only
./run_pipeline_singularity.sh --run_mode_annotation svannasv   # Svanna SV annotation only
./run_pipeline_singularity.sh --run_mode_annotation terp       # TERT promoter analysis only
./run_pipeline_singularity.sh --run_mode_annotation snv        # SNV annotation (Clair3 + ClairS-TO) only
./run_pipeline_singularity.sh --run_mode_annotation rmd        # Markdown report only
```

## Input Requirements

### Sample ID File Format
```
# For annotation pipeline (with tumor content)
sample_id1   0.75    # 75% tumor content
sample_id2          # Auto-calculate with ACE

# For mergebam pipeline (with flowcell)
sample_id1   flowcell_id1
sample_id2   flowcell_id2
```

### Directory Structure

The pipeline uses a standardized directory structure with separate input and output paths:

```
Pipeline directory:
/data/routine_diana/Diana/
├── conf/                         # Configuration files
│   ├── mergebam.config          # Mergebam module config
│   ├── epi2me.config            # Epi2me module config
│   └── annotation.config        # Annotation module config
├── modules/                      # Nextflow modules
├── containers/                   # Singularity container images
├── bin/                         # Helper scripts
├── docs/                        # Documentation
└── smart_sample_monitor_v2.sh  # Automated monitoring script

Pipeline data directory (configured via params.path):
/data/
├── reference/                    # Reference files (GRCh38, BED files, etc.)
└── humandb/                      # Annotation databases

Input data directory (configured via params.input_dir in mergebam.config):
/data/WGS_[DATE]/                # Oxford Nanopore sequencing output
├── SAMPLE_01/                    # Sample directory
│   └── [subdirectory]/          # Any subdirectory structure
│       ├── *.bam                # BAM files from ONT sequencing
│       ├── *.bam.bai            # BAM index files
│       └── final_summary_*_*_*.txt  # Completion marker file
├── SAMPLE_02/
│   └── [subdirectory]/
│       ├── *.bam
│       ├── *.bam.bai
│       └── final_summary_*_*_*.txt
└── ...

Output directory (configured via params.path_output):
routine_diana/
├── sample_ids_bam.txt           # Sample IDs for BAM merging
│
├── routine_bams/                # Processed BAM files (Mergebam module)
│   ├── merge_bams/              # Merged BAM files per sample
│   └── roi_bams/                # Region of interest extracted BAMs
│
├── routine_epi2me/              # Epi2me module results
│   └── [sample_id]/
│       ├── *.wf_mods.bedmethyl.gz     # Methylation calls (modkit)
│       ├── *.sniffles.vcf.gz          # Structural variants (Sniffles2)
│       ├── *_segs.bed                 # CNV segments (QDNAseq)
│       ├── *_bins.bed                 # CNV bins
│       ├── *_copyNumbersCalled.rds    # CNV RDS file for ACE
│       ├── clair3/                    # Germline SNV calling (Clair3)
│       │   └── *.vcf.gz
│       └── clairs-to/                 # Somatic SNV calling (ClairS-TO)
│           └── *.vcf.gz
│
├── routine_annotation/            # Analysis module results (detailed outputs)
│   └── [sample_id]/
│       ├── classifier/          # Tumor classification
│       │   ├── nanodx/         # NanoDx neural network results
│       │   └── sturgeon/       # Sturgeon methylation classifier
│       ├── cnv/                 # CNV analysis
│       │   ├── ace/            # ACE tumor content estimation
│       │   ├── annotatedcnv/   # Annotated CNV calls
│       │   └── *.pdf           # CNV plots (chr7, chr9, full genome)
│       ├── coverage/            # IGV coverage snapshots
│       │   ├── *_egfr_coverage.pdf
│       │   ├── *_idh1_coverage.pdf
│       │   ├── *_idh2_coverage.pdf
│       │   └── *_tertp_coverage.pdf
│       ├── cramino/             # BAM statistics
│       │   └── *_cramino_statistics.txt
│       ├── merge_annot_clair3andclairsto/  # Variant annotation
│       │   └── *_merge_annotation_filter_snvs_allcall.csv
│       ├── methylation/         # MGMT methylation analysis
│       │   └── *_MGMT_results.csv
│       └── structure_variant/   # SV annotation
│           ├── *_circos.pdf    # Circos plot
│           ├── *_fusion_events.tsv  # Fusion events
│           └── *_svanna_annotation.html  # Svanna SV annotation
│
└── routine_results/             # Final published reports (per sample)
    └── [sample_id]/
        ├── [sample_id]_bedmethyl_sturgeon_general.pdf  # Sturgeon classification
        ├── [sample_id]_markdown_pipeline_report.pdf    # Main comprehensive report
        ├── [sample_id]_mnpflex_input.bed               # MNP-Flex input format
        ├── [sample_id]_occ_svanna_annotation.html      # SV annotation HTML
        └── [sample_id]_tsne_plot.html                  # t-SNE visualization
```

## Required Reference Data

### Automated Download (Recommended)

**The `setup_pipeline.sh` script automatically downloads and sets up all required reference files from Zenodo.**

Simply run:
```bash
./setup_pipeline.sh docker    # For Docker users
# or
./setup_pipeline.sh singularity    # For Singularity users
```

The script will:
1. Download reference data from [Zenodo (DOI: 10.5281/zenodo.20761496)](https://doi.org/10.5281/zenodo.20761496)
2. Extract and organize all files into the correct directory structure
3. Set up NanoDx classifier models
4. Configure all required paths

### Manual Setup (Advanced Users Only)

If you prefer manual setup or need to customize the reference files:

**Core reference files** (automatically placed in `data/reference/`):
- `reference_core.tar.gz` - Contains GRCh38 reference genome, BED files, and annotations including:
  - `GRCh38.fa` and `GRCh38.fa.fai` - Human reference genome
  - `EPIC_sites_NEW.bed` - Methylation sites
  - `MGMT_CpG_Island.hg38.bed` - MGMT CpG islands
  - `roi.protein_coding.bed` - Region of interest BED file (protein-coding genes for SNV screening and BAM extraction)
  - `TERTp_variants.bed` - TERT promoter variants
  - `human_GRCh38_trf.bed` - Tandem repeat regions
  - `CNV_genes_tuned.csv` - CNV gene annotations
  - `roi_fusions_genes.txt` - User-defined region of interest gene list for SV/fusion filtering and SNV annotation (one gene per line; can be replaced with any custom gene list)
  - `nanoDx/` - NanoDx neural network classifier (with models from Zenodo)

### SNV Annotation Tool Selection

DIANA supports two SNV annotation backends. The choice is set in `nextflow.config`:

```groovy
params {
    snv_annotator = "annovar"  // default — ANNOVAR
    // snv_annotator = "vep"  // alternative — Ensembl VEP
}
```

Or override on the command line:
```bash
./run_pipeline_singularity.sh --run_mode_order --snv_annotator vep
```

| Annotator | Default | Requirements |
|-----------|---------|--------------|
| **ANNOVAR** | ✓ Yes | ANNOVAR Perl scripts (bundled, subject to ANNOVAR licence) + ANNOVAR databases |
| **VEP** | No | VEP cache + ClinVar VCF + COSMIC VCF (see below) |

**Annotation databases** (automatically placed in `data/humandb/`):
- `humandb.tar.gz` - Contains ANNOVAR-format annotation databases (used by both annotators):
  - `hg38_refGene.txt` - RefGene annotation
  - `hg38_refGeneMrna.fa` - RefGene mRNA sequences
  - `hg38_clinvar_20240611.txt` - ClinVar annotations (ANNOVAR format)
  - `hg38_cosmic100coding2024.txt` - **Placeholder only** — contains empty COSMIC IDs (see COSMIC section below)

### ANNOVAR Scripts (bundled, used by default)

The following ANNOVAR Perl scripts are **bundled** with DIANA in the `bin/` directory:
- `annotate_variation.pl`
- `coding_change.pl`
- `convert2annovar.pl`
- `table_annovar.pl`
- `index_annovar.pl`
- `prepare_annovar_user.pl`

> **License notice:** Usage of ANNOVAR is subject to its own licence terms. See [THIRD_PARTY_LICENSES.md](THIRD_PARTY_LICENSES.md) for details and the required citation.

### COSMIC Database (required, not included)

COSMIC data requires a free institutional registration and **cannot be redistributed**. The `humandb.tar.gz` archive includes a placeholder `hg38_cosmic100coding2024.txt` with empty COSMIC IDs so the pipeline runs without error, but no real COSMIC annotations will appear in reports until you install the full database.

To prepare the COSMIC database (following [ANNOVAR COSMIC instructions](https://annovar.openbioinformatics.org/en/latest/user-guide/filter/#cosmic-annotations)):

```bash
# 1. Register and log in at https://cancer.sanger.ac.uk/cosmic/login
#    then download COSMIC v{version} files from https://cancer.sanger.ac.uk/cosmic/download
#    You need: Cosmic_GenomeScreensMutant and Cosmic_NonCodingVariants (VCF + TSV, GRCh38)
#    Replace {version} with the current COSMIC version (e.g. 100)

# 2. Extract and decompress (replace {version} with your downloaded version)
tar xvf Cosmic_GenomeScreensMutant_Vcf_v{version}_GRCh38.tar
tar xvf Cosmic_GenomeScreensMutant_Tsv_v{version}_GRCh38.tar
gunzip Cosmic_GenomeScreensMutant_v{version}_GRCh38.vcf.gz
gunzip Cosmic_GenomeScreensMutant_v{version}_GRCh38.tsv.gz
tar xvf Cosmic_NonCodingVariants_Tsv_v{version}_GRCh38.tar
tar xvf Cosmic_NonCodingVariants_Vcf_v{version}_GRCh38.tar
gunzip Cosmic_NonCodingVariants_v{version}_GRCh38.vcf.gz
gunzip Cosmic_NonCodingVariants_v{version}_GRCh38.tsv.gz

# 3. Build the ANNOVAR-formatted COSMIC file
echo -e '#Chr\tStart\tEnd\tRef\tAlt\tCOSMIC100' > hg38_cosmic100_raw.txt
prepare_annovar_user.pl -dbtype cosmic \
    Cosmic_GenomeScreensMutant_v{version}_GRCh38.tsv \
    -vcf Cosmic_GenomeScreensMutant_v{version}_GRCh38.vcf >> hg38_cosmic100_raw.txt
prepare_annovar_user.pl -dbtype cosmic \
    Cosmic_NonCodingVariants_v{version}_GRCh38.tsv \
    -vcf Cosmic_NonCodingVariants_v{version}_GRCh38.vcf >> hg38_cosmic100_raw.txt
index_annovar.pl hg38_cosmic100_raw.txt -outfile hg38_cosmic100coding2024.txt

# 4. Copy to DIANA humandb directory (replaces the placeholder)
cp hg38_cosmic100coding2024.txt     /path/to/Diana/data/humandb/
cp hg38_cosmic100coding2024.txt.idx /path/to/Diana/data/humandb/
```

> **Note:** `prepare_annovar_user.pl` and `index_annovar.pl` are part of the ANNOVAR package downloaded in the step above.

### Ensembl VEP Setup (required only when snv_annotator = "vep")

VEP requires three components in `data/humandb/`. The **VEP cache** and **ClinVar VCF** are included in the Zenodo reference archive and extracted automatically by `setup_pipeline.sh`. Only the **COSMIC VCF** must be prepared manually.

#### 1. VEP cache (`homo_sapiens_refseq/`)

Distributed via Zenodo and automatically placed at `data/humandb/homo_sapiens_refseq/` during setup.

To download manually (e.g. to use a newer cache version):
```bash
cd /path/to/Diana/data/humandb
curl -O https://ftp.ensembl.org/pub/release-115/variation/indexed_vep_cache/homo_sapiens_refseq_vep_115_GRCh38.tar.gz
tar xzf homo_sapiens_refseq_vep_115_GRCh38.tar.gz
# Creates: homo_sapiens_refseq/115_GRCh38/
```

Update `conf/annotation.config` if using a different version:
```groovy
vep_cache_version = "115"
```

#### 2. ClinVar VCF (`clinvar.vcf.gz`)

Distributed via Zenodo and automatically placed in `data/humandb/` during setup.

To generate your own from the ANNOVAR ClinVar file (e.g. after a ClinVar update):
```bash
cd /path/to/Diana/data/humandb
python3 /path/to/Diana/bin/annovar_to_vep_vcf.py \
    clinvar hg38_clinvar_20240611.txt clinvar.vcf.gz
# Creates: clinvar.vcf.gz and clinvar.vcf.gz.tbi
```

#### 3. COSMIC VCF (`CosmicCodingMuts.vcf.gz`) — must be prepared manually

After preparing the ANNOVAR COSMIC file (see COSMIC Database section above), convert it to bgzipped VCF for VEP using the provided script:

```bash
cd /path/to/Diana/data/humandb
python3 /path/to/Diana/bin/annovar_to_vep_vcf.py \
    cosmic hg38_cosmic100coding2024.txt CosmicCodingMuts.vcf.gz
# Creates: CosmicCodingMuts.vcf.gz and CosmicCodingMuts.vcf.gz.tbi
```

The final `data/humandb/` structure required for VEP mode:

```
data/humandb/
├── homo_sapiens_refseq/           # VEP cache (Zenodo or manual)
│   └── 115_GRCh38/
├── clinvar.vcf.gz                 # ClinVar VCF for VEP (Zenodo or generated)
├── clinvar.vcf.gz.tbi
├── CosmicCodingMuts.vcf.gz        # COSMIC VCF for VEP (must be generated)
├── CosmicCodingMuts.vcf.gz.tbi
├── hg38_refGene.txt               # ANNOVAR databases (also used by ANNOVAR mode)
├── hg38_clinvar_20240611.txt
└── hg38_cosmic100coding2024.txt
```

**Additional reference files** (automatically extracted to `data/reference/`):
- `general.zip` - Sturgeon classifier model (kept as zip, not extracted) — **see Sturgeon note below**
- `Assembly.zip` - Assembly folder for vcfcircos visualization (automatically extracted)
- `svanna-data.zip` - Svanna structural variant annotation database (optional, automatically extracted)

### Sturgeon Classifier Model (`general.zip`)

The Sturgeon methylation classifier model (`general.zip`) must be present in `data/reference/` for Sturgeon-based classification to run. **If the file is absent the pipeline continues without Sturgeon** — all other analyses (NanoDx, MGMT, SNV, CNV, report) are unaffected.

If you need to download it separately:

> **Download:** [https://www.dropbox.com/s/yzca4exl40x9ukw/general.zip?dl=0](https://www.dropbox.com/s/yzca4exl40x9ukw/general.zip?dl=0)

```bash
# Copy to the reference folder
cp general.zip /path/to/Diana/data/reference/general.zip
```

> The file must remain as a **zip archive** — do not extract it. Sturgeon reads the model directly from the zip.

### Clair3 Basecalling Model

Two model types are provided for R10.4.1 chemistry (400 bps):

| Model folder | Basecalling type | When to use |
|---|---|---|
| `r1041_e82_400bps_sup_v520` | **SUP** (Super Accuracy) | Default — highest SNV calling accuracy, slower basecalling |
| `r1041_e82_400bps_hac_v520` | **HAC** (High Accuracy) | Faster basecalling, slightly lower accuracy |

The default is **SUP**. To switch to HAC, pass `--clair3_model hac` on the command line:

```bash
# Use HAC model
./run_pipeline_singularity.sh --run_mode_order --clair3_model hac

# Use SUP model (default — no flag needed)
./run_pipeline_singularity.sh --run_mode_order
```

You can also change the permanent default in `nextflow.config`:
```groovy
params {
    clair3_model = "sup"   // change to "hac" to make HAC the default
}
```

> **Note:** Both model folders (`r1041_e82_400bps_sup_v520/` and `r1041_e82_400bps_hac_v520/`) must be present in `data/reference/` for the respective model to work (automatically extracted). The models use **PyTorch format** (`.pt` files) and require the PyTorch-based Clair3 container.

**Using a different Clair3 model:** Additional pre-trained models for other chemistries, platforms, or Dorado versions are available at [https://github.com/HKU-BAL/Clair3](https://github.com/HKU-BAL/Clair3). Download the model files, place them in a folder under `data/reference/`, and pass the folder path directly via `--clair3_model_path`:

```bash
./run_pipeline_singularity.sh --run_mode_order \
    --clair3_model_path /path/to/Diana/data/reference/your_custom_model/
```

**Note on roi.protein_coding.bed:** This ROI BED file uses OCC (Onco-Comprehensive-Coverage) genes but can be substituted with any custom ROI BED file. It's used for:
- Extracting regions of interest during BAM merging (mergebam module)
- SNV screening regions for variant calling (ClairS-TO analysis)
- Ensure proper BED format with exactly 10 tab-separated fields per line

**Note on roi_fusions_genes.txt:** Plain-text gene list (one gene symbol per line) used for SV/fusion event filtering and SNV annotation. This file can be replaced with any user custom gene list of interest — for example, a laboratory-specific panel of oncology-relevant genes. The default list contains 204 genes covering common fusion partners and oncogenes.

**Manual download:** If needed, all reference files are available at [Zenodo (DOI: 10.5281/zenodo.20761496)](https://doi.org/10.5281/zenodo.20761496)

### Directory Structure Setup
After downloading the reference files, your directory structure should look like this:

```
data/
├── reference/                    # Reference files
│   ├── GRCh38.fa
│   ├── GRCh38.fa.fai
│   ├── gencode.v48.annotation.gff3
│   ├── Assembly/                # Assembly folder for vcfcircos (from Zenodo)
│   ├── EPIC_sites_NEW.bed
│   ├── MGMT_CpG_Island.hg38.bed
│   ├── roi.protein_coding.bed
│   ├── TERTp_variants.bed
│   ├── human_GRCh38_trf.bed
│   ├── CNV_genes_tuned.csv
│   ├── roi_fusions_genes.txt
│   └── etc
│
└── humandb/                     # Annotation databases
    ├── hg38_refGene.txt
    ├── hg38_refGeneMrna.fa
    ├── hg38_clinvar_20240611.txt
    └── hg38_cosmic100coding2024.txt
```

## ACE Tumor Content Calculation

The pipeline intelligently handles tumor content:
- **Provided value**: Use directly if specified in sample ID file
- **Auto-calculation**: ACE analyzes copy number profiles to estimate tumor cellularity
- **Multiple estimates**: ACE provides several estimates and selects the best fit
- **Results**: Saved in `${sample_id}_ace_results/threshold_value.txt`

## Report Generation

### Standard Report Generation

**PDF reports are automatically generated** when running the pipeline with the following modes:
- `--run_mode_annotation rmd` - Generate reports only
- `--run_mode_order` - Run complete pipeline sequentially and generate reports
- `--run_mode_epiannotation` - Run Epi2me and annotation modules and generate reports

The reports are automatically created in the `routine_results/{sample_id}/` directory with the name `{sample_id}_markdown_pipeline_report.pdf`.

### Additional Report Generation

The `generate_report.sh` script is provided for **additional report generation** in cases where:
- You want to regenerate reports after re-running specific processes
- You need to create reports for samples that were processed separately
- You need to generate reports after the pipeline has already completed


## Configuration

### Path Configuration

The pipeline uses three main path parameters that must be configured:

**1. Pipeline Data Path (`params.path`)** - Reference files and databases
```groovy
// conf/annotation.config, conf/epi2me.config, conf/mergebam.config
params {
    path = "/data/routine_diana/Diana/data"
    // Contains: reference/, humandb/ directories
}
```

**2. Input Data Path (`params.input_dir`)** - ONT sequencing output
```groovy
// conf/mergebam.config
params {
    input_dir = "/data/WGS_27102025"
    // Contains: Sample directories with BAM files
    // Can be overridden via CLI: --input_dir or smart_sample_monitor -d
}
```

**3. Output Path (`params.path_output`)** - Pipeline results
```groovy
// conf/mergebam.config, conf/epi2me.config, conf/annotation.config
params {
    path_output = "/data/routine_diana"
    // Contains: sample_ids_bam.txt, routine_bams/, routine_epi2me/, routine_results/
}
```

**Key Points:**
- `params.path`: Reference data (rarely changes)
- `params.input_dir`: ONT sequencing input (changes per run)
- `params.path_output`: Where all results are stored (consistent location)
- The `input_dir` can be overridden using `--input_dir` flag or `smart_sample_monitor_v2.sh -d`

### SNV Filtering Configuration

The pipeline includes configurable quality thresholds for SNV filtering in the final reports:

```groovy
// conf/annotation.config
params {
    snv_depth_threshold = 10    // Minimum sequencing depth (default: 10)
    snv_gq_threshold = 10       // Minimum Genotype Quality (default: 10)
}
```

**How Filtering Works:**
- **Depth threshold**: Filters out variants with sequencing depth below the threshold
- **GQ threshold**: For variants with multiple GQ values from different callers (e.g., "20,26,41"), keeps the variant if ANY value meets the threshold
- Both filters must pass for a variant to appear in the final report

**Examples:**
```groovy
# Stricter filtering (higher quality variants only)
snv_depth_threshold = 15
snv_gq_threshold = 20

# More permissive filtering (include more variants)
snv_depth_threshold = 5
snv_gq_threshold = 5
```

**Note:** These thresholds only affect the variants shown in the Markdown PDF reports. The raw VCF files contain all called variants regardless of these filters.

### Container Configuration
Choose your preferred container engine and run the unified setup script:

```bash
# For Docker
./setup_pipeline.sh docker

# For Singularity/Apptainer
./setup_pipeline.sh singularity
```

The setup script handles Java, Nextflow, reference files, and container images in one step.

### Work Directory Customization
You can specify a custom temporary work directory using the `-w` flag. This is useful for:
- Managing disk space on different storage locations
- Avoiding permission issues
- Organizing temporary files

**Example:**
```bash
# Docker
./run_pipeline_docker.sh --run_mode_annotation tertp -w /path/to/your/work/dir

# Singularity/Apptainer  
./run_pipeline_singularity.sh --run_mode_annotation tertp -w /home/chbope/extension/trash/tmp
```

**Note:** The `-w` flag sets Nextflow's work directory where temporary files and intermediate results are stored during pipeline execution. By default nextflow create a folder `work` in the working directory.

## Automated Sample Monitoring

The pipeline includes `smart_sample_monitor_v2.sh` for **automated monitoring and processing** of Oxford Nanopore sequencing runs. This intelligent script continuously monitors sample directories and automatically triggers the pipeline when sequencing completes.

### Key Features:

**Monitoring & Execution:**
- **Real-time Monitoring**: Watches for `final_summary_*_*_*.txt` files indicating completed sequencing
- **Automatic Pipeline Triggering**: Starts processing immediately when samples are ready
- **Sequential Processing**: Processes one sample at a time, queuing others
- **Markdown Report Validation**: Verifies successful completion before marking as done

**Version 2 Enhancements:**
- **CLI Data Directory Override**: `--data-dir` takes precedence over `mergebam.config`
- **Resume Control**: Disabled by default for fresh runs; use `-r` to enable caching
- **Symlink Resolution**: Works correctly when installed as global command
- **Portable Execution**: Automatically finds pipeline directory from any location
- **Sample IDs File**: Hardcoded to `/data/routine_diana/sample_ids_bam.txt`

### Basic Usage:

```bash
# Run from pipeline directory with default config (auto-detects Singularity or Docker)
./smart_sample_monitor_v2.sh

# Monitor specific data directory (overrides config)
./smart_sample_monitor_v2.sh -d /data/WGS_27102025

# Enable resume for cached results
./smart_sample_monitor_v2.sh -d /data/WGS_27102025 -r

# Verbose logging
./smart_sample_monitor_v2.sh -d /data/WGS_27102025 -v

# Combination: resume + verbose
./smart_sample_monitor_v2.sh -d /data/WGS_27102025 -r -v

# Force Docker (useful when both Docker and Singularity are available)
./smart_sample_monitor_v2.sh --docker -d /data/WGS_27102025

# Force Singularity/Apptainer
./smart_sample_monitor_v2.sh --singularity -d /data/WGS_27102025

# Explicit engine flag (equivalent to --docker / --singularity)
./smart_sample_monitor_v2.sh -e docker -d /data/WGS_27102025 -r -v
./smart_sample_monitor_v2.sh -e singularity -d /data/WGS_27102025 -r -v
```

### Global Command Installation:

Install the monitor as a global command accessible from any directory:

**User-level installation (Recommended - No sudo required):**
```bash
# Create user bin directory and symbolic link
mkdir -p ~/bin
ln -sf /data/routine_diana/Diana/smart_sample_monitor_v2.sh ~/bin/smart_sample_monitor

# Add ~/bin to PATH (run once)
cat >> ~/.bashrc << 'EOF'

# Add user's bin directory to PATH
if [ -d "$HOME/bin" ]; then
    export PATH="$HOME/bin:$PATH"
fi
EOF

# Activate changes
source ~/.bashrc

# Verify installation
which smart_sample_monitor
```

**System-wide installation (Requires sudo):**
```bash
sudo ln -sf /data/routine_diana/Diana/smart_sample_monitor_v2.sh /usr/local/bin/smart_sample_monitor
```

**Then use from anywhere:**
```bash
# Run from any directory
cd /tmp
smart_sample_monitor -d /data/WGS_27102025 -v

# Monitor with custom work directory
smart_sample_monitor -d /data/WGS_27102025 -w /data/trash -r

# Force Docker from anywhere
smart_sample_monitor --docker -d /data/WGS_27102025 -v
```

### Command-Line Options:

| Option | Long Form | Description | Default |
|--------|-----------|-------------|---------|
| `-d` | `--data-dir` | Base data directory (overrides config) | Auto-detect from config |
| `-p` | `--pipeline` | Pipeline base directory | Auto-detected |
| `-w` | `--workdir` | Nextflow work directory | `/data/trash` |
| `-c` | `--config` | Config file to parse | `conf/mergebam.config` |
| `-i` | `--interval` | Check interval in seconds | 300 (5 min) |
| `-t` | `--timeout` | Maximum wait time in seconds | 432000 (5 days) |
| `-e` | `--engine` | Container engine: `singularity`, `apptainer`, or `docker` | Auto-detect |
| | `--docker` | Shorthand for `--engine docker` | - |
| | `--singularity` | Shorthand for `--engine singularity` | - |
| `-r` | `--resume` | Enable Nextflow resume | Disabled |
| `-v` | `--verbose` | Enable verbose logging | Disabled |
| `-h` | `--help` | Show help message | - |

### Workflow:

1. **Initialize**: Load sample IDs from `/data/routine_diana/sample_ids_bam.txt`
2. **Monitor**: Check each sample directory for `final_summary_*_*_*.txt`
3. **Queue**: Mark ready samples for processing
4. **Execute**: Run `--run_mode_order` for each sample sequentially
5. **Validate**: Check for markdown report generation
6. **Report**: Display final status summary

### Use Case:

This script is essential for **routine ONT sequencing workflows** where:
- Multiple samples complete sequencing at different times
- Immediate processing is desired upon completion
- Manual monitoring would be time-consuming and error-prone
- Consistent processing workflow is required

Instead of manually checking and starting the pipeline for each sample, the monitor **automatically detects completion** and starts processing immediately, **maximizing throughput** and **reducing manual intervention**.

**Important:** Ensure all paths are correctly configured in `conf/mergebam.config`:
- `params.path`: Reference data directory
- `params.input_dir`: Default input directory (can be overridden with `-d`)
- `params.path_output`: Output results directory

**See [docs/GLOBAL_COMMAND_SETUP.md](docs/GLOBAL_COMMAND_SETUP.md) for detailed installation, troubleshooting, and advanced usage.** 

## Troubleshooting

### Common Issues
1. **Container engine conflict**: Ensure only one container system is enabled
2. **Missing reference files**: Download required external files
3. **Permission issues**: Check container and file permissions

### Verification Commands
```bash
# Check containers
docker images | grep vilhelmmagnuslab          # Docker
ls -la containers/*.sif                        # Singularity

# Test pipeline
./test_pipeline_docker.sh                            # Docker
./test_pipeline_singularity.sh                # Singularity
```

## Support

- **Documentation**:
  - [DOCKER_SETUP.md](DOCKER_SETUP.md) - Docker installation and setup
  - [SINGULARITY_SETUP.md](SINGULARITY_SETUP.md) - Singularity/Apptainer setup
  - [docs/GLOBAL_COMMAND_SETUP.md](docs/GLOBAL_COMMAND_SETUP.md) - Global command installation
- **Issues**: [GitHub Issues](https://github.com/VilhelmMagnusLab/Diana/issues)
- **Contact**: 
  - Christian Domilongo Bope (chbope@ous-hf.no / christianbope@gmail.com)
  - Skarphedinn Halldorsson (skahal@ous-hf.no / skabbi@gmail.com)
  - Richard Nagymihaly (ricnag@ous-hf.no)

## Citation

If you use this pipeline in your research, please cite:

> Bope CD, Nagymihaly R, Halldorsson S, et al. *DIANA: Diagnostic Integrated Analytics for Neoplastic Alterations a long-read whole genome sequencing pipeline for molecular neuropathology.* 2026. https://doi.org/10.64898/2026.03.25.714119

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Disclaimer

Diagnostic Integrated Analytics for Neoplastic Alterations pipeline (DIANA) is an investigational research tool that has not undergone full clinical validation. Any clinical use or interpretation of its results is entirely at the discretion and responsibility of the treating physician