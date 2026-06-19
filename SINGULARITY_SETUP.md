# nWGS Pipeline - Singularity/Apptainer Setup Guide

This guide explains how to run the nWGS pipeline using Singularity/Apptainer containers, making it easy to get started in HPC environments or when Docker is not available.

## 🐳 Quick Start with Singularity/Apptainer

The easiest way to run the nWGS pipeline is using Singularity/Apptainer containers. All required images are hosted at [https://hub.docker.com/repositories/vilhelmmagnuslab](https://hub.docker.com/repositories/vilhelmmagnuslab) and will be automatically converted to Singularity format.

### Prerequisites

1. **Singularity or Apptainer**: Install on your system
   - [Apptainer](https://apptainer.org/docs/admin/main/installation.html) (recommended, newer)
   - [Singularity](https://sylabs.io/guides/latest/admin-guide/installation.html) (legacy)

2. **Nextflow**: Will be auto-installed if missing

### One-Command Setup

```bash
chmod +x setup_singularity.sh
./setup_singularity.sh
```

This will:
- Check Singularity/Apptainer installation
- Create necessary directories  
- Pull all required Singularity/Apptainer images
- Create convenient run scripts

### Run the Pipeline

```bash
./run_pipeline_singularity.sh
```

### For More Details
See this guide for comprehensive Singularity/Apptainer setup instructions.

## Manual Setup (Alternative)

If you prefer manual setup:

### 1. Pull Singularity/Apptainer Images

```bash
# Create containers directory
mkdir -p containers/

# Core analysis images
singularity pull --dir containers/ docker://vilhelmmagnuslab/nwgs_default_images:latest
singularity pull --dir containers/ docker://vilhelmmagnuslab/ace_1.24.0:latest
singularity pull --dir containers/ docker://vilhelmmagnuslab/annotcnv_images_27feb1025:latest
singularity pull --dir containers/ docker://vilhelmmagnuslab/clair3_amd64
singularity pull --dir containers/ docker://vilhelmmagnuslab/clairsto_amd64:latest
singularity pull --dir containers/ docker://vilhelmmagnuslab/igv_report_amd64:latest
singularity pull --dir containers/ docker://vilhelmmagnuslab/vcf2circos:latest
singularity pull --dir containers/ docker://vilhelmmagnuslab/nanodx_images_3feb25:latest
singularity pull --dir containers/ docker://vilhelmmagnuslab/markdown_images_28feb2025:latest
singularity pull --dir containers/ docker://vilhelmmagnuslab/mgmt_nanopipe_amd64_18feb2025_cramoni:latest
singularity pull --dir containers/ docker://vilhelmmagnuslab/gviz_amd64ps:latest

# Epi2me images
singularity pull --dir containers/ docker://vilhelmmagnuslab/snifflesv252_update:latest
singularity pull --dir containers/ docker://vilhelmmagnuslab/qdnaseq_amd64:latest
singularity pull --dir containers/ docker://vilhelmmagnuslab/modkit:latest
```

### 2. Create Directory Structure

```bash
mkdir -p data/{reference,humandb,testdata,results}
```

### 3. Run with Nextflow

```bash
nextflow run main.nf -c conf/analysis.config -with-singularity
```

## Configuration

### Singularity/Apptainer Configuration

The pipeline is configured to use Singularity/Apptainer containers in both configuration files:

- `conf/analysis.config` - Main analysis configuration
- `conf/epi2me.config` - Epi2me-specific configuration

Key Singularity/Apptainer settings:
```groovy
apptainer {
    enabled = true
    autoMounts = true
    runOptions = '--bind ${projectDir}/data:${projectDir}/data --bind ${projectDir}/containers:${projectDir}/containers'
}
```

### Available Singularity/Apptainer Images

| Process | Singularity Image | Description |
|---------|------------------|-------------|
| Default | `nwgs_default_images_latest.sif` | General analysis tools |
| ACE TMC | `ace_1.24.0_latest.sif` | ACE copy number analysis |
| AnnotateCNV | `annotcnv_images_27feb1025_latest.sif` | CNV annotation |
| ClairS-TO | `clairs-to_latest.sif` | somatic small variant calling |
| Clair3 | `clair3_latest.sif` | Variant calling |
| IGV Tools | `igv_report_amd64_latest.sif` | IGV report generation |
| Circos Plot | `vcf2circos_latest.sif` | Circos visualization |
| NanoDx | `nanodx_images_3feb25_latest.sif` | NanoDx classification |
| Markdown | `markdown_images_28feb2025_latest.sif` | Report generation |
| Cramino | `mgmt_nanopipe_amd64_18feb2025_cramoni_latest.sif` | Quality assessment |
| Gviz | `gviz_amd64ps_latest.sif` | Genomic visualization |
| Sniffles2 | `snifflesv252_update_latest_latest.sif` | SV calling (Epi2me) |
| QDNAseq | `qdnaseq_amd64_latest_latest.sif` | CNV analysis (Epi2me) |
| Modkit | `modkit_latest_latest.sif` | Modified base calling (Epi2me) |

##  Testing

Run the test script to verify your setup:

```bash
./test_pipeline_singularity.sh
```

This will run a minimal test to ensure all components are working correctly.

##  Directory Structure

After setup, your directory structure should look like:

```
nWGS_pipeline/
├── conf/
│   ├── analysis.config      # Main Singularity/Apptainer configuration
│   └── epi2me.config        # Epi2me Singularity/Apptainer configuration
├── containers/
│   ├── nwgs_default_images_latest.sif
│   ├── ace_1.24.0_latest.sif
│   ├── clair3_amd64_latest.sif
│   └── ... (other .sif files)
├── data/
│   ├── reference/           # Reference files
│   ├── humandb/            # Annotation databases
│   ├── testdata/           # Input data
│   └── results/            # Output results
├── setup_singularity.sh    # Automated setup script
├── run_pipeline_singularity.sh # Pipeline runner script
├── test_pipeline_singularity.sh # Test script
└── SINGULARITY_SETUP.md    # This guide
```


### Common Issues

1. **Singularity/Apptainer not installed**
   ```bash
   # Check installation
   singularity --version
   # or
   apptainer --version
   ```

2. **Permission denied**
   ```bash
   # Make sure you have execute permissions
   chmod +x setup_singularity.sh
   chmod +x run_pipeline_singularity.sh
   ```

3. **Out of disk space**
   ```bash
   # Clean up old images
   rm containers/*.sif
   # Re-run setup
   ./setup_singularity.sh
   ```

4. **Memory issues**
   - Increase memory allocation in configuration files
   - Use fewer parallel processes

### Getting Help

- Check Singularity/Apptainer logs: `singularity exec containers/nwgs_default_images_latest.sif --help`
- Verify image availability: `ls -la containers/*.sif`
- Test individual containers: `singularity exec containers/nwgs_default_images_latest.sif --help`

## Performance Tips

1. **Use local storage** for better I/O performance
2. **Allocate sufficient memory** (8GB+ recommended)
3. **Use multiple cores** by adjusting `cpus` parameter in config
4. **Mount data directories** efficiently using bind mounts

##  Additional Resources

- [Apptainer Documentation](https://apptainer.org/docs/)
- [Singularity Documentation](https://sylabs.io/docs/)
- [Nextflow Documentation](https://www.nextflow.io/docs/latest/)
- [nWGS Pipeline Documentation](README.md)

## Switching Between Container Systems

If you need to switch between Docker and Singularity/Apptainer:

### From Docker to Singularity/Apptainer:
```bash
# Use Singularity/Apptainer setup
./setup_singularity.sh
./run_pipeline_singularity.sh
```

### From Singularity/Apptainer to Docker:
```bash
# Use Docker setup
./setup_docker.sh
./run_pipeline_docker.sh
```

Both setups are compatible and use the same configuration files.

---

**Note**: All Singularity/Apptainer images are automatically pulled from Docker Hub and converted to Singularity format during setup. The pipeline supports both Singularity and Apptainer seamlessly. 