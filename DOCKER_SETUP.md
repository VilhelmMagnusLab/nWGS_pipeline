# nWGS Pipeline - Docker Setup Guide

This guide explains how to run the nWGS pipeline using Docker containers, making it easy to get started without managing local Singularity images.

## 🐳 Quick Start

### Prerequisites

1. **Docker**: Install Docker on your system
   - [Docker Desktop](https://docs.docker.com/desktop/) (Windows/Mac)
   - [Docker Engine](https://docs.docker.com/engine/install/) (Linux)

2. **Nextflow**: The pipeline will automatically install Nextflow if not present

### One-Command Setup

Run the automated setup script:

```bash
chmod +x setup_docker.sh
./setup_docker.sh
```

This script will:
- ✅ Check Docker installation
- 📁 Create necessary directories
- 🐳 Pull all required Docker images from `vilhelmmagnuslab` repository
- 🚀 Create convenient run scripts

### Running the Pipeline

After setup, simply run:

```bash
./run_pipeline_docker.sh
```

Or with custom configuration:

```bash
./run_pipeline_docker.sh conf/analysis.config -profile your_profile
```

## 📋 Manual Setup (Alternative)

If you prefer manual setup:

### 1. Pull Docker Images

```bash
# Core analysis images
docker pull vilhelmmagnuslab/nwgs_default_images:latest
docker pull vilhelmmagnuslab/ace_1.24.0:latest
docker pull vilhelmmagnuslab/annotcnv_images_27feb1025:latest
docker pull vilhelmmagnuslab/clair3:latest
docker pull vilhemmanguslab/clars-to:latest
docker pull vilhelmmagnuslab/igv_report_amd64:latest
docker pull vilhelmmagnuslab/vcf2circos:latest
docker pull vilhelmmagnuslab/nanodx_images_3feb25:latest
docker pull vilhelmmagnuslab/markdown_images_28feb2025:latest
docker pull vilhelmmagnuslab/mgmt_nanopipe_amd64_18feb2025_cramoni:latest
docker pull vilhelmmagnuslab/gviz_amd64ps:latest

# Epi2me images
docker pull vilhelmmagnuslab/snifflesv252_update_latest:latest
docker pull vilhelmmagnuslab/qdnaseq_amd64_latest:latest
docker pull vilhelmmagnuslab/modkit_latest:latest
```

### 2. Create Directory Structure

```bash
mkdir -p data/{reference,humandb,testdata,results}
```

### 3. Run with Nextflow

```bash
nextflow run main.nf -c conf/analysis.config -with-docker
```

## 🔧 Configuration

### Docker Configuration

The pipeline is configured to use Docker containers in multiple configuration files:

- `conf/analysis.config` - Main analysis configuration
- `conf/epi2me.config` - Epi2me-specific configuration  
- `conf/mergebam.config` - Mergebam-specific configuration

Key Docker settings:
```groovy
docker {
    enabled = true
    runOptions = '-v /home/chbope/extension/nWGS_manuscript_data/data:/home/chbope/extension/nWGS_manuscript_data/data'
}
```

### Pipeline Modes

The nWGS pipeline supports three main modes:

1. **Analysis Mode** (default): Complete analysis workflow
   ```bash
   ./run_pipeline_docker.sh --run_mode_analysis
   ```

2. **Epi2me Mode**: Multi-omics analysis with Epi2me tools
   ```bash
   ./run_pipeline_docker.sh --run_mode_epi2me
   ```

3. **Mergebam Mode**: BAM merging and region extraction
   ```bash
   ./run_pipeline_docker.sh --run_mode_mergebam
   ```

### Available Docker Images

| Process | Docker Image | Description |
|---------|-------------|-------------|
| Default | `vilhelmmagnuslab/nwgs_default_images` | General analysis tools |
| ACE TMC | `vilhelmmagnuslab/ace_1.24.0` | ACE copy number analysis |
| AnnotateCNV | `vilhelmmagnuslab/annotcnv_images_27feb1025` | CNV annotation |
| ClairS-TO | `vilhelmmagnuslab/clairsto_amd64` | somatic small variant calling |
| Clair3 | `vilhelmmagnuslab/clair3_amd64` | variant calling |
| IGV Tools | `vilhelmmagnuslab/igv_report_amd64` | IGV report generation |
| Circos Plot | `vilhelmmagnuslab/vcf2circos` | Circos visualization |
| NanoDx | `vilhelmmagnuslab/nanodx_images_3feb25` | NanoDx classification |
| Markdown | `vilhelmmagnuslab/markdown_images_28feb2025` | Report generation |
| Cramino | `vilhelmmagnuslab/mgmt_nanopipe_amd64_18feb2025_cramoni` | Quality assessment |
| Gviz | `vilhelmmagnuslab/gviz_amd64ps` | Genomic visualization |
| Sniffles2 | `vilhelmmagnuslab/snifflesv252_update_latest` | SV calling (Epi2me) |
| QDNAseq | `vilhelmmagnuslab/qdnaseq_amd64_latest` | CNV analysis (Epi2me) |
| Modkit | `vilhelmmagnuslab/modkit_latest` | Modified base calling |

## 🧪 Testing

Run the test script to verify your setup:

```bash
./test_pipeline_docker.sh
```

This will run a minimal test to ensure all components are working correctly.

## 📁 Directory Structure

After setup, your directory structure should look like:

```
nWGS_pipeline/
├── conf/
│   ├── analysis.config      # Main Docker configuration
│   ├── epi2me.config        # Epi2me Docker configuration
│   └── mergebam.config      # Mergebam Docker configuration
├── data/
│   ├── reference/           # Reference files
│   ├── humandb/            # Annotation databases
│   ├── testdata/           # Input data
│   └── results/            # Output results
├── setup_docker.sh         # Automated setup script
├── run_pipeline_docker.sh    # Pipeline runner script (Docker)
├── test_pipeline_docker.sh   # Test script (Docker)
└── DOCKER_SETUP.md         # This guide
```

## 🔍 Troubleshooting

### Common Issues

1. **Docker not running**
   ```bash
   # Start Docker daemon
   sudo systemctl start docker
   ```

2. **Permission denied**
   ```bash
   # Add user to docker group
   sudo usermod -aG docker $USER
   # Log out and back in
   ```

3. **Out of disk space**
   ```bash
   # Clean up Docker images
   docker system prune -a
   ```

4. **Memory issues**
   - Increase Docker memory limit in Docker Desktop settings
   - Or reduce memory allocation in configuration files

5. **Path configuration issues**
   - Ensure all paths in config files point to correct locations
   - Update `params.path` in configuration files to match your data directory

### Getting Help

- Check Docker logs: `docker logs <container_id>`
- Verify image availability: `docker images | grep vilhelmmagnuslab`
- Test individual containers: `docker run --rm vilhelmmagnuslab/nwgs_default_images:latest --help`

## 🚀 Performance Tips

1. **Use SSD storage** for better I/O performance
2. **Allocate sufficient memory** to Docker (8GB+ recommended)
3. **Use multiple cores** by adjusting `cpus` parameter in config
4. **Mount data directories** efficiently using Docker volumes
5. **Use Docker layer caching** by keeping base images updated

## 📋 Configuration Checklist

Before running the pipeline, ensure:

- [ ] Docker is installed and running
- [ ] All Docker images are pulled successfully
- [ ] Reference files are placed in `data/reference/`
- [ ] Input data is placed in `data/testdata/`
- [ ] Paths in configuration files are updated to match your setup
- [ ] Sample IDs file is created (`data/testdata/sample_ids.txt`)

## 📚 Additional Resources

- [Docker Documentation](https://docs.docker.com/)
- [Nextflow Documentation](https://www.nextflow.io/docs/latest/)
- [nWGS Pipeline Documentation](README.md)
- [Singularity Setup Guide](SINGULARITY_SETUP.md)

## 🔄 Migration from Singularity

If you're migrating from Singularity to Docker:

1. Update configuration files to use Docker containers
2. Replace `.sif` container references with Docker image names
3. Update volume mounting syntax for Docker
4. Test with the provided test script

---

**Note**: All Docker images are hosted at [https://hub.docker.com/repositories/vilhelmmagnuslab](https://hub.docker.com/repositories/vilhelmmagnuslab) and are automatically pulled during setup.

**Support**: For issues specific to Docker setup, please check the troubleshooting section above or refer to the main pipeline documentation. 