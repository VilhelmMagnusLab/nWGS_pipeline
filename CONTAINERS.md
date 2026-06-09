# Container Quick Reference

This is a quick reference guide for all Singularity/Docker containers used in the nWGS pipeline.

## Container Images Location

**Singularity:** `${params.nWGS_dir}/containers/`
**Docker Hub:** `vilhelmmagnuslab/[container-name]:latest`

---

## Supported Architectures

DIANA container images are built for **linux/amd64** (also known as **x86_64**) — the standard 64-bit architecture of Intel- and AMD-based servers and workstations. The *amd64* label does not imply that an AMD processor is required; any x86_64 CPU is supported.

The current release is tested exclusively on x86_64 Linux systems. ARM64 platforms (including Apple Silicon Macs running Docker Desktop or Linux ARM servers) are not currently validated and may exhibit reduced performance or unexpected behaviour due to instruction-set emulation.

**Running on non-x86_64 systems:**
- **Docker:** add `--platform linux/amd64` to force emulation (e.g. `docker run --platform linux/amd64 ...`)
- **Singularity/Apptainer:** amd64 images will run under QEMU emulation on ARM hosts if the host kernel supports it, but this is unsupported and not recommended for production use
- For native ARM64 support, rebuild the containers from the provided Dockerfiles targeting `linux/arm64`

---

## Container Index

| Container Name | Purpose | Key Tools | Config File |
|----------------|---------|-----------|-------------|
| **clair3_amd64** | SNV calling | Clair3, samtools 1.15.1 | epi2me.config, annotation.config |
| **clairsto_amd64** | Somatic SNV calling | ClairS-TO, bcftools | epi2me.config, annotation.config |
| **snifflesv252_update** | Structural variants | Sniffles v2.5.2, truvari | epi2me.config |
| **qdnaseq_amd64** | CNV detection | QDNAseq, R 4.3.2 | epi2me.config |
| **ace_1.24.0** | Copy number estimation | ACE, QDNAseq | annotation.config |
| **annotcnv_images_27feb1025** | CNV annotation | cramino 0.16.0, R 4.2.3 | annotation.config |
| **modkit** | Methylation calling | modkit | epi2me.config |
| **sturgeon_amd64_21jan** | Methylation classification | Sturgeon | annotation.config |
| **mgmt_nanopipe_cramino** | BAM statistics | cramino 0.16.0 | epi2me.config |
| **nanodx_images_3feb25** | NN classification | PyTorch 1.12.1 | annotation.config |
| **crossnnumap** | t-SNE/UMAP plots | R t-SNE, UMAP | annotation.config |
| **vcf2circos** | Circos plots | vcf2circos | annotation.config |
| **igv_report_amd64** | IGV.js reports | IGV.js | annotation.config |
| **gviz_amd64** | Gene coverage plots | Gviz | annotation.config |
| **markdown_images_28feb2025** | RMarkdown reports | pandoc, knitr | annotation.config |
| **nwgs_default_images** | General utilities | Various | annotation.config, mergebam.config |

---

## Process to Container Mapping

### Epi2me Workflow (modules/epi2me.nf)

| Process | Container |
|---------|-----------|
| run_epi2me_modkit | modkit |
| run_epi2me_cnv | qdnaseq_amd64 |
| run_epi2me_sv | snifflesv252_update |
| run_clair3 | clair3_amd64 |
| run_clairs_to | clairsto_amd64 |
| cramino_report | mgmt_nanopipe_cramino |

### Annotation Workflow (modules/annotation.nf)

| Process | Container |
|---------|-----------|
| extract_epic | nwgs_default_images |
| sturgeon | sturgeon_amd64_21jan |
| mgmt_promoter | nwgs_default_images |
| nanodx | nanodx_images_3feb25 |
| run_nn_classifier | nwgs_default_images |
| tsne_plot | nwgs_default_images |
| svannasv | nwgs_default_images |
| circosplot | vcf2circos |
| svannasv_fusion_events | nwgs_default_images |
| clair3_annotate | clair3_amd64 |
| clairs_to_annotate | clairsto_amd64 |
| merge_annotation | nwgs_default_images |
| igv_tools | igv_report_amd64 |
| plot_genomic_regions | gviz_amd64 |
| ace_tmc | ace_1.24.0 |
| annotatecnv | annotcnv_images_27feb1025 |
| markdown_report | markdown_images_28feb2025 |
| tertp | nwgs_default_images |
| crossnn_tsne | crossnnumap |

### MergeBam Workflow (modules/mergebam.nf)

| Process | Container |
|---------|-----------|
| mergebam | nwgs_default_images |
| bam_index | nwgs_default_images |

---

## Container Configuration

All containers are configured in the following config files:
- `conf/epi2me.config` - Epi2me workflow containers
- `conf/annotation.config` - Annotation workflow containers
- `conf/mergebam.config` - BAM merging containers

### Switching Between Docker and Singularity

**Singularity (default):**
```groovy
container = "${params.nWGS_dir}/containers/[name].sif"
```

**Docker:**
```groovy
container = 'vilhelmmagnuslab/[name]:latest'
```

To use Docker, comment out the Singularity line and uncomment the Docker line in the respective config file.

---

## Building Containers

### From Docker to Singularity

```bash
# Pull Docker image and convert to Singularity
singularity pull [name].sif docker://vilhelmmagnuslab/[name]:latest

# Or build from Dockerfile
docker build -f dockerfiles/Dockerfile_[name] -t vilhelmmagnuslab/[name]:latest .
singularity build [name].sif docker-daemon://vilhelmmagnuslab/[name]:latest
```

### Container Cache Location

Containers are stored in: `${params.nWGS_dir}/containers/`

Default: `/home/godzilla/nWGS_pipeline/containers/`

---

## Updating Containers

To update all containers from Docker Hub:

```bash
cd ${params.nWGS_dir}/containers/
for img in *.sif; do
    name=$(basename $img .sif | sed 's/_latest$//')
    echo "Pulling $name..."
    singularity pull --force $img docker://vilhelmmagnuslab/$name:latest
done
```

---

## Troubleshooting

### Container Not Found
```
Error: Container image not found
```
**Solution:** Check that the container exists in `${params.nWGS_dir}/containers/`

### Permission Denied
```
Error: Failed to bind mount
```
**Solution:** Check Apptainer runOptions in nextflow.config:
```groovy
runOptions = '--bind /home/godzilla:/home/godzilla --bind /data:/data'
```

### R Package Errors
```
Error: R_LIBS_USER environment variable
```
**Solution:** Some containers need R environment cleared. Check nextflow.config runOptions.

---

## Version Information

For detailed version information about tools inside each container, see:
- **SOFTWARE_VERSIONS.md** - Complete version listing
- **versions.yml** - Machine-readable format
- **dockerfiles/** - Container build definitions

Last updated: 2026-02-26
