#!/bin/bash

# Diana Pipeline Docker Setup Script
# This script helps users set up and run the Diana pipeline using Docker images

set -e

echo "=========================================="
echo "Diana Pipeline Docker Setup"
echo "=========================================="

# Check if Docker is installed
if ! command -v docker &> /dev/null; then
    echo " Docker is not installed. Please install Docker first."
    echo "   Visit: https://docs.docker.com/get-docker/"
    exit 1
fi

# Check if Docker daemon is running
if ! docker info &> /dev/null; then
    echo " Docker daemon is not running. Please start Docker first."
    exit 1
fi

echo "Docker is available and running"

# Ask user for working directory location
echo ""
echo "=========================================="
echo "Working Directory Setup"
echo "=========================================="
echo "The pipeline requires a working directory structure for processing data."
echo "This will create the following directories:"
echo "  - routine_bams/       (Processed BAM files)"
echo "  - routine_epi2me/     (Epi2me module results)"
echo "  - routine_analysis/   (Analysis module results)"
echo "  - routine_results/    (Final results)"
echo ""
read -p "Enter the parent directory path [default: ~/Documents]: " WORK_DIR_PARENT
WORK_DIR_PARENT=${WORK_DIR_PARENT:-~/Documents}

# Expand tilde to home directory
WORK_DIR_PARENT="${WORK_DIR_PARENT/#\~/$HOME}"

# Create the main working directory
WORK_DIR="${WORK_DIR_PARENT}/routine_diana"
echo "Creating working directory structure at: $WORK_DIR"

mkdir -p "$WORK_DIR/routine_bams/merge_bams"
mkdir -p "$WORK_DIR/routine_bams/roi_bams"
mkdir -p "$WORK_DIR/routine_epi2me"
mkdir -p "$WORK_DIR/routine_analysis"
mkdir -p "$WORK_DIR/routine_results"

echo "✓ Working directory structure created successfully"
echo ""

# Create necessary directories in the package
echo "Creating package directories..."
mkdir -p data/reference
mkdir -p data/humandb

# Function to pull Docker image if it doesn't exist
pull_if_not_exists() {
    local image_name=$1
    local image_with_tag="${image_name}:latest"

    if docker images --format "{{.Repository}}:{{.Tag}}" | grep -q "^${image_with_tag}$"; then
        echo "   ✓ $image_with_tag already exists, skipping..."
    else
        echo "   Pulling $image_with_tag..."

        # Temporarily disable Docker authentication for public images
        # by hiding the Docker config file
        local docker_config_backup=""
        if [ -f "$HOME/.docker/config.json" ]; then
            docker_config_backup="$HOME/.docker/config.json.docker_backup_$$"
            mv "$HOME/.docker/config.json" "$docker_config_backup"
        fi

        # Pull the image without authentication (for public images)
        docker pull "$image_with_tag"

        # Restore Docker config if it was backed up
        if [ -n "$docker_config_backup" ] && [ -f "$docker_config_backup" ]; then
            mv "$docker_config_backup" "$HOME/.docker/config.json"
        fi
    fi
}

# Pull Docker images (this will take some time on first run)
echo "🐳 Pulling Docker images from vilhelmmagnuslab repository..."
echo "   This may take several minutes on first run..."

# Core analysis images
echo "Pulling core analysis images..."
pull_if_not_exists "vilhelmmagnuslab/nwgs_default_images"
pull_if_not_exists "vilhelmmagnuslab/ace_1.24.0"
pull_if_not_exists "vilhelmmagnuslab/annotcnv_images_27feb1025"
pull_if_not_exists "vilhelmmagnuslab/clair3_amd64"
pull_if_not_exists "vilhelmmagnuslab/clairsto_amd64"
pull_if_not_exists "vilhelmmagnuslab/igv_report_amd64"
pull_if_not_exists "vilhelmmagnuslab/vcf2circos"
pull_if_not_exists "vilhelmmagnuslab/nanodx_env"
pull_if_not_exists "vilhelmmagnuslab/crossnnumap"
pull_if_not_exists "vilhelmmagnuslab/markdown_images_28feb2025"
pull_if_not_exists "vilhelmmagnuslab/mgmt_nanopipe_amd64_18feb2025_cramoni"
pull_if_not_exists "vilhelmmagnuslab/gviz_amd64"
pull_if_not_exists "vilhelmmagnuslab/sturgeon_amd64_21jan"

# Optional VEP image (only needed when snv_annotator = "vep")
pull_if_not_exists "ensemblorg/ensembl-vep"

# Epi2me images
echo "Pulling Epi2me analysis images..."
pull_if_not_exists "vilhelmmagnuslab/snifflesv252_update"
pull_if_not_exists "vilhelmmagnuslab/qdnaseq_amd64"
pull_if_not_exists "vilhelmmagnuslab/modkit"

echo "✓ All Docker images pulled successfully"

# Create a simple run script
cat > run_pipeline_docker.sh << 'EOF'
#!/bin/bash

# Diana Pipeline Runner Script for Docker
# Usage: ./run_pipeline_docker.sh [run_mode] [other_nextflow_options]

set -e

# Check if Nextflow is installed
if ! command -v nextflow &> /dev/null; then
    echo " Nextflow is not installed. Installing Nextflow..."
    curl -s https://get.nextflow.io | bash
    export PATH="$PWD:$PATH"
fi

# Auto-detect config file based on arguments
CONFIG="conf/analysis.config"  # Default config

# Check if epi2me mode is specified
if [[ "$*" == *"--run_mode_epi2me"* ]]; then
    CONFIG="conf/epi2me.config"
    echo " Using Epi2me configuration: $CONFIG"
elif [[ "$*" == *"--run_mode_mergebam"* ]]; then
    CONFIG="conf/mergebam.config"
    echo " Using Mergebam configuration: $CONFIG"
else
    echo "Using default analysis configuration: $CONFIG"
fi

echo " Starting Diana pipeline with Docker containers..."
echo "   Configuration: $CONFIG"
echo "   Arguments: $@"

# Run the pipeline
nextflow run main.nf \
    -c "$CONFIG" \
    -with-docker \
    "$@"

echo "Pipeline completed successfully!"
EOF

chmod +x run_pipeline_docker.sh

# Create a quick test script
cat > test_pipeline_docker.sh << 'EOF'
#!/bin/bash

# Quick test script for the Diana pipeline (Docker)
# This will run a minimal test to verify everything is working

set -e

echo "🧪 Running quick pipeline test with Docker..."

# Run with test profile (if available)
if [ -f "conf/test.config" ]; then
    ./run_pipeline_docker.sh conf/test.config -profile test
else
    echo " No test configuration found. Running with default config..."
    ./run_pipeline_docker.sh conf/analysis.config -profile test
fi

echo " Test completed!"
EOF

chmod +x test_pipeline_docker.sh

echo ""
echo "=========================================="
echo "Setup Complete!"
echo "=========================================="
echo ""
echo "Working directory created at: $WORK_DIR"
echo "Directory structure:"
echo "  $WORK_DIR/"
echo "  ├── routine_bams/         # Processed BAM files"
echo "  │   ├── merge_bams/       # Merged BAM files per sample"
echo "  │   └── roi_bams/         # Region of interest extracted BAMs"
echo "  ├── routine_epi2me/       # Epi2me module results"
echo "  ├── routine_analysis/     # Analysis module results"
echo "  └── routine_results/      # Final results"
echo ""
echo "Next steps:"
echo "1. Place your reference files in: data/reference/"
echo "2. Update the configuration files with your working directory paths:"
echo "   - conf/analysis.config"
echo "   - conf/epi2me.config"
echo "   - conf/mergebam.config"
echo "3. Run the pipeline with: ./run_pipeline_docker.sh"
echo "4. Test the setup with: ./test_pipeline_docker.sh"
echo ""
echo "Pipeline modes:"
echo "  - Mergebam:    ./run_pipeline_docker.sh --run_mode_mergebam"
echo "  - Epi2me:      ./run_pipeline_docker.sh --run_mode_epi2me"
echo "  - Annotation:  ./run_pipeline_docker.sh --run_mode_annotation"
echo ""
echo "Before starting the pipeline, make sure that all paths for each mode are correctly"
echo "set in the appropriate config files: conf/annotation.config, conf/epi2me.config, and conf/mergebam.config"
echo "For more information, see the README.md file."
echo ""
echo "Happy analyzing! 🧬" 
