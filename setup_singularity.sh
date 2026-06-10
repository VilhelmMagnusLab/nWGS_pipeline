#!/bin/bash

# Diana Pipeline Singularity/Apptainer Setup Script
set -e

echo "=========================================="
echo "Diana Pipeline Singularity/Apptainer Setup"
echo "=========================================="

# Check if Singularity or Apptainer is installed
SINGULARITY_CMD=""
if command -v apptainer &> /dev/null; then
    SINGULARITY_CMD="apptainer"
    echo "✓ Apptainer found"
elif command -v singularity &> /dev/null; then
    SINGULARITY_CMD="singularity"
    echo "✓ Singularity found"
else
    echo "   Neither Singularity nor Apptainer is installed."
    echo "   Please install one of them:"
    echo "   - Apptainer: https://apptainer.org/docs/admin/main/installation.html"
    echo "   - Singularity: https://sylabs.io/guides/latest/admin-guide/installation.html"
    exit 1
fi

echo "Using: $SINGULARITY_CMD"

# Ask user for working directory location
echo ""
echo "=========================================="
echo "Working Directory Setup"
echo "=========================================="
echo "The pipeline requires a working directory structure for processing data."
echo "This will create the following directories:"
echo "  - routine_bams/       (Processed BAM files)"
echo "  - routine_epi2me/     (Epi2me module results)"
echo "  - routine_annotation/   (Analysis module results)"
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
mkdir -p "$WORK_DIR/routine_annotation"
mkdir -p "$WORK_DIR/routine_results"

echo "✓ Working directory structure created successfully"
echo ""

# Create necessary directories in the package
echo "Creating package directories..."
mkdir -p data/reference
mkdir -p data/humandb
mkdir -p containers
mkdir -p .empty_r_overlay

# Pull Singularity/Apptainer images
echo "   Pulling Singularity/Apptainer images..."
echo "   This may take several minutes on first run..."

# Function to pull image if it doesn't exist
pull_if_not_exists() {
    local image_name=$1
    # Extract just the image name from the repository path (e.g., "vilhelmmagnuslab/diana_default_images" -> "diana_default_images")
    local image_basename=$(basename "$image_name")
    local image_file="containers/${image_basename}_latest.sif"
    
    if [ -f "$image_file" ]; then
        echo "   ✓ $image_name already exists, skipping..."
    else
        echo "   Pulling $image_name..."
        $SINGULARITY_CMD pull --dir containers/ docker://$image_name:latest
    fi
}

# Core analysis images
echo "Pulling core analysis images..."
pull_if_not_exists "vilhelmmagnuslab/nwgs_default_images"
pull_if_not_exists "vilhelmmagnuslab/ace_1.24.0"
pull_if_not_exists "vilhelmmagnuslab/annotcnv_images_27feb1025"
pull_if_not_exists "vilhelmmagnuslab/clair3_amd64"
pull_if_not_exists "vilhelmmagnuslab/clairsto_amd64"
pull_if_not_exists "vilhelmmagnuslab/igv_report_amd64"
pull_if_not_exists "vilhelmmagnuslab/vcf2circos"
pull_if_not_exists "vilhelmmagnuslab/nanodx_images_3feb25"
pull_if_not_exists "vilhelmmagnuslab/crossnnumap"
pull_if_not_exists "vilhelmmagnuslab/markdown_images_28feb2025"
pull_if_not_exists "vilhelmmagnuslab/mgmt_nanopipe_amd64_18feb2025_cramoni"
pull_if_not_exists "vilhelmmagnuslab/gviz_amd64ps"
pull_if_not_exists "vilhelmmagnuslab/sturgeon_amd64_21jan"

# Optional VEP image (only needed when snv_annotator = "vep")
pull_if_not_exists "ensemblorg/ensembl-vep"

# Epi2me images
echo "Pulling Epi2me analysis images..."
pull_if_not_exists "vilhelmmagnuslab/snifflesv252_update"
pull_if_not_exists "vilhelmmagnuslab/qdnaseq_amd64"
pull_if_not_exists "vilhelmmagnuslab/modkit"

echo "✓ All Singularity/Apptainer images pulled successfully"

# Create a user-friendly run script
cat > run_pipeline_singularity.sh << 'EOF'
#!/bin/bash

# Diana Pipeline Runner Script for Singularity/Apptainer
set -e

# Check if Nextflow is installed
if ! command -v nextflow &> /dev/null; then
    echo " Nextflow is not installed. Installing Nextflow..."
    curl -s https://get.nextflow.io | bash
    export PATH="$PWD:$PATH"
fi

# Check if Singularity/Apptainer is available
SINGULARITY_CMD=""
if command -v apptainer &> /dev/null; then
    SINGULARITY_CMD="apptainer"
elif command -v singularity &> /dev/null; then
    SINGULARITY_CMD="singularity"
else
    echo " Neither Singularity nor Apptainer is available."
    echo "   Please run setup_singularity.sh first."
    exit 1
fi

echo "Using: $SINGULARITY_CMD"

# Auto-detect config file based on arguments
CONFIG="conf/annotation.config"  # Default config

# Check if epi2me mode is specified
if [[ "$*" == *"--run_mode_epi2me"* ]]; then
    CONFIG="conf/epi2me.config"
    echo " Using Epi2me configuration: $CONFIG"
elif [[ "$*" == *"--run_mode_mergebam"* ]]; then
    CONFIG="conf/mergebam.config"
    echo " Using Mergebam configuration: $CONFIG"
else
    echo " Using default annotation configuration: $CONFIG"
fi

echo " Starting Diana pipeline with Singularity/Apptainer containers..."
echo "   Configuration: $CONFIG"
echo "   Arguments: $@"

# Run the pipeline with Singularity/Apptainer
nextflow run main.nf \
    -c "$CONFIG" \
    -with-singularity \
    "$@"

echo " Pipeline completed successfully!"
EOF

chmod +x run_pipeline_singularity.sh

# Create a quick test script
cat > test_pipeline_singularity.sh << 'EOF'
#!/bin/bash

# Quick test script for the Diana pipeline with Singularity/Apptainer
set -e

echo " Running quick pipeline test with Singularity/Apptainer..."

# Create a minimal test sample file
mkdir -p data/testdata
echo "test_sample" > data/testdata/sample_ids.txt

# Run with test profile (if available)
if [ -f "conf/test.config" ]; then
    ./run_pipeline_singularity.sh -profile test
else
    echo " No test configuration found. Running with default config..."
    ./run_pipeline_singularity.sh -profile test
fi

echo " Test completed!"
EOF

chmod +x test_pipeline_singularity.sh

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
echo "  ├── routine_annotation/     # Analysis module results"
echo "  └── routine_results/      # Final results"
echo ""
echo "Next steps:"
echo "1. Place your reference files in: data/reference/"
echo "2. Update the configuration files with your working directory paths:"
echo "   - conf/annotation.config"
echo "   - conf/epi2me.config"
echo "   - conf/mergebam.config"
echo "3. Run the pipeline with: ./run_pipeline_singularity.sh"
echo "4. Test the setup with: ./test_pipeline_singularity.sh"
echo ""
echo "Pipeline modes:"
echo "  - Mergebam:    ./run_pipeline_singularity.sh --run_mode_mergebam"
echo "  - Epi2me:      ./run_pipeline_singularity.sh --run_mode_epi2me"
echo "  - Annotation:  ./run_pipeline_singularity.sh --run_mode_annotation"
echo ""
echo "Before starting the pipeline, make sure that all paths for each mode are correctly"
echo "set in the appropriate config files: conf/annotation.config, conf/epi2me.config, and conf/mergebam.config"
echo "For more information, see the README.md file."
echo ""
echo "Happy analyzing! 🧬 with Diana pipeline" 
