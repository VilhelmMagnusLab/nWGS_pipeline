#!/bin/bash
################################################################################
# Diana Pipeline - Automated Setup Script
################################################################################
# This script automatically downloads all reference files from Zenodo and
# sets up the pipeline for immediate use.
#
# Zenodo Record: https://doi.org/10.5281/zenodo.20596458
#
# Usage:
#   ./setup_pipeline.sh docker|singularity [OPTIONS]
#
# Arguments:
#   docker          Use Docker containers
#   singularity     Use Singularity/Apptainer containers
#
# Options:
#   --skip-optional          Skip downloading optional large files (~20 GB)
#   --skip-containers        Skip container setup (only download reference files)
#   --skip-reference         Skip reference file download (only setup containers)
#   --config-only            Only create directories and install Nextflow (no images, no Zenodo)
#   --work-dir DIR           Parent directory for routine_diana/ (default: $HOME)
#   --help                   Show this help message
#
# Examples:
#   ./setup_pipeline.sh docker
#   ./setup_pipeline.sh singularity --skip-optional
#   ./setup_pipeline.sh docker --skip-containers
#   ./setup_pipeline.sh singularity --config-only
################################################################################

set -e  # Exit on error

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
NC='\033[0m' # No Color

# Configuration
ZENODO_RECORD="20596458"  # Diana pipeline reference files v3 (DOI: 10.5281/zenodo.20596458)
BASE_URL="https://zenodo.org/record/${ZENODO_RECORD}/files"

# MD5 checksums for Zenodo archives — update these whenever archives are rebuilt
declare -A ARCHIVE_MD5=(
    ["diana_dummy.tar.gz"]="3067976db7e7bd190a415099c5c1cc53"
    ["humandb.tar.gz"]="a63ae92e6245129116b0d28f92d7a512"
    ["reference_core.tar.gz"]="620129cbe4143a3451d83eaf11f1944f"
    ["general.zip"]="d3553e44f1bbab8c820e04168c5b8e59"
    ["Assembly.zip"]="de6bc8bed97e433cbe8be55ed0e1536a"
    ["svanna-data.zip"]="6a5cabf40172cc420553ae8ca6ea3805"
    ["r1041_e82_400bps_sup_v420.zip"]="ab2112991f0cd224eb6c0f43d22acf3d"
)
PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_DIR="${PIPELINE_DIR}/data"
REFERENCE_DIR="${DATA_DIR}/reference"
HUMANDB_DIR="${DATA_DIR}/humandb"

# Options
CONTAINER_SYSTEM=""
SKIP_OPTIONAL=false
SKIP_CONTAINERS=false
SKIP_REFERENCE=false
CONFIG_ONLY=false
WORK_DIR_PARENT=""

################################################################################
# Helper Functions
################################################################################

print_header() {
    echo -e "${BLUE}==========================================${NC}"
    echo -e "${BLUE}$1${NC}"
    echo -e "${BLUE}==========================================${NC}"
}

print_success() {
    echo -e "${GREEN}✅ $1${NC}"
}

print_error() {
    echo -e "${RED}❌ $1${NC}"
}

print_warning() {
    echo -e "${YELLOW}⚠️  $1${NC}"
}

print_info() {
    echo -e "${CYAN}ℹ️  $1${NC}"
}

check_command() {
    if command -v "$1" &> /dev/null; then
        return 0
    else
        return 1
    fi
}

download_file() {
    local filename="$1"
    local destination="$2"
    local url="${BASE_URL}/${filename}"

    echo -e "${CYAN}📥 Downloading ${filename}...${NC}"

    if check_command wget; then
        wget -q --show-progress "${url}" -O "${destination}" 2>&1 | \
            grep --line-buffered "%" | \
            sed -u -e "s,\.,,g" | \
            awk '{printf("\r  Progress: %s", $2+0); fflush()}'
        echo ""
    elif check_command curl; then
        curl -L "${url}" -o "${destination}" --progress-bar
    else
        print_error "Neither wget nor curl found. Please install one of them."
        exit 1
    fi

    if [ $? -eq 0 ]; then
        print_success "Downloaded ${filename}"
    else
        print_error "Failed to download ${filename}"
        print_info "URL: ${url}"
        exit 1
    fi
}

verify_checksum() {
    local file="$1"
    local filename
    filename="$(basename "${file}")"
    local expected="${ARCHIVE_MD5[$filename]}"

    # No checksum registered for this file — skip silently
    if [ -z "${expected}" ]; then
        return 0
    fi

    print_info "Verifying checksum for ${filename}..."

    local actual
    if check_command md5sum; then
        actual=$(md5sum "${file}" | awk '{print $1}')
    elif check_command md5; then
        actual=$(md5 -q "${file}")
    else
        print_warning "md5sum not found — skipping checksum verification"
        return 0
    fi

    if [ "${actual}" = "${expected}" ]; then
        print_success "Checksum OK: ${filename}"
    else
        print_error "Checksum mismatch for ${filename}!"
        print_error "  Expected: ${expected}"
        print_error "  Got:      ${actual}"
        print_error "The file is corrupted (likely an interrupted download)."
        rm -f "${file}"
        print_info "The bad file has been removed. Re-run setup_pipeline.sh to retry."
        exit 1
    fi
}

extract_archive() {
    local archive="$1"
    local destination="$2"

    echo -e "${CYAN}📦 Extracting $(basename ${archive})...${NC}"

    if [[ $archive == *.tar.gz ]] || [[ $archive == *.tgz ]]; then
        tar -xzf "${archive}" -C "${destination}"
    elif [[ $archive == *.tar ]]; then
        tar -xf "${archive}" -C "${destination}"
    elif [[ $archive == *.zip ]]; then
        if check_command unzip; then
            unzip -q "${archive}" -d "${destination}"
        else
            print_error "unzip command not found. Please install unzip."
            exit 1
        fi
    elif [[ $archive == *.gz ]] && [[ $archive != *.tar.gz ]]; then
        gunzip -k "${archive}"
    else
        print_error "Unknown archive format: ${archive}"
        return 1
    fi

    print_success "Extracted $(basename ${archive})"
}

show_usage() {
    grep "^#" "$0" | grep -v "^#!/bin/bash" | sed 's/^# //' | sed 's/^#//'
}

################################################################################
# Parse Command Line Arguments
################################################################################

parse_args() {
    # First argument must be container system
    if [ $# -eq 0 ]; then
        print_error "Missing required argument: container system"
        echo ""
        show_usage
        exit 1
    fi

    case "$1" in
        docker|singularity)
            CONTAINER_SYSTEM="$1"
            shift
            ;;
        --help|-h)
            show_usage
            exit 0
            ;;
        *)
            print_error "Invalid container system: $1"
            echo ""
            print_info "First argument must be either 'docker' or 'singularity'"
            echo ""
            show_usage
            exit 1
            ;;
    esac

    # Parse remaining options
    while [[ $# -gt 0 ]]; do
        case $1 in
            --skip-optional)
                SKIP_OPTIONAL=true
                shift
                ;;
            --skip-containers)
                SKIP_CONTAINERS=true
                shift
                ;;
            --skip-reference)
                SKIP_REFERENCE=true
                shift
                ;;
            --config-only)
                SKIP_REFERENCE=true
                SKIP_CONTAINERS=true
                CONFIG_ONLY=true
                shift
                ;;
            --work-dir)
                WORK_DIR_PARENT="$2"
                shift 2
                ;;
            --help|-h)
                show_usage
                exit 0
                ;;
            *)
                print_error "Unknown option: $1"
                echo "Use --help for usage information"
                exit 1
                ;;
        esac
    done
}

################################################################################
# Main Setup Functions
################################################################################

check_prerequisites() {
    print_header "Checking Prerequisites"

    local missing_deps=()

    # Check for wget or curl
    if ! check_command wget && ! check_command curl; then
        missing_deps+=("wget or curl")
    fi

    # Check for tar
    if ! check_command tar; then
        missing_deps+=("tar")
    fi

    # Check for unzip (needed for Assembly.zip, svanna-data.zip, and model files)
    if ! check_command unzip; then
        print_warning "unzip not found - will be needed for extracting .zip files"
        missing_deps+=("unzip")
    fi

    # Validate container system
    if [ "$SKIP_CONTAINERS" = false ]; then
        if [ "$CONTAINER_SYSTEM" = "docker" ]; then
            if ! check_command docker; then
                print_error "Docker not found but docker mode was specified"
                print_info "Install Docker: https://docs.docker.com/get-docker/"
                exit 1
            fi
            if ! docker info &> /dev/null; then
                print_error "Docker daemon is not running"
                print_info "Please start Docker and try again"
                exit 1
            fi
            print_success "Docker is available and running"
        elif [ "$CONTAINER_SYSTEM" = "singularity" ]; then
            if ! check_command singularity && ! check_command apptainer; then
                print_error "Neither Singularity nor Apptainer found but singularity mode was specified"
                print_info "Install Apptainer: https://apptainer.org/docs/admin/main/installation.html"
                print_info "Install Singularity: https://sylabs.io/guides/latest/admin-guide/installation.html"
                exit 1
            fi
            if check_command apptainer; then
                print_success "Apptainer is available"
            else
                print_success "Singularity is available"
            fi
        fi
    fi

    # Check for nextflow
    if ! check_command nextflow; then
        print_warning "Nextflow not found - will be installed automatically"
    else
        print_success "Nextflow is available"
    fi

    if [ ${#missing_deps[@]} -gt 0 ]; then
        print_error "Missing required dependencies:"
        for dep in "${missing_deps[@]}"; do
            echo "  - $dep"
        done
        echo ""
        exit 1
    fi

    print_success "All prerequisites met!"
    echo ""
}

create_directories() {
    print_header "Creating Directory Structure"

    mkdir -p "${REFERENCE_DIR}"
    mkdir -p "${HUMANDB_DIR}"
    mkdir -p "${PIPELINE_DIR}/logs"
    mkdir -p "${PIPELINE_DIR}/containers"
    mkdir -p "${REFERENCE_DIR}/nanoDx/static"
    mkdir -p "${PIPELINE_DIR}/.empty_r_overlay"

    print_success "Created pipeline directory structure"
    echo ""

    # Create routine_diana working directory
    if [ -z "$WORK_DIR_PARENT" ]; then
        echo "Where should the routine_diana/ working directory be created?"
        read -p "  Parent directory [default: $HOME]: " WORK_DIR_PARENT
        WORK_DIR_PARENT=${WORK_DIR_PARENT:-$HOME}
    fi

    # Expand tilde
    WORK_DIR_PARENT="${WORK_DIR_PARENT/#\~/$HOME}"
    local ROUTINE_DIR="${WORK_DIR_PARENT}/routine_diana"

    print_info "Creating working directory at: ${ROUTINE_DIR}"
    mkdir -p "${ROUTINE_DIR}/routine_bams/merge_bams"
    mkdir -p "${ROUTINE_DIR}/routine_bams/roi_bams"
    mkdir -p "${ROUTINE_DIR}/routine_epi2me"
    mkdir -p "${ROUTINE_DIR}/routine_annotation"
    mkdir -p "${ROUTINE_DIR}/routine_results"

    # Create sample ID files with default demo sample if they don't exist
    if [ ! -f "${ROUTINE_DIR}/sample_ids.txt" ]; then
        echo "diana-001" > "${ROUTINE_DIR}/sample_ids.txt"
    fi
    if [ ! -f "${ROUTINE_DIR}/sample_ids_bam.txt" ]; then
        printf "diana-001\tPBE00000\n" > "${ROUTINE_DIR}/sample_ids_bam.txt"
    fi
    chmod 664 "${ROUTINE_DIR}/sample_ids.txt" "${ROUTINE_DIR}/sample_ids_bam.txt"

    # Save chosen path so other scripts (smart_sample_monitor_v2.sh etc.) can find it
    echo "DIANA_ROUTINE_DIR=${ROUTINE_DIR}" > "${PIPELINE_DIR}/.diana_env"

    print_success "Created routine_diana/ structure at: ${ROUTINE_DIR}"

    # Update config files if the user chose a non-default path
    if [ "${WORK_DIR_PARENT}" != "${HOME}" ]; then
        print_info "Custom path detected — updating config files to use: ${ROUTINE_DIR}"
        local configs=("${PIPELINE_DIR}/conf/annotation.config" \
                       "${PIPELINE_DIR}/conf/epi2me.config" \
                       "${PIPELINE_DIR}/conf/mergebam.config")
        for cfg in "${configs[@]}"; do
            if [ -f "$cfg" ]; then
                sed -i "s|System.getProperty('user.home')|'${WORK_DIR_PARENT}'|g" "$cfg"
                print_success "Updated: $(basename $cfg)"
            fi
        done
    fi
    echo ""
}

verify_md5() {
    local file="$1"
    local md5_file="$2"

    if [ ! -f "$md5_file" ]; then
        print_warning "MD5 checksum file not found: $md5_file"
        return 1
    fi

    print_info "Verifying MD5 checksum for $(basename $file)..."

    if check_command md5sum; then
        local expected_md5=$(cat "$md5_file" | awk '{print $1}')
        local actual_md5=$(md5sum "$file" | awk '{print $1}')

        if [ "$expected_md5" = "$actual_md5" ]; then
            print_success "MD5 checksum verified"
            return 0
        else
            print_error "MD5 checksum mismatch!"
            print_error "Expected: $expected_md5"
            print_error "Got: $actual_md5"
            return 1
        fi
    else
        print_warning "md5sum not found - skipping checksum verification"
        return 0
    fi
}

setup_nanodx_classifier() {
    print_info "Verifying nanoDx classifier setup..."

    # Define paths
    local NANODX_DIR="${REFERENCE_DIR}/nanoDx"
    local NANODX_STATIC="${NANODX_DIR}/static"
    local ZENODO_NANODX="14006255"  # Specific Zenodo record for nanoDx files (fallback only)
    local NANODX_URL="https://zenodo.org/record/${ZENODO_NANODX}/files"

    # Note: nanoDx is included in reference_core.tar.gz, so this function mainly verifies

    # Create directory structure if needed
    mkdir -p "${NANODX_STATIC}"

    # Check if nanoDx folder exists in pipeline root (old location from previous versions)
    if [ -d "${PIPELINE_DIR}/nanoDx" ] && [ ! -d "${NANODX_DIR}" ]; then
        print_info "Moving nanoDx folder from pipeline root to data/reference/..."
        mv "${PIPELINE_DIR}/nanoDx" "${REFERENCE_DIR}/"
        print_success "Moved nanoDx to data/reference/"
    fi

    # Check if nanoDx model files are present (should be from reference_core.tar.gz)
    # Only download from Zenodo as fallback if files are missing
    local files_to_download=(
        "Capper_et_al.h5"
        "Capper_et_al.h5.md5"
        "Capper_et_al_NN.pkl"
    )

    local all_files_present=true
    for file in "${files_to_download[@]}"; do
        if [ ! -f "${NANODX_STATIC}/${file}" ]; then
            all_files_present=false
            break
        fi
    done

    if [ "$all_files_present" = true ]; then
        print_success "nanoDx classifier files verified (from reference_core.tar.gz)"
        echo ""
        return
    fi

    print_warning "nanoDx model files missing - downloading from Zenodo as fallback..."
    print_info "(Note: These should have been included in reference_core.tar.gz)"

    for file in "${files_to_download[@]}"; do
        if [ ! -f "${NANODX_STATIC}/${file}" ]; then
            echo -e "${CYAN}  Downloading ${file}...${NC}"

            if check_command wget; then
                wget -q --show-progress "${NANODX_URL}/${file}" -O "${NANODX_STATIC}/${file}"
            elif check_command curl; then
                curl -L "${NANODX_URL}/${file}" -o "${NANODX_STATIC}/${file}" --progress-bar
            fi

            if [ $? -eq 0 ]; then
                echo -e "  ${GREEN}✓${NC} Downloaded ${file}"
            else
                print_error "Failed to download ${file}"
                exit 1
            fi
        else
            echo -e "  ${GREEN}✓${NC} ${file} already exists"
        fi
    done

    # Verify MD5 checksum for the model file
    if [ -f "${NANODX_STATIC}/Capper_et_al.h5" ] && [ -f "${NANODX_STATIC}/Capper_et_al.h5.md5" ]; then
        verify_md5 "${NANODX_STATIC}/Capper_et_al.h5" "${NANODX_STATIC}/Capper_et_al.h5.md5"
    fi

    print_success "nanoDx classifier setup complete!"
    echo ""
}

setup_vep_cache() {
    local VEP_CACHE_VERSION="115"
    local VEP_CACHE_FILE="homo_sapiens_refseq_vep_${VEP_CACHE_VERSION}_GRCh38.tar.gz"
    local VEP_CACHE_URL="https://ftp.ensembl.org/pub/release-${VEP_CACHE_VERSION}/variation/indexed_vep_cache/${VEP_CACHE_FILE}"
    local VEP_CACHE_DIR="${HUMANDB_DIR}/homo_sapiens_refseq"

    print_info "Setting up Ensembl VEP cache (v${VEP_CACHE_VERSION} GRCh38)..."

    # Skip if already present
    if [ -d "${VEP_CACHE_DIR}" ]; then
        print_success "VEP cache (homo_sapiens_refseq) already present"
        echo ""
        return
    fi

    print_info "Downloading VEP cache from Ensembl (~15 GB, this will take a while)..."
    print_info "URL: ${VEP_CACHE_URL}"

    local download_ok=false
    local tmp_file="${HUMANDB_DIR}/${VEP_CACHE_FILE}"

    if check_command wget; then
        wget -q --show-progress "${VEP_CACHE_URL}" -O "${tmp_file}" 2>&1 | \
            grep --line-buffered "%" | \
            sed -u -e "s,\.,,g" | \
            awk '{printf("\r  Progress: %s", $2+0); fflush()}'
        echo ""
        [ ${PIPESTATUS[0]} -eq 0 ] && download_ok=true
    elif check_command curl; then
        curl -L "${VEP_CACHE_URL}" -o "${tmp_file}" --progress-bar
        [ $? -eq 0 ] && download_ok=true
    else
        print_warning "Neither wget nor curl found — cannot download VEP cache automatically."
    fi

    if [ "$download_ok" = true ] && [ -f "${tmp_file}" ]; then
        print_info "Extracting VEP cache to ${HUMANDB_DIR}..."
        tar -xzf "${tmp_file}" -C "${HUMANDB_DIR}" && \
            rm "${tmp_file}" && \
            print_success "VEP cache installed at ${VEP_CACHE_DIR}"
    else
        # Clean up partial download if it exists
        [ -f "${tmp_file}" ] && rm -f "${tmp_file}"
        print_warning "VEP cache download failed — skipping (VEP annotation will not work)."
        print_warning "To install manually, run:"
        print_warning "  wget ${VEP_CACHE_URL} -O ${HUMANDB_DIR}/${VEP_CACHE_FILE}"
        print_warning "  tar -xzf ${HUMANDB_DIR}/${VEP_CACHE_FILE} -C ${HUMANDB_DIR}"
        print_info "The pipeline will continue and use ANNOVAR as the default annotator."
    fi

    echo ""
}

setup_svanna_database() {
    print_info "Setting up Svanna database..."

    # Check if svanna-data already exists (either as directory or zip)
    if [ -d "${REFERENCE_DIR}/svanna-data" ]; then
        print_success "Svanna database already present"
        echo ""
        return
    fi

    # Extract svanna-data.zip (bundled in reference_core.tar.gz)
    if [ -f "${REFERENCE_DIR}/svanna-data.zip" ]; then
        print_info "Extracting Svanna database..."
        extract_archive "${REFERENCE_DIR}/svanna-data.zip" "${REFERENCE_DIR}"
        print_success "Svanna database setup complete!"
    else
        print_warning "svanna-data.zip not found in reference directory — downloading as fallback..."
        download_file "svanna-data.zip" "${REFERENCE_DIR}/svanna-data.zip"
        verify_checksum "${REFERENCE_DIR}/svanna-data.zip"
        extract_archive "${REFERENCE_DIR}/svanna-data.zip" "${REFERENCE_DIR}"
        print_success "Svanna database setup complete!"
    fi

    echo ""
}

download_reference_files() {
    if [ "$SKIP_REFERENCE" = true ]; then
        print_info "Skipping reference file download (--skip-reference specified)"
        echo ""
        return
    fi

    print_header "Downloading Reference Files from Zenodo"
    print_info "Zenodo Record: https://doi.org/10.5281/zenodo.${ZENODO_RECORD}"
    print_info "This may take 20-40 minutes depending on your internet speed"
    echo ""
    print_warning "Total download size: ~30-35 GB (core) + ~20 GB (optional)"
    echo ""

    # Download demo/test data (diana_dummy)
    if [ ! -d "${DATA_DIR}/diana_dummy" ]; then
        print_info "Downloading demo test data (diana_dummy.tar.gz)..."
        download_file "diana_dummy.tar.gz" "${PIPELINE_DIR}/diana_dummy.tar.gz"
        verify_checksum "${PIPELINE_DIR}/diana_dummy.tar.gz"
        extract_archive "${PIPELINE_DIR}/diana_dummy.tar.gz" "${DATA_DIR}"
        rm "${PIPELINE_DIR}/diana_dummy.tar.gz"
        print_success "Demo data extracted to data/diana_dummy/"
        echo ""
    else
        print_success "Demo data (diana_dummy) already present"
        echo ""
    fi

    # Download core reference bundle (includes nanoDx classifier)
    if [ ! -f "${DATA_DIR}/.reference_core_downloaded" ]; then
        print_info "Downloading core reference files (includes nanoDx classifier)..."
        download_file "reference_core.tar.gz" "${PIPELINE_DIR}/reference_core.tar.gz"
        verify_checksum "${PIPELINE_DIR}/reference_core.tar.gz"
        extract_archive "${PIPELINE_DIR}/reference_core.tar.gz" "${REFERENCE_DIR}"
        rm "${PIPELINE_DIR}/reference_core.tar.gz"
        touch "${DATA_DIR}/.reference_core_downloaded"
        echo ""
    else
        print_success "Core reference files already downloaded"
        print_info "(Remove data/.reference_core_downloaded to re-download)"
        echo ""
    fi

    # Download ANNOVAR humandb
    if [ ! -f "${DATA_DIR}/.humandb_downloaded" ]; then
        print_info "Downloading ANNOVAR databases..."
        download_file "humandb.tar.gz" "${PIPELINE_DIR}/humandb.tar.gz"
        verify_checksum "${PIPELINE_DIR}/humandb.tar.gz"
        extract_archive "${PIPELINE_DIR}/humandb.tar.gz" "${HUMANDB_DIR}"
        rm "${PIPELINE_DIR}/humandb.tar.gz"
        touch "${DATA_DIR}/.humandb_downloaded"
        echo ""
    else
        print_success "ANNOVAR databases already downloaded"
        print_info "(Remove data/.humandb_downloaded to re-download)"
        echo ""
    fi

    # Download VEP cache from Ensembl (optional — non-fatal if download fails)
    setup_vep_cache

    # Download general.zip (Sturgeon classifier) - keep as zip, DO NOT extract
    if [ ! -f "${REFERENCE_DIR}/general.zip" ]; then
        print_info "Downloading Sturgeon classifier (general.zip)..."
        download_file "general.zip" "${REFERENCE_DIR}/general.zip"
        print_success "Sturgeon classifier downloaded (kept as .zip)"
        echo ""
    else
        print_success "Sturgeon classifier (general.zip) already present"
        echo ""
    fi

    # Extract Assembly.zip (bundled in reference_core.tar.gz)
    if [ ! -d "${REFERENCE_DIR}/Assembly" ]; then
        if [ -f "${REFERENCE_DIR}/Assembly.zip" ]; then
            print_info "Extracting Assembly data..."
            extract_archive "${REFERENCE_DIR}/Assembly.zip" "${REFERENCE_DIR}"
            print_success "Assembly data extracted"
        else
            print_warning "Assembly.zip not found in reference directory — downloading as fallback..."
            download_file "Assembly.zip" "${REFERENCE_DIR}/Assembly.zip"
            verify_checksum "${REFERENCE_DIR}/Assembly.zip"
            extract_archive "${REFERENCE_DIR}/Assembly.zip" "${REFERENCE_DIR}"
            print_success "Assembly data extracted"
        fi
        echo ""
    else
        print_success "Assembly data already present"
        echo ""
    fi

    # Download and extract Dorado model
    if [ ! -f "${REFERENCE_DIR}/r1041_e82_400bps_sup_v420.zip" ] && [ ! -d "${REFERENCE_DIR}/r1041_e82_400bps_sup_v420" ]; then
        print_info "Downloading Dorado model (r1041_e82_400bps_sup_v420.zip)..."
        download_file "r1041_e82_400bps_sup_v420.zip" "${REFERENCE_DIR}/r1041_e82_400bps_sup_v420.zip"
        print_info "Extracting Dorado model..."
        extract_archive "${REFERENCE_DIR}/r1041_e82_400bps_sup_v420.zip" "${REFERENCE_DIR}"
        print_success "Dorado model extracted"
        echo ""
    else
        print_success "Dorado model already present"
        echo ""
    fi

    # Verify nanoDx classifier (included in reference_core.tar.gz)
    setup_nanodx_classifier

    # Download optional large files (Svanna database, etc.)
    if [ "$SKIP_OPTIONAL" = false ]; then
        echo ""
        print_warning "Optional files include Svanna database (~15-20 GB)"
        print_info "Downloading in 5 seconds... Press Ctrl+C to skip"
        sleep 5

        # Setup Svanna database (large structural variant annotation database)
        setup_svanna_database

        # Download any other optional files if they exist in a bundle
        if [ ! -f "${DATA_DIR}/.optional_downloaded" ]; then
            # Check if reference_optional.tar.gz exists on Zenodo
            # This is for any additional optional files beyond Svanna
            print_info "Checking for additional optional files..."
            # Uncomment below if you have a reference_optional.tar.gz on Zenodo
            # download_file "reference_optional.tar.gz" "${PIPELINE_DIR}/reference_optional.tar.gz"
            # extract_archive "${PIPELINE_DIR}/reference_optional.tar.gz" "${PIPELINE_DIR}"
            # rm "${PIPELINE_DIR}/reference_optional.tar.gz"
            touch "${DATA_DIR}/.optional_downloaded"
            echo ""
        fi
    else
        print_info "Skipping optional files (use without --skip-optional to download)"
        echo ""
    fi

    print_success "All reference files downloaded and extracted!"
    echo ""
}

setup_containers() {
    if [ "$SKIP_CONTAINERS" = true ]; then
        print_info "Skipping container setup (--skip-containers specified)"
        echo ""
        return
    fi

    print_header "Setting Up ${CONTAINER_SYSTEM^} Containers"
    print_info "Pulling/building container images..."
    print_info "This may take 10-30 minutes on first run..."
    echo ""

    if [ "$CONTAINER_SYSTEM" = "docker" ]; then
        setup_docker_containers
    elif [ "$CONTAINER_SYSTEM" = "singularity" ]; then
        setup_singularity_containers
    fi

    print_success "All container images are ready!"
    echo ""
}

# Docker container setup (based on setup_docker.sh)
setup_docker_containers() {
    print_info "Pulling Docker images from vilhelmmagnuslab repository..."

    # Function to pull Docker image if it doesn't exist
    pull_docker_if_not_exists() {
        local image_name=$1
        local image_with_tag="${image_name}:latest"

        if docker images --format "{{.Repository}}:{{.Tag}}" | grep -q "^${image_with_tag}$"; then
            echo -e "  ${GREEN}✓${NC} $image_with_tag already exists, skipping..."
        else
            echo -e "  ${CYAN}Pulling${NC} $image_with_tag..."
            docker pull "$image_with_tag" 2>&1 | grep -v "Pulling from"
            echo -e "  ${GREEN}✓${NC} Pulled $image_with_tag"
        fi
    }

    # Core analysis images
    echo "Core analysis images:"
    pull_docker_if_not_exists "vilhelmmagnuslab/nwgs_default_images"
    pull_docker_if_not_exists "vilhelmmagnuslab/ace_1.24.0"
    pull_docker_if_not_exists "vilhelmmagnuslab/annotcnv_images_27feb1025"
    pull_docker_if_not_exists "vilhelmmagnuslab/clair3_amd64"
    pull_docker_if_not_exists "vilhelmmagnuslab/clairsto_amd64"
    pull_docker_if_not_exists "vilhelmmagnuslab/igv_report_amd64"
    pull_docker_if_not_exists "vilhelmmagnuslab/vcf2circos"
    pull_docker_if_not_exists "vilhelmmagnuslab/nanodx_env"
    pull_docker_if_not_exists "vilhelmmagnuslab/crossnnumap"
    pull_docker_if_not_exists "vilhelmmagnuslab/markdown_images_28feb2025"
    pull_docker_if_not_exists "vilhelmmagnuslab/mgmt_nanopipe_amd64_18feb2025_cramoni"
    pull_docker_if_not_exists "vilhelmmagnuslab/gviz_amd64"
    pull_docker_if_not_exists "vilhelmmagnuslab/sturgeon_amd64_21jan"
    echo ""

    # Epi2me images
    echo "Epi2me analysis images:"
    pull_docker_if_not_exists "vilhelmmagnuslab/snifflesv252_update"
    pull_docker_if_not_exists "vilhelmmagnuslab/qdnaseq_amd64"
    pull_docker_if_not_exists "vilhelmmagnuslab/modkit"
    echo ""
}

# Singularity container setup (based on setup_singularity.sh)
setup_singularity_containers() {
    # Detect Singularity or Apptainer
    SINGULARITY_CMD=""
    if command -v apptainer &> /dev/null; then
        SINGULARITY_CMD="apptainer"
        print_info "Using Apptainer"
    elif command -v singularity &> /dev/null; then
        SINGULARITY_CMD="singularity"
        print_info "Using Singularity"
    else
        print_error "Neither Singularity nor Apptainer found"
        print_info "This should have been caught by prerequisites check"
        exit 1
    fi

    mkdir -p "${PIPELINE_DIR}/containers"

    print_info "Pulling Singularity/Apptainer images from vilhelmmagnuslab repository..."
    echo ""

    # Function to pull image if it doesn't exist
    pull_singularity_if_not_exists() {
        local image_name=$1
        local image_basename=$(basename "$image_name")
        local image_file="${PIPELINE_DIR}/containers/${image_basename}_latest.sif"

        if [ -f "$image_file" ]; then
            echo -e "  ${GREEN}✓${NC} $image_name already exists, skipping..."
        else
            echo -e "  ${CYAN}Pulling${NC} $image_name..."

            $SINGULARITY_CMD pull --dir "${PIPELINE_DIR}/containers/" "docker://$image_name:latest" 2>&1 | \
                grep -v "INFO:" || true

            # Verify the image was downloaded
            if [ -f "$image_file" ]; then
                echo -e "  ${GREEN}✓${NC} Pulled $image_name"
            else
                echo -e "  ${RED}✗${NC} Failed to pull $image_name"
                print_error "Container image download failed: $image_name"
                print_info "Please check your internet connection and Docker Hub access"
                exit 1
            fi
        fi
    }

    # Core analysis images
    echo "Core analysis images:"
    pull_singularity_if_not_exists "vilhelmmagnuslab/nwgs_default_images"
    pull_singularity_if_not_exists "vilhelmmagnuslab/ace_1.24.0"
    pull_singularity_if_not_exists "vilhelmmagnuslab/annotcnv_images_27feb1025"
    pull_singularity_if_not_exists "vilhelmmagnuslab/clair3_amd64"
    pull_singularity_if_not_exists "vilhelmmagnuslab/clairsto_amd64"
    pull_singularity_if_not_exists "vilhelmmagnuslab/igv_report_amd64"
    pull_singularity_if_not_exists "vilhelmmagnuslab/vcf2circos"
    pull_singularity_if_not_exists "vilhelmmagnuslab/nanodx_images_3feb25"
    pull_singularity_if_not_exists "vilhelmmagnuslab/crossnnumap"
    pull_singularity_if_not_exists "vilhelmmagnuslab/markdown_images_28feb2025"
    pull_singularity_if_not_exists "vilhelmmagnuslab/mgmt_nanopipe_amd64_18feb2025_cramoni"
    pull_singularity_if_not_exists "vilhelmmagnuslab/gviz_amd64ps"
    pull_singularity_if_not_exists "vilhelmmagnuslab/sturgeon_amd64_21jan"
    echo ""

    # Epi2me images
    echo "Epi2me analysis images:"
    pull_singularity_if_not_exists "vilhelmmagnuslab/snifflesv252_update"
    pull_singularity_if_not_exists "vilhelmmagnuslab/qdnaseq_amd64"
    pull_singularity_if_not_exists "vilhelmmagnuslab/modkit"
    echo ""
}

install_java() {
    print_header "Checking Java"

    # Check if java exists and extract major version
    local java_ok=false
    if check_command java; then
        local java_ver
        java_ver=$(java -version 2>&1 | head -1 | sed -E 's/.*version "([0-9]+)[._].*/\1/')
        if [[ "$java_ver" =~ ^[0-9]+$ ]] && [ "$java_ver" -ge 11 ] && [ "$java_ver" -le 21 ]; then
            print_success "Java ${java_ver} is compatible with Nextflow (requires 11–21)"
            echo ""
            return
        elif [[ "$java_ver" =~ ^[0-9]+$ ]]; then
            print_warning "Java ${java_ver} found but Nextflow requires Java 11–21 — installing Java 21..."
        else
            print_warning "Could not determine Java version — installing Java 21..."
        fi
    else
        print_warning "Java not found — installing Java 21 (required by Nextflow)..."
    fi

    # Try package managers in order
    if check_command apt-get; then
        print_info "Installing Java 21 via apt..."
        sudo apt-get update -qq
        sudo apt-get install -y openjdk-21-jdk 2>/dev/null || \
        sudo apt-get install -y openjdk-17-jdk
    elif check_command dnf; then
        print_info "Installing Java 21 via dnf..."
        sudo dnf install -y java-21-openjdk 2>/dev/null || \
        sudo dnf install -y java-17-openjdk
    elif check_command yum; then
        print_info "Installing Java 21 via yum..."
        sudo yum install -y java-21-openjdk 2>/dev/null || \
        sudo yum install -y java-17-openjdk
    else
        # Fallback: download Temurin 21 JRE from Adoptium
        print_info "No package manager found — downloading Temurin 21 JRE from Adoptium..."
        local jre_url="https://api.adoptium.net/v3/binary/latest/21/ga/linux/x64/jre/hotspot/normal/eclipse"
        local jre_tar="${PIPELINE_DIR}/temurin21.tar.gz"
        local jre_dir="${PIPELINE_DIR}/jre"

        if check_command curl; then
            curl -L "$jre_url" -o "$jre_tar" --progress-bar
        else
            wget -q --show-progress "$jre_url" -O "$jre_tar"
        fi

        mkdir -p "$jre_dir"
        tar -xzf "$jre_tar" -C "$jre_dir" --strip-components=1
        rm "$jre_tar"

        export JAVA_HOME="$jre_dir"
        export PATH="$jre_dir/bin:$PATH"
        echo "export JAVA_HOME=${jre_dir}" >> "${PIPELINE_DIR}/.diana_env"
        echo "export PATH=${jre_dir}/bin:\$PATH" >> "${PIPELINE_DIR}/.diana_env"
        print_info "Temurin 21 installed to ${jre_dir} — JAVA_HOME saved to .diana_env"
    fi

    # Verify installation succeeded
    if check_command java; then
        local installed_ver
        installed_ver=$(java -version 2>&1 | head -1 | sed -E 's/.*version "([0-9]+)[._].*/\1/')
        print_success "Java ${installed_ver} installed successfully"
    else
        print_error "Java installation failed — please install Java 11–21 manually"
        print_info "  Ubuntu/Debian: sudo apt-get install openjdk-21-jdk"
        print_info "  RHEL/CentOS:   sudo dnf install java-21-openjdk"
        print_info "  Download:      https://adoptium.net"
        exit 1
    fi
    echo ""
}

install_nextflow() {
    print_header "Installing Nextflow"

    if check_command nextflow; then
        print_success "Nextflow already installed: $(nextflow -version 2>&1 | head -1 | xargs)"
        echo ""
        return
    fi

    cd "${PIPELINE_DIR}"
    curl -s https://get.nextflow.io | bash
    chmod +x nextflow

    # Add to PATH via .diana_env so pipeline runner scripts find it without re-installing
    if ! grep -q "nextflow" "${PIPELINE_DIR}/.diana_env" 2>/dev/null; then
        echo "export PATH=${PIPELINE_DIR}:\$PATH" >> "${PIPELINE_DIR}/.diana_env"
    fi

    print_success "Nextflow installed to ${PIPELINE_DIR}/nextflow"
    print_info "PATH updated in .diana_env — source it or open a new shell to use 'nextflow' directly"
    echo ""
}

validate_setup() {
    print_header "Validating Setup"

    if [ -f "${PIPELINE_DIR}/validate_setup.sh" ]; then
        bash "${PIPELINE_DIR}/validate_setup.sh"
    else
        print_warning "validate_setup.sh not found, skipping validation"

        # Basic validation
        print_info "Performing basic validation..."

        local missing_files=()

        # Check key reference files
        [ ! -f "${REFERENCE_DIR}/GRCh38.fa" ] && missing_files+=("GRCh38.fa")
        [ ! -f "${REFERENCE_DIR}/roi.protein_coding.bed" ] && missing_files+=("roi.protein_coding.bed")
        [ ! -f "${REFERENCE_DIR}/CNV_genes_tuned.csv" ] && missing_files+=("CNV_genes_tuned.csv")

        # Check humandb
        [ ! -f "${HUMANDB_DIR}/hg38_refGene.txt" ] && missing_files+=("hg38_refGene.txt")

        if [ ${#missing_files[@]} -gt 0 ]; then
            print_warning "Some files appear to be missing:"
            for file in "${missing_files[@]}"; do
                echo "  - $file"
            done
            echo ""
            print_info "This might be normal if files are named differently or optional"
        else
            print_success "Key reference files found!"
        fi
    fi

    echo ""
}

create_quick_start_guide() {
    print_header "Creating Quick Start Guide"

    cat > "${PIPELINE_DIR}/QUICKSTART.md" << 'EOF'
# Diana Pipeline - Quick Start Guide

## 🎉 Setup Complete!

Your Diana pipeline is ready to use. All reference files have been downloaded from Zenodo and containers are configured.

---

## Running the Pipeline

### 1. Prepare Your Sample Files

Create a `sample_ids.txt` file with your sample IDs:

```bash
echo "your_sample_id" > data/sample_ids.txt
```

### 2. Choose a Pipeline Mode

#### Option A: Full Pipeline (Order Mode - Recommended)

Runs the complete pipeline sequentially:

```bash
# Docker
./run_pipeline_docker.sh --run_mode_order

# Singularity
./run_pipeline_singularity.sh --run_mode_order
```

#### Option B: Specific Modules

**Epi2me module only** (methylation, CNV, SV calling):
```bash
./run_pipeline_docker.sh --run_mode_epi2me
```

**Annotation module only** (SNV, MGMT, TERTp, etc.):
```bash
./run_pipeline_docker.sh --run_mode_annotation tertp
```

**Combined Epi2me + Annotation** (epiannotation mode):
```bash
./run_pipeline_docker.sh --run_mode_epiannotation
```

**Merge BAM files**:
```bash
./run_pipeline_docker.sh --run_mode_mergebam
```

### 3. Monitor Progress

Logs are automatically saved to `logs/` directory with sample IDs and timestamps:

- `logs/trace-{mode}_{sample_id}_{timestamp}.txt` - Execution trace
- `logs/execution_report-{mode}_{sample_id}_{timestamp}.html` - Interactive report
- `logs/execution_timeline-{mode}_{sample_id}_{timestamp}.html` - Timeline visualization

**View reports:**
```bash
# Open the latest execution report
firefox logs/execution_report-*.html

# Or use your preferred browser
google-chrome logs/execution_report-*.html
```

### 4. View Results

Results are organized in:

```
routine_results/{sample_id}/
├── {sample_id}_markdown_pipeline_report.pdf    # Main PDF report
├── CNV/                                         # Copy number variations
├── SNV/                                         # Single nucleotide variants
├── SV/                                          # Structural variants
├── methylation/                                 # Methylation analysis
└── ...
```

---

## Configuration

Edit configuration files in `conf/` to customize pipeline behavior:

- **`conf/epi2me.config`** - Epi2me module (methylation, CNV, SV)
- **`conf/annotation.config`** - Annotation module (SNV filtering, reports)
- **`conf/mergebam.config`** - BAM merging settings

### Key Parameters

**SNV Filtering** (in `conf/annotation.config`):
```groovy
params {
    snv_depth_threshold = 10    // Minimum sequencing depth
    snv_gq_threshold = 10       // Minimum genotype quality
}
```

**Working Directories** (configure paths for your system):
```groovy
params {
    path = "/path/to/Diana"
    path_output = "/path/to/routine_diana"
}
```

---

## Example: Complete Workflow

```bash
# 1. Setup (already done!)
# ./setup_pipeline.sh docker

# 2. Prepare sample
echo "T25-000" > data/sample_ids.txt

# 3. Run full pipeline
./run_pipeline_docker.sh --run_mode_order

# 4. Wait for completion (monitor logs in real-time)
tail -f logs/trace-order_T25-000_*.txt

# 5. Check results
ls -lh routine_results/T25-000/

# 6. View PDF report
evince routine_results/T25-000/T25-000_markdown_pipeline_report.pdf
```

---

## Advanced Options

### Custom Work Directory

Specify a custom temporary work directory (useful for large datasets):

```bash
./run_pipeline_docker.sh --run_mode_order -w /path/to/work/dir
```

### Resume Failed Runs

Nextflow automatically caches completed processes. Resume from where it stopped:

```bash
./run_pipeline_docker.sh --run_mode_order -resume
```

---

## Troubleshooting

### Pipeline Fails?

1. **Check logs** in `logs/` directory
2. **Run validation**: `./validate_setup.sh`
3. **Check configuration** files in `conf/`
4. **Review error** in trace file: `logs/trace-*.txt`

### Common Issues

**"Container not found"**
- Re-run setup: `./setup_pipeline.sh docker` (or singularity)

**"Reference file not found"**
- Verify files in `data/reference/`
- Re-download if needed: `rm data/.reference_core_downloaded && ./setup_pipeline.sh docker --skip-containers`

**"Out of disk space"**
- Clear work directory: `rm -rf work/`
- Use custom work directory: `-w /path/to/large/disk`

**"Permission denied"**
- Ensure scripts are executable: `chmod +x run_pipeline_*.sh`

---

## Automated Monitoring

For production environments, use the automated monitoring script:

```bash
# Monitor ONT sequencing runs automatically
./smart_sample_monitor_v2.sh
```

This script automatically:
- Detects new samples
- Runs the pipeline
- Monitors progress
- Handles timeouts
- Cleans up on completion

---

## Getting Help

- **Full documentation**: `README.md`
- **Software versions**: `SOFTWARE_VERSIONS.md`
- **Container information**: `CONTAINERS.md`
- **Changelog**: `CHANGELOG.md`
- **GitHub issues**: https://github.com/VilhelmMagnusLab/Diana/issues

---

## Next Steps

1. ✅ Read the full `README.md` for advanced features
2. ✅ Customize configuration files for your needs
3. ✅ Run a test sample to verify everything works
4. ✅ Set up automated monitoring if needed

**Happy analyzing! 🧬**
EOF

    print_success "Created QUICKSTART.md"
    echo ""
}

print_completion_message() {
    echo ""
    echo -e "${GREEN}╔══════════════════════════════════════════════════════════╗${NC}"
    echo -e "${GREEN}║                                                          ║${NC}"
    echo -e "${GREEN}║      🎉  Setup Complete! Pipeline Ready to Use  🎉       ║${NC}"
    echo -e "${GREEN}║                                                          ║${NC}"
    echo -e "${GREEN}╚══════════════════════════════════════════════════════════╝${NC}"
    echo ""

    if [ "$SKIP_REFERENCE" = false ]; then
        echo -e "${BLUE}📂 Reference files installed in:${NC}"
        echo "   └─ ${DATA_DIR}"
        echo ""
    fi

    if [ "$SKIP_CONTAINERS" = false ]; then
        echo -e "${BLUE}🐳 Container system configured:${NC}"
        echo "   └─ ${CONTAINER_SYSTEM^}"
        echo ""
    fi

    echo -e "${BLUE}🚀 Quick Start:${NC}"
    echo ""
    echo "   1. Read the quick start guide:"
    echo -e "      ${GREEN}cat QUICKSTART.md${NC}"
    echo ""
    echo "   2. Prepare your sample:"
    echo -e "      ${GREEN}echo 'your_sample_id' > data/sample_ids.txt${NC}"
    echo ""
    echo "   3. Run the pipeline:"
    if [ "$CONTAINER_SYSTEM" = "docker" ]; then
        echo -e "      ${GREEN}./run_pipeline_docker.sh --run_mode_order${NC}"
    else
        echo -e "      ${GREEN}./run_pipeline_singularity.sh --run_mode_order${NC}"
    fi
    echo ""
    echo -e "${BLUE}📖 Documentation:${NC}"
    echo -e "   ├─ Quick start: ${GREEN}QUICKSTART.md${NC}"
    echo -e "   ├─ Full guide: ${GREEN}README.md${NC}"
    echo -e "   ├─ Software versions: ${GREEN}SOFTWARE_VERSIONS.md${NC}"
    echo -e "   └─ Containers: ${GREEN}CONTAINERS.md${NC}"
    echo ""
    echo -e "${BLUE}🔍 Validation:${NC}"
    echo -e "   └─ Run: ${GREEN}./validate_setup.sh${NC}"
    echo ""
    print_success "You're all set! Happy analyzing! 🧬"
    echo ""
}

################################################################################
# Main Execution
################################################################################

update_options_json() {
    local options_file="${REFERENCE_DIR}/options.json"
    if [ -f "$options_file" ]; then
        sed -i "s|\"Static\": \".*\"|\"Static\": \"${REFERENCE_DIR}\"|" "$options_file"
        print_success "Updated options.json Static path to: ${REFERENCE_DIR}"
    else
        print_info "options.json not found yet — will be updated after reference download"
    fi
}

main() {
    clear

    echo -e "${BLUE}"
    cat << "EOF"
╔═══════════════════════════════════════════════════════════════╗
║                                                               ║
║                    DIANA-Automated Setup                      ║
║    Diagnostic Integrated Analytics for Neoplastic Alterations ║
║                                                               ║
║           Zenodo: 10.5281/zenodo.20596458                     ║
║                                                               ║
╚═══════════════════════════════════════════════════════════════╝
EOF
    echo -e "${NC}"

    print_info "Container system: ${CONTAINER_SYSTEM^}"
    echo ""
    print_info "This script will automatically:"
    echo "  1. Download all reference files from Zenodo (~30-50 GB)"
    echo "  2. Extract and organize files in correct directories"
    echo "  3. Set up ${CONTAINER_SYSTEM^} containers"
    echo "  4. Validate the installation"
    echo "  5. Create quick start guide"
    echo ""
    print_warning "Estimated time: 30-60 minutes (depends on internet speed)"
    echo ""

    read -p "Continue with setup? (y/n) " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        print_info "Setup cancelled by user"
        exit 0
    fi

    echo ""

    # Run setup steps
    check_prerequisites
    create_directories
    install_java
    install_nextflow
    download_reference_files
    update_options_json
    setup_containers
    validate_setup
    create_quick_start_guide
    print_completion_message
}

# Parse arguments and run
parse_args "$@"
main
