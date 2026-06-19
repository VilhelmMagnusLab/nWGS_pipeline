#!/bin/bash
################################################################################
# Package Diana Pipeline Files for Zenodo Upload
################################################################################
# This script packages all reference files for upload to Zenodo record 17589248.
# It creates all necessary archives in the correct format.
#
# Source directories:
#   - /home/godzilla/Diana/data/reference/
#   - /home/godzilla/Diana/data/humandb/
#
# Usage:
#   ./package_for_zenodo.sh [output_directory]
#
# Example:
#   ./package_for_zenodo.sh              # Creates ./zenodo_upload/
#   ./package_for_zenodo.sh ./my_upload  # Creates ./my_upload/
#
################################################################################

set -e

# Color codes
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
NC='\033[0m'

# Configuration
PIPELINE_DIR="/home/godzilla/diana"
DATA_DIR="${PIPELINE_DIR}/data"
REFERENCE_DIR="${DATA_DIR}/reference"
HUMANDB_DIR="${DATA_DIR}/humandb"
OUTPUT_DIR="${1:-${PIPELINE_DIR}/zenodo_upload}"

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

human_readable_size() {
    local bytes=$1
    if [ $bytes -lt 1024 ]; then
        echo "${bytes}B"
    elif [ $bytes -lt 1048576 ]; then
        echo "$(( bytes / 1024 ))KB"
    elif [ $bytes -lt 1073741824 ]; then
        echo "$(( bytes / 1048576 ))MB"
    else
        echo "$(awk "BEGIN {printf \"%.1f\", $bytes/1073741824}")GB"
    fi
}

check_file_exists() {
    local file=$1
    local description=$2

    if [ ! -e "$file" ]; then
        print_error "${description} not found: $file"
        return 1
    fi
    return 0
}

################################################################################
# Main Packaging Functions
################################################################################

create_output_directory() {
    print_header "Creating Output Directory"

    mkdir -p "$OUTPUT_DIR"
    print_success "Output directory: $OUTPUT_DIR"
    echo ""
}

package_reference_core() {
    print_header "Packaging Reference Core Files"

    local output_file="${OUTPUT_DIR}/reference_core.tar.gz"

    # Check if reference directory exists
    if ! check_file_exists "${REFERENCE_DIR}" "Reference directory"; then
        return 1
    fi

    print_info "Creating reference_core.tar.gz from ${REFERENCE_DIR}..."
    print_warning "This may take 5-10 minutes for ~25 GB of files..."

    # Package reference files.
    # Tar from DATA_DIR with 'reference' as the entry (not '.' from inside
    # REFERENCE_DIR) so the archive has a top-level reference/ directory —
    # setup_pipeline.sh extracts this archive into DATA_DIR, expecting that
    # top-level directory to land at data/reference/.
    # nanoDx is INCLUDED in reference_core.tar.gz since it's already in data/reference/nanoDx/
    # Assembly.zip and svanna-data.zip ARE bundled (kept as zips so users can
    # extract on demand); their already-extracted directory counterparts are
    # excluded to avoid shipping the same content twice. general.zip
    # (Sturgeon) is distributed separately, not via Zenodo. The deprecated
    # v420 Dorado model is excluded in favor of the bundled v520 sup/hac models.
    tar -czf "$output_file" \
        -C "${DATA_DIR}" \
        --exclude='reference/general.zip' \
        --exclude='reference/Assembly' \
        --exclude='reference/svanna-data' \
        --exclude='reference/r1041_e82_400bps_sup_v420' \
        --exclude='reference/r1041_e82_400bps_sup_v420.zip' \
        reference

    local size=$(stat -f%z "$output_file" 2>/dev/null || stat -c%s "$output_file")
    local size_human=$(human_readable_size $size)

    print_success "Created reference_core.tar.gz ($size_human)"
    echo ""
}

package_humandb() {
    print_header "Packaging ANNOVAR Databases"

    local output_file="${OUTPUT_DIR}/humandb.tar.gz"

    if ! check_file_exists "${HUMANDB_DIR}" "humandb directory"; then
        return 1
    fi

    print_info "Creating humandb.tar.gz from ${HUMANDB_DIR}..."
    print_warning "This may take 3-5 minutes for ~10 GB of files..."

    # Tar from DATA_DIR with 'humandb' as the entry so the archive has a
    # top-level humandb/ directory (see package_reference_core comment above).
    # homo_sapiens_refseq/ (the VEP cache, ~24 GB) is excluded — it's
    # downloaded directly from Ensembl by setup_pipeline.sh instead of
    # being bundled here.
    tar -czf "$output_file" \
        -C "${DATA_DIR}" \
        --exclude='humandb/homo_sapiens_refseq' \
        humandb

    local size=$(stat -f%z "$output_file" 2>/dev/null || stat -c%s "$output_file")
    local size_human=$(human_readable_size $size)

    print_success "Created humandb.tar.gz ($size_human)"
    echo ""
}

copy_general_zip() {
    print_header "Copying Sturgeon Classifier (general.zip)"

    local source="${REFERENCE_DIR}/general.zip"
    local dest="${OUTPUT_DIR}/general.zip"

    if check_file_exists "$source" "general.zip"; then
        cp "$source" "$dest"

        local size=$(stat -f%z "$dest" 2>/dev/null || stat -c%s "$dest")
        local size_human=$(human_readable_size $size)

        print_success "Copied general.zip ($size_human)"
        print_warning "Remember: general.zip should NOT be extracted (Sturgeon expects zip format)"
    else
        print_warning "general.zip not found - skipping"
        print_info "If you need it, download from previous Zenodo version"
    fi
    echo ""
}

create_assembly_zip() {
    print_header "Handling Assembly.zip"

    local source_zip="${REFERENCE_DIR}/Assembly.zip"
    local source_dir="${REFERENCE_DIR}/Assembly"
    local output_file="${OUTPUT_DIR}/Assembly.zip"

    # Try to copy existing zip first
    if [ -f "$source_zip" ]; then
        cp "$source_zip" "$output_file"

        local size=$(stat -f%z "$output_file" 2>/dev/null || stat -c%s "$output_file")
        local size_human=$(human_readable_size $size)

        print_success "Copied Assembly.zip ($size_human)"
    # Otherwise create from directory
    elif [ -d "$source_dir" ]; then
        print_info "Creating Assembly.zip from ${source_dir}..."

        cd "${REFERENCE_DIR}"
        zip -r "$output_file" Assembly/ -q
        cd "${PIPELINE_DIR}"

        local size=$(stat -f%z "$output_file" 2>/dev/null || stat -c%s "$output_file")
        local size_human=$(human_readable_size $size)

        print_success "Created Assembly.zip ($size_human)"
    else
        print_warning "Assembly not found (neither zip nor directory) - skipping"
        print_info "If you need it, download from previous Zenodo version"
    fi
    echo ""
}

copy_or_create_dorado_model() {
    print_header "Handling Dorado Model"

    local source_zip="${REFERENCE_DIR}/r1041_e82_400bps_sup_v420.zip"
    local source_dir="${REFERENCE_DIR}/r1041_e82_400bps_sup_v420"
    local dest="${OUTPUT_DIR}/r1041_e82_400bps_sup_v420.zip"

    # Try to copy existing zip first
    if [ -f "$source_zip" ]; then
        cp "$source_zip" "$dest"

        local size=$(stat -f%z "$dest" 2>/dev/null || stat -c%s "$dest")
        local size_human=$(human_readable_size $size)

        print_success "Copied r1041_e82_400bps_sup_v420.zip ($size_human)"
    # Otherwise create from directory
    elif [ -d "$source_dir" ]; then
        print_info "Creating r1041_e82_400bps_sup_v420.zip from directory..."

        cd "${REFERENCE_DIR}"
        zip -r "$dest" r1041_e82_400bps_sup_v420/ -q
        cd "${PIPELINE_DIR}"

        local size=$(stat -f%z "$dest" 2>/dev/null || stat -c%s "$dest")
        local size_human=$(human_readable_size $size)

        print_success "Created r1041_e82_400bps_sup_v420.zip ($size_human)"
    else
        print_warning "Dorado model not found (neither zip nor directory) - skipping"
        print_info "If you need it, download from previous Zenodo version"
    fi
    echo ""
}

create_svanna_zip() {
    print_header "Creating Svanna Database (Optional)"

    local source_zip="${REFERENCE_DIR}/svanna-data.zip"
    local source_dir="${REFERENCE_DIR}/svanna-data"
    local dest="${OUTPUT_DIR}/svanna-data.zip"

    # Try to copy existing zip first
    if [ -f "$source_zip" ]; then
        print_info "Copying existing svanna-data.zip..."
        cp "$source_zip" "$dest"

        local size=$(stat -f%z "$dest" 2>/dev/null || stat -c%s "$dest")
        local size_human=$(human_readable_size $size)

        print_success "Copied svanna-data.zip ($size_human)"
    # Otherwise create from directory
    elif [ -d "$source_dir" ]; then
        print_info "Creating svanna-data.zip from directory..."
        print_warning "This is LARGE (~15-20 GB) and may take 10-20 minutes"

        cd "${REFERENCE_DIR}"
        zip -r "$dest" svanna-data/ -q
        cd "${PIPELINE_DIR}"

        local size=$(stat -f%z "$dest" 2>/dev/null || stat -c%s "$dest")
        local size_human=$(human_readable_size $size)

        print_success "Created svanna-data.zip ($size_human)"
    else
        print_warning "svanna-data not found - skipping optional file"
        print_info "This file is optional. Users can skip with --skip-optional"
    fi
    echo ""
}

verify_packages() {
    print_header "Verifying Packaged Files"

    echo "Files in ${OUTPUT_DIR}:"
    echo ""

    local total_size=0
    local file_count=0

    for file in "${OUTPUT_DIR}"/*; do
        if [ -f "$file" ]; then
            local filename=$(basename "$file")
            local size=$(stat -f%z "$file" 2>/dev/null || stat -c%s "$file")
            local size_human=$(human_readable_size $size)

            echo -e "  ${GREEN}✓${NC} $filename ($size_human)"
            total_size=$((total_size + size))
            ((file_count++))
        fi
    done

    echo ""
    local total_human=$(human_readable_size $total_size)
    echo -e "  ${CYAN}Files: $file_count${NC}"
    echo -e "  ${CYAN}Total size: $total_human${NC}"
    echo ""
}

print_next_steps() {
    print_header "Next Steps"

    echo "Files are ready for Zenodo upload to record 17589248!"
    echo ""
    echo -e "${CYAN}Option 1: Automated Upload (Recommended)${NC}"
    echo "  1. Get your Zenodo API token from:"
    echo "     https://zenodo.org/account/settings/applications/tokens/new/"
    echo "     Scopes needed: deposit:write and deposit:actions"
    echo ""
    echo "  2. Run the upload script:"
    echo -e "     ${GREEN}export ZENODO_TOKEN='your_token_here'${NC}"
    echo -e "     ${GREEN}./upload_to_zenodo.sh --token \"\$ZENODO_TOKEN\" --files-dir $OUTPUT_DIR${NC}"
    echo ""
    echo -e "${CYAN}Option 2: Manual Web Upload${NC}"
    echo "  1. Go to: https://zenodo.org/records/17589248"
    echo "  2. Click 'New version'"
    echo "  3. Delete old files"
    echo "  4. Upload files from: $OUTPUT_DIR"
    echo "  5. Update metadata and publish"
    echo ""
    echo -e "${YELLOW}⚠️  Important Notes:${NC}"
    echo "  - general.zip stays as .zip (DO NOT extract)"
    echo "  - Other .zip files will be extracted by setup_pipeline.sh"
    echo "  - Test on Zenodo sandbox first if unsure!"
    echo "  - Upload adds --sandbox flag: ./upload_to_zenodo.sh --sandbox ..."
    echo ""
}

################################################################################
# Main Execution
################################################################################

main() {
    clear

    echo -e "${BLUE}"
    cat << "EOF"
╔═══════════════════════════════════════════════════════════════╗
║                                                               ║
║        Diana Pipeline - Package Files for Zenodo               ║
║          https://zenodo.org/records/17589248                  ║
║                                                               ║
╚═══════════════════════════════════════════════════════════════╝
EOF
    echo -e "${NC}"

    print_info "This script will package all reference files for Zenodo upload"
    echo ""
    print_info "Source directories:"
    echo "  - Reference: ${REFERENCE_DIR}"
    echo "  - ANNOVAR:   ${HUMANDB_DIR}"
    echo ""
    print_info "Output directory: $OUTPUT_DIR"
    echo ""

    read -p "Continue? (y/n) " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        print_info "Cancelled"
        exit 0
    fi

    echo ""

    # Run packaging steps
    create_output_directory
    package_reference_core      # ~25 GB, 5-10 min
    package_humandb             # ~10 GB, 3-5 min
    copy_general_zip            # ~3 GB, instant
    create_assembly_zip         # ~2 GB, 1-2 min
    copy_or_create_dorado_model # ~2 GB, instant or 1 min
    create_svanna_zip           # ~15 GB, 10-20 min (optional)

    # Verify and show results
    verify_packages
    print_next_steps

    print_success "Packaging complete! 🎉"
    echo ""
}

main
