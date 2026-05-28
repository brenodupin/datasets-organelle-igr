#!/usr/bin/env bash

set -e
project_dir="$(dirname "$(dirname "$(realpath "${BASH_SOURCE[0]}")")")"
brms_path="$project_dir/scripts/brms.py"


process_group() {
    local group="$1"
    echo "## --$group"
    
    GROUP_START=$(date +%s)
    
    # Build paths
    local base_path="$project_dir/$group"
    local tsv_file="$base_path/$group.tsv"
    local igs_file="$base_path/summary_igs_intergenic.tsv"

    
    if [ ! -f "$tsv_file" ]; then
        echo "Warning: File $tsv_file does not exist. Skipping $group."
        return
    fi
    
    if [ ! -f "$igs_file" ]; then
        echo "Warning: File $igs_file does not exist. Skipping $group."
        return
    fi
    
    # Execute the brms analysis on the intergenic regions
    python3 "$brms_path" --tsv "$tsv_file" --igs "$igs_file" "${passthrough_flags[@]}"
    
    GROUP_END=$(date +%s)
    echo "$group group completed in $((GROUP_END - GROUP_START)) seconds"
    echo
}

# Parse command-line arguments
passthrough_flags=()
positional_args=()

while [[ $# -gt 0 ]]; do
    case "$1" in
        --2bin|--verbose|--new-tree|--update-taxa|--overwrite)
            passthrough_flags+=("$1")
            ;;
        --output|--an-column)
            passthrough_flags+=("$1" "$2")
            shift
            ;;
        --help)
            echo "This script is a wrapper to automate multiple runs of the BRMS analysis for different groups."
            echo "Printing $brms_path help message, do note that --tsv and --igs are fixed by this wrapper and should not be provided by the user."
            python3 "$brms_path" --help
            exit 0
            ;;
        *)
            positional_args+=("$1")
            ;;
    esac
    shift
done

echo "========================================="
echo "Starting BRMS analysis"
echo "========================================="
echo "Project directory: $project_dir"

[ ${#passthrough_flags[@]} -gt 0 ] && echo "Passthrough flags: ${passthrough_flags[*]}"

# Default groups (all)
default_groups=(
    "fungi_mit"
    "green_algae_mit"
    "protists_mit"
    "plants_mit"
    "green_algae_plt"
    "plants_plt"
    "protists_plt"
    "metazoans_mit"
)

# Use provided arguments or default to all groups
if [ ${#positional_args[@]} -eq 0 ]; then
    groups=("${default_groups[@]}")
    echo "No groups specified. Processing all groups."
else
    groups=("${positional_args[@]}")
    echo "Processing specified groups: ${groups[*]}"
fi

all_start=$(date +%s)

for group in "${groups[@]}"; do
    process_group "$group"
done

all_end=$(date +%s)
total_time=$((all_end - all_start))
hours=$((total_time / 3600))
minutes=$(((total_time % 3600) / 60))
seconds=$((total_time % 60))
echo "Total processing time: ${hours} hr, ${minutes} min, ${seconds} sec"

echo "========================================="
echo "BRMS analysis completed"
echo "========================================="