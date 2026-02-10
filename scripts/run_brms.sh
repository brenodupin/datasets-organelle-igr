#!/usr/bin/env bash

set -e
project_dir="$(dirname "$(dirname "$(realpath "${BASH_SOURCE[0]}")")")"
brms_path="$project_dir/scripts/brms.py"
echo "Project directory: $project_dir"

process_group() {
    local group="$1"
    echo "## --$group"
    
    GROUP_START=$(date +%s)
    
    # Build paths
    local base_path="$project_dir/$group"
    local tsv_file="$base_path/$group.tsv"
    local igs_file="$base_path/summary_igs_intergenic.tsv"
    local output_dir="$base_path/brms_result"

    
    if [ ! -f "$tsv_file" ]; then
        echo "Warning: File $tsv_file does not exist. Skipping $group."
        return
    fi
    
    if [ ! -f "$igs_file" ]; then
        echo "Warning: File $igs_file does not exist. Skipping $group."
        return
    fi
    
    # Execute the brms analysis on the intergenic regions
    python3 "$brms_path" --tsv "$tsv_file" --igs "$igs_file" --output "$output_dir" --overwrite
    
    GROUP_END=$(date +%s)
    echo "$group group completed in $((GROUP_END - GROUP_START)) seconds"
    echo
}

echo "========================================="
echo "Starting BRMS analysis"
echo "========================================="

# Default groups (all)
default_groups=(
    "fungi_mit"
    "green_algae_mit"
    "protists_mit"
    "metazoans_mit"
    "plants_mit"
    "green_algae_plt"
    "plants_plt"
    "protists_plt"
)

# Use provided arguments or default to all groups
if [ $# -eq 0 ]; then
    groups=("${default_groups[@]}")
    echo "No groups specified. Processing all groups."
else
    groups=("$@")
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