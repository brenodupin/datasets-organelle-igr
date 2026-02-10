#!/usr/bin/env bash

set -e
project_dir="$(dirname "$(dirname "$(realpath "${BASH_SOURCE[0]}")")")"
igs_path="$project_dir/scripts/summary_igs.py"
echo "Project directory: $project_dir"

# check if tigre command is available
if ! command -v tigre &> /dev/null
then
    echo "tigre command could not be found. Please install tigre and ensure it is in your PATH."
    echo 'Installation using pip: pip install "tigre[all]"'
    exit 1
fi

process_group() {
    local group="$1"
    echo "## --$group"
    
    GROUP_START=$(date +%s)
    
    # Build paths
    local base_path="$project_dir/$group"
    local gdict_file="$base_path/${group}.gdict"
    local tsv_file="$base_path/$group.tsv"
    local igs_file="$base_path/summary_igs_intergenic.tsv"
    
    # Check if directory exists
    if [ ! -d "$base_path" ]; then
        echo "Warning: Directory $base_path does not exist. Skipping $group."
        return
    fi
    
    # Tigre commands to extract the intergenic regions
    tigre clean multiple -v --log "$base_path/clean_brms.log" --gdict "$gdict_file" --tsv "$tsv_file" --overwrite
    tigre extract multiple -v --log "$base_path/extract_brms.log" --tsv "$tsv_file" --overwrite
    tigre getfasta multiple -v --log "$base_path/getfasta_brms.log" --tsv "$tsv_file" --overwrite

    # Summarize intergenic regions
    python3 "$igs_path" --tsv "$tsv_file" --output "$igs_file" --overwrite
    
    GROUP_END=$(date +%s)
    echo "$group group completed in $((GROUP_END - GROUP_START)) seconds"
    echo
}

echo "========================================="
echo "Starting IGR extraction"
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
echo "IGR extraction completed"
echo "========================================="