#!/usr/bin/env bash
set -eu
project_dir="$(dirname "$(dirname "$(realpath "${BASH_SOURCE[0]}")")")"

###############################################################
# IF YOUR DATA DIRECTORY MOVES, EDIT THIS LINE.
# Every group path is built from it.
base=${project_dir}
# Default points to the dataset-organelle-igr directory.
###############################################################

igs_path=${project_dir}/scripts/summary_igs.py

tnow() {
    date '+%Y-%m-%d %H:%M:%S'
}

# check if tigre command is available
if ! command -v tigre &> /dev/null; then
    echo "tigre command could not be found. Please install tigre and ensure it is in your PATH."
    echo 'Installation using pip: pip install "tigre[all]"'
    exit 1
fi

process_group() {
    local GROUP="$1"

    local base_path="${base}/${GROUP}"
    local gdict_file="${base_path}/${GROUP}.gdict"
    local tsv_file="${base_path}/${GROUP}.tsv"
    local igs_file="${base_path}/${GROUP}_summary_igr.tsv"

    if [ ! -d "$base_path" ]; then
        echo "[$(tnow)] Warning: directory $base_path does not exist. Skipping $GROUP."
        return
    fi

    local start elapsed
    start=$(date +%s)

    echo "[$(tnow)] Starting IGR extraction for $GROUP..."

    tigre clean multiple -v --log "${base_path}/clean_brms.log" --gdict "$gdict_file" --tsv "$tsv_file" --overwrite
    tigre extract multiple -v --log "${base_path}/extract_brms.log" --tsv "$tsv_file" --overwrite
    tigre getfasta multiple -v --log "${base_path}/getfasta_brms.log" --tsv "$tsv_file" --overwrite

    python3 "$igs_path" --tsv "$tsv_file" --output "$igs_file" --overwrite

    elapsed=$(( $(date +%s) - start ))
    echo "[$(tnow)] Finished IGR extraction for $GROUP in $(( elapsed / 60 )) min $(( elapsed % 60 )) sec"
    echo
}

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
    echo "[$(tnow)] No groups specified. Processing all groups."
else
    groups=("$@")
    echo "[$(tnow)] Processing specified groups: ${groups[*]}"
fi

start_time=$(date +%s)
echo "[$(tnow)] Starting IGR extraction for all groups..."
echo

for group in "${groups[@]}"; do
    process_group "$group"
done

elapsed=$(( $(date +%s) - start_time ))
echo "[$(tnow)] Finished all IGR extraction in $(( elapsed / 3600 )) hr $(( (elapsed % 3600) / 60 )) min $(( elapsed % 60 )) sec"
echo "Done!"