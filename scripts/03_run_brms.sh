#!/usr/bin/env bash
set -eu
project_dir="$(dirname "$(dirname "$(realpath "${BASH_SOURCE[0]}")")")"

###############################################################
# IF YOU DATA DIRECTORY MOVES, EDIT THIS LINE.
# Every group path is built from it.
base=${project_dir}
# Default points to the dataset-organelle-igr directory.
###############################################################

R_pol=${project_dir}/scripts/brms_polarity.R
R_type_length=${project_dir}/scripts/brms_type_length.R
R_sigma_icc=${project_dir}/scripts/brms_sigma_icc.R
R_contrast=${project_dir}/scripts/brms_contrast.R 
R_divergence=${project_dir}/scripts/brms_divergence.R

brms_folder_polarity="brms_polarity"
brms_folder_type_length="brms_type_length"
failures=()

tnow() {
    date '+%Y-%m-%d %H:%M:%S'
}

run_brms() {
    local NAME="$1"
    local GROUP="$2"
    local BRMS_FOLDER="$3"
    local R_PATH="$4"
    local RESULTS_APPEND="$5"

    local tree="${base}/${GROUP}/${BRMS_FOLDER}/tree.nwk"
    local igr="${base}/${GROUP}/${BRMS_FOLDER}/filtered_3bin.tsv"
    local models="${base}/${GROUP}/${BRMS_FOLDER}/brms_model.rds"
    local results_txt="${base}/${GROUP}/${BRMS_FOLDER}/brms_${RESULTS_APPEND}.txt"
    local results_tsv="${base}/${GROUP}/${BRMS_FOLDER}/brms_${RESULTS_APPEND}_row.tsv"

    local start elapsed status
    start=$(date +%s)

    echo "[$(tnow)] Starting $NAME for $GROUP..."
    if Rscript "$R_PATH" "$tree" "$igr" "$results_txt" "$models" "$results_tsv" 2>&1; then
        status="Finished"
    else
        status="FAILED"
        failures+=("$GROUP/$NAME")
    fi

    elapsed=$(( $(date +%s) - start ))
    echo "[$(tnow)] $status $NAME for $GROUP in $(( elapsed / 60 )) min $(( elapsed % 60 )) sec"
    sleep 2
}

run_divergence() {
    local GROUP="$1"
    local POLARITY_FOLDER="$2"
    local TYPE_LENGTH_FOLDER="$3"
    local R_PATH="$4"
    local RESULTS_APPEND="$5"

    local polarity_dir="${base}/${GROUP}/${POLARITY_FOLDER}"
    local type_length_dir="${base}/${GROUP}/${TYPE_LENGTH_FOLDER}"
    local polarity_tsv="${polarity_dir}/brms_${RESULTS_APPEND}_row.tsv"
    local type_length_tsv="${type_length_dir}/brms_${RESULTS_APPEND}_row.tsv"

    echo "[$(tnow)] Starting Divergence for $GROUP..."
    if Rscript "$R_PATH" "$polarity_dir" "$polarity_tsv" "$type_length_dir" "$type_length_tsv" ; then
        echo "[$(tnow)] Finished Divergence for $GROUP."
    else
        echo "[$(tnow)] FAILED Divergence for $GROUP."
        failures+=("$GROUP/Divergence")
    fi
}

# Default groups (all)
default_groups=(
    "fungi_mit"
    "green_algae_mit"
    "green_algae_plt"
    "protists_mit"
    "protists_plt"
    "plants_mit"
    "plants_plt"
    "metazoans_mit"
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
echo "[$(tnow)] Starting BRMS runs for all groups..."
echo

for group in "${groups[@]}"; do
    echo "## --$group"
    run_brms "Polarity"    "$group" "$brms_folder_polarity"    "$R_pol"         "polarity"
    run_brms "Type_Length" "$group" "$brms_folder_type_length" "$R_type_length" "type_length"
    run_brms "Sigma_ICC"   "$group" "$brms_folder_type_length" "$R_sigma_icc"   "sigma_icc"
    run_brms "Contrast"    "$group" "$brms_folder_type_length" "$R_contrast"    "contrast"
    run_divergence "$group" "$brms_folder_polarity" "$brms_folder_type_length" "$R_divergence" "divergence"
    echo
done

elapsed=$(( $(date +%s) - start_time ))
echo "[$(tnow)] Finished all BRMS runs in $(( elapsed / 3600 )) hr $(( (elapsed % 3600) / 60 )) min $(( elapsed % 60 )) sec"

if [ ${#failures[@]} -gt 0 ]; then
    echo "${#failures[@]} run(s) failed:"
    printf '  %s\n' "${failures[@]}"
    exit 1
fi
echo "Done!"