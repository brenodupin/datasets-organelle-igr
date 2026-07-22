#!/usr/bin/env bash
set -eu
project_dir="$(dirname "$(dirname "$(realpath "${BASH_SOURCE[0]}")")")"

###############################################################
# IF YOU DATA DIRECTORY MOVES, EDIT THIS LINE.
# Every group path is built from it.
base=${project_dir}
# Default points to the dataset-organelle-igr directory.
###############################################################

brms_folder_polarity="brms_polarity"
brms_folder_type_length="brms_type_length"

output_path="${base}/all_groups_brms_results.txt"

tnow() {
    date +"%Y-%m-%d %H:%M:%S"
}

print_section() {
    local label="$1"
    local file="$2"
    local indent="  "
 
    echo "${indent}>>>>>>> ${label} <<<<<<<<"
    if [ -f "$file" ]; then
        sed "s/^/${indent}${indent}/" "$file"
    else
        echo "${indent}${indent}[missing: $file]"
    fi
    echo
}

extract_text() {
    local GROUP="$1"
    local BRMS_POLARITY="$2"
    local BRMS_TYPE_LENGTH="$3"

    local pol_txt="${base}/${GROUP}/${BRMS_POLARITY}/brms_result.txt"
    local type_length_txt="${base}/${GROUP}/${BRMS_TYPE_LENGTH}/brms_result.txt"
    local sigma_icc_txt="${base}/${GROUP}/${BRMS_TYPE_LENGTH}/brms_sigma_icc_result.txt"
    local contrast_txt="${base}/${GROUP}/${BRMS_TYPE_LENGTH}/brms_contrast_result.txt"

    echo "####### ${GROUP} #######"
    echo

    print_section "${BRMS_POLARITY}/brms_result.txt" "$pol_txt"
    print_section "${BRMS_TYPE_LENGTH}/brms_result.txt" "$type_length_txt"
    print_section "${BRMS_TYPE_LENGTH}/brms_sigma_icc_result.txt" "$sigma_icc_txt"
    print_section "${BRMS_TYPE_LENGTH}/brms_contrast_result.txt" "$contrast_txt"
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

# Build the combined file, one section per group, separated by a blank line.
: > "$output_path"

first=true
for GROUP in "${groups[@]}"; do
    if [ "$first" = true ]; then
        first=false
    else
        printf "\n\n" >> "$output_path"
    fi
    extract_text "$GROUP" "$brms_folder_polarity" "$brms_folder_type_length" >> "$output_path"
done

echo "[$(tnow)] Combined results written to: $output_path"
echo "Done!"