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

sections=(
    "${brms_folder_polarity}:polarity"
    "${brms_folder_type_length}:type_length"
    "${brms_folder_type_length}:sigma_icc"
    "${brms_folder_type_length}:contrast"
)

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

    echo "####### ${GROUP} #######"
    echo

    local section folder append rel_path
    for section in "${sections[@]}"; do
        folder="${section%%:*}"
        append="${section##*:}"
        rel_path="${folder}/brms_${append}.txt"

        print_section "$rel_path" "${base}/${GROUP}/${rel_path}"
    done
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
    extract_text "$GROUP" >> "$output_path"
done

echo "[$(tnow)] Combined results written to: $output_path"
echo "Done!"