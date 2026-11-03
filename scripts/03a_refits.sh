#!/usr/bin/env bash
set -eu
project_dir="$(dirname "$(dirname "$(realpath "${BASH_SOURCE[0]}")")")"

# If your data directory moves, edit this line; every group path is built from it.
base=${project_dir}

R_pol=${project_dir}/scripts/brms_polarity.R
R_type_length=${project_dir}/scripts/brms_type_length.R
R_sigma_icc=${project_dir}/scripts/brms_sigma_icc.R
R_contrast=${project_dir}/scripts/brms_contrast.R
R_divergence=${project_dir}/scripts/brms_divergence.R

# Robustness refits built by 02a_refits.ipynb (R2.1a, R2.1b).
# Every model of a run is fitted inside <group>/brms_refits/<run>/.
# Override with e.g.: RUNS="drop_longest" ./03a_refits.sh metazoans_mit
refits_dir="brms_refits"
read -ra runs <<< "${RUNS:-drop_longest drop_top1pct gene_gene_only}"

default_groups=(
    fungi_mit
    green_algae_mit
    green_algae_plt
    protists_mit
    protists_plt
    plants_mit
    plants_plt
    metazoans_mit
)

failures=()

tnow() { date '+%Y-%m-%d %H:%M:%S'; }

# Output folder of one group/run; every step of it appends to its log.
run_dir() { echo "${base}/$1/${refits_dir}/$2"; }

# Timestamped line to stdout and the run's log.
log_line() {
    local group="$1" run="$2"; shift 2
    echo "[$(tnow)] $*" | tee -a "$(run_dir "$group" "$run")/brms_refits.log"
}

run_brms() {
    local name="$1" group="$2" run="$3" r_path="$4" suffix="$5"
    local dir log start rc status elapsed
    dir="$(run_dir "$group" "$run")"
    log="${dir}/brms_refits.log"
    start=$(date +%s)

    log_line "$group" "$run" "Starting $name for $group/$run..."
    Rscript "$r_path" \
        "${dir}/tree.nwk" \
        "${dir}/brms_input.tsv" \
        "${dir}/brms_${suffix}.txt" \
        "${dir}/brms_model.rds" \
        "${dir}/brms_${suffix}_row.tsv" 2>&1 | tee -a "$log"
    rc=${PIPESTATUS[0]}

    status="Finished"
    if [ "$rc" -ne 0 ]; then
        status="FAILED"
        failures+=("$group/$run/$name")
    fi

    elapsed=$(( $(date +%s) - start ))
    log_line "$group" "$run" "$status $name for $group/$run in $(( elapsed / 60 )) min $(( elapsed % 60 )) sec"
    return "$rc"
}

# Sigma_ICC and Contrast read the Type_Length caches, so they are skipped
# when Type_Length fails.
skip_brms() {
    local name="$1" group="$2" run="$3"
    log_line "$group" "$run" "SKIPPED $name for $group/$run (Type_Length failed)"
    failures+=("$group/$run/$name (skipped)")
}

# All caches of a run live in one folder, so the folder and output file are
# passed for both arguments: both passes check the same caches and write the
# same file.
run_divergence() {
    local group="$1" run="$2"
    local dir rc
    dir="$(run_dir "$group" "$run")"

    log_line "$group" "$run" "Starting Divergence for $group/$run..."
    Rscript "$R_divergence" "$dir" "${dir}/brms_divergence_row.tsv" \
                            "$dir" "${dir}/brms_divergence_row.tsv" 2>&1 | tee -a "${dir}/brms_refits.log"
    rc=${PIPESTATUS[0]}

    if [ "$rc" -eq 0 ]; then
        log_line "$group" "$run" "Finished Divergence for $group/$run."
    else
        log_line "$group" "$run" "FAILED Divergence for $group/$run."
        failures+=("$group/$run/Divergence")
    fi
}

if [ $# -eq 0 ]; then
    groups=("${default_groups[@]}")
    echo "[$(tnow)] No groups specified. Processing all groups."
else
    groups=("$@")
    echo "[$(tnow)] Processing specified groups: ${groups[*]}"
fi
echo "[$(tnow)] Refit runs: ${runs[*]}"

start_time=$(date +%s)
echo "[$(tnow)] Starting BRMS refit runs..."
echo

for group in "${groups[@]}"; do
    for run in "${runs[@]}"; do
        echo "## --$group / $run"
        if [ ! -f "$(run_dir "$group" "$run")/brms_input.tsv" ]; then
            echo "[$(tnow)] MISSING $(run_dir "$group" "$run")/brms_input.tsv; run 02a_refits.ipynb first."
            failures+=("$group/$run (missing input)")
            echo
            continue
        fi

        run_brms Polarity "$group" "$run" "$R_pol" polarity || true

        # gene_gene_only has a single flanking_types level, so only the
        # polarity model can be fitted (R2.1b).
        if [ "$run" != "gene_gene_only" ]; then
            if run_brms Type_Length "$group" "$run" "$R_type_length" type_length; then
                run_brms Sigma_ICC "$group" "$run" "$R_sigma_icc" sigma_icc || true
                run_brms Contrast  "$group" "$run" "$R_contrast"  contrast  || true
            else
                skip_brms Sigma_ICC "$group" "$run"
                skip_brms Contrast  "$group" "$run"
            fi
        fi

        run_divergence "$group" "$run"
        echo
    done
done

elapsed=$(( $(date +%s) - start_time ))
echo "[$(tnow)] Finished all BRMS refit runs in $(( elapsed / 3600 )) hr $(( (elapsed % 3600) / 60 )) min $(( elapsed % 60 )) sec"

if [ ${#failures[@]} -gt 0 ]; then
    echo "${#failures[@]} run(s) failed:"
    printf '  %s\n' "${failures[@]}"
    exit 1
fi
echo "Done!"
