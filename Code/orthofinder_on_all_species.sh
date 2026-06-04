#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# Build combined OrthoFinder input from:
#   funannotate_output/predict_results/*.proteins.fa
# (Does NOT modify any per-species tmp_proteins.)
# Run OrthoFinder into a NEW, non-existent run dir (required by OrthoFinder).
#
# Usage:
#   bash orthofinder_from_funannotate.sh
#   # or override ROOT:
#   ROOT="/datadrive/Species" bash orthofinder_from_funannotate.sh
# ============================================================

ROOT="${ROOT:-../Species}"

echo "> CWD : $(pwd)"
echo "> ROOT: ${ROOT}"

# Threading: auto-detect or use PIPELINE_THREADS override
TOTAL_CORES=$(nproc 2>/dev/null || echo 1)
if [[ -n "${PIPELINE_THREADS:-}" ]]; then
  THREADS="$PIPELINE_THREADS"
else
  RESERVED=4
  THREADS=$(( TOTAL_CORES - RESERVED ))
  if (( THREADS < 1 )); then
    THREADS=1
  fi
fi

echo "Detected $TOTAL_CORES cores, using $THREADS threads for this pipeline"
echo

# --- Prompt for genus and species (space/comma separated) ---
echo "> Enter genus (e.g., Aspergillus):"
read -r genus
[[ -z "${genus}" ]] && { echo "ERROR: no genus provided"; exit 1; }

echo "> Enter species (space/comma separated, e.g., flavus fumigatus niger oryzae):"
read -r species_input
species_input="$(echo "$species_input" | tr '[:upper:],' '[:lower:] ' | xargs)"
[[ -z "${species_input}" ]] && { echo "ERROR: no species provided"; exit 1; }

echo "> Genus  : ${genus}"
echo "> Species: ${species_input}"

# --- abs path helper ---
abspath() {
  if command -v realpath >/dev/null 2>&1; then
    realpath "$1"
  elif command -v readlink >/dev/null 2>&1; then
    readlink -f "$1"
  else
    echo "$1"
  fi
}

# --- Prepare combined input dir (safe to refresh) ---
combo_root="${ROOT}/${genus}/all_combined/orthofinder_output"
combo_input="${combo_root}/input_proteins"
mkdir -p "${combo_input}"
rm -f "${combo_input}/"* 2>/dev/null || true

# --- Collect & link proteins from funannotate predict_results ---
total_linked=0
for species in $species_input; do
  echo
  echo "=== ${genus}/${species} ==="
  fun_dir="${ROOT}/${genus}/${species}/funannotate_output"

  if [[ ! -d "${fun_dir}" ]]; then
    echo "WARN: ${fun_dir} not found; skipping ${species}"
    continue
  fi

  echo "> Searching: ${fun_dir}/predict_results/*.proteins.fa"
  mapfile -t prot_files < <(find -L "${fun_dir}" -type f -path "*/predict_results/*.proteins.fa" | sort)

  if [[ ${#prot_files[@]} -eq 0 ]]; then
    echo "WARN: no *.proteins.fa found for ${species}"
    continue
  fi

  linked=0
  for pf in "${prot_files[@]}"; do
    src="$(abspath "$pf")"
    base="$(basename "$pf")"
    ln -s -f "${src}" "${combo_input}/${species}__${base}"
    ((linked++)) || true
    ((total_linked++)) || true
  done

  echo "OK: linked ${linked} file(s) for ${species}"
done

if [[ ${total_linked} -eq 0 ]]; then
  echo "ERROR: No protein FASTA files were linked into ${combo_input}. Aborting."
  exit 1
fi

echo
echo ">>> Combined ${total_linked} protein files into: ${combo_input}"

## --- Threads: cap at 32 to avoid warnings from DIAMOND/OrthoFinder ---
#NPROC=$(nproc || echo 1)
#CAP=32
#if [[ "$NPROC" -gt "$CAP" ]]; then
#  THREADS="${THREADS:-$CAP}"
#else
#  THREADS="${THREADS:-$NPROC}"
#fi
#[[ "$THREADS" -lt 1 ]] && THREADS=1

# --- Choose a NEW, non-existent run directory for -o (critical) ---
# Do NOT pre-create this path.
run_parent="${combo_root}"
ts="$(date +%Y%m%d_%H%M%S)"
run_dir="${run_parent}/run_${ts}"

echo
echo "> Running OrthoFinder:"
echo "  - input : ${combo_input}"
echo "  - out   : ${run_dir}   (must NOT exist; OrthoFinder will create it + inner Results_*)"
echo "  - cores : ${THREADS}"

# Ensure parent exists, but NOT run_dir itself
mkdir -p "${run_parent}"
if [[ -e "${run_dir}" ]]; then
  echo "ERROR: run_dir already exists (${run_dir}). Please delete or try again."; exit 1
fi

# Launch OrthoFinder (it will create run_dir and then run_dir/Results_*)
orthofinder -f "${combo_input}" -t "${THREADS}" -a "${THREADS}" -S diamond -M msa -A mafft -o "${run_dir}"

# --- Find the actual Results_* directory that OrthoFinder created ---
results_dir="$(ls -dt "${run_dir}"/Results_* 2>/dev/null | head -n1 || true)"
if [[ -z "${results_dir}" ]]; then
  echo "ERROR: Could not find an OrthoFinder Results_* directory under ${run_dir}"
  exit 1
fi

# --- Validate presence of Orthogroups.tsv ---
orthogroups_file="$(find "${results_dir}" -type f -path "*/Orthogroups/Orthogroups.tsv" | head -n 1 || true)"
if [[ -z "${orthogroups_file}" ]]; then
  echo "ERROR: Orthogroups.tsv missing in ${results_dir}"
  exit 1
fi

echo
echo "> OrthoFinder analysis complete"
echo "> Orthogroups.tsv: ${orthogroups_file}"
echo "> Full output    : ${results_dir}"
