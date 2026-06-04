#!/bin/bash
set -uo pipefail  # Do NOT use -e so we can handle failures manually

# =============================================================================
# antiSMASH Secondary Metabolite Detection Script (with side env)
# =============================================================================
# This script:
#  - Creates/uses a small, isolated conda env at ./.envs/antismash8
#  - Runs antiSMASH 8.0.2 from that env so your home env stays unchanged
#  - Downloads/checks antiSMASH databases inside that small env
# =============================================================================

# --- config for the side env ---
ANTISMASH_ENV_PATH="${ANTISMASH_ENV_PATH:-./.envs/antismash8}"
ANTISMASH_PY="${ANTISMASH_PY:-3.11}"
ANTISMASH_VER="${ANTISMASH_VER:-8.0.2}"
BIOPYTHON_VER="${BIOPYTHON_VER:-1.81}"
CONDA_CHANNELS="-c conda-forge -c bioconda"

need_cmd() {
  command -v "$1" >/dev/null 2>&1 || { echo "ERROR: '$1' not found in PATH"; exit 1; }
}

ensure_antismash_env() {
  need_cmd conda

  if conda list -p "$ANTISMASH_ENV_PATH" antismash >/dev/null 2>&1; then
    echo "> Found existing antiSMASH env at $ANTISMASH_ENV_PATH"
  else
    echo "> Creating antiSMASH side env at $ANTISMASH_ENV_PATH ..."
    if command -v mamba >/dev/null 2>&1; then
      mamba create -y -p "$ANTISMASH_ENV_PATH" $CONDA_CHANNELS \
        "python=${ANTISMASH_PY}" "antismash=${ANTISMASH_VER}" "biopython=${BIOPYTHON_VER}"
    else
      conda create -y -p "$ANTISMASH_ENV_PATH" $CONDA_CHANNELS \
        "python=${ANTISMASH_PY}" "antismash=${ANTISMASH_VER}" "biopython=${BIOPYTHON_VER}"
    fi
  fi

  echo "> antiSMASH version check (side env):"
  conda run -p "$ANTISMASH_ENV_PATH" antismash --version || {
    echo "ERROR: antiSMASH not runnable in side env"; exit 1;
  }
}

# Return the antiSMASH package dir inside the side env
get_antismash_pkg_dir() {
  conda run -p "$ANTISMASH_ENV_PATH" python - <<'PY'
import pathlib, antismash
print(pathlib.Path(antismash.__file__).parent)
PY
}

# Ensure databases exist in the side env
ensure_antismash_db() {
  echo "> Checking antiSMASH databases..."
  if conda run -p "$ANTISMASH_ENV_PATH" antismash --check-prereqs; then
    echo "> antiSMASH databases are ready."
  else
    echo "> Downloading missing antiSMASH databases..."
    conda run -p "$ANTISMASH_ENV_PATH" download-antismash-databases || {
      echo "ERROR: database download failed"; exit 1;
    }
  fi
}

# =============================================================================
# Main
# =============================================================================

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

echo "> Enter species name (genus/species format):"
read species

# === SETTINGS ===
fun_dir="../Species/$species/funannotate_output"
tmp_gbk_dir="${fun_dir}/tmp_gbk"
OUTDIR="../Species/$species/antismash_output"

# Prepare side env
ensure_antismash_env
ensure_antismash_db

# Prepare GBK symlink directory
mkdir -p "$tmp_gbk_dir"
rm -f "$tmp_gbk_dir"/*  # Clean up old symlinks

echo "> Searching for *.gbk under $fun_dir..."

# Find and symlink all relevant GBK files
gbk_files=()
while IFS= read -r gbk_file; do
    abs_path=$(realpath "$gbk_file")
    ln -s "$abs_path" "$tmp_gbk_dir/$(basename "$abs_path")"
    gbk_files+=("$abs_path")
done < <(find "$fun_dir" -type f -path "*/predict_results/*.gbk")

if [ ${#gbk_files[@]} -eq 0 ]; then
    echo "ERROR: No .gbk files found under $fun_dir"
    exit 1
fi

echo "> Linked ${#gbk_files[@]} gbk files to $tmp_gbk_dir"
mkdir -p "$OUTDIR"

# Initialize counters
success_count=0
fail_count=0
skipped_count=0

# Run antiSMASH for each GBK using the side env
for gbk in "$tmp_gbk_dir"/*.gbk; do
    base=$(basename "$gbk" .gbk)
    outpath="$OUTDIR/$base"

    # Skip if output already exists and looks complete
    if [[ -f "$outpath/index.html" || -f "$outpath/antismash.html" || -f "$outpath/css/main.css" ]]; then
        echo "Skipping $base (already completed)"
        ((skipped_count++))
        continue
    fi

    # Skip if output directory exists and is non-empty (avoid antiSMASH aborting)
    if [[ -d "$outpath" ]] && [[ -n "$(ls -A "$outpath" 2>/dev/null)" ]]; then
        echo "Skipping $base (output folder already exists and is not empty)"
        ((skipped_count++))
        continue
    fi

    echo "> Running antiSMASH on $base..."

    if conda run -p "$ANTISMASH_ENV_PATH" antismash "$gbk" \
        --taxon fungi \
        --output-dir "$outpath" \
        --genefinding-tool none \
        --asf \
        --cb-knownclusters \
        --cpus "$THREADS" \
        --cb-general \
        --cb-subclusters \
        --pfam2go; then

        if [[ -f "$outpath/index.html" || -f "$outpath/antismash.html" ]]; then
            echo "Finished $base successfully."
            ((success_count++))
        else
            echo "antiSMASH ran but output seems incomplete for $base"
            ((fail_count++))
        fi
    else
        echo "antiSMASH failed on $base"
        ((fail_count++))
    fi
done

# Summary
echo ""
echo "=== antiSMASH Summary for $species ==="
echo "Successful: $success_count"
echo "Skipped:    $skipped_count"
echo "Failed:     $fail_count"

if [[ $success_count -eq 0 && $skipped_count -eq 0 ]]; then
    echo "ERROR: antiSMASH did not produce any usable results."
    exit 1
fi

echo ""
echo "antiSMASH analysis complete for $species"
echo "Results directory: $OUTDIR"
