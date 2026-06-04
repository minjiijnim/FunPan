#!/bin/bash
set -euo pipefail

# =============================================================================
# EggNOG Functional Annotation Script
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
read -r species

# ---- Paths ------------------------------------------------------------------
fun_dir="../Species/$species/funannotate_output"
protein_dir="${fun_dir}/tmp_proteins"
eggnog_dir="../Species/$species/eggnog_output"
data_dir="../Data/EggNOG"
mkdir -p "$protein_dir" "$eggnog_dir" "$data_dir"

# ---- Gather protein files from funannotate predict_results -------------------
echo "> Searching for *.proteins.fa under: $fun_dir"
# Refresh the tmp folder (symlinks only)
find "$protein_dir" -type l -delete || true

mapfile -t prot_sources < <(find "$fun_dir" -type f -path "*/predict_results/*.proteins.fa" | sort)
if [[ ${#prot_sources[@]} -eq 0 ]]; then
  echo "ERROR: No *.proteins.fa found under $fun_dir" >&2
  exit 1
fi

for p in "${prot_sources[@]}"; do
  abs=$(realpath "$p")
  ln -s "$abs" "$protein_dir/$(basename "$abs")"
done

mapfile -t protein_files < <(find "$protein_dir" -type l -name "*.proteins.fa" | sort)
echo "> Linked ${#protein_files[@]} protein file(s) into $protein_dir"

# ---- Ensure eggNOG DB is available ------------------------------------------
have_db=false
# Basic sanity: sqlite DB + at least one DIAMOND file exists
if [[ -f "$data_dir/eggnog.db" ]] && compgen -G "$data_dir/data/"'*_dmnd.dmnd' >/dev/null; then
  have_db=true
fi

if [[ "$have_db" == false ]]; then
  echo "> No eggNOG DB detected in $data_dir"

  # Try to find the official downloader
  DL=""
  if command -v download_eggnog_data.py >/dev/null 2>&1; then
    DL="$(command -v download_eggnog_data.py)"
  else
    # Sometimes it’s installed next to emapper.py
    if command -v emapper.py >/dev/null 2>&1; then
      em_dir="$(dirname "$(command -v emapper.py)")"
      [[ -x "$em_dir/download_eggnog_data.py" ]] && DL="$em_dir/download_eggnog_data.py" || true
    fi
  fi
  
  if [[ -n "$DL" ]]; then
    echo "> Using downloader: $DL (db=euk)"
    # Try non-interactive flag first (some versions support -y)
    if python3 "$DL" --help 2>&1 | grep -qE '\s-y(,| )|--yes|--all'; then
      python3 "$DL" --data_dir "$data_dir" -d euk -y
    else
      # Fallback: auto-answer all prompts with "y"
      yes | python3 "$DL" --data_dir "$data_dir" -d euk
    fi
  else
    echo "ERROR: Could not locate 'download_eggnog_data.py' on PATH or next to emapper.py." >&2
    echo "Please install the DB with one of these commands, then re-run this script:" >&2
    echo "  mamba install -n pipeline_test -c bioconda eggnog-data" >&2
    echo "    # or:" >&2
    echo "  python -m pip install --upgrade eggnog-mapper" >&2
    exit 1
  fi


  # Re-check
  if [[ ! -f "$data_dir/eggnog.db" ]] || ! compgen -G "$data_dir/"'*.dmnd' >/dev/null; then
    echo "ERROR: eggNOG database not found/complete after download into $data_dir." >&2
    exit 1
  fi
fi

# ---- Helper function to check if annotation file is empty -------------------
is_empty_annotation() {
  local annot_file="$1"
  # Check if file contains "0 queries scanned" or lacks a proper header
  if grep -q "^## 0 queries scanned" "$annot_file" 2>/dev/null; then
    return 0  # true - file is empty
  elif ! grep -q "^#query\|^query" "$annot_file" 2>/dev/null; then
    return 0  # true - no header found
  fi
  return 1  # false - file appears valid
}

# ---- Run EggNOG-mapper -------------------------------------------------------
echo "> Running EggNOG-mapper on ${species} protein files"
for protein_file in "${protein_files[@]}"; do
  base_name="$(basename "$protein_file" .proteins.fa)"
  annotation_file="$eggnog_dir/${base_name}.emapper.annotations"
  hits_file="$eggnog_dir/${base_name}.emapper.hits"

  # Check if annotation file exists and is valid
  if [[ -f "$annotation_file" ]]; then
    if is_empty_annotation "$annotation_file"; then
      echo "> Found empty/invalid annotation file for $base_name - will rerun"
      # Remove ALL emapper output files to force fresh run
      rm -f "$eggnog_dir/${base_name}".emapper.*
      echo "> Processing $base_name"
      emapper.py \
        -i "$protein_file" -o "$base_name" \
        --output_dir "$eggnog_dir" --itype proteins \
        --cpu "$THREADS" --data_dir "$data_dir"
    else
      echo "> Skipping $base_name (already annotated)"
      continue
    fi
  elif [[ -f "$hits_file" ]]; then
    echo "> Processing $base_name"
    echo "> Resuming previous run for $base_name"
    emapper.py \
      -i "$protein_file" -o "$base_name" \
      --output_dir "$eggnog_dir" --itype proteins \
      --cpu "$THREADS" --data_dir "$data_dir" --resume
  else
    echo "> Processing $base_name"
    echo "> Running fresh annotation for $base_name"
    emapper.py \
      -i "$protein_file" -o "$base_name" \
      --output_dir "$eggnog_dir" --itype proteins \
      --cpu "$THREADS" --data_dir "$data_dir"
  fi
done

# ---- Verify output -----------------------------------------------------------
annotation_count=$(find "$eggnog_dir" -name "*.emapper.annotations" | wc -l | awk '{print $1}')
if [[ "$annotation_count" -eq 0 ]]; then
  echo "ERROR: No EggNOG annotation files were created in $eggnog_dir" >&2
  exit 1
fi

# Check for empty annotation files
empty_count=0
echo "> Checking for empty annotation files..."
while IFS= read -r annot_file; do
  if is_empty_annotation "$annot_file"; then
    echo "  WARNING: Empty annotation file: $(basename "$annot_file")"
    ((empty_count++))
  fi
done < <(find "$eggnog_dir" -name "*.emapper.annotations")

if [[ $empty_count -gt 0 ]]; then
  echo "WARNING: Found $empty_count empty annotation file(s)." >&2
  echo "These may need manual inspection or rerunning with different parameters." >&2
fi

# ─── Cleanup emapper temp DIAMOND files ─────────────────────────────────────
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
tmp_pattern="${SCRIPT_DIR}/emappertmp_dmdn_*"
if compgen -G "$tmp_pattern" >/dev/null; then
    echo "> Removing emapper temporary files..."
    rm -rf $tmp_pattern
fi

echo "> EggNOG-mapper analysis complete"
echo "> Annotation files: $annotation_count"
echo "> Output directory: $eggnog_dir"


