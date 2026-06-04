#!/bin/bash
set -uo pipefail  # Do NOT use -e so we can handle failures manually

# =============================================================================
# Gubbins Recombination Detection and Phylogenetic Analysis Script (with side env)
# =============================================================================
# This script:
#  - Creates/uses a small, isolated conda env at ./.envs/gubbins
#  - Runs Gubbins from that env so your home env stays unchanged
#  - Takes a whole-genome alignment as input
#  - Detects recombination regions and generates a corrected phylogenetic tree
# =============================================================================

# =============================================================================
# Config
# =============================================================================

# --- config for the side env ---
GUBBINS_ENV_PATH="${GUBBINS_ENV_PATH:-./.envs/gubbins}"
# Convert to absolute path so it works after cd
GUBBINS_ENV_PATH="$(cd "$(dirname "$GUBBINS_ENV_PATH")" 2>/dev/null && pwd)/$(basename "$GUBBINS_ENV_PATH")"
GUBBINS_PY="${GUBBINS_PY:-3.10}"
CONDA_CHANNELS="-c conda-forge -c bioconda"

# Tree builder settings (will try iqtree first, fall back to raxml-ng if it fails)
PRIMARY_BUILDER="iqtree"
FALLBACK_BUILDER="raxmlng"

need_cmd() {
  command -v "$1" >/dev/null 2>&1 || { echo "ERROR: '$1' not found in PATH"; exit 1; }
}

ensure_gubbins_env() {
  need_cmd conda

  if conda list -p "$GUBBINS_ENV_PATH" gubbins >/dev/null 2>&1; then
    echo "> Found existing Gubbins env at $GUBBINS_ENV_PATH"
  else
    echo "> Creating Gubbins side env at $GUBBINS_ENV_PATH ..."
    if command -v mamba >/dev/null 2>&1; then
      mamba create -y -p "$GUBBINS_ENV_PATH" $CONDA_CHANNELS \
        "python=${GUBBINS_PY}" gubbins raxmlng iqtree fasttree harvesttools
    else
      conda create -y -p "$GUBBINS_ENV_PATH" $CONDA_CHANNELS \
        "python=${GUBBINS_PY}" gubbins raxmlng iqtree fasttree harvesttools
    fi
  fi

  echo "> Gubbins version check (side env):"
  conda run -p "$GUBBINS_ENV_PATH" run_gubbins.py --version || {
    echo "ERROR: Gubbins not runnable in side env"; exit 1;
  }
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

# Prepare side env
ensure_gubbins_env
echo

# Prompt for species name
echo "> Enter species name (genus/species format):"
read -r species

# ---- Paths (use absolute paths so they work after cd) ----------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
species_dir="$SCRIPT_DIR/../Species/$species"
gubbins_dir="$SCRIPT_DIR/../Species/$species/gubbins_output"
parsnp_dir="$SCRIPT_DIR/../Species/$species/parsnp_output"

mkdir -p "$gubbins_dir"

# ---- Locate input alignment -------------------------------------------------
echo "> Searching for input alignment..."

input_alignment=""

# Check for Parsnp output first (XMFA format needs conversion)
if [[ -f "$parsnp_dir/parsnp.xmfa" ]]; then
  echo "> Found Parsnp XMFA alignment"
  echo "> Note: XMFA needs conversion to FASTA for Gubbins"

  echo "> Converting XMFA to FASTA using harvesttools (from side env)..."
  converted_fasta="$gubbins_dir/core_genome_alignment.fasta"
  conda run -p "$GUBBINS_ENV_PATH" harvesttools -x "$parsnp_dir/parsnp.xmfa" -M "$converted_fasta" || {
    echo "ERROR: harvesttools conversion failed" >&2
    exit 1
  }
  input_alignment="$converted_fasta"
fi

# Require Parsnp output
if [[ ! -f "$parsnp_dir/parsnp.xmfa" ]]; then
  echo "ERROR: No parsnp.xmfa found. Run run_parsnp.sh first." >&2
  exit 1
fi

echo "> Using alignment: $input_alignment"

# ---- Validate alignment -----------------------------------------------------
echo "> Validating alignment..."

seq_count=$(grep -c "^>" "$input_alignment" 2>/dev/null || echo 0)
echo "> Found $seq_count sequences in alignment"

if [[ $seq_count -lt 3 ]]; then
  echo "ERROR: Gubbins requires at least 3 sequences" >&2
  exit 1
fi

echo "> Using Gubbins from side env: $GUBBINS_ENV_PATH"

# ---- Check for existing output ----------------------------------------------
output_prefix="$gubbins_dir/gubbins"
tree_builder="$PRIMARY_BUILDER"

if [[ -f "${output_prefix}.final_tree.tre" ]]; then
  echo ""
  echo "> Existing Gubbins output found. Overwrite? (y/n):"
  read -r overwrite
  if [[ "$overwrite" != "y" && "$overwrite" != "Y" ]]; then
    echo "> Skipping Gubbins run"
    echo "> Results directory: $gubbins_dir"
    exit 0
  fi
  rm -f "$gubbins_dir"/gubbins.*
fi

# ---- Run Gubbins (with fallback) --------------------------------------------
run_gubbins_with_builder() {
  local builder="$1"
  echo ""
  echo "> Running Gubbins on ${species} alignment with $builder..."
  echo "> This may take a while depending on alignment size..."

  # Change to output directory (Gubbins writes to current directory)
  cd "$gubbins_dir"

  conda run -p "$GUBBINS_ENV_PATH" run_gubbins.py \
    --threads "$THREADS" \
    --tree-builder "$builder" \
    --prefix gubbins \
    --verbose \
    "$input_alignment" \
    2>&1 | tee gubbins.log

  # Return success if final tree was produced
  [[ -f "${output_prefix}.final_tree.tre" ]]
}

# Try primary builder (iqtree), fall back to raxml-ng if it fails
if run_gubbins_with_builder "$PRIMARY_BUILDER"; then
  tree_builder="$PRIMARY_BUILDER"
  echo "> Gubbins succeeded with $PRIMARY_BUILDER"
else
  echo ""
  echo "> WARNING: Gubbins failed with $PRIMARY_BUILDER, trying $FALLBACK_BUILDER..."
  # Clean up failed attempt
  rm -f "$gubbins_dir"/gubbins.*

  if run_gubbins_with_builder "$FALLBACK_BUILDER"; then
    tree_builder="$FALLBACK_BUILDER"
    echo "> Gubbins succeeded with $FALLBACK_BUILDER"
  else
    echo "ERROR: Gubbins failed with both $PRIMARY_BUILDER and $FALLBACK_BUILDER" >&2
    exit 1
  fi
fi

# ---- Summary ----------------------------------------------------------------
recomb_count=0
if [[ -f "${output_prefix}.recombination_predictions.gff" ]]; then
  recomb_count=$(grep -c "^[^#]" "${output_prefix}.recombination_predictions.gff" 2>/dev/null || echo 0)
fi

echo ""
echo "=== Gubbins Summary for $species ==="
echo "Input sequences:      $seq_count"
echo "Tree builder:         $tree_builder"
echo "Recombination events: $recomb_count"
echo ""
echo "Output files:"
echo "  Final tree:         ${output_prefix}.final_tree.tre"
echo "  Recombination GFF:  ${output_prefix}.recombination_predictions.gff"
echo "  Filtered alignment: ${output_prefix}.filtered_polymorphic_sites.fasta"
echo "  Node-labeled tree:  ${output_prefix}.node_labelled.final_tree.tre"

echo ""
echo "> Gubbins analysis complete for $species"
echo "> Results directory: $gubbins_dir"
