#!/bin/bash
set -euo pipefail

# =============================================================================
# IQ-TREE2 Phylogenetic Tree from Gubbins Output
# =============================================================================
# This script:
#  - Takes the recombination-filtered alignment from Gubbins as input
#  - Runs IQ-TREE2 with ModelFinder to select optimal substitution model
#  - Produces a phylogenetic tree with ultrafast bootstrap support
# =============================================================================

# Threading: auto-detect or use PIPELINE_THREADS override
TOTAL_CORES=$(nproc 2>/dev/null || echo 1)
#if [[ -n "${PIPELINE_THREADS:-}" ]]; then
#  THREADS="$PIPELINE_THREADS"
#else
#  RESERVED=4
#  THREADS=$(( TOTAL_CORES - RESERVED ))
#  if (( THREADS < 1 )); then
#    THREADS=1
#  fi
#fi

THREADS=1
# Because IQ-TREE keeps crashing at multi threads

echo "Detected $TOTAL_CORES cores, using $THREADS threads for this pipeline"
echo

# Prompt for species name
echo "> Enter species name (genus/species format):"
read -r species

# ---- Paths ------------------------------------------------------------------
gubbins_dir="../Species/$species/gubbins_output"
iqtree_dir="../Species/$species/iqtree_output_snv"

# ---- Locate Gubbins filtered alignment --------------------------------------
input_alignment="$gubbins_dir/gubbins.filtered_polymorphic_sites.fasta"

if [[ ! -f "$input_alignment" ]]; then
  echo "ERROR: Gubbins filtered alignment not found at $input_alignment" >&2
  echo "Run run_gubbins.sh first." >&2
  exit 1
fi

echo "> Found Gubbins filtered alignment: $input_alignment"

mkdir -p "$iqtree_dir"

# ---- Validate alignment -----------------------------------------------------
echo "> Validating alignment..."

seq_count=$(grep -c "^>" "$input_alignment" 2>/dev/null || echo 0)
echo "> Found $seq_count sequences in alignment"

if [[ $seq_count -lt 3 ]]; then
  echo "ERROR: IQ-TREE requires at least 3 sequences" >&2
  exit 1
fi

# ---- Check IQ-TREE 2 availability -------------------------------------------
# Prefer iqtree2 command (IQ-TREE 2.x) over generic iqtree (might be v3)
IQTREE_CMD=""
if command -v iqtree2 >/dev/null 2>&1; then
  IQTREE_CMD="iqtree2"
elif command -v iqtree >/dev/null 2>&1; then
  # Check if it's version 2.x
  iqtree_version=$(iqtree --version 2>&1 | head -1 || echo "")
  if [[ "$iqtree_version" == *"2."* ]]; then
    IQTREE_CMD="iqtree"
  else
    echo "WARNING: Found iqtree but it appears to be version 3.x" >&2
    echo "This script requires IQ-TREE 2.x" >&2
  fi
fi

if [[ -z "$IQTREE_CMD" ]]; then
  echo "ERROR: IQ-TREE 2 not found in PATH" >&2
  echo "Install with: conda install -c bioconda -c conda-forge 'iqtree>=2,<3'" >&2
  exit 1
fi
echo "> Using '$IQTREE_CMD' command"

# ---- Check for existing output ----------------------------------------------
output_prefix="$iqtree_dir/gubbins_tree"

if [[ -f "${output_prefix}.treefile" ]]; then
  echo ""
  echo "> Existing IQ-TREE output found. Overwrite? (y/n):"
  read -r overwrite
  if [[ "$overwrite" != "y" && "$overwrite" != "Y" ]]; then
    echo "> Skipping IQ-TREE run"
    echo "> Results directory: $iqtree_dir"
    exit 0
  fi
  rm -f "$iqtree_dir"/gubbins_tree.*
fi

# ---- Run IQ-TREE2 -----------------------------------------------------------
echo ""
echo "> Running IQ-TREE2 with ModelFinder on Gubbins filtered alignment..."
echo "> This may take a while depending on alignment size..."

# Run IQ-TREE 2 with ModelFinder
# -m MFP: ModelFinder to select best substitution model
# -bb 1000: 1000 ultrafast bootstrap replicates
# --bnni: reduce risk of overestimating branch supports with UFBoot
$IQTREE_CMD \
  -s "$input_alignment" \
  --prefix "$output_prefix" \
  -m MFP \
  -bb 1000 \
  --bnni \
  -T "$THREADS"

# ---- Verify output -----------------------------------------------------------
if [[ ! -f "${output_prefix}.treefile" ]]; then
  echo "ERROR: IQ-TREE2 did not produce a tree file" >&2
  exit 1
fi

# ---- Extract model info ------------------------------------------------------
best_model="unknown"
if [[ -f "${output_prefix}.iqtree" ]]; then
  best_model=$(grep "Best-fit model" "${output_prefix}.iqtree" | head -1 | sed 's/.*: //' || echo "unknown")
fi

# ---- Summary ----------------------------------------------------------------
echo ""
echo "=== IQ-TREE2 Summary for $species (Gubbins input) ==="
echo "Input alignment:      $input_alignment"
echo "Sequences:            $seq_count"
echo "Best-fit model:       $best_model"
echo "Bootstrap replicates: 1000 (ultrafast)"
echo ""
echo "Output files:"
echo "  Tree file:          ${output_prefix}.treefile"
echo "  Log file:           ${output_prefix}.log"
echo "  IQ-TREE report:     ${output_prefix}.iqtree"

echo ""
echo "> IQ-TREE2 analysis complete for $species"
echo "> Results directory: $iqtree_dir"
