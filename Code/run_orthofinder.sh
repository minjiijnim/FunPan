#!/bin/bash
set -euo pipefail

# =============================================================================
# OrthoFinder Comparative Genomics Script
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

# Prompt for species name
echo "> Enter species name (genus/species format):"
read species

# Directory setup
ortho_dir="../Species/$species/orthofinder_output"
fun_dir="../Species/$species/funannotate_output"
tmp_protein_dir="${fun_dir}/tmp_proteins"

# Create fresh tmp_proteins directory
mkdir -p "$tmp_protein_dir"
rm -f "$tmp_protein_dir"/*  # Clear old symlinks if any

echo "> Searching for *.proteins.fa under $fun_dir"

# Find and symlink all relevant protein FASTA files
prot_files=()
while IFS= read -r prot_file; do
    abs_path=$(realpath "$prot_file")
    ln -s -f "$abs_path" "$tmp_protein_dir/$(basename "$abs_path")"
    prot_files+=("$abs_path")
done < <(find -L "$fun_dir" -type f -path "*/predict_results/*.proteins.fa")

if [ ${#prot_files[@]} -eq 0 ]; then
    echo "ERROR: No .proteins.fa files found under $fun_dir"
    exit 1
fi

echo "> Linked ${#prot_files[@]} protein files to $tmp_protein_dir"


# Run OrthoFinder
echo "> Running OrthoFinder on ${species} proteins"
if ! orthofinder -f "$tmp_protein_dir" -t $THREADS -a $THREADS -S diamond -M msa -A mafft -o "$ortho_dir"; then
    echo "ERROR: OrthoFinder execution failed"
    exit 1
fi

# Validate expected OrthoFinder output
orthogroups_file=$(find "$ortho_dir" -type f -path "*/Orthogroups/Orthogroups.tsv" | head -n 1)

if [[ -z "$orthogroups_file" ]]; then
    echo "ERROR: OrthoFinder output not found (Orthogroups.tsv missing) in $ortho_dir"
    exit 1
fi

echo "> OrthoFinder analysis complete"
echo "> Results are in ${ortho_dir}"
