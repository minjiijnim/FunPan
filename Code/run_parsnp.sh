#!/bin/bash
set -euo pipefail

# =============================================================================
# Parsnp Core Genome Alignment and Phylogenetic Analysis Script
# =============================================================================
# This script:
#  - Locates funannotate clean genome files (*.clean.fa)
#  - Runs Parsnp for rapid core genome alignment
#  - Generates a core genome phylogenetic tree
# =============================================================================

# Threading: auto-detect or use PIPELINE_THREADS override
TOTAL_CORES=$(nproc 2>/dev/null || echo 1)
if [[ -n "${PIPELINE_THREADS:-}" ]]; then
  THREADS="$PIPELINE_THREADS"
else
  RESERVED=2
  THREADS=$(( TOTAL_CORES - RESERVED ))
  if (( THREADS < 1 )); then
    THREADS=1
  fi
fi

echo "Detected $TOTAL_CORES cores, using $THREADS threads for this pipeline"
echo

# Prompt for species name
echo "> Enter species name (genus/species format):"
read -r species

# ---- Paths ------------------------------------------------------------------
funannotate_dir="../Species/$species/funannotate_output"
parsnp_dir="../Species/$species/parsnp_output"
parsnp_input_dir="../Species/$species/parsnp_input"
logfile="../Species/$species/parsnp_$(date +%Y%m%d_%H%M%S).log"

if [[ ! -d "$funannotate_dir" ]]; then
  echo "ERROR: Funannotate output directory not found at $funannotate_dir" >&2
  exit 1
fi

mkdir -p "$parsnp_dir"
mkdir -p "$parsnp_input_dir"

# ---- Gather funannotate clean genome files ----------------------------------
echo "> Scanning for funannotate clean files in: $funannotate_dir"

genome_files=()
while IFS= read -r -d '' f; do
  genome_files+=("$f")
done < <(find -L "$funannotate_dir" -type f -name "*.clean.fa" -print0 2>/dev/null)

if [[ ${#genome_files[@]} -eq 0 ]]; then
  echo "ERROR: No funannotate clean files (*.clean.fa) found in $funannotate_dir" >&2
  exit 1
fi

echo "> Found ${#genome_files[@]} clean genome file(s)"

if [[ ${#genome_files[@]} -lt 3 ]]; then
  echo "ERROR: Parsnp requires at least 3 genomes for analysis" >&2
  exit 1
fi

# ---- Create symlinks for Parsnp input ---------------------------------------
# Parsnp requires all genomes in a single directory
echo "> Creating symlinks in parsnp_input directory..."
rm -f "$parsnp_input_dir"/*.fa 2>/dev/null || true

for f in "${genome_files[@]}"; do
  # Extract strain name from path (directory name)
  strain_dir=$(basename "$(dirname "$f")")
  ln -sf "$(realpath "$f")" "$parsnp_input_dir/${strain_dir}.fa"
done

echo "> Linked ${#genome_files[@]} genomes to $parsnp_input_dir"

# ---- Select best reference genome from metadata -----------------------------
# Priority: 1) Best GCF (RefSeq) genome by quality
#           2) If no GCF, best GCA genome by quality
echo ""
echo "> Selecting best quality reference genome..."

metadata_file="../Species/$species/metadata.csv"
best_reference=""
best_score=-999999
ref_type=""

# Helper function to calculate quality score from metadata
calc_quality_score() {
  local acc="$1"
  local meta="$2"
  awk -F',' -v acc="$acc" '
    NR==1 {
      for (i=1; i<=NF; i++) {
        h = tolower($i)
        gsub(/[ \t\r]+/, "", h)
        if (h == "assemblyaccession") acc_idx=i
        if (h == "buscocomplete") busco_idx=i
        if (h == "buscototal") busco_total_idx=i
        if (h == "assemblystatsnumberofcontigs") contigs_idx=i
        if (h == "assemblystatscontign50" || h == "contign50" || h == "n50") n50_idx=i
      }
      next
    }
    {
      curr_acc = $acc_idx
      gsub(/^[ \t\r]+|[ \t\r]+$/, "", curr_acc)
      if (curr_acc != acc) next

      busco = $busco_idx
      busco_total = $busco_total_idx
      contigs = $contigs_idx
      n50 = $n50_idx

      if (busco_total > 0) {
        busco_pct = (busco / busco_total) * 100
      } else if (busco <= 1) {
        busco_pct = busco * 100
      } else {
        busco_pct = busco
      }

      quality_score = busco_pct - (contigs * 0.01) + (n50 * 0.000001)
      print quality_score
      exit
    }
  ' "$meta"
}

if [[ -f "$metadata_file" ]]; then
  # Step 1: Look for best GCF (RefSeq) genome among clean files
  echo "   Searching for GCF (RefSeq) genomes..."
  best_gcf=""
  best_gcf_score=-999999

  for f in "${genome_files[@]}"; do
    strain_dir=$(basename "$(dirname "$f")")
    if [[ "$strain_dir" == GCF_* ]]; then
      acc=$(echo "$strain_dir" | sed -E 's/^(GCF_[0-9]+\.[0-9]+).*/\1/')
      quality_score=$(calc_quality_score "$acc" "$metadata_file")

      if [[ -n "$quality_score" ]] && (( $(echo "$quality_score > $best_gcf_score" | bc -l) )); then
        best_gcf_score="$quality_score"
        best_gcf="$f"
      fi
    fi
  done

  if [[ -n "$best_gcf" ]]; then
    best_reference="$best_gcf"
    best_score="$best_gcf_score"
    ref_type="GCF (RefSeq)"
    echo "   Found GCF genome: $(basename "$(dirname "$best_gcf")")"
  else
    # Step 2: No GCF found, use best quality GCA genome
    echo "   No GCF genomes found, selecting best quality GCA genome..."

    for f in "${genome_files[@]}"; do
      strain_dir=$(basename "$(dirname "$f")")
      acc=$(echo "$strain_dir" | sed -E 's/^(GC[AF]_[0-9]+\.[0-9]+).*/\1/')
      quality_score=$(calc_quality_score "$acc" "$metadata_file")

      if [[ -n "$quality_score" ]] && (( $(echo "$quality_score > $best_score" | bc -l) )); then
        best_score="$quality_score"
        best_reference="$f"
      fi
    done
    ref_type="GCA (best quality)"
  fi
else
  echo "   WARNING: metadata.csv not found at $metadata_file"
  echo "   Using first genome as reference"
  best_reference="${genome_files[0]}"
  ref_type="first available (no metadata)"
fi

# Build reference option for parsnp
reference_opt=""
if [[ -n "$best_reference" ]]; then
  strain_dir=$(basename "$(dirname "$best_reference")")
  ref_symlink="$parsnp_input_dir/${strain_dir}.fa"
  reference_opt="-r $ref_symlink"
  echo "> Selected reference: $strain_dir"
  echo "   Type: $ref_type"
  echo "   Quality score: $best_score"
else
  echo "> No reference selected, Parsnp will auto-select"
fi

# ---- Check Parsnp availability ----------------------------------------------
if ! command -v parsnp >/dev/null 2>&1; then
  echo "> 'parsnp' not found in PATH, installing via conda..."
  if ! command -v conda >/dev/null 2>&1; then
    echo "ERROR: conda not found. Please install conda first." >&2
    exit 1
  fi
  conda install -y -c bioconda parsnp
  if ! command -v parsnp >/dev/null 2>&1; then
    echo "ERROR: Failed to install parsnp" >&2
    exit 1
  fi
fi
echo "> Using 'parsnp' command"

# ---- Check for existing output ----------------------------------------------
force_overwrite="--force"


#force_overwrite=""
#if [[ -d "$parsnp_dir" && -n "$(ls -A "$parsnp_dir" 2>/dev/null)" ]]; then
#  echo ""
#  echo "> Existing Parsnp output found in $parsnp_dir"
#  echo "> Overwrite? (y/n):"
#  read -r overwrite
#  if [[ "$overwrite" != "y" && "$overwrite" != "Y" ]]; then
#    echo "> Skipping Parsnp run"
#    echo "> Results directory: $parsnp_dir"
#    exit 0
#  fi
#  force_overwrite="--force"
#fi

# ---- Run Parsnp -------------------------------------------------------------
echo ""
echo "> Running Parsnp on ${species} genomes"
echo "> This may take a while depending on genome sizes..."

# Run parsnp with all genomes in the directory
# -c: force inclusion of all input sequences
# -p: number of threads
# -d: input directory (symlinked clean files)
# -o: output directory
# --force-overwrite: overwrite existing output if confirmed

parsnp \
  -d "$parsnp_input_dir" \
  -o "$parsnp_dir" \
  -p "$THREADS" \
  -c \
  $reference_opt \
  $force_overwrite \
  2>&1 | tee "$logfile"

# ---- Verify output -----------------------------------------------------------
if [[ ! -f "$parsnp_dir/parsnp.tree" ]]; then
  echo "ERROR: Parsnp did not produce a tree file" >&2
  exit 1
fi

# If tree and core alignment exist and are non-empty, remove the parsnp input folder
if [[ -s "$parsnp_dir/parsnp.tree" && -s "$parsnp_dir/parsnp.xmfa" ]]; then
  echo "> Parsnp tree and core alignment present and non-empty."
  echo "> Removing parsnp input directory: $parsnp_input_dir"
  rm -rf "$parsnp_input_dir"
else
  echo "> Parsnp output incomplete or empty; keeping parsnp input directory: $parsnp_input_dir"
fi

# ---- Summary ----------------------------------------------------------------
echo ""
echo "=== Parsnp Summary for $species ==="
echo "Input source:  $funannotate_dir (*.clean.fa)"
echo "Input genomes: ${#genome_files[@]}"
if [[ -n "$best_reference" ]]; then
  echo "Reference:     $strain_dir"
  echo "Ref type:      $ref_type"
  echo "Ref score:     $best_score"
fi
echo "Output tree:   $parsnp_dir/parsnp.tree"
echo "Core genome:   $parsnp_dir/parsnp.xmfa"
echo "GGR file:      $parsnp_dir/parsnp.ggr"

echo ""
echo "> Parsnp analysis complete for $species"
echo "> Results directory: $parsnp_dir"
