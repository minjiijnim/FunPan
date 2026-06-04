#!/bin/bash
set -euo pipefail

# -------------------
# Config / threads
# -------------------

# Detect total cores
TOTAL_CORES=$(nproc 2>/dev/null || echo 1)

# If PIPELINE_THREADS is already set, respect it.
# Otherwise, use TOTAL_CORES-4 (but at least 1).
if [[ -n "${PIPELINE_THREADS:-}" ]]; then
  PIPELINE_THREADS="$PIPELINE_THREADS"
else
  # Reserve 4 cores for system/other jobs
  RESERVED=4
  PIPELINE_THREADS=$(( TOTAL_CORES - RESERVED ))
  if (( PIPELINE_THREADS < 1 )); then
    PIPELINE_THREADS=1
  fi
fi

# BUSCO threads: default to PIPELINE_THREADS unless overridden
BUSCO_THREADS=${BUSCO_THREADS:-$PIPELINE_THREADS}

# BUSCO mode: auto / auto_euk / auto_prok
BUSCO_MODE=${BUSCO_MODE:-auto}

# BUSCO offline mode (1 = offline)
BUSCO_OFFLINE=${BUSCO_OFFLINE:-0}

echo "> Total cores: $TOTAL_CORES; using PIPELINE_THREADS=$PIPELINE_THREADS; BUSCO_THREADS=$BUSCO_THREADS"

# Marker dir / inputs file (honor pipeline's MARKER_DIR if exported)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MARKER_DIR="${MARKER_DIR:-${SCRIPT_DIR}/.pipeline_checkpoints}"
PIPELINE_INPUTS="${MARKER_DIR}/pipeline_inputs.txt"
mkdir -p "$MARKER_DIR"

# Upsert helper: upsert_kv KEY VALUE FILE
upsert_kv() {
  local k="$1"; shift
  local v="$1"; shift
  local f="$1"
  if [ -f "$f" ] && grep -qE "^${k}=" "$f"; then
    sed -i.bak -E "s|^(${k}=).*|\1${v}|" "$f"
  else
    printf "%s=%s\n" "$k" "$v" >> "$f"
  fi
}

# --- Make BUSCO + BBTools runnable without naming a conda env ---
# Use the same bin dir as the busco we'll execute
BUSCO_BIN="$(command -v busco)" || { echo "ERROR: busco not found"; exit 127; }
BUSCO_BINDIR="$(dirname "$BUSCO_BIN")"
# Prepend that bin so stats.sh comes from the same place
export PATH="$BUSCO_BINDIR:$PATH"

# Locate BBTools' stats.sh and make sure it's executable
STATS_SH="$(command -v stats.sh)" || { echo "ERROR: stats.sh not found next to $BUSCO_BIN"; exit 127; }
[ -x "$STATS_SH" ] || chmod +x "$STATS_SH"

# Derive the conda prefix hosting busco, then point Java there
# (…/envs/<name>/bin -> …/envs/<name>)
CONDA_PREFIX_FROM_BUSCO="$(cd "$BUSCO_BINDIR"/.. && pwd -P)"
export JAVA_HOME="$CONDA_PREFIX_FROM_BUSCO/lib/jvm"
export PATH="$JAVA_HOME/bin:$PATH"
# Help the loader find the right libjli.so (avoids JLI_StringDup issues)
export LD_LIBRARY_PATH="$JAVA_HOME/lib/jli:$JAVA_HOME/lib"
# Tell BBTools exactly which java to use
export BBMAP_JAVA="$JAVA_HOME/bin/java"

# Final sanity checks (fail early with clear errors)
if ! "$BBMAP_JAVA" -version >/dev/null 2>&1; then
  echo "ERROR: Java in $BBMAP_JAVA cannot run. Check OpenJDK install in $CONDA_PREFIX_FROM_BUSCO."
  exit 127
fi
if ! "$STATS_SH" --version >/dev/null 2>&1; then
  echo "ERROR: 'stats.sh' exists but still can't run (Java/permissions)."
  echo "DEBUG busco=$(command -v busco)"
  echo "DEBUG stats.sh=$STATS_SH"
  echo "DEBUG BBMAP_JAVA=$BBMAP_JAVA"
  exit 127
fi

# Root to hold BUSCO downloads/cache
buscoroot="../Data/BUSCO"
buscoroot=$(realpath "$buscoroot")
download_path="${buscoroot}/busco_downloads"
mkdir -p "$download_path"

echo "> What is the genus name?"
read -r genus
mkdir -p "../Species/$genus"

echo "> What is the species name? $genus __ :"
read -r species_name
species_name=${species_name#$genus }  # strip leading "Genus "

# Save inputs (so the main pipeline can pick them up even if you run this standalone)
upsert_kv "GENUS" "$genus" "$PIPELINE_INPUTS"
upsert_kv "SPECIES_NAME" "$species_name" "$PIPELINE_INPUTS"
upsert_kv "FULL_SPECIES" "$genus/$species_name" "$PIPELINE_INPUTS"

echo "> Downloading genomes for $genus $species_name"
output_dir="../Species/$genus/$species_name"
zip_file="$output_dir/$species_name.zip"
mkdir -p "$output_dir"

# -------------------
# Download & extract
# -------------------
if [ ! -f "$zip_file" ]; then
  datasets download genome taxon "$genus $species_name" \
    --include genome,protein,gbff,gff3 --filename "$zip_file"
  echo "> Downloaded genome data to $zip_file"
else
  echo "> ZIP file already exists. Skipping download."
fi

if [ ! -d "$output_dir/ncbi_dataset" ]; then
  python3 - <<PY
import zipfile, sys
zipfile.ZipFile("$zip_file").extractall("$output_dir")
print("> Extraction complete! Files are in $output_dir")
PY
else
  echo "> Extraction already completed. Skipping."
fi

# -------------------
# Organize files
# -------------------
if [ ! -d "$output_dir/genome" ]; then
  mkdir -p "$output_dir/genome" "$output_dir/protein" "$output_dir/genbank" "$output_dir/gff3"

  if [ -f "$output_dir/ncbi_dataset/data/assembly_data_report.jsonl" ]; then
    mv "$output_dir/ncbi_dataset/data/assembly_data_report.jsonl" "$output_dir/"
    echo "> Moved assembly_data_report.jsonl to $output_dir"
  fi

  for asm_dir in "$output_dir/ncbi_dataset/data"/*; do
    [ -d "$asm_dir" ] || continue
    acc=$(basename "$asm_dir")  # e.g. GCF_000184455.2 or GCA_000269785.2

    # find the genomic fna and keep its original full name
    fna=$(find "$asm_dir" -maxdepth 1 -type f -name "*_genomic.fna" -print -quit)
    if [ -n "$fna" ]; then
      fna_name=$(basename "$fna")
      cp "$fna" "$output_dir/genome/${fna_name}"

      # Derive base name from genome filename (remove _genomic.fna)
      base_name="${fna_name%_genomic.fna}"

      # protein - use same naming as genome file
      prot="$asm_dir/protein.faa"
      [ -f "$prot" ] && cp "$prot" "$output_dir/protein/${base_name}.faa"

      # genbank - use full name
      gbk="$asm_dir/genomic.gbff"
      [ -f "$gbk" ] && cp "$gbk" "$output_dir/genbank/${base_name}.gbk"

      # gff3 - use full name
      gff="$asm_dir/genomic.gff"
      [ -f "$gff" ] && cp "$gff" "$output_dir/gff3/${base_name}.gff3"
    fi
  done

  rm -rf "$output_dir/ncbi_dataset"
  echo "> Moved FNA->genome, FAA->protein, GBK->genbank, GFF3->gff3 and removed ncbi_dataset"

  # -------------------
  # Auto RNA-seq discovery using assembly_data_report.jsonl
  # -------------------
  genome_dir="${output_dir}/genome"

  jsonl="$output_dir/assembly_data_report.jsonl"
  RNA_ROOT="$output_dir/rna"
  mkdir -p "$RNA_ROOT"   # always make top-level rna folder

  if [[ -f "$jsonl" ]]; then
      echo "> Extracting BioSamples from assembly_data_report.jsonl"

      # Extract AssemblyAccession <TAB> BioSample into biosamples.tsv
      biosample_tsv="$output_dir/biosamples.tsv"
      python3 - <<PY
import json

jsonl_file = "${jsonl}"
out_file = "${biosample_tsv}"

with open(jsonl_file) as f, open(out_file, "w") as out:
    for line in f:
        if not line.strip():
            continue
        rec = json.loads(line)
        acc = rec.get("accession")
        bios = rec.get("assemblyInfo", {}).get("biosample", {}).get("accession")
        if acc and bios:
            out.write(f"{acc}\t{bios}\n")
PY
      echo "> Wrote BioSample map to $biosample_tsv"

      # Check for Entrez Direct + download tool (curl/wget). No SRA Toolkit needed.
      if ! command -v esearch >/dev/null; then
          echo "WARNING: 'esearch' not found (entrez-direct). Skipping RNA-auto-discovery."
      elif ! command -v curl >/dev/null && ! command -v wget >/dev/null; then
          echo "WARNING: Neither 'curl' nor 'wget' found. Only writing run lists (no downloads)."
          DOWNLOAD_FROM_ENA=0
      else
          DOWNLOAD_FROM_ENA=1
          echo "> Searching SRA for RNA-Seq linked to BioSamples (downloading from ENA)..."
      fi

      while IFS=$'\t' read -r asm biosample <&3; do
          [[ -z "$biosample" ]] && continue

          # Find the full genome filename to get consistent naming
          genome_file=$(compgen -G "$genome_dir/${asm}*_genomic.fna" | head -1)
          if [[ -n "$genome_file" ]]; then
              full_name=$(basename "$genome_file")
              base_name="${full_name%_genomic.fna}"
          else
              base_name="$asm"
          fi

          echo "  [$base_name] BioSample=$biosample: querying RNA-Seq..."

          # Temporary runlist (not yet inside rna/$base_name)
          temp_runlist="$(mktemp)"

          # Query SRA (via Entrez) for this BioSample & filter only RNA-Seq runs
          esearch -db biosample -query "$biosample" \
            | elink -target sra \
            | efetch -format runinfo \
            | awk -F',' '
                  NR==1 {
                      for (i=1;i<=NF;i++) {
                          if ($i=="Run") r=i;
                          if ($i=="LibraryStrategy") ls=i;
                      }
                      next
                  }
                  $ls=="RNA-Seq" {print $r}
            ' > "$temp_runlist"

          # If no RNA-Seq runs, cleanup temp file and DO NOT create a subfolder
          if [[ ! -s "$temp_runlist" ]]; then
              echo "    -> No RNA-Seq found."
              rm -f "$temp_runlist"
              continue
          fi

          # Create subfolder only for genomes that actually have RNA-Seq
          strain_rna_dir="$RNA_ROOT/$base_name"
          mkdir -p "$strain_rna_dir"

          # Move temp runlist into the strain folder
          runlist="$strain_rna_dir/runs.txt"
          mv "$temp_runlist" "$runlist"

          echo "    -> Found $(wc -l < "$runlist") RNA-Seq run(s)."

          # Optionally download each run from ENA
          if [[ "${DOWNLOAD_FROM_ENA:-0}" -eq 1 ]]; then
              while read -r run; do
                  [[ -z "$run" ]] && continue

                  echo "      Downloading $run from ENA..."
                  ena_url="https://www.ebi.ac.uk/ena/browser/api/fastq/${run}?download=txt"

                  if command -v curl >/dev/null 2>&1; then
                      if curl -L "$ena_url" -o "$strain_rna_dir/${run}.fastq.gz"; then
                          echo "        -> Downloaded ${run}.fastq.gz"
                      else
                          echo "        ERROR: ENA download failed for $run, skipping..."
                      fi
                  elif command -v wget >/dev/null 2>&1; then
                      if wget -O "$strain_rna_dir/${run}.fastq.gz" "$ena_url"; then
                          echo "        -> Downloaded ${run}.fastq.gz"
                      else
                          echo "        ERROR: ENA download failed for $run, skipping..."
                      fi
                  fi
              done < "$runlist"
          else
              echo "    -> Skipping download (no curl/wget)."
          fi

          echo "    -> Finished RNA-seq for $asm"

      done 3< "$biosample_tsv"
  fi


else
  echo "> Genome file already exists. Skipping file categorization. (If needed, delete $output_dir/genome to redo)"
fi


# -------------------
# Filter GCA/GCF duplicates: Keep GCF, remove GCA
# -------------------
echo "> Checking for GCA/GCF duplicate genomes..."
genome_dir="${output_dir}/genome"
protein_dir="${output_dir}/protein"
genbank_dir="${output_dir}/genbank"
gff3_dir="${output_dir}/gff3"

# Build a map of base accession numbers (without version) to files
declare -A accession_map
declare -A gca_files
declare -A gcf_files

for fna in "$genome_dir"/*.fna; do
  [[ ! -f "$fna" ]] && continue

  filename=$(basename "$fna")
  # Extract accession: GCA_123456789 or GCF_123456789 (without version)
  if [[ "$filename" =~ ^(GC[AF])_([0-9]+) ]]; then
    prefix="${BASH_REMATCH[1]}"  # GCA or GCF
    number="${BASH_REMATCH[2]}"  # numeric part
    base_acc="${number}"         # unique identifier

    if [[ "$prefix" == "GCA" ]]; then
      gca_files["$base_acc"]="$filename"
    elif [[ "$prefix" == "GCF" ]]; then
      gcf_files["$base_acc"]="$filename"
    fi
  fi
done

# Remove GCA files where GCF exists
removed_count=0
for base_acc in "${!gca_files[@]}"; do
  if [[ -n "${gcf_files[$base_acc]:-}" ]]; then
    gca_file="${gca_files[$base_acc]}"
    gcf_file="${gcf_files[$base_acc]}"

    echo "  Found duplicate: GCA vs GCF for accession $base_acc"
    echo "    Keeping:  $gcf_file (RefSeq)"
    echo "    Removing: $gca_file (GenBank)"

    # Remove the base name without _genomic suffix for finding related files
    gca_base="${gca_file%_genomic.fna}"

    # Remove genome file
    rm -f "$genome_dir/$gca_file"

    # Remove associated protein, genbank, and gff3 files
    rm -f "$protein_dir/${gca_base}.faa"
    rm -f "$genbank_dir/${gca_base}.gbk"
    rm -f "$gff3_dir/${gca_base}.gff3"

    removed_count=$((removed_count + 1))
  fi
done

if [[ $removed_count -gt 0 ]]; then
  echo "> Removed $removed_count GCA genome(s) that have GCF equivalents"

  # Clean up orphaned BUSCO directories for removed GCA files
  busco_output_dir="${output_dir}/busco_output"
  if [[ -d "$busco_output_dir" ]]; then
    echo "> Cleaning up orphaned BUSCO directories..."
    orphaned_count=0
    for busco_dir in "$busco_output_dir"/*; do
      [[ ! -d "$busco_dir" ]] && continue
      busco_name=$(basename "$busco_dir")
      # Extract accession from BUSCO directory name (e.g., GCA_009684875.1 from GCA_009684875.1_AoryTK5_1.0_genomic.fna)
      busco_acc=$(echo "$busco_name" | sed -E 's/^(GC[AF]_[0-9]+\.[0-9]+).*/\1/')
      # Check if any genome file with this accession exists
      if ! compgen -G "$genome_dir/${busco_acc}*" >/dev/null; then
        echo "  Removing orphaned BUSCO results: $busco_name"
        rm -rf "$busco_dir" || echo "WARNING: Failed to remove $busco_dir"
        orphaned_count=$((orphaned_count + 1))
      fi
    done
    if [[ $orphaned_count -gt 0 ]]; then
      echo "> Removed $orphaned_count orphaned BUSCO directory(ies)"
    fi
  fi
else
  echo "> No GCA/GCF duplicates found"
fi

echo "> BUSCO root set to $buscoroot"

# -------------------
# Helpers for BUSCO parsing/cleanup
# -------------------
# Extract lineage from the *filename* of short_summary.specific.* if present,
# otherwise fall back to any short_summary.* (still from filename), and *finally*
# fall back to scanning file contents.
get_best_lineage() {
  local run_dir="$1"
  local f bn lineage

  # a) specific summary wins
  if compgen -G "$run_dir/short_summary.specific.*.txt" >/dev/null; then
    f=$(ls -1 "$run_dir"/short_summary.specific.*.txt | tail -n1)
    bn=$(basename "$f")
    # short_summary.specific.<LINEAGE>.<GENOME>.txt
    lineage=$(sed -E 's/^short_summary\.specific\.([^.]+)\..*$/\1/' <<<"$bn")
    echo "$lineage"
    return 0
  fi

  # b) otherwise try any short_summary.* from filename
  if compgen -G "$run_dir/short_summary.*.txt" >/dev/null; then
    f=$(ls -1 "$run_dir"/short_summary.*.txt | tail -n1)
    bn=$(basename "$f")
    # short_summary.(generic|specific).<LINEAGE>.<GENOME>.txt  OR  short_summary.<LINEAGE>.<GENOME>.txt
    lineage=$(sed -E 's/^short_summary\.(generic|specific)\.([^.]+)\..*$/\2/; t; s/^short_summary\.([^.]+)\..*$/\1/' <<<"$bn")
    echo "$lineage"
    return 0
  fi

  # c) last resort: grep inside (rarely needed)
  if compgen -G "$run_dir/short_summary*.txt" >/dev/null; then
    f=$(ls -1 "$run_dir"/short_summary*.txt | tail -n1)
    grep -oE '[A-Za-z0-9_]+_odb[0-9]+' "$f" | tail -n1 || true
  fi
}

# After a BUSCO run, keep only the most-specific results
keep_only_most_specific_busco() {
  local run_dir="$1"
  local best_lineage
  best_lineage="$(get_best_lineage "$run_dir")"

  if [ -z "$best_lineage" ]; then
    echo "> [cleanup] Could not determine best lineage in $run_dir — skipping cleanup."
    return 0
  fi

  echo "> [cleanup] Keeping lineage: $best_lineage"

  # 1) remove run_* directories that aren't the best lineage
  shopt -s nullglob
  for d in "$run_dir"/run_*; do
    [ -d "$d" ] || continue
    local bn; bn="$(basename "$d")"
    [[ "$bn" == "run_${best_lineage}" ]] || rm -rf "$d"
  done

  # 2) remove summary files that aren't the best lineage
  for f in "$run_dir"/short_summary.*.txt "$run_dir"/short_summary.*.json; do
    [[ "$f" =~ $best_lineage ]] || rm -f "$f"
  done
  shopt -u nullglob
}

# -------------------
# BUSCO (detect once, reuse lineage, save to pipeline_inputs.txt)
# -------------------
busco_dir="${output_dir}/busco_output"
mkdir -p "$busco_dir"

# Optional: allow forcing a lineage (e.g. BUSCO_LINEAGE=penicillium_odb12)
BUSCO_LINEAGE="${BUSCO_LINEAGE:-}"

# Per-species cache file
lineage_lock="${busco_dir}/lineage.lock"

# Common BUSCO args (drop --limit; ignored in genome modes)
common_busco_args=(-m genome -c "$BUSCO_THREADS" -f --download_path "$download_path")
[ "${BUSCO_OFFLINE}" = "1" ] && common_busco_args+=(--offline)

echo "> Running BUSCO on genomes (detect once, reuse, and save lineage)..."
fasta_dir="${output_dir}/genome"
shopt -s nullglob
mapfile -t genomes < <(ls -1 "$fasta_dir"/*.fna 2>/dev/null | sort)
shopt -u nullglob

if [ ${#genomes[@]} -eq 0 ]; then
  echo "> No genomes found in $fasta_dir — nothing to do."
else
  # Load cached or forced lineage first
  if [ -z "$BUSCO_LINEAGE" ] && [ -s "$lineage_lock" ]; then
    BUSCO_LINEAGE="$(cat "$lineage_lock")"
    echo "> Loaded cached lineage: $BUSCO_LINEAGE"
  fi
  if [ -z "$BUSCO_LINEAGE" ] && grep -q '^BUSCO_LINEAGE=' "$PIPELINE_INPUTS" 2>/dev/null; then
    # If launched standalone, reuse previously saved lineage
    # shellcheck disable=SC1090
    source "$PIPELINE_INPUTS"
  fi

  if [ -z "${BUSCO_LINEAGE:-}" ]; then
    first_genome="${genomes[0]}"
    first_name=$(basename "$first_genome")
    first_run_dir="${busco_dir}/${first_name}"

    if compgen -G "${first_run_dir}/short_summary*.txt" >/dev/null; then
      # Prefer specific summary if present, else the last short_summary
      if compgen -G "${first_run_dir}/short_summary.specific.*.txt" >/dev/null; then
        summary_file=$(ls -1 "${first_run_dir}"/short_summary.specific.*.txt | tail -n1)
      else
        summary_file=$(ls -1 "${first_run_dir}"/short_summary*.txt | tail -n1)
      fi
      BUSCO_LINEAGE="$(get_best_lineage "$first_run_dir")"
      keep_only_most_specific_busco "$first_run_dir"
    else
      case "${BUSCO_MODE:-auto}" in
        auto)       mode_flag="--auto-lineage" ;;
        auto_euk)   mode_flag="--auto-lineage-euk" ;;
        auto_prok)  mode_flag="--auto-lineage-prok" ;;
        *)          mode_flag="--auto-lineage" ;;
      esac
      echo "> Detecting lineage from first genome: $first_name (mode: ${BUSCO_MODE:-auto})"
      mkdir -p "$first_run_dir"
      busco -i "$first_genome" -o "$first_run_dir" $mode_flag "${common_busco_args[@]}"
      # Prefer specific summary; fallback to last summary found
      if compgen -G "${first_run_dir}/short_summary.specific.*.txt" >/dev/null; then
        summary_file=$(ls -1 "${first_run_dir}"/short_summary.specific.*.txt | tail -n1)
      else
        summary_file=$(ls -1 "${first_run_dir}"/short_summary*.txt | tail -n1 || true)
      fi
      [ -n "${summary_file:-}" ] && BUSCO_LINEAGE="$(get_best_lineage "$first_run_dir")"
      keep_only_most_specific_busco "$first_run_dir"
    fi

    if [ -n "${BUSCO_LINEAGE:-}" ]; then
      echo "$BUSCO_LINEAGE" > "$lineage_lock"
      upsert_kv "BUSCO_LINEAGE" "$BUSCO_LINEAGE" "$PIPELINE_INPUTS"
      echo "> Chosen lineage cached and saved to pipeline_inputs.txt: $BUSCO_LINEAGE"
    else
      echo "> WARNING: Could not parse lineage; will fall back to auto-lineage per genome."
    fi
  else
    echo "$BUSCO_LINEAGE" > "$lineage_lock"
    upsert_kv "BUSCO_LINEAGE" "$BUSCO_LINEAGE" "$PIPELINE_INPUTS"
    echo "> Using forced lineage and saved to pipeline_inputs.txt: $BUSCO_LINEAGE"
  fi

  for genome_file in "${genomes[@]}"; do
    genome_name=$(basename "$genome_file")
    run_dir="${busco_dir}/${genome_name}"

    if compgen -G "${run_dir}/short_summary*.txt" >/dev/null; then
      echo "> BUSCO output already exists for $genome_name. Ensuring cleanup..."
      keep_only_most_specific_busco "$run_dir"
      continue
    fi

    echo "> Running BUSCO for $genome_name..."
    mkdir -p "$run_dir"

    if [ -n "${BUSCO_LINEAGE:-}" ]; then
      busco -i "$genome_file" -o "$run_dir" -l "$BUSCO_LINEAGE" "${common_busco_args[@]}"
    else
      case "${BUSCO_MODE:-auto}" in
        auto)       mode_flag="--auto-lineage" ;;
        auto_euk)   mode_flag="--auto-lineage-euk" ;;
        auto_prok)  mode_flag="--auto-lineage-prok" ;;
        *)          mode_flag="--auto-lineage" ;;
      esac
      busco -i "$genome_file" -o "$run_dir" $mode_flag "${common_busco_args[@]}"
    fi

    # After each run, keep only the most specific lineage material
    keep_only_most_specific_busco "$run_dir"
    echo "> BUSCO analysis complete for $genome_name!"
  done
fi

# Remove any busco*.log files next to the script
find "$(dirname "$0")" -maxdepth 1 -type f -name "busco*.log" -delete
echo "> BUSCO analysis complete! Results are in $busco_dir"

# -------------------
# Merge BUSCO -> JSONL
# -------------------
echo "> Merging BUSCO results into JSONL..."
jsonl_file="${output_dir}/assembly_data_report.jsonl"
updated_jsonl="${output_dir}/busco_updated.jsonl"

declare -A busco_scores
declare -A busco_lineages

shopt -s nullglob
for summary_file in "${busco_dir}"/*/short_summary*.txt; do
  [ -f "$summary_file" ] || continue

  run_dir=$(dirname "$summary_file")
  genome_name=$(basename "$run_dir")               # e.g. GCA_..._genomic.fna
  accession_id=$(echo "$genome_name" | sed -E 's/^([^_]+_[^_]+).*/\1/')

  complete=$(grep -E "^[[:space:]]*[0-9]+[[:space:]]+Complete BUSCOs" "$summary_file" | awk '{print $1}')
  fragmented=$(grep -E "^[[:space:]]*[0-9]+[[:space:]]+Fragmented BUSCOs" "$summary_file" | awk '{print $1}')
  missing=$(grep -E "^[[:space:]]*[0-9]+[[:space:]]+Missing BUSCOs" "$summary_file" | awk '{print $1}')
  total=$(grep -E "^[[:space:]]*[0-9]+[[:space:]]+Total BUSCO groups" "$summary_file" | awk '{print $1}')
  lineage=$(grep -oE '[A-Za-z0-9_]+_odb[0-9]+' "$summary_file" | tail -n1)

  busco_scores["$accession_id"]="$complete,$fragmented,$missing,$total"
  [ -n "$lineage" ] && busco_lineages["$accession_id"]="$lineage"
  printf "> Processed BUSCO summary for %s (accession: %s, lineage: %s)\n" "$genome_name" "$accession_id" "${lineage:-unknown}"
done
shopt -u nullglob

echo "> BUSCO scores collected."

updated_jsonl=$(realpath "$updated_jsonl")
export updated_jsonl output_dir

python3 - <<'PY'
import json, os, re, pandas as pd
from glob import glob
from pathlib import Path

jsonl_file   = os.path.join(os.environ["output_dir"], "assembly_data_report.jsonl")
updated_jsonl= os.environ["updated_jsonl"]
csv_output   = os.path.join(os.environ["output_dir"], "metadata.csv")
busco_dir    = Path(os.environ["output_dir"]) / "busco_output"

# ---------------- BUSCO parsing ----------------
def parse_summary(sf):
    d = {"complete":None,"fragmented":None,"missing":None,"total":None,"lineage":None}
    with open(sf) as f:
        for line in f:
            line=line.strip()
            if "Complete BUSCOs" in line and line.split():
                d["complete"]=int(line.split()[0])
            elif "Fragmented BUSCOs" in line and line.split():
                d["fragmented"]=int(line.split()[0])
            elif "Missing BUSCOs" in line and line.split():
                d["missing"]=int(line.split()[0])
            elif "Total BUSCO groups" in line and line.split():
                d["total"]=int(line.split()[0])
            elif ("The lineage dataset is:" in line
                  or "Results from dataset" in line
                  or "Results from generic domain" in line):
                m = re.search(r'([A-Za-z0-9_]+_odb\d+)', line)
                if m:
                    d["lineage"] = m.group(1)
    return d

busco_by_accession = {}
for sf in glob(str(busco_dir / "*" / "short_summary*.txt")):
    run_dir = Path(sf).parent
    genome_name = run_dir.name                          # e.g. GCA_..._genomic.fna
    accession_id = "_".join(genome_name.split("_")[:2]) # GCA_xxx.x
    busco_by_accession[accession_id] = parse_summary(sf)

# ---------------- helpers ----------------
MISSING_TOKENS = {"missing","not applicable","na","none","unknown","nan",""}

def normkey(s: str) -> str:
    return re.sub(r'[\s_\-]+','', s.strip().lower()) if isinstance(s,str) else s

def attrs_to_map(attrs):
    m = {}
    for item in attrs or []:
        name = item.get("name")
        if not name: 
            continue
        m[normkey(name)] = item.get("value")
    return m

def is_meaningful(v):
    if v is None:
        return False
    s = str(v).strip()
    return s.lower() not in MISSING_TOKENS

def unique_join(values):
    """Order-preserving de-dup + join with comma+space."""
    seen = set()
    out = []
    for v in values:
        if not is_meaningful(v):
            continue
        s = str(v).strip()
        key = s.lower()
        if key in seen:
            continue
        seen.add(key)
        out.append(s)
    return ", ".join(out) if out else None

def first_non_empty(*vals):
    for v in vals:
        if is_meaningful(v):
            return str(v).strip()
    return None

def clean_comment_text(s):
    if not isinstance(s, str):
        return None
    s = s.strip().strip('"').strip("'")
    # NEW: collapse any runs of whitespace (incl. \n, \r, \t) to a single space
    s = re.sub(r'\s+', ' ', s)
    return s if s else None

def find_all_comments(obj):
    """
    Recursively collect string values from keys named 'comment' or 'comments'.
    - If value is a string -> capture
    - If value is a list -> capture any strings in it
    - If value is a dict -> recurse
    """
    out = []

    def walk(x):
        if isinstance(x, dict):
            for k, v in x.items():
                kl = k.lower()
                if kl in ("comment", "comments"):
                    if isinstance(v, str):
                        t = clean_comment_text(v)
                        if t:
                            out.append(t)
                    elif isinstance(v, list):
                        for vv in v:
                            if isinstance(vv, str):
                                t = clean_comment_text(vv)
                                if t:
                                    out.append(t)
                    # still walk nested objects
                walk(v)
        elif isinstance(x, list):
            for el in x:
                walk(el)

    walk(obj)
    return out

# ---------------- target columns ----------------
columns = [
    "Assembly Accession","Assembly Name","Organism Name","Organism Taxonomic ID",
    "ANI Check status","Organism Infraspecific Names Breed","Organism Infraspecific Names Strain",
    "Organism Infraspecific Names Cultivar","Organism Infraspecific Names Ecotype",
    "Organism Infraspecific Names Isolate","Organism Infraspecific Names Sex","Annotation Name",
    "Assembly Stats Total Sequence Length","Assembly Stats Total Number of Chromosomes",
    "Assembly Level","Assembly Type","Assembly Release Date","WGS project accession",
    "Assembly Stats Contig N50","Assembly Stats Scaffold N50",
    "Assembly Stats Number of Scaffolds","Assembly Stats Number of Contigs",
    "BUSCO Complete","BUSCO Fragmented","BUSCO Missing","BUSCO Total",
    "BUSCO Single Copy","BUSCO Duplicated","BUSCO Lineage",
    "Assembly Sequencing Tech","Assembly Submitter",
    "Assembly BioProject Accession","Assembly BioSample Accession",
    "Annotation Count Gene Total","Annotation Count Gene Protein-coding","Annotation Count Gene Pseudogene",
    "Type Material Display Text","CheckM Completeness","CheckM Contamination",
    "Isolation Source","Geo Location","Collection Date","Lat/Lon",
    "BioSample Description Title",   # <-- NEW COLUMN
    "Comment"
]

# ---------------- assemble rows ----------------
rows = []
with open(jsonl_file, "r") as f:
    for line in f:
        if not line.strip():
            continue
        rec = json.loads(line)
        acc = rec.get("accession")
        b   = busco_by_accession.get(acc, {})

        asm  = rec.get("assemblyInfo", {}) or {}
        org  = rec.get("organism", {}) or {}
        ann  = rec.get("annotationInfo", {}) or {}
        aSt  = rec.get("assemblyStats", {}) or {}
        wgs  = rec.get("wgsInfo", {}) or {}
        bios = asm.get("biosample", {}) or {}
        attrs= attrs_to_map(bios.get("attributes"))

        # ---- collect Isolation Source candidates (ALL), de-dup & join ----
        iso_candidates = [
            bios.get("isolationSource"),
            attrs.get("isolationsource"),
            attrs.get("envmedium"),
            attrs.get("envlocalscale"),
            attrs.get("envbroadscale")
        ]
        isolation_source_joined = unique_join(iso_candidates)

        # ---- Geo Location (first non-empty) ----
        geo_loc = first_non_empty(bios.get("geoLocName"), attrs.get("geolocname"))

        # ---- Collection Date (first non-empty) ----
        collection_date = first_non_empty(bios.get("collectionDate"), attrs.get("collectiondate"))

        # ---- Lat/Lon: collect ALL variants, de-dup & join (as strings) ----
        latlon_joined = unique_join([bios.get("latLon"), attrs.get("latlon")])

        # ---- Comments: recursively collect all strings from keys 'comment'/'comments', de-dup & join ----
        comments_joined = unique_join(find_all_comments(rec))

        # ---- BioSample description title ----
        biosample_title = (
            ((bios.get("description") or {}).get("title"))
            if isinstance(bios.get("description"), dict)
            else None
        )

        # Build row
        rows.append([
            rec.get("accession"),
            asm.get("assemblyName"),
            org.get("organismName"),
            org.get("taxId"),
            None,                                           # ANI Check status
            None,                                           # Breed
            (org.get("infraspecificNames") or {}).get("strain"),
            (org.get("infraspecificNames") or {}).get("cultivar"),
            (org.get("infraspecificNames") or {}).get("ecotype"),
            (org.get("infraspecificNames") or {}).get("isolate"),
            (org.get("infraspecificNames") or {}).get("sex"),
            ann.get("name"),
            aSt.get("totalSequenceLength"),
            aSt.get("totalNumberOfChromosomes"),
            asm.get("assemblyLevel"),
            asm.get("assemblyType"),
            ann.get("releaseDate"),
            wgs.get("wgsProjectAccession"),
            aSt.get("contigN50"),
            aSt.get("scaffoldN50"),
            aSt.get("numberOfScaffolds"),
            aSt.get("numberOfContigs"),
            b.get("complete"),
            b.get("fragmented"),
            b.get("missing"),
            b.get("total"),
            None, None,                                     # BUSCO Single/Duplicated (not in short summary)
            b.get("lineage"),
            asm.get("sequencingTech"),
            asm.get("submitter"),
            asm.get("bioprojectAccession"),
            bios.get("accession"),
            (ann.get("stats", {}) or {}).get("geneCounts", {}).get("total"),
            (ann.get("stats", {}) or {}).get("geneCounts", {}).get("proteinCoding"),
            (ann.get("stats", {}) or {}).get("geneCounts", {}).get("pseudogene"),
            None,                                           # Type Material Display Text
            None,                                           # CheckM Completeness
            None,                                           # CheckM Contamination
            isolation_source_joined,                        # Isolation Source (comma-joined uniques)
            geo_loc,                                        # Geo Location
            collection_date,                                # Collection Date
            latlon_joined,                                  # Lat/Lon (comma-joined uniques)
            biosample_title,                                # <-- NEW FIELD
            comments_joined                                 # Comments
        ])

df = pd.DataFrame(rows, columns=columns)
df.to_json(updated_jsonl, orient="records", lines=True)
df.to_csv(csv_output, index=False)
print(f"> BUSCO updated metadata saved to {updated_jsonl} and {csv_output}")
PY



# -------------------
# Cleanup zip if CSV exists
# -------------------
echo "> Cleaning up..."
csv_output="${output_dir}/metadata.csv"
if [ -f "$csv_output" ]; then
  rm -f "$zip_file"
  echo "> Deleted $zip_file"
else
  echo "> Keeping $genus $species_name zip file for retry"
fi

echo "> Process complete"
