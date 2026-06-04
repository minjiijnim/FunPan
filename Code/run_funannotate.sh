#!/usr/bin/env bash
set -euo pipefail

# ─── Prompt ─────────────────────────────────────────────────────────────────
echo -n "> Enter species name (genus/species format): "
read species

# ─── Settings ───────────────────────────────────────────────────────────────
GENOME_DIR=$(realpath "../Species/$species/filtered_genome")
PROTEIN_DIR=$(realpath "../Species/$species/filtered_protein")
GFF3_DIR=$(realpath "../Species/$species/gff3")
RNA_DIR=$(realpath "../Species/$species/rna")
OUTDIR=$(realpath "../Species/$species/funannotate_output")
FUNANNOTATE_DB=$(realpath ../Data/Funannotate)

# Threading: match download_genome_and_BUSCO.sh logic
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

export FUNANNOTATE_DB

mkdir -p "$FUNANNOTATE_DB"

# Get the directory where this script is located
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ─── Funannotate DB check ───────────────────────────────────────────────────
if [ -z "$(ls -A "$FUNANNOTATE_DB" 2>/dev/null)" ]; then
    echo "No Funannotate database found at $FUNANNOTATE_DB. Setting up database..."
    funannotate setup -d "$FUNANNOTATE_DB" -i all
    
    # Cleanup setup log if it exists
    setup_log="${SCRIPT_DIR}/funannotate-setup.log"
    if [ -f "$setup_log" ]; then
        rm -f "$setup_log"
        echo "Removed $setup_log"
    fi
fi

# ─── Ensure funannotate aux scripts are executable ──────────────────────────
if [[ -n "${CONDA_PREFIX:-}" ]]; then
    aux_dir="$CONDA_PREFIX/lib/python*/site-packages/funannotate/aux_scripts"
    for d in $aux_dir; do
        if [[ -d "$d" ]]; then
            echo "Fixing permissions in $d"
            find "$d" -type f -name "*.py" ! -perm -111 -exec chmod 755 {} \;
        fi
    done
else
    echo "WARNING: No active conda environment detected; skipping aux script permission fix." >&2
fi

# ─── Evidence Modeler (EVM) discovery ───────────────────────────────────────
if [[ -z "${CONDA_PREFIX:-}" ]]; then
  echo "ERROR: No conda environment detected. Please activate your env before running." >&2
  exit 1
fi

evm_dir="$(find "$CONDA_PREFIX/opt" -maxdepth 1 -type d -name 'evidencemodeler-*' | head -n1 || true)"
if [[ -z "$evm_dir" ]]; then
  echo "ERROR: Could not locate EVM under $CONDA_PREFIX/opt" >&2
  exit 1
fi

export EVM_HOME="$evm_dir"
export PATH="$EVM_HOME:$EVM_HOME/EVM:$EVM_HOME/EvmUtils:$EVM_HOME/bin:$PATH"

if ! command -v evidence_modeler.pl >/dev/null 2>&1; then
  echo "ERROR: evidence_modeler.pl not found even after setting PATH." >&2
  exit 1
fi
echo "Using EVM_HOME=${EVM_HOME}"

mkdir -p "$OUTDIR"

# ─── Step 1: Rename genome FASTA headers and update GFF3 ────────────────────
echo " Renaming FASTA headers to scaffold_N and updating GFF3…"

for f in "$GENOME_DIR"/*.fna; do
  [[ "$f" == *_renamed.fna ]] && continue

  base_orig=$(basename "${f%.fna}")
  core=${base_orig%_genomic}
  mapfile="${f%.fna}.map"
  renamed_fna="${f%.fna}_renamed.fna"

  rm -f "$mapfile"

  awk '
    BEGIN { sc=0 }
    /^>/ {
      sc++
      orig = substr($0, 2)
      split(orig, a, " ")
      print a[1] "\tscaffold_" sc > "'"$mapfile"'"
      print ">scaffold_" sc
      next
    }
    { print }
  ' "$f" > "$renamed_fna"
  echo "   ⇨ Created renamed FASTA: $renamed_fna"

  mkdir -p "$GFF3_DIR/mismatched_gffs" "$GFF3_DIR/non_evm_format"

  raw_gff="$GFF3_DIR/${core}.gff3"
  renamed_gff="$GFF3_DIR/${core}_renamed.gff3"

  if [[ -f "$raw_gff" ]]; then
    echo "   Validating GFF3 for: $core"
    gff_ids=$(awk '$1 !~ /^#/ { print $1 }' "$raw_gff" | sort -u)
    map_ids=$(cut -f1 "$mapfile" | sort -u)
    mismatches=$(comm -23 <(echo "$gff_ids") <(echo "$map_ids"))

    if [[ -n "$mismatches" ]]; then
      echo "   GFF3 contigs not in genome. Moving to mismatched_gffs/"
      mv "$raw_gff" "$GFF3_DIR/mismatched_gffs/"
      continue
    fi

    if ! grep -qP "\tmRNA\t" "$raw_gff"; then
      echo "   GFF3 missing mRNA features. Moving to non_evm_format/"
      mv "$raw_gff" "$GFF3_DIR/non_evm_format/"
      continue
    fi

    dup_ids=$(awk -F'\t' '
      $0 !~ /^#/ {
        match($9, /ID=([^;]+)/, m)
        if (m[1] != "") {
          if (++id[m[1]] > 1) print m[1]
        }
      }
    ' "$raw_gff" | sort -u)

    if [[ -n "$dup_ids" ]]; then
      echo "   GFF3 has duplicate feature IDs. Moving to non_evm_format/"
      mv "$raw_gff" "$GFF3_DIR/non_evm_format/"
      continue
    fi

    awk -F $'\t' -v OFS=$'\t' -v mapfile="$mapfile" '
      BEGIN {
        while ((getline < mapfile) > 0) {
          split($0, a, "\t"); map[a[1]] = a[2]
        }
      }
      /^#/ { print; next }
      { if ($1 in map) $1 = map[$1]; print }
    ' "$raw_gff" > "$renamed_gff"

    echo "   Renamed GFF3 created: $renamed_gff"
  else
    echo "   No GFF3 found for $raw_gff"
  fi
done

# ─── Step 2: Clean contigs BEFORE masking (funannotate clean) ───────────────
# produce a cleaned FASTA per genome; skip if already present
echo "Cleaning assemblies (funannotate clean)…"
for genome in "$GENOME_DIR"/*_renamed.fna; do
  base=$(basename "$genome" .fna)
  core=${base%_genomic_renamed}
  output_dir="$OUTDIR/$base"
  mkdir -p "$output_dir"

  cleaned_fa="$output_dir/${base}.clean.fa"     # [ADDED] cleaned output path

  if [[ -s "$cleaned_fa" ]]; then
    echo "   [$core] Found existing cleaned FASTA: $cleaned_fa"
  else
    echo "   [$core] funannotate clean -i $genome -o $cleaned_fa"
    funannotate clean -i "$genome" -o "$cleaned_fa"
  fi
done


# ─── Step 3: Mask repeats WITHOUT Dfam (RepeatModeler → RepeatMasker; shared lib) ────────
echo "Running de novo masking via RepeatModeler + RepeatMasker (shared library)…"

# 0) Tool sanity-check first
for t in RepeatModeler RepeatMasker rmblastn trf; do
  command -v "$t" >/dev/null || { echo "ERROR: $t not found. Please activate the conda environment with RepeatModeler installed."; exit 1; }
done

# Resolve BuildDatabase even if it isn't directly on PATH
RM_BIN="$(dirname "$(command -v RepeatModeler)")"
BUILDDB_BIN="${RM_BIN}/BuildDatabase"
if ! command -v BuildDatabase >/dev/null 2>&1; then
  [[ -x "$BUILDDB_BIN" ]] || { echo "ERROR: BuildDatabase not found (checked PATH and $BUILDDB_BIN)"; exit 1; }
fi

# Detect RepeatModeler threading flag (v2: -threads; legacy: -pa)
# Prefer feature-detection over version parsing to be robust across builds
if RepeatModeler -h 2>&1 | grep -q -- '-threads'; then
  RM_HAS_THREADS=1
else
  RM_HAS_THREADS=0
fi

# One shared RepeatModeler library per species
SHARED_REP_DIR="$FUNANNOTATE_DB/repeats"
mkdir -p "$SHARED_REP_DIR"
SHARED_LIB="$SHARED_REP_DIR/${species//\//_}.RepeatModeler.lib.fa"

# 1) Build shared library once (from the first genome; prefer cleaned)
built_shared=0
if [[ ! -s "$SHARED_LIB" ]]; then
  first_genome=""
  for g in "$GENOME_DIR"/*_renamed.fna; do
    base_g="$(basename "$g" .fna)"
    cleaned="$OUTDIR/$base_g/${base_g}.clean.fa"
    if [[ -s "$cleaned" ]]; then first_genome="$cleaned"; else first_genome="$g"; fi
    [[ -n "$first_genome" ]] && break
  done
  [[ -n "$first_genome" ]] || { echo "ERROR: No genomes found to train RepeatModeler library." >&2; exit 1; }

  # sample/outdir name (avoid dirs named like '*.fa')
  if [[ "$first_genome" == *".clean.fa" ]]; then
    sample="$(basename "$(dirname "$first_genome")")"
    outdir="$OUTDIR/$sample"
  else
    fbase="$(basename "$first_genome")"
    sample="${fbase%.fna}"; sample="${sample%.fa}"; sample="${sample%.fasta}"
    outdir="$OUTDIR/$sample"
  fi
  mkdir -p "$outdir/logfiles"

  echo "   [$sample] Building shared RepeatModeler library at: $SHARED_LIB"
  (
    cd "$outdir" || exit 1

    # clean up any partial RM_* from previous failed attempts
    rm -rf RM_* 2>/dev/null || true

    RMDB="${sample}_rmdb"
    "${BUILDDB_BIN:-BuildDatabase}" -name "$RMDB" "$first_genome" |& tee builddatabase.log

    # RepeatModeler 2.x: use -threads (total threads). Do NOT use -pa.
    RM_TOTAL_THREADS=${THREADS:-1}
    (( RM_TOTAL_THREADS < 1 )) && RM_TOTAL_THREADS=1

    echo "   [$sample] RepeatModeler -database ${RMDB} -threads ${RM_TOTAL_THREADS} -LTRStruct"
    RepeatModeler -database "$RMDB" -threads "$RM_TOTAL_THREADS" -LTRStruct |& tee repeatmodeler.log

    found_lib=$(ls RM_*/consensi.fa.classified 2>/dev/null | tail -n1 || true)
    [[ -n "$found_lib" ]] || { echo "ERROR: RepeatModeler produced no consensi.fa.classified"; exit 1; }
    cp -f "$found_lib" "$SHARED_LIB"
  )

  built_shared=1
else
  echo "   Using existing shared RepeatModeler library: $SHARED_LIB"
fi

# 2) Mask all genomes with the shared library (no FamDB involvement)
for genome in "$GENOME_DIR"/*_renamed.fna; do
  base="$(basename "$genome" .fna)"
  core="${base%_genomic_renamed}"
  output_dir="$OUTDIR/$base"
  mkdir -p "$output_dir/logfiles"

  cleaned_fa="$output_dir/${base}.clean.fa"
  mask_input="$genome"; [[ -s "$cleaned_fa" ]] && mask_input="$cleaned_fa"

  masked_fa="$output_dir/${base}.softmasked.fa"
  if [[ -s "$masked_fa" ]]; then
    echo "   [$core] Masked FASTA exists"
    continue
  fi

  echo "   [$core] RepeatMasker -lib (shared): $(basename "$SHARED_LIB")"
  (
    cd "$output_dir" || exit 1
    in_base="$(basename "$mask_input")"
    # RepeatMasker parallelism; avoid oversubscribing tiny boxes
    RM_PA=$(( THREADS>4 ? THREADS/2 : 1 ))
    # -no_is prevents IS-finding/DB use; keeps us strictly local to -lib
    RepeatMasker -lib "$SHARED_LIB" -e rmblast -xsmall -pa "$RM_PA" -no_is "$mask_input"
    [[ -s "${in_base}.masked" ]] && cp -f "${in_base}.masked" "$masked_fa"
  )

  # optional summary
  rm_tbl=$(ls "$output_dir"/*.tbl 2>/dev/null | head -n1 || true)
  if [[ -n "$rm_tbl" ]]; then
    echo "   [$core] Repeat summary ($(basename "$rm_tbl")):"
    grep -E "bases masked|Low complexity|Simple repeats|DNA|LINE|SINE|LTR|RC|Unknown" "$rm_tbl" || true
  fi
done

[[ $built_shared -eq 1 ]] && echo "Shared RepeatModeler library created: $SHARED_LIB"



# ─── Step 3.5: Find best quality proteome for protein evidence ─────
echo "Finding best quality proteome for protein evidence…"
best_reference_proteome=""
best_proteome_score=-999999
metadata_file="../Species/$species/metadata.csv"

# First, try to find best GCF (RefSeq) proteome
for faa in "$PROTEIN_DIR"/GCF_*.faa; do
  [[ ! -f "$faa" ]] && continue

  filename=$(basename "$faa")
  acc=$(echo "$filename" | sed -E 's/^(GCF_[0-9]+\.[0-9]+).*/\1/')

  if [[ -f "$metadata_file" ]]; then
    quality_score=$(awk -F',' -v acc="$acc" '
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
    ' "$metadata_file")

    if [[ -n "$quality_score" ]] && (( $(echo "$quality_score > $best_proteome_score" | bc -l) )); then
      best_proteome_score="$quality_score"
      best_reference_proteome="$faa"
    fi
  else
    best_reference_proteome="$faa"
    break
  fi
done

# If no GCF found, find best quality genome from any available proteome
if [[ -z "$best_reference_proteome" ]]; then
  echo "   No GCF proteome found. Searching for best quality genome in filtered set…"

  for faa in "$PROTEIN_DIR"/*.faa; do
    [[ ! -f "$faa" ]] && continue

    filename=$(basename "$faa" .faa)
    # Extract accession (works for both GCA_xxx.x_name and GCA_xxx.x formats)
    acc=$(echo "$filename" | sed -E 's/^(GC[AF]_[0-9]+\.[0-9]+).*/\1/')

    if [[ -f "$metadata_file" ]]; then
      quality_score=$(awk -F',' -v acc="$acc" '
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
      ' "$metadata_file")

      if [[ -n "$quality_score" ]] && (( $(echo "$quality_score > $best_proteome_score" | bc -l) )); then
        best_proteome_score="$quality_score"
        best_reference_proteome="$faa"
      fi
    else
      # If no metadata, use first available proteome
      best_reference_proteome="$faa"
      break
    fi
  done
fi

if [[ -n "$best_reference_proteome" ]]; then
  if [[ "$best_reference_proteome" == *"GCF_"* ]]; then
    echo "   Selected best GCF (RefSeq) proteome: $(basename "$best_reference_proteome")"
  else
    echo "   Selected best quality proteome: $(basename "$best_reference_proteome")"
  fi
else
  echo "   WARNING: No proteome found in $PROTEIN_DIR"
fi

# ─── Step 4: Run Funannotate predict (prefer soft-masked; else cleaned) ─────
echo "Starting Funannotate annotation…"
for genome in "$GENOME_DIR"/*_renamed.fna; do
  base=$(basename "$genome" .fna)
  core=${base%_genomic_renamed}
  faa="$PROTEIN_DIR/$core.faa"
  renamed_gff="$GFF3_DIR/${core}_renamed.gff3"
  filtered_gff="$GFF3_DIR/${core}_renamed.filtered.gff3"
  output_dir="$OUTDIR/$base"

  echo; echo "Processing $core"

  if [[ -d "$output_dir/predict_results" && -n "$(ls -A "$output_dir/predict_results" 2>/dev/null)" ]]; then
    echo "   Skipping $core — already completed"
    continue
  fi
  mkdir -p "$output_dir"

  masked_fa="$output_dir/${base}.softmasked.fa"
  cleaned_fa="$output_dir/${base}.clean.fa"
  input_genome="$genome"

  if   [[ -s "$masked_fa" ]]; then
    input_genome="$masked_fa"; echo "   Using soft-masked genome: $input_genome"
  elif [[ -s "$cleaned_fa" ]]; then
    input_genome="$cleaned_fa"; echo "   Using cleaned genome: $input_genome"
  else
    echo "   Using renamed genome: $input_genome"
  fi

  flags=( -i "$input_genome" -o "$output_dir" --cpus "$THREADS" --species "$core" --EVM_HOME "$EVM_HOME" )

  # Add protein evidence: collect all proteomes as space-separated list
  protein_evidence_files=()

  # Always use best reference proteome if available (GCF or best quality)
  if [[ -n "$best_reference_proteome" ]]; then
    if [[ "$best_reference_proteome" == *"GCF_"* ]]; then
      echo "   Using best GCF (RefSeq) protein evidence: $(basename "$best_reference_proteome")"
    else
      echo "   Using best quality protein evidence: $(basename "$best_reference_proteome")"
    fi
    protein_evidence_files+=( "$best_reference_proteome" )
  fi

  # Also add strain-specific proteome if it exists and is different from reference
  if [[ -f "$faa" && "$faa" != "$best_reference_proteome" ]]; then
    echo "   Adding strain-specific protein evidence: $(basename "$faa")"
    protein_evidence_files+=( "$faa" )
  fi

  # Add all protein evidence files as single --protein_evidence argument
  if [[ ${#protein_evidence_files[@]} -gt 0 ]]; then
    flags+=( --protein_evidence "${protein_evidence_files[@]}" )
  fi

  # Add RNA-seq evidence if available for this strain
  # Try multiple possible RNA directory names (with and without _genomic suffix)
  strain_rna_dir=""
  for possible_dir in "$RNA_DIR/${core}_genomic" "$RNA_DIR/$core" "$RNA_DIR/$base"; do
    if [[ -d "$possible_dir" ]]; then
      strain_rna_dir="$possible_dir"
      break
    fi
  done

  if [[ -n "$strain_rna_dir" ]]; then
    # Find all fastq.gz files in the RNA directory
    rna_files=()
    for rna_file in "$strain_rna_dir"/*.fastq.gz; do
      [[ -f "$rna_file" ]] && rna_files+=("$rna_file")
    done

    if [[ ${#rna_files[@]} -gt 0 ]]; then
      echo "   Found ${#rna_files[@]} RNA-seq file(s) for $core in $(basename "$strain_rna_dir")"
      # Combine all RNA-seq files as a comma-separated list for --rna_bam
      rna_list=""
      for rna_file in "${rna_files[@]}"; do
        echo "      - $(basename "$rna_file")"
        if [[ -z "$rna_list" ]]; then
          rna_list="$rna_file"
        else
          rna_list="$rna_list,$rna_file"
        fi
      done
      flags+=( --rna_bam "$rna_list" )
    fi
  fi

  echo "   Filtering renamed GFF3: $renamed_gff"
  if [[ -s "$renamed_gff" ]]; then
    { grep -P "\t(?:gene|mRNA|CDS)\t" "$renamed_gff" || true; } \
      | sed -E \
          -e 's/gene_id=([^;]+)/ID=\1/' \
          -e 's/;transcript_id=([^;]+)/;Parent=\1;ID=\1/' \
          -e 's/;protein_id=([^;]+)/;Parent=\1/' \
      > "$filtered_gff"
    if [[ -s "$filtered_gff" ]]; then
      echo "     Created filtered GFF: $filtered_gff"
      flags+=( --other_gff "$filtered_gff":10 )
    else
      echo "     No gene/mRNA/CDS lines after filtering; skipping --other_gff"
    fi
  else
    echo "     Renamed GFF not found; skipping --other_gff"
  fi

  flags+=( --force )
  echo "   Running: funannotate predict ${flags[*]}"
  funannotate predict "${flags[@]}"

  run_log="$output_dir/logfiles/funannotate-predict.log"
  if [[ -f "$run_log" ]] && grep -q "bad contigs" "$run_log"; then
    echo "[$core] Funannotate skipped bad contigs:" >> "$OUTDIR/skipped_bad_contigs.log"
    grep "Found .* bad contigs" "$run_log" >> "$OUTDIR/skipped_bad_contigs.log"
    { grep -A 5 "bad contigs" "$run_log" || true; } \
    | { grep -v 'bad contigs' || true; } >> "$OUTDIR/skipped_bad_contigs.log"
    echo "----------------------------------------" >> "$OUTDIR/skipped_bad_contigs.log"
  fi

  echo "Completed $core"
done


echo
echo "All genomes cleaned, masked, and annotated with Funannotate!"
