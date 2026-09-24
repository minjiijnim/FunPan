#!/usr/bin/env bash
set -euo pipefail

# =============================================================================
# deduplicate_and_filter_genome_QC.sh
#
# Extended version of filter_genome_QC.sh that adds:
#   1. ANI species verification (fastANI against GCF reference)
#   2. Mash-based deduplication (remove near-identical genomes)
#   3. Original QC filters (contigs, N50, BUSCO)
#
# Usage: bash deduplicate_and_filter_genome_QC.sh Aspergillus/flavus [ANI_THRESHOLD] [MASH_THRESHOLD]
#   ANI_THRESHOLD:  minimum ANI% to species reference (default: 95)
#   MASH_THRESHOLD: maximum Mash distance for duplicate detection (default: 0.0001 ≈ 99.99% identity)
# =============================================================================

# Use only comma as IFS for CSV parsing
IFS=,

# ── Arguments ────────────────────────────────────────────────────────────────
if [[ $# -ge 1 ]]; then
    species_name="$1"
else
    echo -n "> Enter species name: "
    read species_name
fi

ANI_THRESHOLD="${2:-95}"
MASH_THRESHOLD="${3:-0.0001}"

echo "==========================================================================="
echo "  Deduplicate & Filter Genome QC"
echo "  Species: $species_name"
echo "  ANI threshold: $ANI_THRESHOLD%"
echo "  Mash dedup threshold: $MASH_THRESHOLD"
echo "==========================================================================="

metadata_file="../Species/${species_name}/metadata.csv"
base_dir="../Species/${species_name}"
genome_dir="${base_dir}/genome"
protein_dir="${base_dir}/protein"
busco_dir="${base_dir}/busco_output"
out_genome_dir="${base_dir}/filtered_genome"
out_protein_dir="${base_dir}/filtered_protein"

# Working directory for intermediate files
work_dir="${base_dir}/qc_work"
mkdir -p "$work_dir" "$out_genome_dir" "$out_protein_dir"

# ── Dependency check (auto-install via conda if missing) ────────────────────
ensure_tool() {
    local tool="$1"
    local pkg="$2"
    if command -v "$tool" &>/dev/null; then
        return 0
    fi
    echo "> $tool not found — installing $pkg via conda..."
    command -v conda >/dev/null 2>&1 || { echo "ERROR: conda not found in PATH" >&2; exit 1; }
    if command -v mamba >/dev/null 2>&1; then
        mamba install -y -c conda-forge -c bioconda "$pkg"
    else
        conda install -y -c conda-forge -c bioconda "$pkg"
    fi
    if ! command -v "$tool" &>/dev/null; then
        echo "ERROR: $tool still not available after install" >&2
        exit 1
    fi
}

ensure_tool fastANI fastani
ensure_tool mash    mash

# =============================================================================
# STEP 0: Parse metadata columns (same as original)
# =============================================================================
read -r -a headers < "$metadata_file"
for i in "${!headers[@]}"; do
    h_norm="$(echo "${headers[i]}" | tr '[:upper:]' '[:lower:]' | tr -d ' ')"
    case "$h_norm" in
        assemblyaccession)                       acc_idx=$((i+1)) ;;
        assemblystatsnumberofcontigs)            contigs_idx=$((i+1)) ;;
        buscocomplete)                           busco_idx=$((i+1)) ;;
        buscototal)                              busco_total_idx=$((i+1)) ;;
        assemblystatscontign50|contign50|n50)    n50_idx=$((i+1)) ;;
    esac
done

if [[ -z "${contigs_idx:-}" || -z "${busco_idx:-}" || -z "${busco_total_idx:-}" || -z "${acc_idx:-}" || -z "${n50_idx:-}" ]]; then
    echo "Could not find required columns in $metadata_file" >&2
    exit 1
fi

# =============================================================================
# STEP 1: FIND REFERENCE GENOME (GCF preferred)
# =============================================================================
echo ""
echo "─── Step 1: Finding reference genome for ANI check ───"

# Look for GCF genome — prioritize filtered_genome (same reference used by funannotate)
ref_genome=""

# Try filtered_genome first (consistent with funannotate reference selection)
gcf_file=$(find "$out_genome_dir" -name "GCF_*.fna" -type f 2>/dev/null | head -1)
if [[ -n "$gcf_file" ]]; then
    ref_genome="$gcf_file"
fi

# Try genome dir
if [[ -z "$ref_genome" ]]; then
    gcf_file=$(find "$genome_dir" -name "GCF_*.fna" -type f 2>/dev/null | head -1)
    if [[ -n "$gcf_file" ]]; then
        ref_genome="$gcf_file"
    fi
fi

# Try busco_output if not found
if [[ -z "$ref_genome" ]]; then
    gcf_file=$(find "$busco_dir" -maxdepth 1 -name "GCF_*.fna" -type f 2>/dev/null | head -1)
    if [[ -n "$gcf_file" ]]; then
        ref_genome="$gcf_file"
    fi
fi

if [[ -z "$ref_genome" ]]; then
    echo "WARNING: No GCF reference genome found. ANI check will be skipped." >&2
    ANI_CHECK=false
else
    ANI_CHECK=true
    ref_acc=$(basename "$ref_genome" | grep -oP 'GCF_\d+\.\d+')
    echo "  Reference genome: $(basename "$ref_genome")"
    echo "  Reference accession: $ref_acc"
fi

# =============================================================================
# STEP 2: ANI CHECK — Compare all genomes against reference
# =============================================================================
echo ""
echo "─── Step 2: ANI species verification (fastANI) ───"

# Build list of query genomes
find "$genome_dir" -name "*.fna" -type f > "$work_dir/all_genomes.txt" 2>/dev/null
# Also check busco_output for genomes
find "$busco_dir" -maxdepth 1 -name "*.fna" -type f >> "$work_dir/all_genomes.txt" 2>/dev/null
# Deduplicate list
sort -u "$work_dir/all_genomes.txt" > "$work_dir/all_genomes_uniq.txt"
mv "$work_dir/all_genomes_uniq.txt" "$work_dir/all_genomes.txt"

n_genomes=$(wc -l < "$work_dir/all_genomes.txt")
echo "  Total genomes to check: $n_genomes"

ani_passed_file="$work_dir/ani_passed.txt"
ani_failed_file="$work_dir/ani_failed.txt"
> "$ani_passed_file"
> "$ani_failed_file"

if [[ "$ANI_CHECK" == true && $n_genomes -gt 0 ]]; then
    ani_results="$work_dir/fastani_results.txt"

    echo "  Running fastANI (all vs reference)..."
    fastANI \
        --ql "$work_dir/all_genomes.txt" \
        --ref "$ref_genome" \
        -o "$ani_results" \
        -t 8 \
        2>/dev/null

    # Parse results
    n_pass=0
    n_fail=0
    while IFS=$'\t' read -r query ref ani frags total; do
        acc=$(basename "$query" | grep -oP 'GC[AF]_\d+\.\d+' || basename "$query" | grep -oP 'GC[AF]_\d+')
        if (( $(echo "$ani >= $ANI_THRESHOLD" | bc -l) )); then
            echo "$acc" >> "$ani_passed_file"
            ((n_pass++))
        else
            echo "$acc	$ani" >> "$ani_failed_file"
            ((n_fail++))
            echo "  FAILED ANI: $acc (ANI=$ani% < $ANI_THRESHOLD%)"
        fi
    done < "$ani_results"

    # Check for genomes not in fastANI output (too divergent to compute)
    while read -r genome_path; do
        acc=$(basename "$genome_path" | grep -oP 'GC[AF]_\d+\.\d+' || basename "$genome_path" | grep -oP 'GC[AF]_\d+')
        if ! grep -q "$acc" "$ani_passed_file" && ! grep -q "$acc" "$ani_failed_file"; then
            echo "$acc	NO_ANI" >> "$ani_failed_file"
            ((n_fail++))
            echo "  FAILED ANI: $acc (no ANI computed — too divergent)"
        fi
    done < "$work_dir/all_genomes.txt"

    echo "  ANI check: $n_pass passed, $n_fail failed (threshold: $ANI_THRESHOLD%)"
else
    echo "  ANI check skipped (no reference genome)"
    # Pass all genomes
    while read -r genome_path; do
        acc=$(basename "$genome_path" | grep -oP 'GC[AF]_\d+\.\d+' || basename "$genome_path" | grep -oP 'GC[AF]_\d+')
        echo "$acc" >> "$ani_passed_file"
    done < "$work_dir/all_genomes.txt"
fi

echo ""
if [[ -s "$ani_failed_file" ]]; then
    echo "  ANI-failed genomes:"
    cat "$ani_failed_file" | while IFS=$'\t' read -r acc ani; do
        echo "    $acc: ANI=$ani"
    done
fi

# =============================================================================
# STEP 3: MASH DEDUPLICATION
# =============================================================================
echo ""
echo "─── Step 3: Mash deduplication ───"

# Build list of ANI-passed genome files
> "$work_dir/ani_passed_genomes.txt"
while read -r acc; do
    # Find the genome file
    gfile=$(find "$genome_dir" "$busco_dir" -maxdepth 1 -name "${acc}*_genomic.fna" -type f 2>/dev/null | head -1)
    if [[ -n "$gfile" ]]; then
        echo "$gfile" >> "$work_dir/ani_passed_genomes.txt"
    fi
done < "$ani_passed_file"

n_passed=$(wc -l < "$work_dir/ani_passed_genomes.txt")
echo "  Genomes passing ANI: $n_passed"

if [[ $n_passed -gt 1 ]]; then
    echo "  Sketching genomes with Mash..."
    mash sketch -l "$work_dir/ani_passed_genomes.txt" -o "$work_dir/all_sketches" -s 10000 2>/dev/null

    echo "  Computing all-vs-all Mash distances..."
    mash dist "$work_dir/all_sketches.msh" "$work_dir/all_sketches.msh" -t > "$work_dir/mash_distances.tsv" 2>/dev/null

    # Cluster near-identical genomes and keep best representative
    echo "  Clustering duplicates (threshold: $MASH_THRESHOLD)..."

    python3 - "$work_dir/mash_distances.tsv" "$MASH_THRESHOLD" "$metadata_file" "$acc_idx" "$busco_idx" "$busco_total_idx" "$contigs_idx" "$n50_idx" << 'PYSCRIPT'
import sys, csv, re
from collections import defaultdict

dist_file = sys.argv[1]
threshold = float(sys.argv[2])
metadata_file = sys.argv[3]
acc_idx = int(sys.argv[4]) - 1
busco_idx = int(sys.argv[5]) - 1
busco_total_idx = int(sys.argv[6]) - 1
contigs_idx = int(sys.argv[7]) - 1
n50_idx = int(sys.argv[8]) - 1

def get_acc(path):
    m = re.search(r'(GC[AF]_\d+\.\d+)', path)
    if m: return m.group(1)
    m = re.search(r'(GC[AF]_\d+)', path)
    return m.group(1) if m else path

# Load quality scores from metadata
quality = {}
with open(metadata_file) as f:
    reader = csv.reader(f)
    next(reader)  # skip header
    for row in reader:
        if len(row) <= max(acc_idx, busco_idx, busco_total_idx, contigs_idx, n50_idx):
            continue
        acc = row[acc_idx].strip()
        try:
            busco = float(row[busco_idx])
            btotal = float(row[busco_total_idx]) if row[busco_total_idx] else 1
            contigs = float(row[contigs_idx]) if row[contigs_idx] else 0
            n50 = float(row[n50_idx]) if row[n50_idx] else 0
            busco_pct = (busco / btotal * 100) if btotal > 0 else busco
            quality[acc] = busco_pct - (contigs * 0.01) + (n50 * 0.000001)
        except (ValueError, IndexError):
            pass

# Parse distance matrix
with open(dist_file) as f:
    header = f.readline().strip().split('\t')
    genomes = [get_acc(h) for h in header[1:]]  # skip first empty col or label

    # If header starts with #query or similar, adjust
    if header[0].startswith('#') or header[0].startswith('/'):
        genomes = [get_acc(h) for h in header]

    dist_matrix = {}
    for line in f:
        parts = line.strip().split('\t')
        row_acc = get_acc(parts[0])
        for j, val in enumerate(parts[1:]):
            col_acc = genomes[j] if j < len(genomes) else f"col{j}"
            try:
                d = float(val)
                if row_acc != col_acc and d <= threshold:
                    dist_matrix.setdefault(row_acc, set()).add(col_acc)
            except ValueError:
                pass

# Greedy clustering: pick best quality as representative
removed = set()
clusters = []
for acc in sorted(dist_matrix.keys()):
    if acc in removed:
        continue
    neighbors = dist_matrix.get(acc, set()) - removed
    cluster = {acc} | neighbors
    # Pick best quality as representative
    best = max(cluster, key=lambda a: quality.get(a, 0))
    dups = cluster - {best}
    removed |= dups
    if dups:
        clusters.append((best, dups))

# Report
if clusters:
    print(f"  Found {len(clusters)} duplicate clusters, removing {len(removed)} redundant genomes:")
    for rep, dups in clusters:
        q_rep = quality.get(rep, 0)
        dup_list = ', '.join(sorted(dups))
        print(f"    Keep {rep} (q={q_rep:.1f}), remove: {dup_list}")
else:
    print("  No duplicates found at threshold")

# Write keep list
all_accs = set()
with open(metadata_file) as f:
    reader = csv.reader(f)
    next(reader)
    for row in reader:
        if len(row) > acc_idx:
            all_accs.add(row[acc_idx].strip())

kept = all_accs - removed
with open(sys.argv[1].replace('mash_distances.tsv', 'dedup_kept.txt'), 'w') as f:
    for acc in sorted(kept):
        f.write(acc + '\n')

with open(sys.argv[1].replace('mash_distances.tsv', 'dedup_removed.txt'), 'w') as f:
    for acc in sorted(removed):
        f.write(acc + '\n')

print(f"  Total: {len(all_accs)} → {len(kept)} after deduplication ({len(removed)} removed)")
PYSCRIPT

else
    echo "  Only $n_passed genomes — skipping deduplication"
    cp "$ani_passed_file" "$work_dir/dedup_kept.txt"
    > "$work_dir/dedup_removed.txt"
fi

# =============================================================================
# STEP 4: ORIGINAL QC FILTERS (contigs, N50, BUSCO) — same as filter_genome_QC.sh
# =============================================================================
echo ""
echo "─── Step 4: QC filters (contigs, N50, BUSCO) ───"

# Calculate thresholds
median_contigs=$(
  awk -F',' -v idx="$contigs_idx" '
    NR>1 && $idx!="" { v[++n]=$idx+0 }
    END{
      if (n==0){print 0; exit}
      asort(v)
      if (n%2)  print v[(n+1)/2]
      else      print int((v[n/2]+v[n/2+1])/2)
    }' "$metadata_file"
)

p75_n50=$(
  awk -F',' -v idx="$n50_idx" '
    NR>1 && $idx!="" { v[++n]=$idx+0 }
    END{
      if (n==0){print 0; exit}
      asort(v)
      r = int(0.75*n); if (r<1) r=1; if (r>n) r=n
      print v[r]
    }' "$metadata_file"
)

busco_cutoff=$(
  awk -F',' -v idx="$busco_idx" '
    NR>1 && $idx!="" && $idx+0>max { max=$idx+0 }
    END{ if (max<=1) print 0.95; else print 95 }' "$metadata_file"
)

echo "  Thresholds: contigs≤$median_contigs, N50≥$p75_n50, BUSCO>$busco_cutoff%"

# Load ANI-passed and dedup-kept sets
declare -A ani_ok dedup_ok
while read -r acc; do ani_ok["$acc"]=1; done < "$ani_passed_file"
if [[ -f "$work_dir/dedup_kept.txt" ]]; then
    while read -r acc; do dedup_ok["$acc"]=1; done < "$work_dir/dedup_kept.txt"
fi

# Clear output dirs
rm -f "$out_genome_dir"/*.fna "$out_genome_dir"/*.map "$out_genome_dir"/*_renamed.fna
rm -f "$out_protein_dir"/*.faa

copied_genomes=0
copied_proteins=0
skipped_ani=0
skipped_dedup=0
skipped_qc=0
declare -a missing_genomes=() missing_proteins=()

set +e
while read -a fields; do
    contigs="${fields[$((contigs_idx-1))]}"
    n50="${fields[$((n50_idx-1))]}"
    busco="${fields[$((busco_idx-1))]}"
    busco_total="${fields[$((busco_total_idx-1))]}"
    acc="${fields[$((acc_idx-1))]}"
    acc="$(echo "$acc" | xargs | tr -d '\r\t')"

    [[ -z "$acc" ]] && continue

    # Check ANI
    acc_base=$(echo "$acc" | grep -oP 'GC[AF]_\d+')
    if [[ "$ANI_CHECK" == true ]] && [[ -z "${ani_ok[$acc]:-}" ]] && [[ -z "${ani_ok[$acc_base]:-}" ]]; then
        # Check with version stripped too
        found=false
        for key in "${!ani_ok[@]}"; do
            if [[ "$key" == "$acc_base"* ]]; then
                found=true
                break
            fi
        done
        if [[ "$found" == false ]]; then
            ((skipped_ani++))
            continue
        fi
    fi

    # Check dedup
    if [[ -f "$work_dir/dedup_removed.txt" ]] && grep -q "$acc" "$work_dir/dedup_removed.txt"; then
        ((skipped_dedup++))
        continue
    fi

    # QC filter (same as original)
    if awk -v c="$contigs" -v n="$n50" -v b="$busco" -v btot="$busco_total" \
          -v med="$median_contigs" -v p75="$p75_n50" -v bcut="$busco_cutoff" '
         BEGIN{
           if (btot > 0) { busco_pct = (b / btot) * 100 }
           else if (b <= 1) { busco_pct = b * 100 }
           else { busco_pct = b }
           pass = (c+0) <= (med+0) && (n+0) >= (p75+0) && busco_pct > (bcut+0)
           exit(pass?0:1)
         }'; then

        # Copy genome
        if compgen -G "${genome_dir}/${acc}_*_genomic.fna" > /dev/null; then
            for src in ${genome_dir}/${acc}_*_genomic.fna; do
                cp -p "$src" "$out_genome_dir/"
                ((copied_genomes++))
            done
        else
            missing_genomes+=("$acc")
        fi

        # Copy protein
        if compgen -G "${protein_dir}/${acc}*.faa" > /dev/null; then
            for src in ${protein_dir}/${acc}*.faa; do
                cp -p "$src" "$out_protein_dir/"
                ((copied_proteins++))
            done
        else
            missing_proteins+=("$acc")
        fi
    else
        ((skipped_qc++))
    fi
done < <(tail -n +2 "$metadata_file")
set -e

# =============================================================================
# STEP 5: Ensure GCF reference is included (same as original)
# =============================================================================
gcf_count=$(find "$out_genome_dir" -type f -name "GCF_*" | wc -l)
if [[ $gcf_count -eq 0 && "$ANI_CHECK" == true ]]; then
    echo "  Adding GCF reference genome to filtered set..."
    if compgen -G "${genome_dir}/${ref_acc}_*_genomic.fna" > /dev/null; then
        for src in ${genome_dir}/${ref_acc}_*_genomic.fna; do
            cp -p "$src" "$out_genome_dir/"
            ((copied_genomes++))
        done
    fi
    if compgen -G "${protein_dir}/${ref_acc}*.faa" > /dev/null; then
        for src in ${protein_dir}/${ref_acc}*.faa; do
            cp -p "$src" "$out_protein_dir/"
            ((copied_proteins++))
        done
    fi
fi

# =============================================================================
# STEP 6: Summary
# =============================================================================
echo ""
echo "==========================================================================="
echo "  SUMMARY"
echo "==========================================================================="
echo "  Genomes in metadata:     $n_genomes"
echo "  Failed ANI check:        $skipped_ani"
echo "  Removed as duplicates:   $skipped_dedup"
echo "  Failed QC filters:       $skipped_qc"
echo "  Copied to filtered dir:  $copied_genomes genomes, $copied_proteins proteins"
echo ""

if [[ -s "$ani_failed_file" ]]; then
    echo "  ANI-failed samples:"
    cat "$ani_failed_file" | while IFS=$'\t' read -r acc ani; do
        echo "    $acc (ANI=$ani)"
    done
fi

if [[ -s "$work_dir/dedup_removed.txt" ]]; then
    echo ""
    echo "  Deduplicated (removed) samples:"
    cat "$work_dir/dedup_removed.txt" | while read -r acc; do
        echo "    $acc"
    done
fi

echo ""
echo "  Output: $out_genome_dir"
echo "          $out_protein_dir"
echo "  Work:   $work_dir"
