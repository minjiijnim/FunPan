#!/usr/bin/env bash
set -euo pipefail

# dbCAN CAZyme Annotation (official CLI pattern, no extra state files)
# - DB bootstrap via: run_dbcan database --db_dir <dir>
# - Per-file run via:  run_dbcan CAZyme_annotation ... --mode protein
# - Skips samples that already have "overview.tsv"
# - Uses native --resume if available (detected), else just skips completed

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

METHODS=${DBCAN_METHODS:-}          # optional: e.g. "hmm,diamond,dbCANsub" or "hmm" (see below)

echo "> Enter species name (genus/species format):"
read -r species

# Layout
fun_dir="../Species/$species/funannotate_output"
protein_dir="${fun_dir}/tmp_proteins"
out_root="../Species/$species/dbcan_output"
db_dir="../Data/dbCAN"

mkdir -p "$protein_dir" "$out_root" "$db_dir"

# ---- locate runner -----------------------------------------------------------
if command -v run_dbcan >/dev/null 2>&1; then
  RUN_DB="$(command -v run_dbcan)"
elif command -v run_dbcan.py >/dev/null 2>&1; then
  RUN_DB="$(command -v run_dbcan.py)"
else
  echo "ERROR: run_dbcan not found on PATH. Activate the conda env with dbCAN." >&2
  exit 1
fi
echo "> Using runner: $RUN_DB"

# ---- ensure DB (official command; only if missing) ---------------------------
need_db=false
[[ -f "$db_dir/CAZy.dmnd" ]]     || need_db=true
[[ -f "$db_dir/dbCAN-sub.hmm" ]] || need_db=true

if $need_db; then
  echo "> Downloading dbCAN database (official): $db_dir"
  "$RUN_DB" database --db_dir "$db_dir"
fi

# minimal sanity
[[ -f "$db_dir/CAZy.dmnd" ]]     || { echo "ERROR: $db_dir/CAZy.dmnd missing after download." >&2; exit 1; }
[[ -f "$db_dir/dbCAN-sub.hmm" ]] || { echo "ERROR: $db_dir/dbCAN-sub.hmm missing after download." >&2; exit 1; }

if ! compgen -G "$db_dir/dbCAN-sub.hmm.h3?" >/dev/null; then
  echo "> Pressing HMM database"
  hmmpress "$db_dir/dbCAN-sub.hmm"
fi

# ---- gather protein files (same pattern as your EggNOG flow) -----------------
echo "> Searching for *.proteins.fa under: $fun_dir"
find "$protein_dir" -type l -delete || true
mapfile -t prot_sources < <(find "$fun_dir" -type f -path "*/predict_results/*.proteins.fa" | sort)
[[ ${#prot_sources[@]} -gt 0 ]] || { echo "ERROR: No *.proteins.fa under $fun_dir" >&2; exit 1; }

for p in "${prot_sources[@]}"; do
  ln -s "$(realpath "$p")" "$protein_dir/$(basename "$p")"
done

mapfile -t protein_files < <(find "$protein_dir" -type l -name "*.proteins.fa" | sort)
echo "> Linked ${#protein_files[@]} protein file(s) into $protein_dir"

# ---- feature detection (threads, resume, methods) ----------------------------
help_txt="$("$RUN_DB" CAZyme_annotation --help 2>/dev/null || true)"

THREAD_FLAG=()
grep -q -- "--threads" <<<"$help_txt" && THREAD_FLAG=(--threads "$THREADS")

RESUME_FLAG=()
grep -q -- "--resume"  <<<"$help_txt" && RESUME_FLAG=(--resume)

# Methods: official docs allow repeating --methods; we’ll map comma or space list into repeats
declare -a METHOD_FLAGS=()
if [[ -n "$METHODS" ]]; then
  IFS=',' read -r -a _methods <<<"${METHODS// /,}"
  for m in "${_methods[@]}"; do
    [[ -n "$m" ]] && METHOD_FLAGS+=(--methods "$m")
  done
fi

# ---- run dbCAN per file (skip if already has results) ------------------------
echo "> Running dbCAN CAZyme_annotation on ${species} protein files"
for protein_file in "${protein_files[@]}"; do
  base="$(basename "$protein_file" .proteins.fa)"
  out_dir="$out_root/$base"
  sentinel="$out_dir/overview.tsv"
  mkdir -p "$out_dir"

  if [[ -f "$sentinel" ]]; then
    echo "> Skipping $base (results present: overview.tsv)"
    continue
  fi

  echo "> Processing $base"
  # Official syntax:
  # run_dbcan CAZyme_annotation --input_raw_data <INPUT> --out_dir <OUT> --db_dir <DB> --mode protein
  # + optional: --threads, --methods (hmm/diamond/dbCANsub), --resume (if supported)
  "$RUN_DB" CAZyme_annotation \
    --input_raw_data "$protein_file" \
    --output_dir "$out_dir" \
    --db_dir "$db_dir" \
    --mode protein \
    "${THREAD_FLAG[@]}" \
    "${METHOD_FLAGS[@]}" \
    "${RESUME_FLAG[@]}"
done

# ---- quick summary -----------------------------------------------------------
count_done=$(find "$out_root" -mindepth 1 -maxdepth 1 -type d -exec test -f "{}/overview.tsv" ';' -print | wc -l | awk '{print $1}')
echo "> Completed: $count_done run(s)"
echo "> Output root: $out_root"
