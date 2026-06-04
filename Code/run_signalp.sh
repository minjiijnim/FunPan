#!/usr/bin/env bash
set -euo pipefail

# SignalP Signal Peptide Prediction (auto-parallel CPU)
# - Predicts signal peptides in protein sequences
# - Skips samples that already have output files
# - Processes all protein files from funannotate output
# - Auto-tunes parallel jobs + torch threads based on CPU cores
# - Use TEST_MODE=1 to run on a single file for testing
#
# Overrides (optional):
#   PIPELINE_CORES=96            # force total cores seen/used
#   PIPELINE_RESERVED_CORES=4    # keep headroom
#   SIGNALP_TORCH_THREADS=3      # torch threads per SignalP process
#   SIGNALP_JOBS=30              # number of parallel SignalP processes
#   SIGNALP_ORG=euk              # euk, gram+, gram-, arch
#   SIGNALP_FORMAT=txt           # v6: txt, png, eps, all, none
#   TEST_MODE=1                  # only first proteome

ORGANISM=${SIGNALP_ORG:-euk}
FORMAT=${SIGNALP_FORMAT:-txt}
TEST_MODE=${TEST_MODE:-0}

# ---- CPU / parallelism auto-tuning ------------------------------------------
TOTAL_CORES=$(nproc 2>/dev/null || echo 1)
TOTAL_CORES=${PIPELINE_CORES:-$TOTAL_CORES}

RESERVED=${PIPELINE_RESERVED_CORES:-4}
USE_CORES=$(( TOTAL_CORES - RESERVED ))
(( USE_CORES < 1 )) && USE_CORES=1

# Good default for PyTorch CPU inference: small threads per proc, many procs
THREADS_PER_JOB=${SIGNALP_TORCH_THREADS:-3}

PARALLEL_JOBS=${SIGNALP_JOBS:-$(( USE_CORES / THREADS_PER_JOB ))}
(( PARALLEL_JOBS < 1 )) && PARALLEL_JOBS=1

echo "Detected cores:          $TOTAL_CORES"
echo "Reserved cores:          $RESERVED"
echo "Cores available to use:  $USE_CORES"
echo "Torch threads per job:   $THREADS_PER_JOB"
echo "Parallel SignalP jobs:   $PARALLEL_JOBS"
echo

echo -n "> Enter species name (genus/species format): "
read -r species

# Layout
fun_dir=$(realpath "../Species/$species/funannotate_output")
out_root=$(realpath "../Species/$species/signalp_output")
protein_dir="$out_root/tmp_proteins"
mkdir -p "$out_root" "$protein_dir"

# ---- locate/install runner ---------------------------------------------------
if command -v signalp >/dev/null 2>&1; then
  SIGNALP_BIN="signalp"
elif command -v signalp6 >/dev/null 2>&1; then
  SIGNALP_BIN="signalp6"
  echo "> Note: Using SignalP6 (newer version)"
else
  echo "ERROR: signalp/signalp6 not found in PATH" >&2
  exit 1
fi
echo "> Using runner: $SIGNALP_BIN"

# ---- version detect ----------------------------------------------------------
VERSION_CHECK=$($SIGNALP_BIN --version 2>&1 || $SIGNALP_BIN -V 2>&1 || echo "unknown")
IS_V6=false
if [[ "$VERSION_CHECK" =~ 6\. ]] || [[ "$SIGNALP_BIN" == "signalp6" ]]; then
  IS_V6=true
  echo "> Detected SignalP 6.x"
else
  echo "> Using SignalP legacy (<6)"
fi
echo

# ---- Ensure SignalP6 model weights exist in site-packages (optional) ---------
if $IS_V6; then
  SIGNALP_PACKAGE=$(python -c "import signalp, os; print(os.path.dirname(signalp.__file__))" 2>/dev/null || echo "")
  if [[ -n "$SIGNALP_PACKAGE" ]]; then
    MODEL_WEIGHTS_DIR="$SIGNALP_PACKAGE/model_weights"
    EXPECTED_MODEL="$MODEL_WEIGHTS_DIR/distilled_model_signalp6.pt"
    # Adjust this source path if your model file is elsewhere:
    SOURCE_MODEL="/datadrive/Data/SignalP/models/distilled_model_signalp6.pt"

    if [[ ! -f "$EXPECTED_MODEL" ]]; then
      echo "> Model file missing in site-packages: $EXPECTED_MODEL"
      if [[ -f "$SOURCE_MODEL" ]]; then
        echo "> Copying model file into site-packages (one-time setup)..."
        mkdir -p "$MODEL_WEIGHTS_DIR"
        cp "$SOURCE_MODEL" "$EXPECTED_MODEL" || {
          echo "WARNING: Could not copy model file. SignalP may fail." >&2
        }
      else
        echo "WARNING: Source model file not found at: $SOURCE_MODEL" >&2
        echo "         SignalP may fail without model weights." >&2
      fi
    else
      echo "> Model file present in site-packages"
    fi
  else
    echo "WARNING: Could not locate Python 'signalp' package path; skipping model check." >&2
  fi
  echo
fi

# ---- Gather protein files from funannotate predict_results -------------------
echo "> Searching for *.proteins.fa under: $fun_dir"
find "$protein_dir" -type l -delete 2>/dev/null || true

mapfile -t prot_sources < <(find "$fun_dir" -type f -path "*/predict_results/*.proteins.fa" | sort)
if [[ ${#prot_sources[@]} -eq 0 ]]; then
  echo "ERROR: No *.proteins.fa found under $fun_dir" >&2
  exit 1
fi

for p in "${prot_sources[@]}"; do
  abs=$(realpath "$p")
  ln -sf "$abs" "$protein_dir/$(basename "$abs")"
done

mapfile -t protein_files < <(find "$protein_dir" -type l -name "*.proteins.fa" | sort)
echo "> Linked ${#protein_files[@]} protein file(s) into $protein_dir"

if [[ $TEST_MODE -eq 1 ]]; then
  echo "> TEST MODE: Processing only the first protein file"
  protein_files=("${protein_files[0]}")
fi
echo

# Clamp PARALLEL_JOBS to number of files
if (( PARALLEL_JOBS > ${#protein_files[@]} )); then
  PARALLEL_JOBS=${#protein_files[@]}
  (( PARALLEL_JOBS < 1 )) && PARALLEL_JOBS=1
fi
echo "> Final parallel job count: $PARALLEL_JOBS"
echo

# ---- job limiter -------------------------------------------------------------
wait_for_slot () {
  while (( $(jobs -rp | wc -l) >= PARALLEL_JOBS )); do
    sleep 2
  done
}

# ---- run SignalP per file in parallel ---------------------------------------
echo "> Running SignalP on ${species} protein files"
for protein_file in "${protein_files[@]}"; do
  base="$(basename "$protein_file" .proteins.fa)"
  out_dir="$out_root/$base"
  mkdir -p "$out_dir"

  if $IS_V6; then
    sentinel="$out_dir/prediction_results.txt"
    if [[ -f "$sentinel" ]]; then
      echo "> Skipping $base (results present)"
      continue
    fi

    # Map organism types: euk -> eukarya for SignalP6
    SP6_ORG="$ORGANISM"
    [[ "$ORGANISM" == "euk" ]] && SP6_ORG="eukarya"

    # Map legacy-ish format strings if user passes them
    SP6_FORMAT="$FORMAT"
    [[ "$FORMAT" == "short" ]] && SP6_FORMAT="txt"
    [[ "$FORMAT" == "long"  ]] && SP6_FORMAT="all"

    echo "> Queueing $base (SignalP6)"
    wait_for_slot

    (
      export OMP_NUM_THREADS="$THREADS_PER_JOB"
      export MKL_NUM_THREADS="$THREADS_PER_JOB"
      export OPENBLAS_NUM_THREADS="$THREADS_PER_JOB"

      $SIGNALP_BIN \
        --fastafile "$protein_file" \
        --output_dir "$out_dir" \
        --organism "$SP6_ORG" \
        --format "$SP6_FORMAT" \
        --torch_num_threads "$THREADS_PER_JOB"
    ) >"$out_dir/signalp.log" 2>&1 &

  else
    sentinel="$out_dir/${base}.signalp.out"
    if [[ -f "$sentinel" ]]; then
      echo "> Skipping $base (results present)"
      continue
    fi

    echo "> Queueing $base (legacy)"
    wait_for_slot

    (
      $SIGNALP_BIN \
        -t "$ORGANISM" \
        -f "$FORMAT" \
        "$protein_file" \
        > "$sentinel"
    ) >"$out_dir/signalp.log" 2>&1 &
  fi
done

wait

# ---- quick summary -----------------------------------------------------------
if $IS_V6; then
  count_done=$(find "$out_root" -mindepth 1 -maxdepth 1 -type d -exec test -f "{}/prediction_results.txt" ';' -print | wc -l | awk '{print $1}')
else
  count_done=$(find "$out_root" -mindepth 1 -maxdepth 1 -type d -exec bash -lc 'ls "{}"/*.signalp.out >/dev/null 2>&1' ';' -print | wc -l | awk '{print $1}')
fi

echo
echo "==================================="
echo "SignalP Pipeline Complete"
echo "==================================="
echo "Completed: $count_done run(s)"
echo "Output root: $out_root"
echo "SignalP version: $(if $IS_V6; then echo "6.x"; else echo "legacy (<6)"; fi)"
echo
echo "Example output files:"
for sample_dir in "$out_root"/*/; do
  [[ -d "$sample_dir" ]] || continue
  sample=$(basename "$sample_dir")
  if $IS_V6; then
    [[ -f "$sample_dir/prediction_results.txt" ]] && echo "  - $sample_dir/prediction_results.txt" && break
  else
    [[ -f "$sample_dir/${sample}.signalp.out" ]] && echo "  - $sample_dir/${sample}.signalp.out" && break
  fi
done
