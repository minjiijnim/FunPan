#!/usr/bin/env bash
set -euo pipefail

# InterProScan Protein Annotation
# - Scans protein sequences against InterPro signatures
# - Skips samples that already have output TSV files
# - Processes all protein files from funannotate output
# - Use TEST_MODE=1 to run on a single file for testing

# Default: Skip Hamap (bacterial), ProSitePatterns (slow), PRINTS (outdated)
# This runs 15 tools instead of 18, saving ~30-40% time
# Override with: INTERPRO_APPS="Pfam,SMART,..." to specify custom tools
APPLICATIONS=${INTERPRO_APPS:-}      # optional: e.g. "Pfam,SMART,CDD" (comma-separated)
TEST_MODE=${TEST_MODE:-0}            # set TEST_MODE=1 to process only first file

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

echo -n "> Enter species name (genus/species format): "
read -r species

# Layout
fun_dir=$(realpath "../Species/$species/funannotate_output")
out_root=$(realpath "../Species/$species/interproscan_output")
protein_dir="$out_root/tmp_proteins"

# Create output directories
mkdir -p "$out_root" "$protein_dir"

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


# ---- locate/install runner ---------------------------------------------------
INTERPROSCAN_DIR=$(realpath "../Data/InterProScan")
INTERPROSCAN_VERSION="5.76-107.0"
INTERPROSCAN_BIN_PATH="${INTERPROSCAN_DIR}/interproscan-${INTERPROSCAN_VERSION}/interproscan.sh"

if command -v interproscan.sh >/dev/null 2>&1; then
  INTERPRO_BIN="interproscan.sh"
  echo "> Found interproscan.sh in PATH"
elif command -v interproscan >/dev/null 2>&1; then
  INTERPRO_BIN="interproscan"
  echo "> Found interproscan in PATH"
elif [[ -x "$INTERPROSCAN_BIN_PATH" ]]; then
  INTERPRO_BIN="$INTERPROSCAN_BIN_PATH"
  echo "> Found local InterProScan installation"
else
  echo "> InterProScan not found. Starting automatic installation..."
  echo "> This is a one-time setup that will download ~2.5 GB"
  echo

  # Save current directory
  ORIGINAL_DIR=$(pwd)

  # Create installation directory
  mkdir -p "$INTERPROSCAN_DIR"
  cd "$INTERPROSCAN_DIR"

  # Download InterProScan
  TARBALL="interproscan-${INTERPROSCAN_VERSION}-64-bit.tar.gz"
  DOWNLOAD_URL="https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/${INTERPROSCAN_VERSION}/${TARBALL}"
  MD5_URL="${DOWNLOAD_URL}.md5"

  echo "> Downloading InterProScan ${INTERPROSCAN_VERSION}..."
  echo "> URL: $DOWNLOAD_URL"

  if [[ ! -f "$TARBALL" ]]; then
    wget -q --show-progress "$DOWNLOAD_URL" || {
      echo "ERROR: Failed to download InterProScan" >&2
      exit 1
    }

    wget -q "$MD5_URL" || {
      echo "WARNING: Could not download MD5 checksum" >&2
    }

    # Verify checksum if MD5 file was downloaded
    if [[ -f "${TARBALL}.md5" ]]; then
      echo "> Verifying checksum..."
      if md5sum -c "${TARBALL}.md5"; then
        echo "> Checksum verified successfully"
      else
        echo "ERROR: Checksum verification failed. Please try downloading again." >&2
        rm -f "$TARBALL" "${TARBALL}.md5"
        exit 1
      fi
    fi
  else
    echo "> Using existing downloaded file: $TARBALL"
  fi

  # Extract if not already extracted
  if [[ ! -d "interproscan-${INTERPROSCAN_VERSION}" ]]; then
    echo "> Extracting InterProScan..."
    tar -pxzf "$TARBALL" || {
      echo "ERROR: Failed to extract InterProScan" >&2
      exit 1
    }
    echo "> Extraction complete"
  else
    echo "> InterProScan already extracted"
  fi

  # Index HMM models
  cd "interproscan-${INTERPROSCAN_VERSION}"

  if [[ ! -f ".setup_complete" ]]; then
    echo "> Indexing HMM models (this may take several minutes)..."
    python3 setup.py -f interproscan.properties || {
      echo "ERROR: Failed to index HMM models" >&2
      exit 1
    }
    touch .setup_complete
    echo "> HMM model indexing complete"
  else
    echo "> HMM models already indexed"
  fi

  # Return to original directory
  cd "$ORIGINAL_DIR"

  # Set the binary path using absolute path
  INTERPRO_BIN="$INTERPROSCAN_BIN_PATH"

  echo
  echo "==================================="
  echo "InterProScan Installation Complete"
  echo "==================================="
  echo "Installation location: ${INTERPROSCAN_DIR}/interproscan-${INTERPROSCAN_VERSION}"
  echo "Binary: $INTERPRO_BIN"
  echo
fi

# Verify the binary is executable
if [[ ! -x "$INTERPRO_BIN" ]]; then
  echo "ERROR: InterProScan binary is not executable: $INTERPRO_BIN" >&2
  echo "Attempting to fix permissions..."
  chmod +x "$INTERPRO_BIN" 2>/dev/null && echo "> Permissions fixed" || {
    echo "ERROR: Could not fix permissions. Please run: chmod +x $INTERPRO_BIN" >&2
    exit 1
  }
fi

# Check for Java (InterProScan requires Java 11+)
if ! command -v java >/dev/null 2>&1; then
  echo "ERROR: Java not found. InterProScan requires Java 11 or higher." >&2
  echo "" >&2
  echo "Install Java 11 using conda:" >&2
  echo "  conda install -c conda-forge openjdk=11 -y" >&2
  echo "" >&2
  echo "Or install system-wide:" >&2
  echo "  sudo apt-get install openjdk-11-jdk  # Ubuntu/Debian" >&2
  echo "  sudo yum install java-11-openjdk     # CentOS/RHEL" >&2
  exit 1
fi

# Check Java version
JAVA_VERSION=$(java -version 2>&1 | head -n1 | awk -F'"' '{print $2}' | awk -F'.' '{print $1}')
if [[ "$JAVA_VERSION" -lt 11 ]] 2>/dev/null; then
  echo "WARNING: Java version may be too old. InterProScan requires Java 11+" >&2
  echo "Current Java version: $(java -version 2>&1 | head -n1)" >&2
fi

echo "> Using runner: $INTERPRO_BIN"
echo "> Java version: $(java -version 2>&1 | head -n1 | awk -F'"' '{print $2}')"

# ---- Gather protein files from funannotate predict_results -------------------
echo "> Searching for *.proteins.fa under: $fun_dir"
# Refresh the tmp folder (remove old symlinks)
find "$protein_dir" -type l -delete 2>/dev/null || true

mapfile -t prot_sources < <(find "$fun_dir" -type f -path "*/predict_results/*.proteins.fa" | sort)
if [[ ${#prot_sources[@]} -eq 0 ]]; then
  echo "ERROR: No *.proteins.fa found under $fun_dir" >&2
  exit 1
fi

# Create symlinks to protein files
for p in "${prot_sources[@]}"; do
  abs=$(realpath "$p")
  ln -s "$abs" "$protein_dir/$(basename "$abs")"
done

mapfile -t protein_files < <(find "$protein_dir" -type l -name "*.proteins.fa" | sort)
echo "> Linked ${#protein_files[@]} protein file(s) into $protein_dir"

# Test mode: only process first file
if [[ $TEST_MODE -eq 1 ]]; then
  echo "> TEST MODE: Processing only the first protein file"
  protein_files=("${protein_files[0]}")
fi

# ---- build command flags -----------------------------------------------------
THREAD_FLAG=()
THREAD_FLAG=(-cpu "$THREADS")

APP_FLAG=()
if [[ -n "$APPLICATIONS" ]]; then
  APP_FLAG=(-appl "$APPLICATIONS")
else
  # Default: exclude Hamap (bacterial), ProSitePatterns (slow), PRINTS (outdated)
  # This saves ~30-40% runtime while keeping all relevant fungal annotations
  DEFAULT_APPS="AntiFam,CDD,Coils,FunFam,Gene3D,MobiDBLite,NCBIfam,PANTHER,Pfam,PIRSF,PIRSR,ProSiteProfiles,SFLD,SMART,SUPERFAMILY"
  APP_FLAG=(-appl "$DEFAULT_APPS")
  echo "> Using optimized tool set (excludes: Hamap, ProSitePatterns, PRINTS)"
fi

# ---- run InterProScan per file (skip if already has results) -----------------
echo "> Running InterProScan on ${species} protein files"
for protein_file in "${protein_files[@]}"; do
  base="$(basename "$protein_file" .proteins.fa)"
  out_dir="$out_root/$base"
  sentinel="$out_dir/${base}.tsv"
  mkdir -p "$out_dir"

  if [[ -f "$sentinel" ]]; then
    echo "> Skipping $base (results present: ${base}.tsv)"
    continue
  fi

  echo "> Processing $base"
  # InterProScan syntax:
  # interproscan.sh -i <INPUT> -b <BASENAME> -f TSV,GFF3,XML -cpu <N> [-appl <APPS>]
  # Note: -b and -d are mutually exclusive, use -b to specify output basename
  $INTERPRO_BIN \
    -i "$protein_file" \
    -b "${out_dir}/${base}" \
    -f TSV,GFF3,XML \
    "${THREAD_FLAG[@]}" \
    "${APP_FLAG[@]}"
done

# ---- quick summary -----------------------------------------------------------
count_done=$(find "$out_root" -mindepth 1 -maxdepth 1 -type d -exec test -f "{}/*.tsv" ';' -print | wc -l | awk '{print $1}')
echo
echo "==================================="
echo "InterProScan Pipeline Complete"
echo "==================================="
echo "Completed: $count_done run(s)"
echo "Output root: $out_root"
echo
echo "Example output files:"
for sample_dir in "$out_root"/*/; do
  [[ -d "$sample_dir" ]] || continue
  sample=$(basename "$sample_dir")
  if [[ -f "$sample_dir/${sample}.tsv" ]]; then
    echo "  - ${sample}.tsv"
    break
  fi
done
