#!/bin/bash
set -euo pipefail
trap 'echo "ERROR at line $LINENO"; exit 1' ERR

# =============================================================================
# BiG-SCAPE BGC Similarity Network Analysis (auto-installs from local clone)
# =============================================================================

echo "> Enter species name (genus/species format):"
read species

# ---------------------------
# Config / Paths
# ---------------------------
ANTISMASH_DIR="../Species/$species/antismash_output"
RENAMED_DIR="$ANTISMASH_DIR/antismash_gbk_unique"
BIGSCAPE_DIR="../Species/$species/bigscape_output"
PFAM_PATH="../Data/BiG-SCAPE/Pfam-A.hmm"                  # auto-download/index if missing
BIGSCAPE_SRC="${BIGSCAPE_SRC:-../Data/BiG-SCAPE}"        # local BiG-SCAPE clone

# ---------------------------
# Helpers
# ---------------------------
ensure_hmmer() {
  if ! command -v hmmpress >/dev/null 2>&1; then
    echo "ERROR: hmmpress (HMMER) not found. Install: conda install -c bioconda hmmer"
    exit 1
  fi
}

ensure_pfam() {
  local pfam_path="$1"
  local pfam_dir
  pfam_dir="$(dirname "$pfam_path")"
  mkdir -p "$pfam_dir"

  if [[ ! -s "$pfam_path" ]]; then
    echo "PFAM not found at $pfam_path — downloading..."
    local gz="${pfam_path}.gz"
    local url_primary="https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz"
    local url_backup="https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam35.0/Pfam-A.hmm.gz"
    if command -v wget >/dev/null 2>&1; then
      wget -O "$gz" "$url_primary" || wget -O "$gz" "$url_backup"
    elif command -v curl >/dev/null 2>&1; then
      curl -L "$url_primary" -o "$gz" || curl -L "$url_backup" -o "$gz"
    else
      echo "ERROR: Need wget or curl to download Pfam."
      exit 1
    fi
    gunzip -f "$gz"
  fi

  # Index if any HMMER index is missing
  if [[ ! -s "${pfam_path}.h3m" || ! -s "${pfam_path}.h3i" || ! -s "${pfam_path}.h3f" || ! -s "${pfam_path}.h3p" ]]; then
    echo "Indexing Pfam with hmmpress..."
    hmmpress "$pfam_path"
  fi
  echo "Pfam ready at $pfam_path"
}


ensure_bigscape_src() {
  local target_dir="$BIGSCAPE_SRC"
  if [[ -d "$target_dir" && -f "$target_dir/bigscape.py" ]]; then
    echo "BiG-SCAPE source found at $target_dir"
    return 0
  fi
  echo "BiG-SCAPE source not found — cloning into $target_dir..."
  mkdir -p "$(dirname "$target_dir")"
  if command -v git >/dev/null 2>&1; then
    git clone https://github.com/medema-group/BiG-SCAPE.git "$target_dir"
  else
    echo "ERROR: git not installed — cannot clone BiG-SCAPE"
    exit 1
  fi
}


ensure_bigscape_cli() {
  if command -v bigscape >/dev/null 2>&1; then
    return 0
  fi
  echo "bigscape CLI not found — installing from local clone: $BIGSCAPE_SRC"

  if [[ ! -d "$BIGSCAPE_SRC" ]]; then
    echo "ERROR: BIGSCAPE_SRC not found at $BIGSCAPE_SRC"
    echo "       Set BIGSCAPE_SRC=/path/to/BiG-SCAPE or clone the repo there."
    exit 1
  fi

  if ! python -m pip --version >/dev/null 2>&1; then
    echo "ERROR: pip not available in this environment."
    exit 1
  fi

  # Try standard install; fall back to editable
  if (cd "$BIGSCAPE_SRC" && python -m pip install .); then
    :
  else
    echo "Standard install failed, trying editable install..."
    (cd "$BIGSCAPE_SRC" && python -m pip install -e .)
  fi

  hash -r
  if ! command -v bigscape >/dev/null 2>&1; then
    echo "ERROR: bigscape still not available after installation."
    exit 1
  fi
  echo "bigscape installed and available"
}

# ---------------------------
# Pre-checks
# ---------------------------
if [[ ! -d "$ANTISMASH_DIR" ]]; then
  echo "ERROR: antiSMASH output directory not found: $ANTISMASH_DIR"
  exit 1
fi

# Collect region GBKs, but **exclude** RENAMED_DIR to avoid self-copies
mapfile -t REGION_GBKS < <(
  find -L "$ANTISMASH_DIR" -path "$RENAMED_DIR" -prune -o -type f -name "*.region*.gbk" -print | sort
)

if [[ ${#REGION_GBKS[@]} -eq 0 ]]; then
  echo "ERROR: No antiSMASH region files (*.region*.gbk) found in $ANTISMASH_DIR"
  exit 1
fi

echo "Found region files in these genome dirs:"
printf "%s\n" "${REGION_GBKS[@]}" | awk -F/ '{NF--; print $0}' | sort -u

# ---------------------------
# Prepare renamed BGCs
# ---------------------------
mkdir -p "$RENAMED_DIR"
# Do NOT nuke existing renamed files blindly; BiG-SCAPE can be re-run incrementally.

echo "Renaming antiSMASH BGC files with genome prefix..."
rename_count=0

for gbk_file in "${REGION_GBKS[@]}"; do
  genome_dir="$(dirname "$gbk_file")"
  genome="$(basename "$genome_dir")"

  # accession from genome dir (fallback to dir name)
  accession="$genome"
  if [[ "$genome" =~ (G[CA]F?_[0-9]+\.[0-9]+) ]]; then
    accession="${BASH_REMATCH[1]}"
  fi

  gbk_base="$(basename "$gbk_file")"
  region="${gbk_base%.gbk}"
  if [[ "$gbk_base" =~ (scaffold[^/]*\.region[0-9]+) ]]; then
    region="${BASH_REMATCH[1]}"
  fi

  dest="$RENAMED_DIR/${accession}_${region}.gbk"

  # Skip if dest is literally the same file (paranoia) or already exists
  if [[ -e "$dest" ]]; then
    # If it's the same inode (same file), skip; else keep existing to avoid overwrite
    if [[ "$gbk_file" -ef "$dest" ]]; then
      echo "Skipping self-copy: $dest"
    else
      echo "Skipping existing: $dest"
    fi
    continue
  fi

  cp "$gbk_file" "$dest"
  rename_count=$((rename_count + 1))
done

if [[ "$rename_count" -eq 0 ]]; then
  echo "No new BGC region files needed renaming (all up to date)."
else
  echo "Renamed $rename_count new .gbk files into $RENAMED_DIR"
fi

# ---------------------------
# Ensure deps & DBs
# ---------------------------
ensure_hmmer
ensure_pfam "$PFAM_PATH"
ensure_bigscape_src
ensure_bigscape_cli

# ---------------------------
# Run BiG-SCAPE
# ---------------------------
echo "Running BiG-SCAPE on renamed BGCs..."
mkdir -p "$BIGSCAPE_DIR"

if ! bigscape cluster -i "$RENAMED_DIR" -o "$BIGSCAPE_DIR" \
  --pfam-path "$PFAM_PATH" --mibig-version 4.0 --gcf-cutoffs 0.3; then
  echo "ERROR: BiG-SCAPE failed during execution."
  exit 1
fi

# ---------------------------
# Output checks / Post-process
# ---------------------------
#if [[ ! -d "$BIGSCAPE_DIR/network_files" ]]; then
#  echo "ERROR: BiG-SCAPE did not produce expected network_files output."
#  exit 1
#fi

if [[ -f "$BIGSCAPE_DIR/bigscape_output.db" ]]; then
  mv "$BIGSCAPE_DIR/bigscape_output.db" "$BIGSCAPE_DIR/data_sqlite.db"
  echo "Renamed bigscape_output.db to data_sqlite.db"
fi

if [[ -f "$BIGSCAPE_DIR/index.html" ]]; then
  sed -i 's|/data_sqlite.db|data_sqlite.db|g' "$BIGSCAPE_DIR/index.html"
  echo "Patched index.html to fix database reference"
fi

echo ""
echo "BiG-SCAPE completed successfully"
echo "Output directory: $BIGSCAPE_DIR"
