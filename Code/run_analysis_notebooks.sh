#!/bin/bash
set -euo pipefail

# =============================================================================
# Run all FunPan analysis notebooks (NB0-NB5) in order
# =============================================================================
#
# Usage: bash run_analysis_notebooks.sh [--from NB1] [--only NB3]
#
# Runs notebooks with jupyter nbconvert --execute, saving outputs in-place.
# Results appear in the notebooks themselves AND in NB*_Results/ folders.
# =============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ANALYSIS_DIR="${SCRIPT_DIR}/../Analysis"
cd "$ANALYSIS_DIR"

# Color definitions
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
CYAN='\033[0;36m'
WHITE='\033[1;37m'
NC='\033[0m'

# Notebook list (in dependency order)
declare -a NOTEBOOKS=(
    "NB0_DataPrep.ipynb"
    "NB1_Pangenome.ipynb"
    "NB2_PanGWAS.ipynb"
    "NB3_Convergence.ipynb"
    "NB4_RareGenome.ipynb"
    "NB5_CoreGenome.ipynb"
)

declare -a NB_DESCRIPTIONS=(
    "QC & Phenotype Classification"
    "Pangenome Architecture & Closure"
    "Pan-GWAS & Functional Enrichment"
    "Cross-Species Convergence Testing"
    "Rare Genome & Lifestyle Prediction"
    "Core Genome Cross-Species Conservation"
)

TOTAL=${#NOTEBOOKS[@]}
START_FROM=0
ONLY=""

# Parse arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        --from)
            shift
            for i in $(seq 0 $((TOTAL-1))); do
                if [[ "${NOTEBOOKS[$i]}" == *"$1"* ]]; then
                    START_FROM=$i
                    break
                fi
            done
            ;;
        --only)
            shift
            ONLY="$1"
            ;;
        -h|--help)
            echo "Usage: bash run_analysis_notebooks.sh [--from NB1] [--only NB3]"
            echo ""
            echo "Options:"
            echo "  --from NBx   Start from notebook NBx (skip earlier ones)"
            echo "  --only NBx   Run only notebook NBx"
            echo ""
            echo "Notebooks (in order):"
            for i in $(seq 0 $((TOTAL-1))); do
                printf "  %s  %s\n" "${NOTEBOOKS[$i]}" "${NB_DESCRIPTIONS[$i]}"
            done
            exit 0
            ;;
    esac
    shift
done

# Check jupyter is available
if ! command -v jupyter &>/dev/null; then
    echo -e "${RED}ERROR: jupyter not found. Activate the funpan conda env first:${NC}"
    echo "  conda activate funpan && pip install jupyter nbconvert"
    exit 1
fi

# Header
echo ""
echo -e "${WHITE}================================================================================${NC}"
echo -e "${WHITE}  FunPan Analysis Pipeline — Running Jupyter Notebooks                          ${NC}"
echo -e "${WHITE}================================================================================${NC}"
echo ""

# Status overview
echo -e "${CYAN}Notebooks to run:${NC}"
for i in $(seq 0 $((TOTAL-1))); do
    if [[ -n "$ONLY" ]]; then
        if [[ "${NOTEBOOKS[$i]}" == *"$ONLY"* ]]; then
            echo -e "  ${YELLOW}►${NC} ${NOTEBOOKS[$i]}  ${NB_DESCRIPTIONS[$i]}"
        else
            echo -e "  ${WHITE}○${NC} ${NOTEBOOKS[$i]}  (skipped)"
        fi
    elif [[ $i -ge $START_FROM ]]; then
        echo -e "  ${YELLOW}○${NC} ${NOTEBOOKS[$i]}  ${NB_DESCRIPTIONS[$i]}"
    else
        echo -e "  ${GREEN}✓${NC} ${NOTEBOOKS[$i]}  (skipped)"
    fi
done
echo ""

PASSED=0
FAILED=0
SKIPPED=0

for i in $(seq 0 $((TOTAL-1))); do
    NB="${NOTEBOOKS[$i]}"
    DESC="${NB_DESCRIPTIONS[$i]}"

    # Skip logic
    if [[ -n "$ONLY" && "$NB" != *"$ONLY"* ]]; then
        ((SKIPPED++))
        continue
    fi
    if [[ $i -lt $START_FROM ]]; then
        ((SKIPPED++))
        continue
    fi

    echo -e "${WHITE}────────────────────────────────────────────────────────────────────────────────${NC}"
    echo -e "${YELLOW}► [$((i+1))/${TOTAL}] Running: ${NB}${NC}"
    echo -e "${CYAN}  ${DESC}${NC}"
    echo ""

    START_TIME=$(date +%s)

    if jupyter nbconvert \
        --to notebook \
        --execute \
        --inplace \
        --ExecutePreprocessor.timeout=3600 \
        --ExecutePreprocessor.kernel_name=funpan \
        "$NB" 2>&1 | while IFS= read -r line; do
            # Show progress lines from the notebook
            if [[ "$line" == *"==="* || "$line" == *"Saved"* || "$line" == *"Loading"* || \
                  "$line" == *"WARNING"* || "$line" == *"genomes"* || "$line" == *"OGs"* || \
                  "$line" == *"Complete"* || "$line" == *"Summary"* ]]; then
                echo -e "  ${WHITE}${line}${NC}"
            fi
        done; then

        END_TIME=$(date +%s)
        ELAPSED=$(( END_TIME - START_TIME ))
        MINS=$(( ELAPSED / 60 ))
        SECS=$(( ELAPSED % 60 ))
        echo ""
        echo -e "  ${GREEN}✓ Completed: ${NB} (${MINS}m ${SECS}s)${NC}"
        ((PASSED++))
    else
        END_TIME=$(date +%s)
        ELAPSED=$(( END_TIME - START_TIME ))
        echo ""
        echo -e "  ${RED}✗ FAILED: ${NB} (after ${ELAPSED}s)${NC}"
        echo -e "  ${RED}  Check the notebook for error details.${NC}"
        ((FAILED++))

        # Ask whether to continue
        if [[ $i -lt $((TOTAL-1)) ]]; then
            echo ""
            echo -e "${YELLOW}  Continue with remaining notebooks? (y/n)${NC}"
            read -r CONT
            if [[ "$CONT" != "y" && "$CONT" != "Y" ]]; then
                echo -e "${RED}  Pipeline stopped.${NC}"
                break
            fi
        fi
    fi
done

# Summary
echo ""
echo -e "${WHITE}================================================================================${NC}"
echo -e "${WHITE}  SUMMARY                                                                       ${NC}"
echo -e "${WHITE}================================================================================${NC}"
echo -e "  ${GREEN}Passed:  ${PASSED}${NC}"
echo -e "  ${RED}Failed:  ${FAILED}${NC}"
echo -e "  ${WHITE}Skipped: ${SKIPPED}${NC}"
echo ""

if [[ $FAILED -eq 0 ]]; then
    echo -e "${GREEN}All notebooks completed successfully.${NC}"
    echo ""
    echo -e "${CYAN}Results:${NC}"
    for d in NB0_Results NB1_Results NB2_Results NB3_Results NB4_Results NB5_Results; do
        if [[ -d "$d" ]]; then
            n_files=$(find "$d" -type f | wc -l)
            echo -e "  ${WHITE}${d}/${NC}  (${n_files} files)"
        fi
    done
else
    echo -e "${RED}Some notebooks failed. Fix errors and re-run with --from NBx.${NC}"
fi
echo ""
