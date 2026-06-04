#!/bin/bash
set -euo pipefail

# =============================================================================
# Complete Genome Analysis Pipeline
# =============================================================================
# This script combines all analysis steps in the correct order:
#  1. Download genomes and run BUSCO analysis
#  2. ANI species verification, Mash deduplication, and QC filtering
#  3. Run Funannotate annotation
#  4. Run EggNOG functional annotation
#  5. Run OrthoFinder for comparative genomics
#  6. Run antiSMASH for secondary metabolite detection
#  7. Run InterProScan for protein domain annotation
#  8. Run SignalP for signal peptide prediction
#  9. Run dbCAN for CAZyme annotation
# 10. Run Parsnp for core genome alignment
# 11. Run BiG-SCAPE for biosynthetic gene cluster analysis
# 12. Run Gubbins for recombination detection
# 13. Run IQ-TREE on Gubbins output (SNV-based tree)
#
# Dependency structure:
#   Steps 1 -> 2 -> 3 are sequential
#   Steps 4-10 all depend on Step 3 (independent of each other)
#   Step 11 depends on Step 6 (antiSMASH)
#   Step 12 depends on Step 10 (Parsnp)
#   Step 13 depends on Step 12 (Gubbins)
# =============================================================================


SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Create a folder for all pipeline logs
LOG_DIR="${SCRIPT_DIR}/.pipeline_logfiles"
mkdir -p "$LOG_DIR"

# Timestamped log file inside the log folder
LOG_FILE="${LOG_DIR}/pipeline_$(date +%Y%m%d_%H%M%S).log"

# Pipeline configuration
TOTAL_STEPS=13
CURRENT_STEP=0

# Dynamic CPU/Thread Detection
CPU_CORES=$(nproc)
CPU_THREADS=$(lscpu | grep '^Thread(s) per core:' | awk '{print $4}')
TOTAL_CPUS=$((CPU_CORES * CPU_THREADS))

# Use 75% of available CPUs to leave room for system processes
PIPELINE_THREADS=$((TOTAL_CPUS * 3 / 4))
# Ensure at least 1 thread
if [ $PIPELINE_THREADS -lt 1 ]; then
    PIPELINE_THREADS=1
fi

echo "Detected $CPU_CORES cores with $CPU_THREADS threads per core (Total: $TOTAL_CPUS CPUs)"
echo "Pipeline will use $PIPELINE_THREADS threads for computations"

# Checkpoint marker directory
MARKER_DIR="${SCRIPT_DIR}/.pipeline_checkpoints"

# Step progress tracking
STEP_PROGRESS_CURRENT=0
STEP_PROGRESS_TOTAL=0
STEP_PROGRESS_LABEL=""

# Color definitions
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
PURPLE='\033[0;35m'
CYAN='\033[0;36m'
WHITE='\033[1;37m'
NC='\033[0m' # No Color

# =============================================================================
# Pre-flight: check required tools & Python packages in the current environment
# =============================================================================
check_environment() {
    local repo_root="${SCRIPT_DIR}/.."
    local env_yml="${repo_root}/environment.yml"

    echo ""
    echo -e "${CYAN}Checking required tools and Python packages...${NC}"

    # Core CLI tools the pipeline runs from the active (funpan) environment
    local req_cmds=(python datasets fastANI mash busco funannotate emapper.py \
                    orthofinder run_dbcan parsnp harvesttools iqtree2 mafft diamond blastp jupyter)
    # Core Python packages used by the analysis notebooks/modules
    local req_py=(pandas numpy scipy matplotlib seaborn statsmodels sklearn Bio requests)
    # Tools that live in their OWN envs/installs (warn only, do not block)
    local ext_cmds=(antismash interproscan.sh signalp6 bigscape run_gubbins.py)

    local missing_cmd=() missing_py=() missing_ext=()
    local c m
    for c in "${req_cmds[@]}"; do command -v "$c" >/dev/null 2>&1 || missing_cmd+=("$c"); done
    for m in "${req_py[@]}";  do python -c "import $m" >/dev/null 2>&1 || missing_py+=("$m"); done
    for c in "${ext_cmds[@]}"; do command -v "$c" >/dev/null 2>&1 || missing_ext+=("$c"); done
    python -c "import nbconvert" >/dev/null 2>&1 || missing_ext+=("nbconvert(py)")

    if [ ${#missing_ext[@]} -gt 0 ]; then
        echo -e "${YELLOW}  Note: optional/externally-managed tools not on PATH: ${missing_ext[*]}${NC}"
        echo -e "${YELLOW}        (these run in their own envs/installs: antismash8, gubbins, InterProScan, SignalP, BiG-SCAPE)${NC}"
    fi

    if [ ${#missing_cmd[@]} -eq 0 ] && [ ${#missing_py[@]} -eq 0 ]; then
        echo -e "${GREEN}  All core tools and Python packages found.${NC}"
        echo ""
        return 0
    fi

    echo -e "${RED}  Missing core requirements:${NC}"
    [ ${#missing_cmd[@]} -gt 0 ] && echo -e "${RED}    Tools:          ${missing_cmd[*]}${NC}"
    [ ${#missing_py[@]}  -gt 0 ] && echo -e "${RED}    Python modules: ${missing_py[*]}${NC}"
    echo ""
    echo -e "${YELLOW}  Recommended: create and activate the 'funpan' conda environment:${NC}"
    if [ -f "$env_yml" ]; then
        echo -e "${WHITE}    conda env create -f ${env_yml}${NC}"
    else
        echo -e "${WHITE}    conda env create -f environment.yml${NC}"
    fi
    echo -e "${WHITE}    conda activate funpan${NC}"
    if [ "${CONDA_DEFAULT_ENV:-}" = "funpan" ]; then
        echo ""
        echo -e "${YELLOW}  (You appear to be in the 'funpan' env but items are still missing — the env may be incomplete; try re-creating it.)${NC}"
    fi
    echo ""
    exit 1
}

check_environment

# Step names array for display
declare -a STEP_NAMES=(
    ""  # index 0 unused
    "Download genomes and run BUSCO analysis"
    "ANI verify, deduplicate, and QC-filter genomes"
    "Run Funannotate annotation"
    "Run EggNOG functional annotation"
    "Run OrthoFinder for comparative genomics"
    "Run antiSMASH for secondary metabolite detection"
    "Run InterProScan for protein domain annotation"
    "Run SignalP for signal peptide prediction"
    "Run dbCAN for CAZyme annotation"
    "Run Parsnp for core genome alignment"
    "Run BiG-SCAPE for BGC network analysis"
    "Run Gubbins for recombination detection"
    "Run IQ-TREE on Gubbins output (SNV-based tree)"
)

declare -a STEP_DESCRIPTIONS=(
    ""  # index 0 unused
    "Genome Download and Quality Assessment"
    "ANI Species Check, Mash Deduplication & QC Filtering"
    "Genome Annotation with Funannotate"
    "Functional Annotation with EggNOG"
    "Comparative Genomics with OrthoFinder"
    "Secondary Metabolite Detection with antiSMASH"
    "Protein Domain Annotation with InterProScan"
    "Signal Peptide Prediction with SignalP"
    "CAZyme Annotation with dbCAN"
    "Core Genome Alignment with Parsnp"
    "BGC Network Analysis with BiG-SCAPE"
    "Recombination Detection with Gubbins"
    "SNV-based Phylogeny with IQ-TREE"
)

# Function to print a progress bar
print_progress_bar() {
    local current=$1
    local total=$2
    local width=50
    local percentage=$((current * 100 / total))
    local filled=$((current * width / total))
    local empty=$((width - filled))

    printf "${CYAN}Progress: [${NC}"
    printf "%*s" $filled | tr ' ' '#'
    printf "%*s" $empty | tr ' ' '-'
    printf "${CYAN}] %d%% (%d/%d)${NC}" $percentage $current $total
}

# Function to update step progress
update_step_progress() {
    local current=$1
    local total=$2
    local label="$3"

    STEP_PROGRESS_CURRENT=$current
    STEP_PROGRESS_TOTAL=$total
    STEP_PROGRESS_LABEL="$label"

    # Refresh the display
    show_fixed_header
}

# Function to reset step progress
reset_step_progress() {
    STEP_PROGRESS_CURRENT=0
    STEP_PROGRESS_TOTAL=0
    STEP_PROGRESS_LABEL=""
}

# Function to clear screen and show fixed header
show_fixed_header() {
    clear
    echo -e "${WHITE}################################################################################${NC}"
    echo -e "${WHITE}#                           PAN-GENOME ANALYSIS PIPELINE                       #${NC}"
    echo -e "${WHITE}################################################################################${NC}"
    echo ""
    echo -e "${CYAN}Pipeline Steps Overview:${NC}"

    # Show all steps with status indicators
    for i in $(seq 1 $TOTAL_STEPS); do
        local status_icon
        local color
        if [ $i -lt $CURRENT_STEP ]; then
            status_icon="✓"
            color="$GREEN"
        elif [ $i -eq $CURRENT_STEP ]; then
            status_icon="►"
            color="$YELLOW"
        else
            status_icon="○"
            color="$WHITE"
        fi

        printf "  ${color}${status_icon} Step %2d: %s${NC}\n" "$i" "${STEP_NAMES[$i]}"
    done

    echo ""
    # Show current step with progress if available
    if [ $CURRENT_STEP -gt 0 ] && [ $CURRENT_STEP -le $TOTAL_STEPS ]; then
        echo -e "${CYAN}Step ${CURRENT_STEP}: ${STEP_DESCRIPTIONS[$CURRENT_STEP]}${NC}"

        # Show progress bar if step is in progress
        if [ $STEP_PROGRESS_TOTAL -gt 0 ]; then
            print_progress_bar $STEP_PROGRESS_CURRENT $STEP_PROGRESS_TOTAL
            echo ""
        fi

        echo -e "${WHITE}=================================================================================${NC}"
        echo ""
    else
        echo ""
    fi
}

# Function to log messages with enhanced formatting
log_message() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a "$LOG_FILE"
}

# Function to display step header with progress
show_step_header() {
    local step_num=$1
    local step_name="$2"
    local step_description="$3"

    CURRENT_STEP=$step_num

    # Show fixed header with updated current step
    show_fixed_header

    echo -e "${YELLOW}${step_description}${NC}"
    echo ""
}

# Function to show step completion
show_step_completion() {
    local step_name="$1"

    # Update the fixed header to show completion
    show_fixed_header

    echo -e "${GREEN}[COMPLETED] Step ${CURRENT_STEP}/${TOTAL_STEPS}: ${step_name}${NC}"
    echo ""
    sleep 2  # Brief pause to show completion before next step
}

# Function to show step failure
show_step_failure() {
    local step_name=$1
    echo ""
    echo -e "${RED}[FAILED] Step ${CURRENT_STEP}/${TOTAL_STEPS}: ${step_name}${NC}"
    echo -e "${RED}Pipeline execution stopped.${NC}"
    echo ""
}

# Function to show substep progress
show_substep() {
    local substep_name="$1"
    echo -e "${YELLOW}> ${substep_name}...${NC}"
}

# Function to show substep with progress update
show_substep_with_progress() {
    local substep_name="$1"
    local current=$2
    local total=$3
    local label="$4"

    update_step_progress $current $total "$label"
    echo -e "${YELLOW}> ${substep_name}...${NC}"
}

# Function to create checkpoint marker
create_checkpoint() {
    local step_num=$1
    local step_name="$2"
    mkdir -p "$MARKER_DIR"
    echo "$(date '+%Y-%m-%d %H:%M:%S')" > "${MARKER_DIR}/step_${step_num}.completed"
    log_message "Checkpoint created for Step ${step_num}: ${step_name}"
}

# Function to check if step is completed
is_step_completed() {
    local step_num=$1
    [[ -f "${MARKER_DIR}/step_${step_num}.completed" ]]
}

# Function to find first incomplete step
find_first_incomplete_step() {
    for i in $(seq 1 $TOTAL_STEPS); do
        if ! is_step_completed $i; then
            echo $i
            return
        fi
    done
    echo $((TOTAL_STEPS + 1))  # All steps completed
}

# Function to clear all checkpoints (for fresh start)
clear_checkpoints() {
    if [[ -d "$MARKER_DIR" ]]; then
        rm -rf "$MARKER_DIR"
        log_message "All checkpoints cleared"
    fi
}

# Function to run a script and handle errors
run_script() {
    local script_name="$1"
    local step_name="$2"
    local script_path="${SCRIPT_DIR}/${script_name}"

    if [[ ! -f "$script_path" ]]; then
        log_message "ERROR: Script not found: $script_path"
        show_step_failure "$step_name"
        exit 1
    fi

    # Show step progress phases
    show_substep_with_progress "Initializing $script_name" 1 3 "$step_name"
    log_message "Starting: $script_name"

    show_substep_with_progress "Executing $script_name" 2 3 "$step_name"
    if bash "$script_path"; then
        show_substep_with_progress "Completing $script_name" 3 3 "$step_name"
        log_message "Completed successfully: $script_name"
        reset_step_progress
        show_step_completion "$step_name"
        create_checkpoint "$CURRENT_STEP" "$step_name"
    else
        log_message "ERROR: Failed to complete: $script_name"
        reset_step_progress
        show_step_failure "$step_name"
        exit 1
    fi
}

# Function to run script with genus and species inputs
run_script_with_genus_species() {
    local script_name="$1"
    local step_name="$2"
    local script_path="${SCRIPT_DIR}/${script_name}"

    if [[ ! -f "$script_path" ]]; then
        log_message "ERROR: Script not found: $script_path"
        show_step_failure "$step_name"
        exit 1
    fi

    # Show step progress phases
    show_substep_with_progress "Preparing inputs for $script_name" 1 4 "$step_name"
    log_message "Starting: $script_name (genus: $GENUS, species: $SPECIES_NAME)"

    show_substep_with_progress "Processing genus: $GENUS" 2 4 "$step_name"
    show_substep_with_progress "Processing species: $SPECIES_NAME" 3 4 "$step_name"

    show_substep_with_progress "Executing $script_name" 4 4 "$step_name"
    # Pipe genus and species inputs to the script
    if printf "%s\n%s\n" "$GENUS" "$SPECIES_NAME" | bash "$script_path"; then
        log_message "Completed successfully: $script_name"
        reset_step_progress
        show_step_completion "$step_name"
        create_checkpoint "$CURRENT_STEP" "$step_name"
    else
        log_message "ERROR: Failed to complete: $script_name"
        reset_step_progress
        show_step_failure "$step_name"
        exit 1
    fi
}

# Function to run script with species input (genus/species format)
run_script_with_species() {
    local script_name="$1"
    local step_name="$2"
    local script_path="${SCRIPT_DIR}/${script_name}"

    if [[ ! -f "$script_path" ]]; then
        log_message "ERROR: Script not found: $script_path"
        show_step_failure "$step_name"
        exit 1
    fi

    # Show step progress phases
    show_substep_with_progress "Preparing species data: $FULL_SPECIES" 1 3 "$step_name"
    log_message "Starting: $script_name (species: $FULL_SPECIES)"

    show_substep_with_progress "Executing $script_name" 2 3 "$step_name"

    # Handle different script input methods
    case "$script_name" in
        "ani_and_filter_genome_QC.sh")
            # This script accepts command line arguments but may exit with non-zero even on success
            # Run without strict error handling to prevent pipeline exit
            set +e
            bash "$script_path" "$FULL_SPECIES" 2>&1
            exit_code=$?
            set -e
            log_message "Filter script exited with code $exit_code"
            # Check if files were actually created (success indicator)
            if [[ -d "/datadrive/Species/$GENUS/$SPECIES_NAME/filtered_genome" ]] &&
               [[ $(ls -1 "/datadrive/Species/$GENUS/$SPECIES_NAME/filtered_genome" 2>/dev/null | wc -l) -gt 0 ]]; then
                success=true
                log_message "Filter script produced output files - treating as successful"
            else
                success=false
                log_message "Filter script did not produce expected output files"
            fi
            ;;
        *)
            # Most other scripts expect interactive input via stdin
            if printf "%s\n" "$FULL_SPECIES" | bash "$script_path"; then
                success=true
            else
                success=false
            fi
            ;;
    esac

    if [[ "$success" == "true" ]]; then
        show_substep_with_progress "Finalizing $script_name results" 3 3 "$step_name"
        log_message "Completed successfully: $script_name"
        reset_step_progress
        show_step_completion "$step_name"
        create_checkpoint "$CURRENT_STEP" "$step_name"
    else
        log_message "ERROR: Failed to complete: $script_name"
        reset_step_progress
        show_step_failure "$step_name"
        exit 1
    fi
}

# =============================================================================
# Main Pipeline Execution
# =============================================================================

log_message "Starting complete genome analysis pipeline"
log_message "Log file: $LOG_FILE"

# Initialize current step to 0 and show initial header
CURRENT_STEP=0
show_fixed_header

# Check for existing checkpoints and present options to user
FIRST_INCOMPLETE_STEP=$(find_first_incomplete_step)

if [[ $FIRST_INCOMPLETE_STEP -le $TOTAL_STEPS ]]; then
    # Load previous species name if available
    PREVIOUS_SPECIES=""
    if [[ -f "${MARKER_DIR}/pipeline_inputs.txt" ]]; then
        source "${MARKER_DIR}/pipeline_inputs.txt"
        PREVIOUS_SPECIES=" for $GENUS $SPECIES_NAME"
    fi

    echo -e "${YELLOW}Previous pipeline run detected${PREVIOUS_SPECIES}. Checkpoints found for steps 1-$((FIRST_INCOMPLETE_STEP-1)).${NC}"
    echo -e "${CYAN}Options:${NC}"
    echo "  1) Resume from step $FIRST_INCOMPLETE_STEP (recommended)"
    echo "  2) Start fresh (clear all checkpoints)"
    echo "  3) Start from a specific step"
    echo ""
    echo -e "${YELLOW}> Choose option (1, 2, or 3):${NC}"
    read RESUME_CHOICE

    case "$RESUME_CHOICE" in
        1)
            echo -e "${GREEN}Resuming pipeline from step $FIRST_INCOMPLETE_STEP.${NC}"
            ;;
        2)
            clear_checkpoints
            FIRST_INCOMPLETE_STEP=1
            echo -e "${GREEN}All checkpoints cleared. Starting fresh pipeline.${NC}"
            ;;
        3)
            echo -e "${CYAN}Select step to start from:${NC}"
            for i in $(seq 1 $TOTAL_STEPS); do
                printf "  %2d) %s\n" "$i" "${STEP_NAMES[$i]}"
            done
            echo ""
            echo -e "${YELLOW}> Enter step number (1-${TOTAL_STEPS}):${NC}"
            read START_STEP

            if [[ "$START_STEP" =~ ^[0-9]+$ ]] && [ "$START_STEP" -ge 1 ] && [ "$START_STEP" -le $TOTAL_STEPS ]; then
                FIRST_INCOMPLETE_STEP=$START_STEP
                echo -e "${GREEN}Starting pipeline from step $START_STEP.${NC}"
            else
                echo -e "${RED}Invalid step number. Starting from step 1.${NC}"
                FIRST_INCOMPLETE_STEP=1
            fi
            ;;
        *)
            echo -e "${YELLOW}Invalid option. Resuming from step $FIRST_INCOMPLETE_STEP.${NC}"
            ;;
    esac
    echo ""
else
    # No checkpoints found, show start options
    echo -e "${CYAN}No previous pipeline run detected.${NC}"
    echo -e "${CYAN}Options:${NC}"
    echo "  1) Start from the beginning"
    echo "  2) Start from a specific step"
    echo ""
    echo -e "${YELLOW}> Choose option (1 or 2):${NC}"
    read START_CHOICE

    case "$START_CHOICE" in
        1)
            FIRST_INCOMPLETE_STEP=1
            echo -e "${GREEN}Starting pipeline from the beginning.${NC}"
            ;;
        2)
            echo -e "${CYAN}Select step to start from:${NC}"
            for i in $(seq 1 $TOTAL_STEPS); do
                printf "  %2d) %s\n" "$i" "${STEP_NAMES[$i]}"
            done
            echo ""
            echo -e "${YELLOW}> Enter step number (1-${TOTAL_STEPS}):${NC}"
            read START_STEP

            if [[ "$START_STEP" =~ ^[0-9]+$ ]] && [ "$START_STEP" -ge 1 ] && [ "$START_STEP" -le $TOTAL_STEPS ]; then
                FIRST_INCOMPLETE_STEP=$START_STEP
                echo -e "${GREEN}Starting pipeline from step $START_STEP.${NC}"
            else
                echo -e "${RED}Invalid step number. Starting from step 1.${NC}"
                FIRST_INCOMPLETE_STEP=1
            fi
            ;;
        *)
            echo -e "${YELLOW}Invalid option. Starting from step 1.${NC}"
            FIRST_INCOMPLETE_STEP=1
            ;;
    esac
    echo ""
fi

# Collect all required inputs upfront
if [[ $FIRST_INCOMPLETE_STEP -eq 1 ]] || [[ -z "${GENUS:-}" ]]; then
    echo "Please provide the following information for the pipeline:"
    echo ""
    echo "> What is the genus name?"
    read GENUS
    echo "> What is the species name? $GENUS __ :"
    read SPECIES_NAME
    # Remove genus + _ if it's included in the input
    SPECIES_NAME=${SPECIES_NAME#$GENUS }
    FULL_SPECIES="$GENUS/$SPECIES_NAME"

    echo ""
    echo "Pipeline will analyze: $GENUS $SPECIES_NAME"
    echo "Species path format: $FULL_SPECIES"
    echo ""

    # Save inputs to a temporary file for resume
    mkdir -p "$MARKER_DIR"
    echo "GENUS=$GENUS" > "${MARKER_DIR}/pipeline_inputs.txt"
    echo "SPECIES_NAME=$SPECIES_NAME" >> "${MARKER_DIR}/pipeline_inputs.txt"
    echo "FULL_SPECIES=$FULL_SPECIES" >> "${MARKER_DIR}/pipeline_inputs.txt"
else
    # Load saved inputs from previous run
    if [[ -f "${MARKER_DIR}/pipeline_inputs.txt" ]]; then
        source "${MARKER_DIR}/pipeline_inputs.txt"
        echo "Loaded previous inputs: $GENUS $SPECIES_NAME"
        echo ""
    else
        echo "ERROR: Cannot resume - previous inputs not found. Please start fresh."
        exit 1
    fi
fi

# Export variables for use in subscripts
export GENUS
export SPECIES_NAME
export FULL_SPECIES
export PIPELINE_THREADS

log_message "Pipeline configured for species: $GENUS $SPECIES_NAME"

# =============================================================================
# Step 1: Download genomes and BUSCO
# =============================================================================
if [[ $FIRST_INCOMPLETE_STEP -le 1 ]]; then
    show_step_header 1 "${STEP_DESCRIPTIONS[1]}" "Downloading genomes from NCBI and running BUSCO quality analysis"
    run_script_with_genus_species "download_genome_and_BUSCO.sh" "${STEP_DESCRIPTIONS[1]}"
else
    log_message "Skipping Step 1: Already completed (checkpoint found)"
fi

# =============================================================================
# Step 2: Filter genomes
# =============================================================================
if [[ $FIRST_INCOMPLETE_STEP -le 2 ]]; then
    show_step_header 2 "${STEP_DESCRIPTIONS[2]}" "Verifying species identity via ANI, removing near-duplicates with Mash, and applying QC thresholds"
    run_script_with_species "ani_and_filter_genome_QC.sh" "${STEP_DESCRIPTIONS[2]}"
else
    log_message "Skipping Step 2: Already completed (checkpoint found)"
fi

# =============================================================================
# Step 3: Funannotate annotation
# =============================================================================
if [[ $FIRST_INCOMPLETE_STEP -le 3 ]]; then
    show_step_header 3 "${STEP_DESCRIPTIONS[3]}" "Performing structural and functional genome annotation"
    run_script_with_species "run_funannotate.sh" "${STEP_DESCRIPTIONS[3]}"
else
    log_message "Skipping Step 3: Already completed (checkpoint found)"
fi

# =============================================================================
# Steps 4-10: Independent annotation and analysis steps (all depend on Step 3)
# =============================================================================

# Step 4: EggNOG annotation
if [[ $FIRST_INCOMPLETE_STEP -le 4 ]]; then
    show_step_header 4 "${STEP_DESCRIPTIONS[4]}" "Adding orthology and functional annotations to predicted genes"
    run_script_with_species "run_eggnog.sh" "${STEP_DESCRIPTIONS[4]}"
else
    log_message "Skipping Step 4: Already completed (checkpoint found)"
fi

# Step 5: OrthoFinder
if [[ $FIRST_INCOMPLETE_STEP -le 5 ]]; then
    show_step_header 5 "${STEP_DESCRIPTIONS[5]}" "Identifying orthologous gene groups across genomes"
    run_script_with_species "run_orthofinder.sh" "${STEP_DESCRIPTIONS[5]}"
else
    log_message "Skipping Step 5: Already completed (checkpoint found)"
fi

# Step 6: antiSMASH
if [[ $FIRST_INCOMPLETE_STEP -le 6 ]]; then
    show_step_header 6 "${STEP_DESCRIPTIONS[6]}" "Detecting and analyzing biosynthetic gene clusters"
    run_script_with_species "run_antismash.sh" "${STEP_DESCRIPTIONS[6]}"
else
    log_message "Skipping Step 6: Already completed (checkpoint found)"
fi

# Step 7: InterProScan
if [[ $FIRST_INCOMPLETE_STEP -le 7 ]]; then
    show_step_header 7 "${STEP_DESCRIPTIONS[7]}" "Scanning proteins against InterPro signature databases"
    run_script_with_species "run_interproscan.sh" "${STEP_DESCRIPTIONS[7]}"
else
    log_message "Skipping Step 7: Already completed (checkpoint found)"
fi

# Step 8: SignalP
if [[ $FIRST_INCOMPLETE_STEP -le 8 ]]; then
    show_step_header 8 "${STEP_DESCRIPTIONS[8]}" "Predicting signal peptides in protein sequences"
    run_script_with_species "run_signalp.sh" "${STEP_DESCRIPTIONS[8]}"
else
    log_message "Skipping Step 8: Already completed (checkpoint found)"
fi

# Step 9: dbCAN
if [[ $FIRST_INCOMPLETE_STEP -le 9 ]]; then
    show_step_header 9 "${STEP_DESCRIPTIONS[9]}" "Annotating carbohydrate-active enzymes (CAZymes)"
    run_script_with_species "run_dbcan.sh" "${STEP_DESCRIPTIONS[9]}"
else
    log_message "Skipping Step 9: Already completed (checkpoint found)"
fi

# Step 10: Parsnp
if [[ $FIRST_INCOMPLETE_STEP -le 10 ]]; then
    show_step_header 10 "${STEP_DESCRIPTIONS[10]}" "Building core genome alignment for phylogenetic analysis"
    run_script_with_species "run_parsnp.sh" "${STEP_DESCRIPTIONS[10]}"
else
    log_message "Skipping Step 10: Already completed (checkpoint found)"
fi

# =============================================================================
# Steps 11-13: Downstream analyses with specific dependencies
# =============================================================================

# Step 11: BiG-SCAPE (depends on Step 6: antiSMASH)
if [[ $FIRST_INCOMPLETE_STEP -le 11 ]]; then
    show_step_header 11 "${STEP_DESCRIPTIONS[11]}" "Comparing biosynthetic gene clusters and building similarity networks"
    run_script_with_species "run_bigscape.sh" "${STEP_DESCRIPTIONS[11]}"
else
    log_message "Skipping Step 11: Already completed (checkpoint found)"
fi

# Step 12: Gubbins (depends on Step 10: Parsnp)
if [[ $FIRST_INCOMPLETE_STEP -le 12 ]]; then
    show_step_header 12 "${STEP_DESCRIPTIONS[12]}" "Detecting and removing recombination from core genome alignment"
    run_script_with_species "run_gubbins.sh" "${STEP_DESCRIPTIONS[12]}"
else
    log_message "Skipping Step 12: Already completed (checkpoint found)"
fi

# Step 13: IQ-TREE on Gubbins output (depends on Step 12: Gubbins)
if [[ $FIRST_INCOMPLETE_STEP -le 13 ]]; then
    show_step_header 13 "${STEP_DESCRIPTIONS[13]}" "Building phylogenetic tree from recombination-filtered SNVs"
    run_script_with_species "run_iqtree_gubbins.sh" "${STEP_DESCRIPTIONS[13]}"
else
    log_message "Skipping Step 13: Already completed (checkpoint found)"
fi

# =============================================================================
# Pipeline completion
# =============================================================================
CURRENT_STEP=$TOTAL_STEPS
show_fixed_header
echo -e "${WHITE}################################################################################${NC}"
echo -e "${WHITE}#                           PIPELINE COMPLETED SUCCESSFULLY!                  #${NC}"
echo -e "${WHITE}################################################################################${NC}"
echo ""
echo -e "${GREEN}[SUCCESS] All ${TOTAL_STEPS} steps completed successfully for ${GENUS} ${SPECIES_NAME}${NC}"
echo ""
log_message "Complete pipeline finished successfully"
log_message "Check individual output directories for results"

echo -e "${CYAN}Pipeline execution log: ${WHITE}$LOG_FILE${NC}"
echo ""
echo -e "${YELLOW}Results can be found in the following directories:${NC}"
echo -e "${WHITE}  Genome data:     ${NC}../Species/${GENUS}/${SPECIES_NAME}/genome"
echo -e "${WHITE}  Funannotate:     ${NC}../Species/${GENUS}/${SPECIES_NAME}/funannotate_output/"
echo -e "${WHITE}  EggNOG:          ${NC}../Species/${GENUS}/${SPECIES_NAME}/eggnog_output/"
echo -e "${WHITE}  OrthoFinder:     ${NC}../Species/${GENUS}/${SPECIES_NAME}/orthofinder_output/"
echo -e "${WHITE}  antiSMASH:       ${NC}../Species/${GENUS}/${SPECIES_NAME}/antismash_output/"
echo -e "${WHITE}  InterProScan:    ${NC}../Species/${GENUS}/${SPECIES_NAME}/interproscan_output/"
echo -e "${WHITE}  SignalP:         ${NC}../Species/${GENUS}/${SPECIES_NAME}/signalp_output/"
echo -e "${WHITE}  dbCAN:           ${NC}../Species/${GENUS}/${SPECIES_NAME}/dbcan_output/"
echo -e "${WHITE}  Parsnp:          ${NC}../Species/${GENUS}/${SPECIES_NAME}/parsnp_output/"
echo -e "${WHITE}  BiG-SCAPE:       ${NC}../Species/${GENUS}/${SPECIES_NAME}/bigscape_output/"
echo -e "${WHITE}  Gubbins:         ${NC}../Species/${GENUS}/${SPECIES_NAME}/gubbins_output/"
echo -e "${WHITE}  IQ-TREE (SNV):   ${NC}../Species/${GENUS}/${SPECIES_NAME}/iqtree_output_snv/"
echo ""

# =============================================================================
# Optional: Cross-species OrthoFinder (all species combined)
# =============================================================================
# Not part of the numbered per-species pipeline above. Builds a combined
# OrthoFinder run across all annotated species (prompts for genus/species list).
echo -e "${YELLOW}Run OrthoFinder across all species combined? (y/n):${NC}"
read RUN_ORTHO_ALL
if [[ "$RUN_ORTHO_ALL" == "y" || "$RUN_ORTHO_ALL" == "Y" ]]; then
    ORTHO_ALL_SCRIPT="${SCRIPT_DIR}/orthofinder_on_all_species.sh"
    if [[ -f "$ORTHO_ALL_SCRIPT" ]]; then
        log_message "Running cross-species OrthoFinder (all species combined)..."
        bash "$ORTHO_ALL_SCRIPT"
    else
        echo -e "${RED}OrthoFinder (all species) script not found at: $ORTHO_ALL_SCRIPT${NC}"
    fi
fi
echo ""

# =============================================================================
# Optional: Run analysis notebooks
# =============================================================================
echo -e "${YELLOW}Run analysis notebooks? (y/n):${NC}"
read RUN_NB
if [[ "$RUN_NB" == "y" || "$RUN_NB" == "Y" ]]; then
    ANALYSIS_SCRIPT="${SCRIPT_DIR}/run_analysis_notebooks.sh"
    if [[ -f "$ANALYSIS_SCRIPT" ]]; then
        log_message "Running analysis notebooks..."
        bash "$ANALYSIS_SCRIPT"
    else
        echo -e "${RED}Analysis script not found at: $ANALYSIS_SCRIPT${NC}"
    fi
fi
