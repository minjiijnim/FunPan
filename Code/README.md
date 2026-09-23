# FunPan Pipeline Scripts

Reference for the 13-step shell pipeline in `Code/`. Run it end to end with
`run_complete_pipeline.sh`, or run any step on its own.

For an overview of the project, environment setup and quick start, see the
[main README](../README.md). For the downstream analysis notebooks, see
[Analysis/README.md](../Analysis/README.md).

---

## Table of Contents

1. [Tools & Scripts Documentation](#tools--scripts-documentation)
2. [Tool-Specific Parameters](#tool-specific-parameters)
3. [Additional Standalone Scripts](#additional-standalone-scripts)
4. [Log File Locations](#log-file-locations)

---

## Tools & Scripts Documentation

### 1. `run_complete_pipeline.sh` - Master Pipeline Controller

**Purpose:** Orchestrates all 13 pipeline steps with progress tracking and error handling.

**Key Features:**
- Interactive step selection (resume, fresh start, specific step)
- Real-time progress bar with colored output
- Automatic CPU detection (uses 75% of available cores)
- Checkpoint-based resumption
- Comprehensive logging to `.pipeline_logfiles/`

**Environment Variables:**
- `PIPELINE_THREADS`: Number of threads to use (default: auto-detected)

**Usage:**
```bash
bash run_complete_pipeline.sh
```

**Checkpoint Files:**
```
.pipeline_checkpoints/
├── step_1.completed
├── step_2.completed
├── ...
└── pipeline_inputs.txt  # Saved genus, species, lineage
```

---

### 2. `download_genome_and_BUSCO.sh` - Genome Acquisition & QC

**Purpose:** Download genomes from NCBI and assess quality with BUSCO.

**Parameters:**

| Parameter | Default | Description |
|-----------|---------|-------------|
| `PIPELINE_THREADS` | auto (total-4) | CPU threads for BUSCO |
| `BUSCO_THREADS` | `PIPELINE_THREADS` | BUSCO-specific threads |
| `BUSCO_MODE` | `auto` | Lineage detection: `auto`, `auto_euk`, `auto_prok` |
| `BUSCO_OFFLINE` | `0` | Set to `1` for offline mode |
| `BUSCO_LINEAGE` | auto-detected | Force specific lineage (e.g., `eurotiomycetes_odb10`) |

**Process Flow:**
1. **Download:** Uses NCBI Datasets CLI to fetch:
   - Genome assemblies (`.fna`)
   - Protein sequences (`.faa`)
   - GenBank files (`.gbk`)
   - GFF3 annotations (`.gff3`)

2. **RNA-seq Discovery:**
   - Extracts BioSample IDs from assembly metadata
   - Queries SRA for linked RNA-seq runs
   - Downloads FASTQ files from ENA
   - Organizes by strain: `rna/{strain_name}/`

3. **Duplicate Handling:**
   - Identifies GCA/GCF pairs (same genome, GenBank vs RefSeq)
   - Keeps GCF (RefSeq), removes GCA
   - Cleans up orphaned BUSCO results

4. **BUSCO Analysis:**
   - Auto-detects lineage from first genome
   - Reuses lineage for all subsequent genomes
   - Saves lineage to `pipeline_inputs.txt`
   - Cleans up redundant lineage results

5. **Metadata Generation:**
   - Merges BUSCO scores into assembly metadata
   - Extracts isolation source, location, collection date
   - Outputs `metadata.csv` and `busco_updated.jsonl`

**Output:**
```
Species/{genus}/{species}/
├── genome/                      # Genome FASTA files
├── protein/                     # Protein FASTA files
├── genbank/                     # GenBank format files
├── gff3/                        # GFF3 annotation files
├── rna/{strain_name}/           # RNA-seq FASTQ files
│   ├── runs.txt                # List of SRA run IDs
│   └── {SRR}.fastq.gz          # Downloaded reads
├── busco_output/{genome_name}/  # BUSCO results per genome
├── metadata.csv                 # Comprehensive metadata table
├── busco_updated.jsonl          # JSONL with BUSCO scores
└── biosamples.tsv               # Assembly-BioSample mapping
```

---

### 3. `ani_and_filter_genome_QC.sh` - ANI Species Check, Deduplication & QC Filtering

**Purpose:** Verify species identity, remove near-identical duplicates, and filter genomes by assembly quality metrics.

**Usage:**
```bash
bash ani_and_filter_genome_QC.sh Aspergillus/flavus [ANI_THRESHOLD] [MASH_THRESHOLD]
```

**Process Flow:**

1. **ANI Species Verification (fastANI):**
   - Compares all genomes against the GCF reference genome
   - Removes genomes below the ANI threshold (default: 95%)
   - Genomes too divergent to compute ANI are also removed

2. **Mash Deduplication:**
   - Sketches ANI-passed genomes with Mash (sketch size 10,000)
   - Computes all-vs-all Mash distances
   - Clusters near-identical genomes (default threshold: 0.0001 = 99.99% identity)
   - Keeps the highest-quality representative per cluster

3. **QC Filtering (same criteria as before):**

| Metric | Threshold | Rationale |
|--------|-----------|-----------|
| **BUSCO Complete** | > 95% | Ensures genome completeness |
| **Contig Count** | ≤ median | Prefers more contiguous assemblies |
| **Contig N50** | ≥ 75th percentile | Prefers longer contigs |
| **GCF Requirement** | ≥ 1 | Always includes at least one RefSeq genome |

**Parameters:**

| Parameter | Default | Description |
|-----------|---------|-------------|
| `ANI_THRESHOLD` | 95 | Minimum ANI% to GCF reference |
| `MASH_THRESHOLD` | 0.0001 | Maximum Mash distance for duplicate detection (~99.99% identity) |

**Dependencies:** `fastANI`, `mash` — auto-installed via conda/mamba if not found on PATH.

**Input:**
- `metadata.csv` from Step 1
- Genome files in `genome/`
- Protein files in `protein/`

**Output:**
```
filtered_genome/     # High-quality genome assemblies
filtered_protein/    # Corresponding protein sequences
qc_work/             # Intermediate files (ANI results, Mash distances, dedup lists)
```

---

### 4. `run_funannotate.sh` - Structural & Functional Annotation

**Purpose:** Predict genes and annotate genomes using ab initio and evidence-based methods.

**Process Flow:**

1. **FASTA Header Renaming:** Converts contig names to `scaffold_1`, `scaffold_2`, etc.
2. **GFF3 Validation:** Checks for contig mismatches, mRNA features, duplicate IDs
3. **Assembly Cleaning:** Removes short contigs (`funannotate clean`)
4. **Repeat Masking:** De novo RepeatModeler library + RepeatMasker soft-masking
5. **Protein Evidence Selection:** Best GCF proteome as reference + strain-specific
6. **RNA-seq Evidence:** Auto-detects RNA-seq files per strain
7. **Gene Prediction:** EvidenceModeler integration with ab initio predictors

**Output:**
```
funannotate_output/{genome}/
├── {genome}.clean.fa            # Cleaned assembly
├── {genome}.softmasked.fa       # Repeat-masked assembly
├── predict_results/
│   ├── {genome}.gbk            # GenBank annotation
│   ├── {genome}.gff3           # GFF3 annotation
│   ├── {genome}.proteins.fa    # Predicted proteins
│   └── {genome}.transcripts.fa # Predicted transcripts
└── logfiles/
```

---

### 5. `run_eggnog.sh` - Functional Annotation

**Purpose:** Assign functional annotations (GO terms, KEGG pathways, domains) via orthology.

**Input:** `*.proteins.fa` from Funannotate predict results

**Output:**
```
eggnog_output/
├── {genome}.emapper.annotations  # Tab-delimited functional annotations
├── {genome}.emapper.hits         # DIAMOND hits
└── {genome}.emapper.seed_orthologs
```

**Annotation Fields:** GO terms, KEGG KO numbers and pathways, COG functional categories

---

### 6. `run_orthofinder.sh` - Comparative Genomics

**Purpose:** Identify orthogroups and infer evolutionary relationships.

**Input:** `*.proteins.fa` from Funannotate predict results

**Options:** `-S diamond -M msa -A mafft`

**Output:**
```
orthofinder_output/Results_{date}/
├── Orthogroups/
│   ├── Orthogroups.tsv              # Gene-to-orthogroup mapping
│   ├── Orthogroups.GeneCount.tsv    # Counts per species
│   └── Orthogroups_SingleCopyOrthologues.txt
├── MultipleSequenceAlignments/      # MSA per orthogroup
├── Species_Tree/
│   └── SpeciesTree_rooted.txt       # Newick format
└── Gene_Trees/                      # Individual gene trees
```

---

### 7. `run_antismash.sh` - Secondary Metabolite Detection

**Purpose:** Identify and characterize biosynthetic gene clusters (BGCs).

**Isolated Environment:** Creates separate conda env at `.envs/antismash8` to avoid conflicts.

**antiSMASH Options:**
- `--taxon fungi`: Fungal-specific gene cluster detection
- `--genefinding-tool none`: Uses existing Funannotate annotations
- `--asf`, `--cb-knownclusters`, `--cb-general`, `--cb-subclusters`, `--pfam2go`

**Input:** GenBank files (`.gbk`) from Funannotate predict results

**Output:**
```
antismash_output/{genome}/
├── index.html                   # Interactive report
├── {genome}.gbk                 # Annotated GenBank with BGCs
├── *.region001.gbk              # Individual BGC GenBank files
└── knownclusterblast/           # MIBiG comparisons
```

---

### 8. `run_interproscan.sh` - Protein Domain Annotation

**Purpose:** Scan protein sequences against InterPro signature databases (Pfam, SMART, CDD, PANTHER, etc.).

**Features:**
- Auto-installs InterProScan if not found (downloads ~2.5 GB)
- Optimized tool set by default (excludes Hamap, ProSitePatterns, PRINTS for ~30-40% speedup)
- `TEST_MODE=1` to process only first file for testing

**Input:** `*.proteins.fa` from Funannotate predict results

**Output:**
```
interproscan_output/{genome}/
├── {genome}.tsv                 # Tab-separated annotations
├── {genome}.gff3                # GFF3 format annotations
└── {genome}.xml                 # XML format annotations
```

---

### 9. `run_signalp.sh` - Signal Peptide Prediction

**Purpose:** Predict signal peptides in protein sequences using SignalP 6.

**Features:**
- Auto-tunes parallel jobs and torch threads based on CPU cores
- Supports both SignalP 6.x and legacy versions
- `TEST_MODE=1` to process only first file for testing

**Input:** `*.proteins.fa` from Funannotate predict results

**Output:**
```
signalp_output/{genome}/
├── prediction_results.txt       # Signal peptide predictions
└── signalp.log                  # Run log
```

---

### 10. `run_dbcan.sh` - CAZyme Annotation

**Purpose:** Annotate carbohydrate-active enzymes (CAZymes) using the dbCAN database.

**Features:**
- Auto-downloads dbCAN database if not present
- Supports multiple detection methods (HMM, DIAMOND, dbCANsub)
- Skips completed samples automatically

**Input:** `*.proteins.fa` from Funannotate predict results

**Output:**
```
dbcan_output/{genome}/
└── overview.tsv                 # CAZyme annotations summary
```

---

### 11. `run_parsnp.sh` - Core Genome Alignment

**Purpose:** Perform rapid core genome alignment using Parsnp.

**Features:**
- Auto-selects best quality reference genome (prefers GCF/RefSeq)
- Requires at least 3 genomes
- Uses funannotate cleaned genomes (`*.clean.fa`)

**Input:** `*.clean.fa` from Funannotate output directories

**Output:**
```
parsnp_output/
├── parsnp.tree                  # Newick phylogenetic tree
├── parsnp.xmfa                  # Core genome alignment (XMFA)
└── parsnp.ggr                   # Gingr visualization file
```

---

### 12. `run_bigscape.sh` - BGC Network Analysis

**Purpose:** Compare BGCs across genomes and group into Gene Cluster Families (GCFs).

**Dependencies:** Auto-clones BiG-SCAPE from GitHub, downloads Pfam-A.hmm

**Process Flow:**
1. Collect `*.region*.gbk` from antiSMASH output
2. Rename with genome accession prefix into `antismash_gbk_unique/`
3. Run BiG-SCAPE clustering (GCF cutoff: 0.3)

**Output:**
```
bigscape_output/
├── index.html                   # Interactive network viewer
├── data_sqlite.db               # SQLite database of results
└── network_files/               # Cytoscape-compatible networks
```

---

### 13. `run_gubbins.sh` - Recombination Detection

**Purpose:** Detect and remove recombination regions from core genome alignment.

**Isolated Environment:** Creates separate conda env at `.envs/gubbins` with Gubbins, IQ-TREE, RAxML-NG, FastTree, and harvesttools.

**Features:**
- Converts Parsnp XMFA to FASTA via harvesttools
- Tries IQ-TREE first, falls back to RAxML-NG if it fails
- Requires at least 3 sequences

**Input:** `parsnp_output/parsnp.xmfa` from Step 10

**Output:**
```
gubbins_output/
├── gubbins.final_tree.tre                       # Recombination-corrected tree
├── gubbins.recombination_predictions.gff        # Detected recombination regions
├── gubbins.filtered_polymorphic_sites.fasta     # Filtered alignment
└── gubbins.node_labelled.final_tree.tre         # Node-labeled tree
```

---

### 14. `run_iqtree_gubbins.sh` - SNV-based Phylogeny

**Purpose:** Build phylogenetic tree from recombination-filtered SNVs using IQ-TREE 2.

**Features:**
- ModelFinder for optimal substitution model selection
- 1000 ultrafast bootstrap replicates with BNNI correction
- Requires IQ-TREE 2.x specifically

**Input:** `gubbins_output/gubbins.filtered_polymorphic_sites.fasta` from Step 12

**Output:**
```
iqtree_output_snv/
├── gubbins_tree.treefile        # Maximum likelihood tree
├── gubbins_tree.iqtree          # Full IQ-TREE report
└── gubbins_tree.log             # Run log
```

---

## Tool-Specific Parameters

**Funannotate:**
- Augustus training: automatic per species
- Protein evidence weight: best quality proteome + strain-specific
- GFF3 evidence weight: 10 (high confidence)
- Minimum protein length: 50 aa (Funannotate default)

**OrthoFinder:**
- Inflation parameter: 1.5 (MCL default)
- E-value threshold: 1e-3
- MSA tool: MAFFT
- Tree inference: FastTree (default) or IQ-TREE

**antiSMASH:**
- Detection strictness: relaxed (default)
- Minimum cluster size: follows antiSMASH defaults
- All detection modules enabled

**BiG-SCAPE:**
- Cutoff: 0.3 (balanced)
  - 0.1-0.2: very strict (many small GCFs)
  - 0.3-0.5: moderate (recommended)
  - 0.6-0.9: loose (few large GCFs)

**InterProScan:**
- Default applications: AntiFam, CDD, Coils, FunFam, Gene3D, MobiDBLite, NCBIfam, PANTHER, Pfam, PIRSF, PIRSR, ProSiteProfiles, SFLD, SMART, SUPERFAMILY
- Override with `INTERPRO_APPS="Pfam,SMART,CDD"`

**SignalP:**
- Organism: `euk` (eukaryote, default)
- Torch threads per job: 3 (default)
- Override parallel jobs with `SIGNALP_JOBS=30`

**dbCAN:**
- Methods: HMM, DIAMOND, dbCANsub (all by default)
- Override with `DBCAN_METHODS="hmm,diamond"`

---

## Additional Standalone Scripts

### `orthofinder_on_all_species.sh` - Cross-Species OrthoFinder

Runs OrthoFinder on protein files from multiple species simultaneously. Required before NB3 (convergence analysis).

```bash
bash orthofinder_on_all_species.sh
# Prompts for: genus, space-separated species list
# Output: Species/{genus}/all_combined/orthofinder_output/
```

---

## Log File Locations

**Pipeline Logs:**
```bash
/datadrive/Code/.pipeline_logfiles/pipeline_YYYYMMDD_HHMMSS.log
```

**Tool-Specific Logs:**
```bash
# Funannotate
Species/{genus}/{species}/funannotate_output/{genome}/logfiles/

# BUSCO
Species/{genus}/{species}/busco_output/{genome}/logs/

# Gubbins
Species/{genus}/{species}/gubbins_output/gubbins.log
```

