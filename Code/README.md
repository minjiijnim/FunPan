# FunPan Pipeline Scripts

Thirteen steps, one script each, run per species. Every path below is relative to
`Species/{genus}/{species}/`.

See also: [project overview](../README.md), [analysis notebooks](../Analysis/README.md).

## Running

```bash
bash run_complete_pipeline.sh
```

`run_complete_pipeline.sh` first offers a fresh start, a resume, or one step, then
asks for genus and species. It sets `PIPELINE_THREADS` to 75% of the cores and
exports it, writes a `.pipeline_checkpoints/step_N.completed` marker per step and
logs to `.pipeline_logfiles/pipeline_*.log`. Saved answers, including the BUSCO
lineage, go to `.pipeline_checkpoints/pipeline_inputs.txt`.

Any step also runs on its own. Step 1 prompts for genus and species on separate
lines, step 2 takes `genus/species` as an argument, and steps 3 to 13 prompt for
`genus/species`:

```bash
echo Aspergillus/flavus | bash run_funannotate.sh
bash ani_and_filter_genome_QC.sh Aspergillus/flavus
```

## Steps

| # | Script | Does | Input | Output |
|---|---|---|---|---|
| 1 | `download_genome_and_BUSCO.sh` | NCBI Datasets download, SRA/ENA RNA-seq fetch, GCA/GCF deduplication, BUSCO | genus + species | `genome/`, `protein/`, `genbank/`, `gff3/`, `rna/{strain}/`, `busco_output/`, `metadata.csv` |
| 2 | `ani_and_filter_genome_QC.sh` | fastANI species check, Mash deduplication, quality filter | `metadata.csv`, `genome/`, `protein/` | `filtered_genome/`, `filtered_protein/`, `qc_work/` |
| 3 | `run_funannotate.sh` | Contig renaming, GFF3 validation, RepeatModeler/RepeatMasker, EvidenceModeler gene prediction | `filtered_genome/`, `gff3/`, `rna/` | `funannotate_output/{genome}/predict_results/` |
| 4 | `run_eggnog.sh` | eggNOG-mapper: GO, KEGG KO and pathway, COG category | `*.proteins.fa` | `eggnog_output/{genome}.emapper.annotations` |
| 5 | `run_orthofinder.sh` | Orthogroups and species tree (`-S diamond -M msa -A mafft`) | `*.proteins.fa` | `orthofinder_output/Results_*/` |
| 6 | `run_antismash.sh` | Biosynthetic gene clusters, fungal taxon | `*.gbk` | `antismash_output/{genome}/` |
| 7 | `run_interproscan.sh` | InterPro signatures: Pfam, SMART, CDD, PANTHER and others | `*.proteins.fa` | `interproscan_output/{genome}/{genome}.tsv` |
| 8 | `run_signalp.sh` | SignalP 6 signal peptides | `*.proteins.fa` | `signalp_output/{genome}/prediction_results.txt` |
| 9 | `run_dbcan.sh` | dbCAN CAZymes | `*.proteins.fa` | `dbcan_output/{genome}/overview.tsv` |
| 10 | `run_parsnp.sh` | Core genome alignment, reference auto-selected | `*.clean.fa` | `parsnp_output/parsnp.{xmfa,tree,ggr}` |
| 11 | `run_bigscape.sh` | BiG-SCAPE gene cluster families, cutoff 0.3 | `*.region*.gbk` | `bigscape_output/` |
| 12 | `run_gubbins.sh` | Recombination detection on the core alignment | `parsnp_output/parsnp.xmfa` | `gubbins_output/gubbins.filtered_polymorphic_sites.fasta` |
| 13 | `run_iqtree_gubbins.sh` | IQ-TREE 2 phylogeny, ModelFinder, 1000 UFBoot | `gubbins_output/gubbins.filtered_polymorphic_sites.fasta` | `iqtree_output_snv/gubbins_tree.treefile` |

Steps 6 and 12 build their own conda environments under `.envs/` to keep
antiSMASH and Gubbins off the main environment. Steps 2, 7, 9 and 11 install
their tools or databases on first run if they are missing.

## Parameters

Set any of these in the environment before running a step.

| Step | Variable | Default | Effect |
|---|---|---|---|
| all | `PIPELINE_THREADS` | 75% of cores under the controller, cores minus 4 standalone | Thread count |
| 1 | `BUSCO_THREADS` | `PIPELINE_THREADS` | Threads for BUSCO only |
| 1 | `BUSCO_MODE` | `auto` | `auto`, `auto_euk` or `auto_prok` |
| 1 | `BUSCO_LINEAGE` | auto-detected | Forces a lineage, for example `eurotiomycetes_odb10` |
| 1 | `BUSCO_OFFLINE` | `0` | `1` uses local lineage files only |
| 7, 8 | `TEST_MODE` | `0` | `1` processes the first genome only |
| 7 | `INTERPRO_APPS` | 15 applications | Comma-separated application list |
| 8 | `SIGNALP_JOBS` | auto | Parallel jobs, 3 torch threads each |
| 9 | `DBCAN_METHODS` | `hmm,diamond,dbcansub` | Detection methods |

Step 2 also takes two optional arguments after the species:
`ANI_THRESHOLD` (default 95, minimum ANI percent to the RefSeq reference) and
`MASH_THRESHOLD` (default 0.0001, the distance below which two assemblies count
as duplicates).

Step 2 quality thresholds: BUSCO complete above 95%, contig count at or below the
median, N50 at or above the 75th percentile, at least one RefSeq assembly kept.

Fixed tool settings worth knowing: OrthoFinder inflation 1.5 and e-value 1e-3;
Funannotate minimum protein length 50 aa and GFF3 evidence weight 10; SignalP
organism `euk`; antiSMASH relaxed strictness with all detection modules on;
BiG-SCAPE cutoff 0.3, where 0.1 to 0.2 is strict and 0.6 to 0.9 is loose.

## Standalone scripts

| Script | Does |
|---|---|
| `orthofinder_on_all_species.sh` | OrthoFinder across several species at once, into `Species/{genus}/all_combined/orthofinder_output/`. NB3 needs it. Prompts for genus and a species list |
| `run_analysis_notebooks.sh` | Executes the six analysis notebooks. See [Analysis/README.md](../Analysis/README.md) |

## Logs

```
Code/.pipeline_logfiles/pipeline_YYYYMMDD_HHMMSS.log
Species/{genus}/{species}/funannotate_output/{genome}/logfiles/
Species/{genus}/{species}/busco_output/{genome}/logs/
Species/{genus}/{species}/gubbins_output/gubbins.log
```
