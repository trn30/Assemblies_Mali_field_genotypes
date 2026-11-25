# Improved long-read sequencing technology reveals full genome architecture of Malian *Plasmodium falciparum* field isolates

Code for PacBio HiFi genome assembly, correction, annotation and quality assessment of eight clonal *Plasmodium falciparum* field isolates from Mali, including coverage experiments and downstream analyses of frameshifts, telomeres and multigene-family motifs.

## About this repository

This repository contains the analysis code associated with the manuscript **"Improved long-read sequencing technology reveals full genome architecture of Malian *Plasmodium falciparum* field isolates"**.

The workflow was developed for PacBio HiFi sequencing of eight field-derived clonal *Plasmodium falciparum* isolates from Mali. It includes de novo genome assembly using Canu, PacBio-only ILRA processing, optional Pilon polishing, Illumina-assisted correction, coverage experiments, and downstream analyses of genome structure, frameshifts, telomeric repeats and multigene-family motifs.

## Abstract

Understanding genome plasticity in malaria parasites and the diversity of their repetitive gene families requires high-quality reference genomes derived from parasites recently isolated from patients. However, generating such genomes in Plasmodium falciparum remains challenging owing to the parasite's extreme AT content, extensive repetitive regions, and the frequent occurrence of multiclonal infections in endemic settings. As a result, most available reference genomes are derived from laboratory-adapted strains, that have experienced prolonged in vitro maintenance, and that incompletely represent natural parasite diversity.

We developed and optimized a de novo assembly strategy based on PacBio HiFi long-read sequencing and applied it to eight recently-isolated and cloned P. falciparum parasites from Mali. HiFi sequencing enabled highly accurate genome reconstruction, yielding assemblies with a  mean of 18 contigs, 301 pseudogenes, and a frameshift rate of 0.04% following Liftoff annotation. Assemblies showed high structural integrity and near-complete chromosomes, including improved recovery of subtelomeric regions and multigene families, without requiring short-read correction. Despite overall conservation of chromosome organization, substantial variation was observed in chromosome-end architecture, subtelomeric structure, and telomeric repeat content among parasite genomes. 

Our approach simplifies the generation of high-quality P. falciparum genomes and reduces sequencing requirements fourfold. The resulting assemblies enabled characterization of natural variation in telomere plasticity and multicopy gene family diversity in field isolates. These findings establish PacBio HiFi sequencing as an effective and scalable strategy for producing highly contiguous and accurate de novo assemblies of P. falciparum, supporting the development of geographically-diverse, field-relevant malaria parasite reference genomes.

## Workflow overview

The repository contains a main PacBio HiFi assembly workflow and independent downstream analysis branches.

### Main PacBio HiFi workflow

```text
PacBio HiFi reads (FASTQ)
          |
          v
   BAM-to-FASTQ (if needed)
          |
          v
       FastQC
          |
          v
     Canu assembly
          |
          v
     ILRA processing
      (PacBio only)
          |
          v
    ILRA assemblies
          |
          +----------------------+
          |                      |
          v                      v
 Optional Pilon             Downstream
 polishing                  analyses
 (04_pilon_main.sh)         (telomeres,
                            frameshifts,
                            motifs)
```

### Illumina-assisted correction

Illumina-assisted correction is implemented as a separate workflow:

```text
Canu assembly
     |
     v
ILRA + Illumina reads
     |
     v
Corrected assemblies
(x3 and x5 output directories)
```

This branch is configured for isolates with Illumina data available in `config/samples.tsv`.

### Additional analyses

The repository also includes:

- Coverage/downsampling experiments.
- A separate 3D7 control experiment.
- Frameshift analysis.
- Telomeric repeat identification and visualization.
- LARS, Rifin and Stevor motif detection.
- Known resistance mutation identification.

## Repository structure

```text
Assemblies_Mali_field_genotypes/
├── LICENSE
├── README.md
├── config/
│   ├── config.sh
│   ├── samples.tsv
│   ├── downsampling.tsv
│   ├── frameshift_samples.tsv
│   └── canu_hifi.spec
├── data/
│   └── motifs/
│       ├── VAR.patter.LARS.fasta
│       ├── Rifin.motif.fasta
│       └── Stevor.motif.fasta
├── docs/
│   └── annotation_and_qc.md
└── scripts/
    ├── 00_check_dependencies.sh
    ├── 01_prepare_reads.sh
    ├── 02_canu_main.sh
    ├── 03_ilra_main.sh
    ├── 04_pilon_main.sh
    ├── 05_illumina_correction.sh
    ├── 06_coverage_experiment.sh
    ├── 07_3d7_control.sh
    ├── 08_resistant_mutation_identification.R
    ├── count_frameshifts.R
    ├── count_motifs.sh
    ├── run_motif_counts.sh
    ├── telomere_analysis.sh
    ├── plot_telomeres_repetitions.R
    └── run_main.sh
```

## Samples

The repository uses biological isolate names throughout the analysis. Original sequencing barcodes are retained in raw FASTQ filenames where needed to locate input files.

The eight isolates included in `config/samples.tsv` are:

- BG1905
- BG0801
- KK2504
- KK2202
- BG1701
- KK2403
- KK2401
- KK2502

Illumina read prefixes are configured for BG1701 and KK2403. The remaining isolates are processed using PacBio data only in the main assembly workflow.


## Configuration

No personal or institution-specific filesystem paths are stored in the repository.

Paths and computational parameters are supplied through environment variables, with generic defaults in `config/config.sh`.

### Main configuration variables

| Variable | Default or purpose |
|---|---|
| `DATA_ROOT` | Root directory containing sequencing data and reference files |
| `WORK_ROOT` | Root directory for generated results |
| `REFERENCE_FASTA` | 3D7 reference genome FASTA |
| `REFERENCE_GFF` | 3D7 reference genome GFF |
| `THREADS` | Number of computational threads |
| `CANU_GENOME_SIZE` | Estimated genome size |
| `CANU_READ_COVERAGE` | Canu read coverage parameter |
| `CANU_USE_GRID` | Whether Canu uses a grid |
| `ASSEMBLY_ENV` | Conda environment for assembly |
| `ILRA_ENV` | Conda environment for ILRA |
| `ILRA_SH` | ILRA executable |
| `ILRA_PATH` | Optional ILRA initialization script |
| `PILON_JAR` | Pilon JAR file |
| `PILON_MEMORY` | Java memory allocation for Pilon |

The default data and results directories are:

```text
DATA_ROOT = ./data/
WORK_ROOT = ./results/
```

### Example configuration

```bash
export DATA_ROOT=/path/to/project_data
export WORK_ROOT=/path/to/project_results

export REFERENCE_FASTA=/path/to/PlasmoDB63-3D7.fasta
export REFERENCE_GFF=/path/to/PlasmoDB63-3D7.gff

export ILRA_SH=ILRA.sh
export ILRA_PATH=/path/to/ILRA/path_to_source

export PILON_JAR=/path/to/pilon/pilon.jar
export PILON_MEMORY=100G

export THREADS=30
```

The exact raw FASTQ filenames and Illumina prefixes are specified in `config/samples.tsv`.

Analysis-specific annotation inputs are specified in `config/frameshift_samples.tsv`.

Large sequencing files, BAM files, annotations and intermediate results should normally be stored outside GitHub.

## Software requirements

## Dependency checking

The script `scripts/00_check_dependencies.sh` checks the software and input files required for each workflow.

Run the default check:

```bash
bash scripts/00_check_dependencies.sh
```

Check a specific stage:

```bash
bash scripts/00_check_dependencies.sh --stage STAGE
```

Available stages:

| Stage | Analysis |
|---|---|
| `main` | PacBio HiFi assembly |
| `pilon` | Pilon polishing |
| `illumina` | Illumina-assisted correction |
| `coverage` | Coverage experiments |
| `3d7` | 3D7 control |
| `telomeres` | Telomere analysis |
| `motifs` | Motif analysis |
| `frameshifts` | Frameshift analysis |
| `resistance` | Resistance mutations |
| `all` | All dependencies |

The script reports missing requirements and exits with an error if required dependencies or files are unavailable.

### Downstream analysis

Additional software is required depending on the analysis:

| Analysis | Software |
|---|---|
| Frameshift analysis | R |
| Telomere analysis | TIDK, samtools, R, ggplot2 |
| Motif analysis | BLAST+ (`makeblastdb`, `tblastn`) |
| Resistance mutation identification | R, Biostrings, argparse, stringr, dplyr, tidyr |

The dependency-checking script reports the availability of the main command-line tools. Additional software for downstream analyses should be checked separately.

## How to run the main assembly workflow

Before running the analysis:

1. Review `config/config.sh`.
2. Check all input paths in `config/samples.tsv`.
3. Verify that the reference FASTA and GFF files are available.
4. Configure the required Conda environments and software.
5. Check that the PacBio FASTQ files exist under `DATA_ROOT`.

### Main workflow

Run the complete default workflow:

```bash
bash scripts/run_main.sh
```

Alternatively, execute each stage independently:

```bash
bash scripts/00_check_dependencies.sh

bash scripts/01_prepare_reads.sh

bash scripts/02_canu_main.sh

bash scripts/03_ilra_main.sh
```

The main stages are described below.

### 00. Dependency checking

Script:

`scripts/00_check_dependencies.sh`

Checks the availability of the main command-line dependencies, ILRA and Pilon configuration.

Run:

```bash
bash scripts/00_check_dependencies.sh
```

Review the output before proceeding, especially any missing dependencies or configuration warnings.

### 01. Read preparation

Script:

`scripts/01_prepare_reads.sh`

This script:

- Searches for BAM files in the configured shipment RAW directories.
- Converts available BAM files to FASTQ using `bam2fastq`.
- Runs FastQC on FASTQ files.
- Stores FastQC results in the corresponding RAW directory.

Run:

```bash
bash scripts/01_prepare_reads.sh
```

The main Canu workflow subsequently uses the FASTQ paths specified in `config/samples.tsv`.

### 02. Canu HiFi assembly

Script:

`scripts/02_canu_main.sh`

Performs de novo genome assembly using PacBio HiFi reads and Canu.

The assembly parameters are configured in:

```text
config/canu_hifi.spec
config/config.sh
```

Run:

```bash
bash scripts/02_canu_main.sh
```

Expected output:

```text
results/canu/
├── BG1905/
│   └── BG1905.contigs.fasta
├── BG0801/
│   └── BG0801.contigs.fasta
└── ...
```

The exact results directory can be changed using `WORK_ROOT`.

### 03. PacBio-only ILRA processing

Script:

`scripts/03_ilra_main.sh`

Processes Canu assemblies using ILRA.

This branch uses PacBio long-read information without Illumina-assisted correction.

The current configuration includes:

- `-C no`
- `-L pb`
- `-p no`

Run:

```bash
bash scripts/03_ilra_main.sh
```

Expected output:

```text
results/ilra/
├── BG1905/
├── BG0801/
├── KK2504/
└── ...
```

The resulting ILRA FASTA assemblies are also used by the downstream telomere analysis.

### 04. Optional Pilon polishing

Script:

`scripts/04_pilon_main.sh`

Performs an optional Pilon polishing run using the raw long-read FASTQ files and the corresponding ILRA assembly.

Run:

```bash
bash scripts/04_pilon_main.sh
```

The script aligns reads with BWA-MEM and passes the resulting BAM to Pilon.

It performs one Pilon invocation per isolate when the expected output FASTA does not already exist.

Expected output:

```text
results/pilon/
└── <ISOLATE>/
    ├── alignment/
    │   └── <ISOLATE>.sorted.bam
    ├── <ISOLATE>_pilon.fasta
    ├── <ISOLATE>_pilon.changes
    └── <ISOLATE>_pilon.vcf
```

The Pilon step is kept separate from the main HiFi assembly and ILRA workflow.

### 05. Illumina-assisted correction

Script:

`scripts/05_illumina_correction.sh`

Performs separate ILRA runs using Illumina reads for isolates with Illumina data specified in `config/samples.tsv`.

The script:

- Skips samples without an Illumina prefix.
- Uses the original Canu assembly as input.
- Supplies the corresponding Illumina read prefix to ILRA.
- Runs two configurations, writing results into `x3` and `x5` output directories.

Run:

```bash
bash scripts/05_illumina_correction.sh
```

Expected output:

```text
results/illumina_correction/
├── BG1701/
│   ├── x3/
│   └── x5/
└── KK2403/
    ├── x3/
    └── x5/
```

The `x5` run explicitly uses five iterations through the ILRA `-i 5` option. The `x3` directory identifies the first configured run; its iteration setting follows the ILRA default in the current script.

Verify the Illumina prefixes and paired-end FASTQ files before running this step.

## Coverage experiment

Script:

`scripts/06_coverage_experiment.sh`

The coverage experiment evaluates the effect of reduced PacBio read coverage on Canu assemblies.

The sampling fractions and random seeds are stored in:

```text
config/downsampling.tsv
```

The downsampling table includes:

- 5X
- 10X
- 20X
- 25X

The script also runs Canu using the full read sets with target working coverage of:

- 50X
- 100X

The coverage experiment is configured for six isolates:

- BG1905
- BG0801
- KK2504
- KK2202
- BG1701
- KK2403

Run:

```bash
bash scripts/06_coverage_experiment.sh
```

Expected output:

```text
results/coverage_experiment/
├── reads/
└── canu/
```

## 3D7 control experiment

Script:

`scripts/07_3d7_control.sh`

A separate control experiment is implemented using the 3D7 reference strain.

The PacBio FASTQ input and Illumina prefix are configured through:

```bash
THREE_D7_PACBIO_FASTQ
THREE_D7_ILLUMINA_PREFIX
```

The experiment evaluates:

- 25X
- 50X
- 100X
- 200X
- 450X

Illumina-assisted ILRA correction is performed for the 200X and 450X assemblies.

Run:

```bash
bash scripts/07_3d7_control.sh
```

Expected output:

```text
results/3D7_control/
├── canu/
│   ├── 25X/
│   ├── 50X/
│   ├── 100X/
│   ├── 200X/
│   └── 450X/
└── ilra/
    ├── 200X/
    └── 450X/
```

## Downstream genomic analyses

The downstream analyses are implemented as independent scripts and should be run after the required assemblies and annotation files are available.

### 1. Frameshift analysis

Script:

`scripts/count_frameshifts.R`

The script compares the number of CDS features between each isolate and its corresponding 3D7 one-to-one ortholog.

The analysis:

1. Reads the Companion one-to-one ortholog tables.
2. Counts CDS features in the isolate GFF3 annotation.
3. Counts CDS features in the 3D7 reference GFF annotation.
4. Identifies differences in CDS feature counts.
5. Evaluates pseudogenic exon features when isolate genes have no CDS features.
6. Generates per-isolate tables and a summary table.

Configuration:

```text
config/frameshift_samples.tsv
```

The table contains three columns:

```text
isolate    orthologs    gff
```

The annotation files are external inputs and are not included in the repository.

Example:

```text
annotation/
├── reference/
│   └── PlasmoDB-3D7.gff
├── BG0801/
│   ├── orthologs.tsv
│   └── BG0801.gff3
├── BG1701/
│   ├── orthologs.tsv
│   └── BG1701.gff3
└── ...
```

The ortholog tables should contain two tab-separated columns:

```text
isolate_gene    PF3D7_gene
```

Run from the repository root:

```bash
Rscript scripts/count_frameshifts.R \
    --sample-table config/frameshift_samples.tsv \
    --reference-gff annotation/reference/PlasmoDB-3D7.gff \
    --output-dir results/frameshifts
```

The script resolves input paths relative to `DATA_ROOT` when that variable is set; otherwise, it resolves them relative to the repository root. The output directory is relative to the repository root.

Expected output:

```text
results/frameshifts/
├── <ISOLATE>_frameshift_all.tsv
├── <ISOLATE>_frameshift_non_pseudogene.tsv
├── <ISOLATE>_frameshift_pseudogene.tsv
└── frameshift_summary.tsv
```

Additional details are provided in:

`docs/annotation_and_qc.md`

### 2. Telomere analysis

Scripts:

- `scripts/telomere_analysis.sh`
- `scripts/plot_telomeres_repetitions.R`

The telomere analysis identifies and visualizes telomeric repeat sequences in the 3D7 reference and the ILRA assemblies.

The workflow uses TIDK to search for the telomeric repeat:

```text
TTTAGGG
```

The analysis uses a 10,000 bp window.

The Bash script:

1. Checks the availability of TIDK, samtools and Rscript.
2. Processes the 3D7 reference FASTA.
3. Searches for ILRA assemblies for each isolate.
4. Generates FASTA indexes using `samtools faidx`.
5. Runs TIDK on the reference and isolate assemblies.
6. Creates an input table containing the FASTA and TIDK output paths.
7. Calls the R plotting script.

The R script generates chromosome-level visualizations for chromosomes 1-14, using the chromosome information in the FASTA index files and the TIDK repeat counts.

Windows with more than ten total telomeric repeats are considered for visualization.

The categories used in the plots are:

- 10-30 repeats
- 30-100 repeats
- More than 100 repeats

Run:

```bash
bash scripts/telomere_analysis.sh
```

Required software:

- TIDK
- samtools
- Rscript
- R package ggplot2

Required input:

- A valid 3D7 reference FASTA specified by `REFERENCE_FASTA`.
- Completed ILRA assemblies under `WORK_ROOT/ilra/`.
- One unambiguous `*ILRA.fasta` file per isolate.

The analysis stops if an expected isolate assembly is missing or if more than one matching ILRA FASTA is found for an isolate.

Expected output:

```text
results/telomeres/
├── telomere_inputs.tsv
├── tidk/
│   ├── 3D7_telomeric_repeat_windows.tsv
│   ├── BG1905_telomeric_repeat_windows.tsv
│   └── ...
└── figures/
    ├── telomeric_repeats.pdf
    └── telomeric_repeats.png
```

The output root is controlled by `WORK_ROOT`.

The chromosome-level ideogram includes chromosome names that can be identified as chromosomes 1-14 from the FASTA index. Non-chromosomal contigs are excluded from this visualization.

### 3. Multigene-family motif analysis

Scripts:

- `scripts/count_motifs.sh`
- `scripts/run_motif_counts.sh`

This analysis identifies LARS, Rifin and Stevor motifs in genome assemblies using BLAST+.

The motif reference sequences are stored in:

```text
data/motifs/
├── VAR.patter.LARS.fasta
├── Rifin.motif.fasta
└── Stevor.motif.fasta
```

The script uses:

- `makeblastdb`
- `tblastn`

The BLAST database is generated in a temporary directory for each assembly and is removed after the analysis unless `--keep-db` is specified.

The original motif-counting thresholds are retained:

**LARS**

- Percentage identity greater than 90%.
- Alignment length greater than 10 amino acids.

**Stevor**

- Percentage identity greater than 90%.

**Rifin**

- All BLAST hits are retained.
- Motif 2 hits are counted separately using the query identifier `rifin.motif2`.

#### Batch analysis

The batch script reads the isolate information from:

```text
config/samples.tsv
```

It searches for the corresponding Pilon assembly:

```text
WORK_ROOT/pilon/<ISOLATE>/<ISOLATE>_pilon.fasta
```

Samples without an available Pilon FASTA are skipped with a warning.

Run:

```bash
bash scripts/run_motif_counts.sh
```

Expected output:

```text
results/motifs/
├── BG1701.tsv
├── KK2403.tsv
└── <ISOLATE>.tsv
```

The output directory is controlled by `WORK_ROOT`.

#### Individual assembly analysis

Individual assemblies can also be analysed independently:

```bash
bash scripts/count_motifs.sh \
    --genome /path/to/assembly.fasta \
    --output /path/to/results.tsv \
    --motif-dir data/motifs \
    --threads 8
```

The script accepts both absolute and relative paths for the genome FASTA.

Relative input paths are resolved using the current directory, `DATA_ROOT` or the repository root.

Relative output paths are written under the repository root; absolute output paths are preserved.

Additional implementation details are provided in:

`docs/annotation_and_qc.md`

### 4. Known resistance mutation identification

Script:

`scripts/08_resistant_mutation_identification.R`

This script compares reference and isolate protein FASTA sequences and identifies amino acid differences.

It also compares observed amino acids with a table of known resistance-associated mutations.

The script additionally detects specific microsatellite motifs in the:

```text
PF3D7_1303500.1-p1
```

protein sequence.

Required inputs:

- Reference protein FASTA.
- Assembly protein FASTA.
- Known mutation CSV file.

The mutation CSV should include the gene and mutation information expected by the script (`Gene` and `mutation` columns).

Run:

```bash
Rscript scripts/08_resistant_mutation_identification.R \
    -r reference_proteins.fasta \
    -a assembly_proteins.fasta \
    -m known_mutations.csv \
    -o output_prefix
```

The output prefix is used to generate:

```text
output_prefix.all_mutations.tsv
output_prefix.known_mutations.tsv
output_prefix.microsatellite.tsv
```

The complete mutation and microsatellite outputs depend on the presence of the corresponding sequences and detected differences.

Required R packages:

- Biostrings
- argparse
- stringr
- dplyr
- tidyr

## Annotation and manual genome inspection

Manual genome inspection using Artemis and ACT is described separately in:

```text
docs/annotation_and_qc.md
```

This documentation covers the annotation inputs, quality assessment and manual inspection procedures that are not fully automated by the main pipeline.


## Data availability

To adhere to code repository best practices, no raw data, intermediate BAM files or generated result tables are hosted in this repository.

- Raw reads and assembly files associated with this study are deposited in ENA under accession **PRJEB98187**.
- Reference genomes were obtained from [PlasmoDB](https://plasmodb.org/).

Large raw sequencing files and annotation resources should be downloaded separately and placed under the paths specified in the configuration files.

## Citation

If you use the code, pipelines or methodological approaches provided in this repository for your research, please cite the associated publication:

**Díaz-Terenti, B., Diakite, M., Randford-Cartwright, L. C., Pradham, S., Martos-Contreras, M. C., Otto, T. D., & Gómez-Díaz, E.** *Improved long-read sequencing technology reveals full genome architecture of Malian Plasmodium falciparum field isolates*. (2026).

**BibTeX:**

```bibtex
@article{diazterenti_malian_pfalciparum,
  title = {Improved long-read sequencing technology reveals full genome architecture of Malian Plasmodium falciparum field isolates},
  author = {D{\'i}az-Terenti, B. and Diakite, M. and Randford-Cartwright, L. C. and Pradham, S. and Otto, T. D. and G{\'o}mez-D{\'i}az, E.},
  year = {2026},
  doi = {PENDING}
}
```


