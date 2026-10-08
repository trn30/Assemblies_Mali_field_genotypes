# Improved long-read sequencing technology reveals full genome architecture of Malian Plasmodium falciparum field isolates

Code for PacBio HiFi genome assembly, correction, annotation and quality assessment of eight clonal Plasmodium falciparum field isolates from Mali, including coverage experiments and downstream analyses of frameshifts and multigene-family motifs.

## About this repository

This repository contains the analysis code associated with the manuscript **“Improved long-read sequencing technology reveals full genome architecture of Malian Plasmodium falciparum field isolates”**.

The workflow was developed for PacBio HiFi sequencing of eight clonal *Plasmodium falciparum* field isolates from Mali and covers de novo genome assembly, ILRA processing, read alignment and polishing, coverage experiments, and downstream analyses of genome structure, frameshifts and multigene-family motifs.

## Abstract

Understanding genome plasticity in malaria parasites and the diversity of their repetitive gene families requires high-quality reference genomes derived from field isolates. However, generating such genomes in Plasmodium falciparum remains challenging owing to the parasite's extreme AT content, extensive repetitive regions, and the frequent occurrence of multiclonal infections in endemic settings. As a result, most available reference genomes are derived from laboratory-adapted strains that incompletely represent natural parasite diversity.

We developed and optimized a de novo assembly strategy based on PacBio HiFi long-read sequencing and applied it to eight field-derived clonal P. falciparum lines from Mali. HiFi sequencing enabled highly accurate genome reconstruction, yielding assemblies with an average of 18 contigs, 301 pseudogenes, and a frameshift rate of 0.04% following Liftoff annotation. Assemblies showed high structural integrity and near-complete chromosomes, including improved recovery of subtelomeric regions and multigene families, without requiring short-read correction. Despite overall conservation of chromosome organization, substantial variation was observed in chromosome-end architecture, subtelomeric structure, and telomeric repeat content among isolates.

Our approach simplifies the generation of high-quality P. falciparum genomes and reduces sequencing requirements by up to fourfold. The resulting assemblies enabled characterization of natural variation in telomere plasticity and multicopy gene family diversity in field isolates. These findings establish PacBio HiFi sequencing as an effective and scalable strategy for producing highly contiguous and accurate de novo assemblies of P. falciparum, supporting the development of geographically diverse, field-relevant malaria parasite reference genomes.

## Workflow

The main analysis is organized as:

```text
PacBio BAM
   │
   ├── bam2fastq
   └── FastQC
        │
        ▼
   Canu assembly
        │
        ▼
      ILRA
        │
        ▼
 BWA-MEM alignment
        │
        ▼
     Pilon
```

Additional analyses are kept as separate scripts:

- Illumina-assisted ILRA correction for BG1701 and KK2403.
- Coverage/downsampling experiments.
- The 3D7 control experiment.

## Repository structure

```text
assemblies_Mali_isolates/
├── README.md
├── CITATION.cff
├── .gitignore
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
├── scripts/
│   ├── 00_check_dependencies.sh
│   ├── 01_prepare_reads.sh
│   ├── 02_canu_main.sh
│   ├── 03_ilra_main.sh
│   ├── 04_pilon_main.sh
│   ├── 05_illumina_correction.sh
│   ├── 06_coverage_experiment.sh
│   ├── 07_3d7_control.sh
│   ├── count_frameshifts.R
│   ├── count_motifs.sh
│   └── run_main.sh
└── docs/
    ├── annotation_and_qc.md

```

## Samples

The repository uses the biological isolate names throughout the analysis:

The repository uses the biological isolate names throughout the analysis.
Original sequencing barcodes are retained only inside raw FASTQ filenames
where needed to locate the input files.

## Configuration

No personal or institution-specific filesystem path is stored in the
repository.

Paths are supplied at runtime through environment variables, with generic
relative defaults in `config/config.sh`.

For example:

```bash
export DATA_ROOT=/path/to/project_data
export REFERENCE_FASTA=/path/to/PlasmoDB63-3D7.fasta
export REFERENCE_GFF=/path/to/PlasmoDB63-3D7.gff
export ILRA_SH=ILRA.sh
export ILRA_PATH=/path/to/ILRA/path_to_source
export PILON_JAR=pilon.jar
```

Pilon is configured by default as:

```bash
export PILON_JAR="pilon.jar"
```

The exact raw FASTQ names and Illumina prefixes are kept in
`config/samples.tsv`, while analysis-specific annotation inputs are listed in
`config/frameshift_samples.tsv`.

## How to run the analysis

The main workflow is divided into computational stages. Before running it, review `config/config.sh` and `config/samples.tsv` and make sure the required software is available.

```bash
bash scripts/00_check_dependencies.sh
bash scripts/01_prepare_reads.sh
bash scripts/02_canu_main.sh
bash scripts/03_ilra_main.sh
bash scripts/04_pilon_main.sh
```

Alternatively:

```bash
bash scripts/run_main.sh
```

Additional analyses are run separately:

```bash
bash scripts/05_illumina_correction.sh
bash scripts/06_coverage_experiment.sh
bash scripts/07_3d7_control.sh
```

Frameshift analysis:

```bash
Rscript scripts/count_frameshifts.R \
  --sample-table config/frameshift_samples.tsv \
  --reference-gff annotation/reference/PlasmoDB-3D7.gff \
  --output-dir results/frameshifts
```

Motif analysis for the final assemblies:

```bash
bash scripts/run_motif_counts.sh
```

## Main analysis

Run:

```bash
bash scripts/00_check_dependencies.sh
bash scripts/01_prepare_reads.sh
bash scripts/02_canu_main.sh
bash scripts/03_ilra_main.sh
bash scripts/04_pilon_main.sh
```

Or:

```bash
bash scripts/run_main.sh
```

These jobs are intentionally separate because assembly and polishing can be
resource-intensive on an HPC cluster.

## Illumina-assisted correction

The original analysis contains Illumina-assisted ILRA corrections for BG1701
and KK2403. Run:

```bash
bash scripts/05_illumina_correction.sh
```

Verify the Illumina prefixes in `config/samples.tsv` before running this step.

## Coverage experiment

The coverage analysis evaluates 5X, 10X, 20X, 25X, 50X and 100X.

The sampling fractions and random seeds retained from the original analysis
are in `config/downsampling.tsv`.

Run:

```bash
bash scripts/06_coverage_experiment.sh
```

## 3D7 control

The original analysis also contains a separate 3D7 control using 25X, 50X,
100X, 200X and 450X coverage, with Illumina-assisted correction for 200X and
450X.

Run:

```bash
bash scripts/07_3d7_control.sh
```

## Annotation, frameshifts and motif analysis

The downstream analyses are implemented in:

```text
scripts/count_frameshifts.R
scripts/count_motifs.sh
```

The frameshift analysis uses Companion one-to-one ortholog tables and GFF3
annotations. The motif analysis uses the motif FASTA files in
`data/motifs/`.

The original Artemis/ACT manual inspection is still documented separately in
`docs/annotation_and_qc.md`, because it is interactive rather than a fully
automated pipeline.


