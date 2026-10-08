# Downstream annotation and motif analyses

## 1. Frameshift analysis

`count_frameshifts.R` replaces the original six copy/pasted sample blocks with
one parameterized script. It retains the original analysis logic:

1. Read the one-to-one Companion ortholog table.
2. Count `CDS` features for each isolate gene and its PF3D7 ortholog.
3. Flag genes where the CDS counts differ.
4. For isolate genes with zero CDS features, count `pseudogenic_exon` features
   and compare those counts with the PF3D7 CDS count.
5. Write per-isolate tables and one summary table.

The original script performed this calculation separately for each field isolate; the repository now uses the biological isolate names consistently.

### Required input organization

The repository expects the files listed in
`config/frameshift_samples.tsv`. These are not included because the GFF and
Companion ortholog tables are large project data.

Example:

```text
annotation/
├── BG0801/
│   ├── orthologs.tsv
│   └── BG0801.gff3
├── BG1701/
│   ├── orthologs.tsv
│   └── BG1701.gff3
...
```

The ortholog tables should have two tab-separated columns:

```text
isolate_gene    PF3D7_gene
```

The 3D7 GFF is provided separately to the script with `--reference-gff`.

Run:

```bash
Rscript scripts/count_frameshifts.R \
  --sample-table config/frameshift_samples.tsv \
  --reference-gff annotation/reference/PlasmoDB-3D7.gff \
  --output-dir results/frameshifts
```

### Outputs

For every isolate:

```text
<ISOLATE>_frameshift_all.tsv
<ISOLATE>_frameshift_non_pseudogene.tsv
<ISOLATE>_frameshift_pseudogene.tsv
```

and:

```text
frameshift_summary.tsv
```

## 2. Motif analysis

`count_motifs.sh` replaces both the original motif-counting script and the
manual test script for individual assemblies.

The three motif resources are stored in:

```text
data/motifs/
├── VAR.patter.LARS.fasta
├── Rifin.motif.fasta
└── Stevor.motif.fasta
```

The original thresholds are retained:

- LARS: `pident > 90` and alignment length `> 10`.
- Stevor: `pident > 90`.
- Rifin: all `tblastn` hits are retained; motif 2 hits are counted separately.

The BLAST database is built in a temporary directory and is never written into
the genome directory.

Run:

```bash
bash scripts/count_motifs.sh \
  --genome path/to/assembly.fasta \
  --output results/motifs/BG0801.tsv
```

For an assembly outside the repository, set `DATA_ROOT`:

```bash
DATA_ROOT=/path/to/data \
bash scripts/count_motifs.sh \
  --genome assemblies/BG0801.fasta \
  --output results/motifs/BG0801.tsv
```

Only the runtime location is supplied this way; no personal server path is
stored in the repository.

## 3. Why `script_count_VAR_genes.sh` is not included

The supplied `script_count_VAR_genes.sh` is a manual test script containing
four hard-coded assemblies and repeating the LARS search. Its functionality is
covered by `count_motifs.sh`, which accepts any assembly as input, so the
manual test script is not needed in the publication repository.

## 4. Files to verify before publication

Check that:

- `config/frameshift_samples.tsv` points to the final Companion ortholog
  tables and isolate GFF files.
- The reference GFF version is the same version used for each sample analysis.
- The Illumina/annotation paths in the main configuration are correct.
- The motif FASTA files are the exact motif definitions used for the paper.


## Software required for these analyses

The frameshift script uses base R only; it does not require `ape`, `dplyr` or
`stringr`. The original script loaded those packages but did not use them.

The motif script requires BLAST+ (`makeblastdb` and `tblastn`). No BLAST database
files are stored in GitHub; they are generated temporarily for each assembly.

## Where to put the annotation inputs

The table `config/frameshift_samples.tsv` is deliberately written with
repository-relative paths. Place the Companion ortholog tables and isolate GFF3
files under the directories shown there, or edit the table to match the final
data layout.

For example:

```text
annotation/
├── reference/
│   └── PlasmoDB-3D7.gff
├── BG0801/
│   ├── orthologs.tsv
│   └── BG0801.gff3
├── BG1905/
│   ├── orthologs.tsv
│   └── BG1905.gff3
├── KK2504/
│   ├── orthologs.tsv
│   └── KK2504.gff3
├── KK2202/
│   ├── orthologs.tsv
│   └── KK2202.gff3
├── BG1701/
│   ├── orthologs.tsv
│   └── BG1701.gff3
└── KK2403/
    ├── orthologs.tsv
    └── KK2403.gff3
```

The files do not have to be committed to GitHub if they are large or are
already deposited elsewhere; the configuration only needs stable paths or
accession/documentation sufficient to retrieve them.


## Runtime paths

The repository contains no machine-specific absolute paths. Set `DATA_ROOT`,
`REFERENCE_FASTA`, `REFERENCE_GFF`, `ILRA_SH`, `ILRA_PATH` and `PILON_JAR` in
the shell environment when running the analysis.

The default Pilon value is `pilon.jar`.
