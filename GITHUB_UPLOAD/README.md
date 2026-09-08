# LIFE703 — Pan-gene Set Analysis of *Candida auris*

MSc Bioinformatics research project (LIFE703), University of Liverpool.

**Student:** Faiza Khan (201936907)
**Supervisor:** Prof. Andy Jones · **Co-supervisor:** Dr Evelina Basenko

This repository holds all code, scripts and documentation for the project. Large
data files and job outputs are **not** tracked (see `.gitignore`); they live on
the HPC cluster and are described in `data/README.md`.

---

## Project overview

*Candida auris* is an emerging, multidrug-resistant fungal pathogen that has
spread as several distinct global clades. The project builds and analyses the
**pan-gene set** of *C. auris*: producing consistent annotations across many
genomes, clustering genes, and classifying them as **core / shell / cloud** to
study genomic diversity, clade-specific genes and resistance-associated genes.

### Phases

1. **Reference comparison** — compare the B8441 protein annotations
   (FungiDB v68, GenBank v2 & v3) to establish a reliable reference set.
2. **Annotation** — transfer the GenBank v3 annotation onto target genomes with
   **Liftoff**; re-annotate genes that do not lift over with **Funannotate**
   (BRAKER3 / HELIXER as alternatives). Extract proteins using NCBI
   **Translation Table 12** (CTG→Ser, the *C. auris* code).
3. **Pan-gene construction** — cluster genes and build a presence–absence matrix
   with the **GET_PANGENES** pipeline, scaling to tens of *C. auris* genomes.
4. **Analysis & interpretation** — classify **core / shell / cloud** genes,
   compute gene statistics, add functional annotation with **InterProScan**, and
   produce visualisations (gene-per-genome, composition, heatmaps).

> Pipeline note: the pan-gene step uses **GET_PANGENES** (Contreras-Moreira
> *et al.*), as confirmed by the supervisor — not PANPASU.

---

## Repository structure

```
.
├── README.md                # this file
├── WORKFLOW.md              # full end-to-end pipeline with exact commands
├── requirements.txt         # Python dependencies
├── .gitignore               # excludes data, outputs, logs
├── scripts/                 # analysis code
│   ├── task1_compare_proteins.py   # Phase 1: protein-set comparison + Venn
│   ├── extract_proteins.py         # Phase 2: CDS to protein (Table 12)
│   ├── liftoff_compare.py          # Phase 2: lift-over efficiency
│   ├── ref_vs_lifted_scatter.py    # Phase 2: transfer-fidelity scatter
│   ├── find_missing_genes.py       # Phase 2: genes that failed to lift
│   ├── gene_stats.py               # Phase 4: gene/exon/intron stats
│   ├── pangenome_classify.py       # Phase 4: core/shell/cloud classification
│   └── README.md
├── slurm/                   # HPC job scripts (submit with sbatch)
│   ├── run_liftoff.sh
│   ├── run_extract_proteins.sh
│   ├── run_gene_stats.sh
│   ├── run_get_pangenes.sh
│   └── README.md
├── funannotate/             # RNA-seq prep for de novo re-annotation
│   ├── rnaseq_accessions.txt       # the 17 SRA runs used
│   ├── download_rnaseq.sh          # run on the LOGIN NODE (no internet on compute)
│   ├── funannotate_fastqc.SLURM.sh
│   └── funannotate_fix_headers.slurm.sh
├── pipeline/
│   └── run_cauris_plots.sh         # pan-genome stats + 6 figures
├── data/                    # inputs (NOT tracked) - see data/README.md
│   └── README.md
├── results/                 # outputs (NOT tracked)
│   └── README.md
└── docs/                    # notes, guides, plan
    ├── hpc_cheatsheet.md
    ├── meeting_notes_2026-06-10.md
    └── project_plan.md
```

---

## Quick start

```bash
# 1. Install Python dependencies (laptop or HPC virtualenv)
pip install -r requirements.txt

# 2. Phase 1 — compare B8441 protein sets
python3 scripts/task1_compare_proteins.py \
  --fungidb data/FungiDB-68_CaurisB8441_AnnotatedProteins.fasta \
  --genbank_v2 data/GCA_002759435.2_proteins.fasta \
  --genbank_v3 data/GCA_002759435.3_proteins.fasta \
  --out results/venn_diagram.png

# 3. Phase 4 — gene statistics from GFF3 annotation(s)
python3 scripts/gene_stats.py data/*.gff3 --outdir results/gene_stats_out --plot

# On the HPC, run the heavy steps via SLURM (never on the login node):
sbatch slurm/run_liftoff.sh
sbatch slurm/run_gene_stats.sh data/B8441.gff3 data/strain1.gff3 ...
```

## Key accessions (B8441 reference)

| Dataset      | Source      | Accession / filename                            |
|--------------|-------------|-------------------------------------------------|
| FungiDB v68  | fungidb.org | FungiDB-68_CaurisB8441_AnnotatedProteins.fasta  |
| GenBank v2   | NCBI        | GCA_002759435.2                                 |
| GenBank v3   | NCBI        | GCA_002759435.3                                 |

## Status — project complete

| Phase | Item | Status | Key result |
|-------|------|--------|------------|
| 1 | B8441 protein-set comparison | Complete | 5,360 proteins shared across all three annotations |
| 2 | Liftoff annotation transfer | Complete | 98.6-99.9% of ~5,600 genes mapped; r = 0.996-0.999 |
| 2 | Funannotate de novo (strain 6684) | Complete | 5,780 genes, 5,606 mRNAs, 175 tRNAs |
| 3 | GET_PANGENES pan-gene construction | Complete | 5,848 pan-gene clusters across 6 genomes |
| 4 | Core/shell classification + figures | Complete | 4,982 core (85.2%), 866 accessory (14.8%) |

### Headline findings

- **Pan-gene set:** 5,848 clusters across six genomes.
- **Composition:** 4,982 core (85.2%), 866 accessory (14.8%) - a largely conserved, relatively closed pan-genome.
- **Per-genome membership:** 5,217-5,553 clusters; the de novo re-annotated 6684 contributed 5,553, in line with the others.
- **Between-genome similarity:** POCS 94.6-99.0% across all pairs.
- **Annotation transfer fidelity:** Pearson r = 0.996-0.999 between reference and lifted gene lengths.

See `WORKFLOW.md` for the full end-to-end pipeline with exact commands and parameters.

## References

See the preliminary report for the full Harvard reference list (Liftoff,
Funannotate, BRAKER, InterProScan, FungiDB, GET_PANGENES, Vernikos *et al.*).

## Updates (Aug 2026)

New analysis scripts in `scripts/`:
- `liftoff_compare.py` — lift-over efficiency (mapped vs unmapped) across parameter sets.
- `ref_vs_lifted_scatter.py` — reference vs lifted gene-length scatter (sanity check of the transfer).
- `find_missing_genes.py` — genes present in reference but not lifted.

New `funannotate/` folder — RNA-seq preparation for de novo re-annotation:
- `rnaseq_accessions.txt` — 17 C. auris RNA-seq runs (SRA).
- `download_rnaseq.sh` — download commands (run on the LOGIN NODE in a screen session; compute nodes have no internet).
- `funannotate_fastqc.SLURM.sh`, `funannotate_fix_headers.slurm.sh` — QC and header-fix jobs.
