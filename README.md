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
LIFE703-cauris-pangenome/
├── README.md               # this file
├── requirements.txt        # Python dependencies
├── .gitignore              # excludes data, outputs, logs
├── scripts/                # analysis code
│   ├── task1_compare_proteins.py   # Phase 1: protein-set comparison + Venn
│   ├── extract_proteins.py         # Phase 2: CDS → protein (Table 12)
│   ├── gene_stats.py               # Phase 4: gene/exon/intron stats
│   ├── pangenome_classify.py       # Phase 4: core/shell/cloud classification
│   └── README.md
├── slurm/                  # HPC job scripts (submit with sbatch)
│   ├── run_extract_proteins.sh
│   ├── run_gene_stats.sh
│   ├── run_liftoff.sh
│   ├── run_get_pangenes.sh
│   └── README.md
├── data/                   # inputs (NOT tracked) — see data/README.md
│   └── README.md
├── results/                # outputs (NOT tracked) — figures/tables kept locally
│   └── README.md
└── docs/                   # notes, guides, plan
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

## Status

| Phase | Item | Status |
|-------|------|--------|
| 1 | B8441 protein comparison | ✅ Done (5,360 shared proteins) |
| 2 | Liftoff across target genomes | ⏳ In progress |
| 2 | Funannotate on un-lifted genes | ⏳ Pending |
| 3 | GET_PANGENES pan-gene construction | ⏳ Pending |
| 4 | Core/shell/cloud + stats + InterProScan | ⏳ Pending |

## References

See the preliminary report for the full Harvard reference list (Liftoff,
Funannotate, BRAKER, InterProScan, FungiDB, GET_PANGENES, Vernikos *et al.*).
