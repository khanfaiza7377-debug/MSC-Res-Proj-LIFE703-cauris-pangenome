# Pan-gene Set Analysis of *Candida auris*

MSc Bioinformatics research project (LIFE703), University of Liverpool.

**Student:** Faiza Khan (201936907)
**Supervisor:** Prof. Andy Jones · **Co-supervisor:** Dr Evelina Basenko

Code, job scripts and summary results for the construction and analysis of a
six-genome pan-gene set of *Candida auris*. Genome sequences, read data and large
intermediates are not tracked — they live on the University of Liverpool HPC and are
described in [`data/README.md`](data/README.md).

---

## What this project did

*C. auris* genomes have been annotated by different pipelines at different times, so
gene models are not directly comparable. Clustering such heterogeneous annotations
inflates the apparent accessory genome: a gene that is genuinely present can appear
absent simply because it was never annotated. This project treats **annotation
consistency as the controlled variable**.

1. **Reference comparison** — the three existing B8441 protein sets (FungiDB, GenBank
   v2, GenBank v3) were compared to establish a reliable reference. 5,360 proteins were
   common to all three.
2. **Annotation transfer** — the GenBank v3 B8441 annotation was lifted onto five target
   genomes with Liftoff, under two parameter sets, as a test of annotation consistency.
   98.6–99.9% of genes mapped; reference-versus-lifted gene lengths correlated at
   Pearson *r* = 0.996–0.999.
3. **De novo re-annotation** — strain 6684, whose existing annotation was fragmented, was
   rebuilt from 17 paired-end RNA-seq runs with Funannotate, yielding 5,780 gene models.
4. **Pan-gene construction** — six consistent annotations were clustered with
   GET_PANGENES into 5,848 pan-gene clusters: 4,982 (85.2%) core, 866 (14.8%) accessory.

## Results at a glance

| Output | File |
|---|---|
| Pan-genome composition and per-genome membership | [`results/pangenome/stats_summary.txt`](results/pangenome/stats_summary.txt) |
| Occupancy, POCS, multi-copy and UpSet figures | `results/pangenome/Fig*.png` |
| Liftoff efficiency, both parameter sets | [`results/liftoff/liftoff_efficiency.tsv`](results/liftoff/liftoff_efficiency.tsv) |
| Reference vs lifted gene-length concordance | [`results/liftoff/ref_vs_lifted_summary.tsv`](results/liftoff/ref_vs_lifted_summary.tsv) |
| Per-genome gene statistics | [`results/liftoff/genome_summary.tsv`](results/liftoff/genome_summary.tsv) |
| Three-way B8441 annotation overlap | `results/reference_comparison/venn_diagram.png` |
| RNA-seq run selection and rationale | [`results/funannotate/sra_selection_summary.tsv`](results/funannotate/sra_selection_summary.tsv) |

## Software versions

| Tool | Version | Notes |
|---|---|---|
| Liftoff | 1.6.3 | `module load liftoff/1.6.3` |
| Funannotate | 1.8.17 | Singularity container `funannotate.sif` |
| GET_PANGENES | not versioned | group installation; see the local usage guide |
| FastQC / MultiQC | cluster module set | read QC |
| Python / Biopython | cluster module set | see `requirements.txt` |

### Liftoff parameter sets

Neither set is the Liftoff default — both enable `-polish` and `-copies`.

| | `-sc` | `-flank` | `-mismatch` | `-gap_open` | `-gap_extend` |
|---|---|---|---|---|---|
| **Set A** | 0.90 | 0.3 | 1 | 1 | 1 |
| **Set B** | 0.85 | 0.1 | 1 | 2 | 1 |

Both returned identical mapping percentages and identical unmapped counts for every
genome, so the transfer is insensitive to the differences tested.

## Repository layout

```
scripts/      Python analysis and plotting
slurm/        SLURM job scripts (Liftoff, GET_PANGENES, protein extraction, gene stats)
funannotate/  De novo re-annotation of strain 6684, RNA-seq download and QC
pipeline/     Convenience wrapper for the plotting steps
results/      Summary outputs and figures (small files only)
docs/         Project plan, HPC notes, meeting notes
data/         Placeholder — describes the inputs held on the HPC
```

Full step-by-step commands are in [`WORKFLOW.md`](WORKFLOW.md).

## Reproducing the analysis

```bash
pip install -r requirements.txt

# 1. Three-way reference comparison
python3 scripts/task1_compare_proteins.py

# 2. Annotation transfer (HPC) — run once per parameter set
sbatch slurm/run_liftoff.sh A
sbatch slurm/run_liftoff.sh B
python3 scripts/liftoff_compare.py
python3 scripts/ref_vs_lifted_scatter.py

# 3. De novo re-annotation of strain 6684 (HPC)
bash funannotate/download_rnaseq.sh
sbatch funannotate/funannotate_fastqc.SLURM.sh
sbatch funannotate/run_funannotate_6684.sh

# 4. Pan-gene construction and classification (HPC)
sbatch slurm/run_get_pangenes.sh
python3 scripts/pangenome_classify.py

# 5. Pan-gene properties by occupancy
python3 scripts/pangene_occupancy_analysis.py \
    --matrix  data/pangene_matrix.tab \
    --summary data/Cauris_ref_genes_Summary_project.tsv \
    --outdir  results/pangenome/
```

Step 5 reproduces the occupancy-versus-gene-property analysis of
Contreras-Moreira *et al.* (2026), *Genome Research* 36(1):226–238, using a FungiDB
gene table for the B8441 reference.

## Data sources

- **B8441 reference** — GenBank `GCA_002759435.3` (v3); RefSeq `GCF_002759435.1`
- **Target genomes** — FungiDB release 69 (B11220, B11221, B11243, B11245, 6684)
- **Pan-gene clustering inputs** — FungiDB release 68 annotations, `GCF_002759435.1`,
  and the Funannotate annotation of 6684 produced here
- **RNA-seq** — 17 paired-end SRA runs, listed in
  [`funannotate/rnaseq_accessions.txt`](funannotate/rnaseq_accessions.txt)

## Acknowledgements

Liftoff parameters from Dr Adriana Ludwig; GET_PANGENES guidance from
Dr Bruno Contreras-Moreira; FungiDB support from Dr Evelina Basenko. Jobs were run in
the VEuPathDB annotation group's HPC allocation.
