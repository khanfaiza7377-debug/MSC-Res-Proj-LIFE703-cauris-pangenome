# scripts/

Analysis code for the project. All scripts print usage if run without arguments.

| Script | Phase | What it does |
|--------|-------|--------------|
| `task1_compare_proteins.py` | 1 | Compares the three B8441 protein sets and draws a 3-way Venn diagram. Needs `biopython`, `matplotlib-venn`. |
| `extract_proteins.py` | 2 | Translates CDS (GFF3 + genome FASTA) → protein FASTA using NCBI Translation Table 12 (CTG→Ser). Standard library only. |
| `gene_stats.py` | 4 | Per-genome gene counts, gene lengths, exon and intron counts from GFF3; writes per-gene tables, a cross-genome summary, and an optional gene-length box plot. |
| `pangenome_classify.py` | 4 | Classifies pan-gene clusters as core / soft-core / shell / cloud from a presence-absence matrix; writes counts, per-cluster table, and an optional composition bar chart. |

Run heavy jobs through SLURM (see `../slurm/`), never on the login node.
