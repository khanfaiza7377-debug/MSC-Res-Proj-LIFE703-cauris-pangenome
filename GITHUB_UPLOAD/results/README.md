# results/

Generated outputs (figures, tables, lifted annotations, pan-genome matrices).
These are **not** tracked in git except this README — regenerate them by running
the scripts/SLURM jobs. Keep a copy of key figures and the final tables in your
OneNote record book for the report and for examiners.

Typical contents once the pipeline runs:

```
results/
├── venn_diagram.png                 # Phase 1
├── liftoff/<label>.liftoff.gff3     # Phase 2
├── gene_stats_out/                  # Phase 4 (genome_summary.tsv, box plot)
└── pangenome/                       # Phase 3-4 (presence_absence, composition)
```
