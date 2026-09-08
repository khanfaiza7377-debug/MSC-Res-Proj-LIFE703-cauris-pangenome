# Meeting notes — 10 June 2026

**Present:** Faiza Khan · Prof. Andy Jones · Dr Evelina Basenko

## Summary
Task 1 (protein-set comparison) complete — the three B8441 annotations are
highly concordant (5,360 shared proteins). The project now moves into the
annotation-transfer phase using GenBank v3 as the reference.

## Next steps
- **Task 2 — Liftoff:** transfer the GenBank v3 annotation onto the other
  *C. auris* genomes (one run per target genome). Genes that fail to lift over
  carry forward to Task 3.
- **Task 3 — Funannotate:** re-annotate un-lifted genes; write a gene-statistics
  script (gene counts, gene lengths, intron/exon counts per gene).
- **Later — InterProScan:** functional/domain annotation once annotations are done.

## Analysis & plotting
- Eve to send a SLURM script + instructions for protein-length descriptive stats.
- Produce a box plot of protein length (always plot numerical data to sense-check).

## Follow-ups / admin
- IT access: HPC login on uni PCs, Cisco VPN on Mac, Teams group-chat access.
- Read the Liftoff paper (Shumate & Salzberg, 2021).

## Update (15 June 2026 — supervisor feedback on preliminary report)
- Pipeline name correction: the pan-gene pipeline is **GET_PANGENES** (not PANPASU).
- Plan was a little under-ambitious: Liftoff + Funannotate should finish within a
  few weeks (by end of June); InterProScan is quick afterwards.
- Scale up: annotate a **large set** of *C. auris* genomes (tens) and call
  **core / shell / cloud** across them.
- Add extra analysis steps modelled on the group's rice pan-gene paper.
