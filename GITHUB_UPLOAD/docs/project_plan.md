# Project plan (summary)

The formal preliminary report (project plan with aims, task/milestone/deliverable
table and Gantt chart) is the assessed Word document submitted on CANVAS. This
file is a working summary kept with the code, updated after supervisor feedback.

## Aim
Construct, analyse and interpret the pan-gene set of *Candida auris* across a
**large set of genomes (tens)** — producing consistent annotations, building a
presence–absence matrix with **GET_PANGENES** (Contreras-Moreira *et al.*, 2023),
and classifying genes as **core / shell / cloud** to study genomic diversity,
clade-specific genes and antifungal-resistance determinants.

## Phases & key deliverables
1. **Reference comparison** — B8441 protein-set comparison (✅ done; 5,360 shared).
2. **Annotation at scale** — Liftoff transfer + Funannotate (BRAKER3/HELIXER) on
   remaining genes across tens of genomes; proteins via Table 12.
   **Target: complete by ~4 July.**
3. **Functional annotation** — InterProScan on all proteins (**by ~18 July**).
4. **Pan-gene construction** — GET_PANGENES across all genomes → presence–absence
   matrix (**by ~1 Aug**).
5. **Occupancy & PAV** — core/shell/cloud classification; genes absent in B8441 but
   present in others; gene statistics (**by ~15 Aug**).
6. **Enrichment & interpretation** — domain/functional enrichment of accessory and
   clade-specific genes (resistance/adaptation); visualisations (gene-per-genome,
   composition, heatmap, pan-genome growth curve), following the group's rice
   PanOryza analysis (Contreras-Moreira *et al.*, 2025) (**by ~29 Aug**).

## Key dates (LIFE703 handbook)
- Preliminary report (CANVAS): **29 June 2026**
- Oral presentation: w/c 6 July 2026
- Draft of final report to supervisor: by 19 August 2026
- Final report + Project Work (A&T): **10 September 2026**

## Supervisor feedback incorporated (15 June 2026)
- Pipeline is **GET_PANGENES**, not PANPASU.
- Tighter timeline — Liftoff/Funannotate by end of June, InterProScan soon after.
- **Scale up** to tens of genomes; call **core/shell/cloud**.
- Added analyses modelled on the group's rice pan-gene paper (PAV, occupancy,
  functional enrichment, pan-genome growth curve).
