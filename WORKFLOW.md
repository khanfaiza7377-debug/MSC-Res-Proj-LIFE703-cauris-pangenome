# Full analysis workflow — Pan-gene Set Analysis of *Candida auris*

End-to-end pipeline, in order, with the commands to reproduce each stage. Heavy
steps run on the University of Liverpool HPC via SLURM. Reference genome: B8441
(GCA_002759435.3). Target genomes: 6684, B11220, B11221, B11243, B11245.

## 1. Reference annotation comparison
Compare the three B8441 protein sets (FungiDB v68, GenBank v2, v3).
```
python3 scripts/task1_compare_proteins.py \
  --fungidb FungiDB-68_CaurisB8441_AnnotatedProteins.fasta \
  --genbank_v2 GCA_002759435.2_proteins.fasta \
  --genbank_v3 GCA_002759435.3_proteins.fasta \
  --out results/venn_diagram.png
```
Result: 5,360 proteins shared by all three.

## 2. Annotation transfer (Liftoff)
Transfer the GenBank v3 (B8441) annotation onto each target genome. Two parameter
sets were compared.
Default (Adriana):
```
liftoff -g $REF_GFF -o ${base}_liftoff.gff -u ${base}_unmapped.txt -p 16 \
  -polish -copies -sc 0.9 -flank 0.3 -mismatch 1 -gap_open 1 -gap_extend 1
```
Tuned for C. auris:
```
liftoff -g $REF_GFF -o ${base}_liftoff.gff -u ${base}_unmapped.txt -p 16 \
  -polish -copies -sc 0.85 -flank 0.1 -mismatch 1 -gap_open 2 -gap_extend 1
```
Result: 98.6-99.9% of ~5,600 genes mapped; the two sets were near-identical.

## 3. Protein extraction (Translation Table 12)
```
python3 scripts/extract_proteins.py <genome>.gff <genome>.fasta <genome>.proteins.fasta
# on the HPC: sbatch slurm/run_extract_proteins.sh <gff> <fasta> <out>
```

## 4. Liftoff QC and gene statistics
```
python3 scripts/liftoff_compare.py faiza_liftoff_adriana_param faiza_liftoff_new_param --outdir results --plot
python3 scripts/ref_vs_lifted_scatter.py B8441_reference.gff lifted_*.gff --outdir results
python3 scripts/gene_stats.py *.gff --outdir results --plot
```
Result: mapped/unmapped per genome; reference-vs-lifted r = 0.996-0.999; gene stats.

## 5. RNA-seq preparation (for Funannotate)
The 17 selected runs (funannotate/rnaseq_accessions.txt). Download on the LOGIN
NODE (compute nodes have no internet), in a screen session:
```
screen -S sra_download
module load sra-tools/3.0.3
prefetch --option-file rnaseq_accessions.txt -O sra_files --max-size 100G
mkdir -p fastq
while read -r acc; do fasterq-dump --split-files --threads 4 -O fastq "sra_files/${acc}/${acc}.sra"; done < rnaseq_accessions.txt
gzip fastq/*.fastq
# then, as SLURM jobs:
cat fastq/*_1.fastq.gz > combined_left_trinity.fastq.gz
cat fastq/*_2.fastq.gz > combined_right_trinity.fastq.gz
sbatch funannotate/funannotate_fastqc.SLURM.sh
sbatch funannotate/funannotate_fix_headers.slurm.sh
```

## 6. De novo annotation of 6684 (Funannotate)
Strain 6684 was re-annotated de novo because its existing annotation was poor and
fragmented. Run in the funannotate Singularity container on the HPC:
`funannotate train -> predict -> update -> annotate`
(exact container bind-mount commands and library fixes are in the team guides:
Funannotate_HPC_Complete_Guide.pdf). Summarise the run with:
```
sbatch funannotate_summary_job.sh   # uses funannotate_summary.py
```
Result: 6684 annotation = 5,780 genes, 5,606 mRNAs, 175 tRNAs; Pfam 4,497, CAZyme 121, MEROPS 193.
Note: InterProScan/GO/BUSCO were not run (not needed for the pan-genome).

## 7. Pan-gene construction (GET_PANGENES)
Run across the six genomes (Funannotate 6684 + four FungiDB genomes + B8441),
deployed with `pangene_job.sh`. Produces the presence-absence matrix
(pangene_matrix.tab), gene-ID matrix (pangene_matrix_genes.tab) and similarity
matrix (POCS.matrix.tab).

## 8. Pan-genome stats and figures
Using the prepared micromamba `dataviz` environment:
```
micromamba activate /mnt/hc-storage/groups/arjones_student_projects/faiza/envs/dataviz
cd Cauris/analysis
sbatch pipeline/run_cauris_plots.sh   # runs cauris_pangenome_stats_and_plots.py
```
Outputs: stats_summary.txt + 6 figures (composition, per-genome clusters,
occupancy histogram, UpSet, POCS heatmap, single/multi-copy).
