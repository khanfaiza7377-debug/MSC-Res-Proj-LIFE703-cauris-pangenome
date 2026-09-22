#!/bin/bash
#SBATCH --job-name=funannotate_6684
#SBATCH --partition=low
#SBATCH --cpus-per-task=48
#SBATCH --mem=128G
#SBATCH --time=200:00:00

# ============================================================
# run_funannotate_6684.sh
# De novo re-annotation of C. auris strain 6684 from 17 paired-end SRA runs.
# Funannotate v1.8.17, run from a Singularity container (funannotate.sif).
#
# Trimmomatic and Trinity k-mer normalisation both fail on the full paired-end
# set (Trimmomatic exhausts the Java heap), so reads are normalised beforehand
# and passed in with --no_trimmomatic --no_normalize_reads.
#
# Reconstructed from Funannotate_6884_Annotation_Guide_20AUG2026.
# ============================================================

set -euo pipefail

HC="${HC:?Set HC to your hc-storage root}"
SPECIESDIR="${SPECIESDIR:?Set SPECIESDIR to the 6684 working directory}"
SIF="${HC}/containers/funannotate.sif"
OUT="funannotate_6684_training"
GENOME="${SPECIESDIR}/repeatmask/softmasked-genome.sorted.fasta"
CPUS=48

run() { singularity exec -B "${SPECIESDIR}:${SPECIESDIR}" -B "${HC}:${HC}" --pwd "${SPECIESDIR}" "${SIF}" "$@"; }

# ---- 1. train: assemble reads, train the ab initio predictors ----
run funannotate train \
    -i "${GENOME}" \
    -o "${OUT}" \
    --left  "${OUT}/training/normalize/tmp_normalized_reads/left.norm.fq" \
    --right "${OUT}/training/normalize/tmp_normalized_reads/right.norm.fq" \
    --no_normalize_reads \
    --no_trimmomatic \
    --species Candida_auris \
    --strain 6684 \
    --cpus ${CPUS} \
    --max_intronlen 5000

# ---- 2. predict: call gene models ----
run funannotate predict \
    -i "${GENOME}" \
    -o "${OUT}" \
    -s "Candida_auris" \
    --strain "6684" \
    --cpus ${CPUS} \
    --organism fungus \
    --rna_bam "${OUT}/training/funannotate_train.coordSorted.bam" \
    --transcript_evidence "${OUT}/training/trinity.fasta.clean" \
    --stringtie "${OUT}/training/funannotate_train.stringtie.gtf" \
    --pasa_gff "${OUT}/training/funannotate_train.pasa.gff3" \
    --max_intronlen 5000

# ---- 3. annotate: functional assignment (Pfam / CAZyme / MEROPS) ----
run funannotate annotate \
    --gff   "${OUT}/predict_results/Candida_auris_6684.gff3" \
    --fasta "${OUT}/predict_results/Candida_auris_6684.scaffolds.fa" \
    --species "Candida auris" \
    --strain 6684 \
    --busco_db saccharomycetes \
    --out "${OUT}" \
    --cpus ${CPUS}

# ---- 4. update: add UTRs from the RNA-seq evidence ----
run funannotate update \
    --input "${OUT}/annotate_results/Candida_auris_6684.gbk" \
    --left  ${SPECIESDIR}/raw_data/*_1.fastq.gz \
    --right ${SPECIESDIR}/raw_data/*_2.fastq.gz \
    --out "${OUT}" \
    --cpus ${CPUS} \
    --max_intronlen 5000

echo "Funannotate 6684 pipeline complete."
