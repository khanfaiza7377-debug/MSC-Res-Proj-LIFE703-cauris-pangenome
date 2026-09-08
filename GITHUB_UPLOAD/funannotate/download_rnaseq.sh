#!/usr/bin/env bash
# ============================================================
# download_rnaseq.sh — download & extract RNA-seq for Funannotate
# Candida auris pan-gene project — Faiza Khan (201936907)
# 17 RNA-seq runs (see rnaseq_accessions.txt). Adapted from A. Ludwig's workflow.
#
# This download is LARGE (several very deep runs). Run it in an interactive
# session with network access, NOT as a heavy job on the login node:
#     srun --pty --cpus-per-task=8 --mem=32G --time=08:00:00 bash
#     bash download_rnaseq.sh
# ============================================================
set -euo pipefail

# >>> EDIT this to your own working directory on hc-storage <<<
WORKDIR="/mnt/hc-storage/groups/arjones_student_projects/faiza/RNAseq"
cd "$WORKDIR"

module load sra-tools/3.0.3

# 1. Download .sra files (--max-size 100G so the deep runs aren't capped)
prefetch --option-file rnaseq_accessions.txt -O sra_files --max-size 100G

# 2. Extract to FASTQ (paired-end auto-splits into _1 / _2)
mkdir -p fastq
while read -r acc; do
    [ -z "$acc" ] && continue
    fasterq-dump --split-files --threads 8 -O fastq "sra_files/${acc}/${acc}.sra"
done < rnaseq_accessions.txt

# 3. Compress
gzip fastq/*.fastq

echo "Done. FASTQ files are in ${WORKDIR}/fastq/"
