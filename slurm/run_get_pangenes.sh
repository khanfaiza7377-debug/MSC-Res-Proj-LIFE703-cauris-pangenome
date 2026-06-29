#!/bin/bash
#SBATCH --job-name=get_pangenes
#SBATCH --output=logs/get_pangenes_%j.out
#SBATCH --error=logs/get_pangenes_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --partition=batch

# ============================================================
# run_get_pangenes.sh  (TEMPLATE — confirm exact options with the supervisor)
# Builds the C. auris pan-gene set with the GET_PANGENES pipeline across all
# annotated genomes, producing a gene presence-absence matrix that feeds
# scripts/pangenome_classify.py (core / soft-core / shell / cloud).
#
# Place each genome's GFF + matching FASTA in data/pangenome_input/ following
# the GET_PANGENES input convention, then submit from the repo root:
#   sbatch slurm/run_get_pangenes.sh
# NEVER run directly on the login node.
# ============================================================

set -euo pipefail
mkdir -p logs results/pangenome

INPUT_DIR="data/pangenome_input"     # one subfolder/pair per genome
OUT_DIR="results/pangenome"

# --- Environment (set to match the cluster installation) ---
# module load get_pangenes
# source activate get_pangenes
# export PATH=$PATH:/path/to/get_pangenes/

echo "Job ${SLURM_JOB_ID} on ${SLURMD_NODENAME} — started $(date)"
echo "Input : ${INPUT_DIR}"
echo "Output: ${OUT_DIR}"

# Example invocation — replace with the exact GET_PANGENES command/flags:
# get_pangenes.pl -d "${INPUT_DIR}" -o "${OUT_DIR}" -t "${SLURM_CPUS_PER_TASK}"

echo "TODO: insert the GET_PANGENES command confirmed with Andy."
echo "Then classify with: python3 scripts/pangenome_classify.py ${OUT_DIR}/presence_absence.tsv --outdir ${OUT_DIR} --plot"
echo "Finished $(date)"
