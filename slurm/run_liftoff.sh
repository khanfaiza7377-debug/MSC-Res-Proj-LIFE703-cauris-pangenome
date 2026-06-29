#!/bin/bash
#SBATCH --job-name=liftoff
#SBATCH --output=logs/liftoff_%j.out
#SBATCH --error=logs/liftoff_%j.err
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --partition=batch

# ============================================================
# run_liftoff.sh  (TEMPLATE — edit the paths for your files)
# Transfers the GenBank v3 (B8441) annotation onto a target genome.
# Run one job per target genome (or use a SLURM array — see below).
#
# Submit from the repo root:
#   sbatch slurm/run_liftoff.sh data/target_strain1.fasta strain1
# NEVER run directly on the login node.
# ============================================================

set -euo pipefail
mkdir -p logs results/liftoff

# --- Reference (B8441, GenBank v3) ---
REF_FASTA="data/B8441_GCA_002759435.3.fasta"     # reference genome
REF_GFF="data/B8441_GCA_002759435.3.gff3"        # reference annotation to lift over

# --- Target (passed as arguments) ---
TARGET_FASTA="${1:?Usage: sbatch slurm/run_liftoff.sh <target.fasta> <label>}"
LABEL="${2:?Provide an output label, e.g. strain1}"
OUT_GFF="results/liftoff/${LABEL}.liftoff.gff3"
UNMAPPED="results/liftoff/${LABEL}.unmapped.txt"

# --- Environment (choose ONE that matches your cluster) ---
# module load liftoff
# source activate liftoff
# module load python/3.10.0 && pip install --user liftoff

echo "Job ${SLURM_JOB_ID} on ${SLURMD_NODENAME} — started $(date)"
echo "Reference: ${REF_FASTA} / ${REF_GFF}"
echo "Target   : ${TARGET_FASTA}  (label ${LABEL})"

liftoff \
  -g "${REF_GFF}" \
  -o "${OUT_GFF}" \
  -u "${UNMAPPED}" \
  -p "${SLURM_CPUS_PER_TASK}" \
  "${TARGET_FASTA}" \
  "${REF_FASTA}"

echo "Finished $(date)"
echo "Lifted annotation: ${OUT_GFF}"
echo "Unmapped genes    : $(wc -l < "${UNMAPPED}") (carry forward to Funannotate)"
