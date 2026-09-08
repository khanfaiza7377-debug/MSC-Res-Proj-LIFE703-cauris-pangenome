#!/bin/bash
#SBATCH --job-name=gene_stats
#SBATCH --output=logs/gene_stats_%j.out
#SBATCH --error=logs/gene_stats_%j.err
#SBATCH --time=00:20:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --partition=batch

# Run gene_stats.py on one or more GFF3 annotations.
# Submit from the repo root:  sbatch slurm/run_gene_stats.sh data/*.gff3
# NEVER run directly on the login node.

set -euo pipefail
mkdir -p logs

module load python/3.10.0   # adjust to what `module avail` shows

if [[ $# -lt 1 ]]; then
    echo "Usage: sbatch slurm/run_gene_stats.sh <gff1> [gff2 ...]"
    exit 1
fi

echo "Job ${SLURM_JOB_ID} on ${SLURMD_NODENAME} — started $(date)"
python3 scripts/gene_stats.py "$@" --outdir results/gene_stats_out --plot
echo "Finished $(date) — outputs in results/gene_stats_out/"
