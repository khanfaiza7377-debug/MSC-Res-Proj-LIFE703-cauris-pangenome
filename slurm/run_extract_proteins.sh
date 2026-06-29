#!/bin/bash
#SBATCH --job-name=extract_proteins
#SBATCH --output=logs/extract_proteins_%j.out
#SBATCH --error=logs/extract_proteins_%j.err
#SBATCH --time=00:30:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1

# Translate C. auris CDS (GFF3 + genome FASTA) -> protein FASTA (Table 12).
# Submit from the repo root:
#   sbatch slurm/run_extract_proteins.sh data/B8441.gff3 data/B8441.fasta results/B8441.proteins.fasta
# NEVER run directly on the login node.

set -euo pipefail
mkdir -p logs

if [[ $# -ne 3 ]]; then
    echo "Usage: sbatch slurm/run_extract_proteins.sh <gff> <genome.fasta> <output.fasta>"
    exit 1
fi
GFF="$1"; GENOME="$2"; OUTPUT="$3"

for f in "$GFF" "$GENOME"; do
    [[ -f "$f" ]] || { echo "ERROR: input not found: $f"; exit 1; }
done

# module load python/3.10.0   # uncomment if needed

echo "Job ${SLURM_JOB_ID} on ${SLURMD_NODENAME} — started $(date)"
python3 scripts/extract_proteins.py "$GFF" "$GENOME" "$OUTPUT"
echo "Finished $(date)"
echo "Protein count: $(grep -c '^>' "$OUTPUT")"
