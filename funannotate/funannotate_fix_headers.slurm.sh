#!/bin/bash
#SBATCH --partition=low
#SBATCH --job-name=fix_headers
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=fix_headers_%j.out
#SBATCH --error=fix_headers_%j.err
# Fix FASTQ headers for Trinity/funannotate. Submit with: sbatch funannotate_fix_headers.slurm.sh
# BEFORE running: copy fix_headers.py into this folder and edit the 4 paths + the
# '+' line fix (see the guide).
set -euo pipefail

module load python/3.11.9

# >>> EDIT to your own working directory <<<
WORKDIR="/mnt/hc-storage/groups/arjones_student_projects/faiza/RNAseq"
cd "$WORKDIR"

python3 fix_headers.py
