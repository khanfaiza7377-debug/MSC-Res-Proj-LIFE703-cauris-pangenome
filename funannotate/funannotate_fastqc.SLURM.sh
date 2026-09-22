#!/usr/bin/env bash
#SBATCH --partition=low
#SBATCH --job-name=fastqc
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=fastqc_%j.out
#SBATCH --error=fastqc_%j.err
# Quality-check the combined RNA-seq reads. Submit with: sbatch funannotate_fastqc.SLURM.sh
set -eo pipefail

# >>> EDIT to your own working directory <<<
WORKDIR="/mnt/hc-storage/groups/arjones_student_projects/faiza/RNAseq"
cd "$WORKDIR"

module load java/17.0.10+7
module load fastqc/0.11.9
module load multiqc/1.14

mkdir -p fastqc_combined
fastqc -t 8 -o fastqc_combined \
    combined_left_trinity.fastq.gz \
    combined_right_trinity.fastq.gz

multiqc fastqc_combined -o multiqc_combined
