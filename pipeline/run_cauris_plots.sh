#!/bin/bash
#SBATCH --job-name=CaurisPlots_Faiza
#SBATCH --output=log_files/plots_%j.log
#SBATCH --error=err_files/plots_%j.err
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=02:30:00
# Pangenome stats & figures from a completed get_pangenes run.
# Requires the prepared micromamba 'dataviz' env. Submit from the analysis dir.
module purge
eval "$(micromamba shell hook --shell bash)"
micromamba activate /mnt/hc-storage/groups/arjones_student_projects/faiza/envs/dataviz
WORK_DIR=/mnt/hc-storage/groups/arjones_student_projects/faiza/Cauris/analysis
DATA_DIR=/mnt/hc-storage/groups/arjones_student_projects/faiza/Cauris/Cauris_pangenes_output
cd "$WORK_DIR"
python3 cauris_pangenome_stats_and_plots.py \
  --matrix "$DATA_DIR/pangene_matrix.tab" \
  --genes-matrix "$DATA_DIR/pangene_matrix_genes.tab" \
  --pocs "$DATA_DIR/POCS.matrix.tab" \
  --outdir figures
