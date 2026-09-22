#!/bin/bash
#SBATCH --job-name=liftoff_cauris_B8441
#SBATCH --partition=low
#SBATCH --output=log_files/liftoff_%j.out
#SBATCH --error=log_files/liftoff_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=200:00:00

# ============================================================
# run_liftoff.sh
# Transfers the GenBank v3 B8441 annotation (GCA_002759435.3) onto the five
# target genomes, under the two parameter sets compared in the report.
#
#   Set A : -sc 0.90  -flank 0.3  -mismatch 1  -gap_open 1  -gap_extend 1
#   Set B : -sc 0.85  -flank 0.1  -mismatch 1  -gap_open 2  -gap_extend 1
#
# NOTE: neither set is the Liftoff default — both enable -polish and -copies
# and change the flanking, mismatch and gap penalties.
#
# Usage:  sbatch slurm/run_liftoff.sh A
#         sbatch slurm/run_liftoff.sh B
# ============================================================

set -euo pipefail

SET="${1:?Usage: sbatch slurm/run_liftoff.sh <A|B>}"

case "$SET" in
  A) SC=0.9  ; FLANK=0.3 ; GAP_OPEN=1 ;;
  B) SC=0.85 ; FLANK=0.1 ; GAP_OPEN=2 ;;
  *) echo "Parameter set must be A or B" >&2; exit 1 ;;
esac

INPUT_DIR="${INPUT_DIR:-data/liftoff}"
OUTPUT_DIR="${OUTPUT_DIR:-results/liftoff_set_${SET}}"

module purge
module load liftoff/1.6.3

mkdir -p "${OUTPUT_DIR}/liftoff_output" log_files err_files

# Reference genome and annotation (B8441 version 3)
REF_FASTA="${INPUT_DIR}/GCA_002759435.3_Cand_auris_B8441_V3_genomic.fna"
REF_GFF="${INPUT_DIR}/GCA_002759435.3_Cand_auris_B8441_V3_genomic.gff"

TARGET_GENOMES=(
    "${INPUT_DIR}/FungiDB-69_Cauris6684_Genome.fasta"
    "${INPUT_DIR}/FungiDB-69_CaurisB11220_Genome.fasta"
    "${INPUT_DIR}/FungiDB-69_CaurisB11221_Genome.fasta"
    "${INPUT_DIR}/FungiDB-69_CaurisB11243_Genome.fasta"
    "${INPUT_DIR}/FungiDB-69_CaurisB11245_Genome.fasta"
)

echo "Parameter set ${SET}: -sc ${SC} -flank ${FLANK} -gap_open ${GAP_OPEN}"

for target_genome in "${TARGET_GENOMES[@]}"; do
    base=$(basename "$target_genome" _Genome.fasta)
    echo ">>> Processing: $base"

    liftoff \
        -g "$REF_GFF" \
        -o "${OUTPUT_DIR}/liftoff_output/${base}_liftoff.gff" \
        -u "${OUTPUT_DIR}/liftoff_output/${base}_unmapped.txt" \
        -p 16 \
        -polish \
        -copies -sc "${SC}" \
        -flank "${FLANK}" \
        -mismatch 1 \
        -gap_open "${GAP_OPEN}" \
        -gap_extend 1 \
        "$target_genome" \
        "$REF_FASTA"

    echo ">>> Done: $base"
done

echo "All liftoff jobs completed for parameter set ${SET}."
