# slurm/

HPC batch scripts. **Always submit with `sbatch` from the repository root** so
that the relative paths (`scripts/`, `data/`, `results/`, `logs/`) resolve.
Never run analyses directly on the login node.

| Script | Submits |
|--------|---------|
| `run_extract_proteins.sh` | `sbatch slurm/run_extract_proteins.sh <gff> <genome.fasta> <out.fasta>` |
| `run_liftoff.sh` | `sbatch slurm/run_liftoff.sh <target.fasta> <label>` (one per target genome) |
| `run_get_pangenes.sh` | `sbatch slurm/run_get_pangenes.sh` (after editing input dir + command) |
| `run_gene_stats.sh` | `sbatch slurm/run_gene_stats.sh data/*.gff3` |

Before submitting:
1. `mkdir -p logs` (scripts also do this) — error logs are `logs/*_%j.err`.
2. Set the correct `module load` / `conda activate` line for each tool.
3. Check jobs with `squeue -u <username>`; cancel with `scancel <jobid>`.

`run_liftoff.sh` and `run_get_pangenes.sh` are **templates** — edit the input
paths and confirm exact tool options/flags with the supervisor before running.
