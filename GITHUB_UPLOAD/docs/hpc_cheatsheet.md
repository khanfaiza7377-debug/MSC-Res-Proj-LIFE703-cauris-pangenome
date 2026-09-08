# HPC cheat sheet (CBF cluster)

Source: CBF HPC Cluster Guide (Dr E. Basenko) + personal notes.

## The golden rule
**Never run jobs on the login node** (`login01.pgb.liv.ac.uk`). It is shared by
everyone — heavy jobs there can crash the system and get your account suspended.
Everything heavy goes through **SLURM** (`sbatch`).

| OK on login node | Must use SLURM |
|------------------|----------------|
| Editing files, transferring files | Liftoff, Funannotate, BLAST, InterProScan |
| Loading modules, checking job status | Any multi-core / multi-threaded work, interactive analysis |

## Log in
```bash
ssh <username>@login01.pgb.liv.ac.uk      # off campus: connect to Cisco VPN first
```

## Storage
- `/home/<username>` — scripts/config/small results (small quota, backed up). Do **not** run pipelines here.
- `/hc-storage/` — large data and outputs. Use the `hc` folder; make your own directory inside it.

## Modules
```bash
module avail
module load python/3.10.0
module list
```

## SLURM
```bash
sbatch job.sh                 # submit
squeue -u <username>          # your jobs
scancel <jobid>              # cancel
```
Keep a `logs/` folder; check `*_%j.err` first if something fails. Liftoff can run
for hours — submit and check the next day.

## Transfer files
```bash
scp myfile.fasta <username>@login01.pgb.liv.ac.uk:/hc-storage/<folder>/   # upload
scp <username>@login01.pgb.liv.ac.uk:/hc-storage/<folder>/results.txt .   # download
wget <url>                                                                # straight onto HPC
```
