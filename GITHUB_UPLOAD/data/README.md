# data/

Input data is **not** tracked in git (FASTA/GFF files are large — see the
root `.gitignore`). Keep the actual files on the HPC under `hc-storage`, and
download them straight onto the cluster with `wget` rather than via your laptop.

## Files needed

| File | Source | Notes |
|------|--------|-------|
| `FungiDB-68_CaurisB8441_AnnotatedProteins.fasta` | fungidb.org | Phase 1 reference proteins |
| `GCA_002759435.2_proteins.fasta` | NCBI | GenBank v2 proteins |
| `GCA_002759435.3_proteins.fasta` | NCBI | GenBank v3 proteins |
| `B8441_GCA_002759435.3.fasta` + `.gff3` | NCBI / FungiDB | Liftoff reference (genome + annotation) |
| target genome FASTAs | NCBI / FungiDB | the *C. auris* genomes to annotate (scale to tens) |
| `pangenome_input/` | derived | per-genome GFF+FASTA for GET_PANGENES |

## Example download (on the HPC)

```bash
cd /hc-storage/<hc>/<your_folder>/data
wget <fungidb_or_ncbi_url> -O B8441_GCA_002759435.3.fasta
```

Do not commit anything from this folder except this README.
