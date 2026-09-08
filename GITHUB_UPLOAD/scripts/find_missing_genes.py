#!/usr/bin/env python3
"""
find_missing_genes.py — find genes present in an existing annotation (e.g.
FungiDB) that have NO overlap with any gene in a Liftoff annotation, so they
can be merged back into the Liftoff set to complete it.

This implements the task from Andy's email (30 June 2026):
  "find any genes from current gene models (downloaded from FungiDB) that do
   not overlap ones in the LiftOff output. We could then add these into the
   LiftOff set, to complete the annotation with any genes that were missed
   by LiftOff."

Primary method: bedtools intersect -v (report 'existing' features with zero
overlap in 'liftoff'). Requires bedtools on PATH (`module load bedtools` on
the HPC). Optionally also runs gffcompare as a cross-check, since it reports
partial-overlap/structural differences that intersect -v does not surface —
"these genes compare different" from the notes.

Outputs (in --outdir)
  missing_genes.gff3      full records (gene + all descendant features) for
                           existing-annotation genes with no Liftoff overlap
  missing_genes_summary.txt   counts: existing genes, liftoff genes, missing
  gffcompare/               (only if --gffcompare and the binary is found)

Usage
  python3 find_missing_genes.py existing.gff3 liftoff.gff3 \
      --outdir results/missing_genes --type gene --gffcompare

Author: Faiza Khan (201936907) — MSc Bioinformatics, University of Liverpool
"""

import argparse
import os
import shutil
import subprocess
import sys


def gff_to_bed(gff_path, feature_types, bed_path):
    """Write a BED of matching features; return {id: line_index} count."""
    import re
    n = 0
    with open(gff_path) as fh, open(bed_path, "w") as out:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 9 or c[2] not in feature_types:
                continue
            m = re.search(r"ID=([^;]+)", c[8])
            gid = m.group(1) if m else f"feature_{n}"
            start0 = int(c[3]) - 1
            strand = c[6] if c[6] in ("+", "-") else "."
            out.write(f"{c[0]}\t{start0}\t{c[4]}\t{gid}\t.\t{strand}\n")
            n += 1
    return n


def run_bedtools_intersect_v(a_bed, b_bed, out_bed):
    if shutil.which("bedtools") is None:
        sys.exit("ERROR: bedtools not found on PATH. On the HPC: module load bedtools")
    with open(out_bed, "w") as out:
        subprocess.run(["bedtools", "intersect", "-v", "-a", a_bed, "-b", b_bed],
                        stdout=out, check=True)


def missing_ids_from_bed(bed_path):
    ids = set()
    with open(bed_path) as fh:
        for line in fh:
            c = line.rstrip("\n").split("\t")
            if len(c) >= 4:
                ids.add(c[3])
    return ids


def build_parent_map(gff_path):
    """Return {feature_id: parent_id_or_None} and {feature_id: feature_type}."""
    import re
    parent_of, type_of = {}, {}
    with open(gff_path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 9:
                continue
            attrs = c[8]
            idm = re.search(r"ID=([^;]+)", attrs)
            if not idm:
                continue
            fid = idm.group(1)
            type_of[fid] = c[2]
            pm = re.search(r"Parent=([^;]+)", attrs)
            parent_of[fid] = pm.group(1).split(",")[0] if pm else None
    return parent_of, type_of


def resolve_top_gene(fid, parent_of, type_of, gene_types):
    """Walk the Parent chain up to the top-level gene feature's ID."""
    seen = set()
    cur = fid
    top = fid if type_of.get(fid) in gene_types else None
    while cur is not None and cur not in seen:
        seen.add(cur)
        if type_of.get(cur) in gene_types:
            top = cur
        cur = parent_of.get(cur)
    return top


def extract_gene_records(gff_path, missing_gene_ids, gene_types, out_gff):
    import re
    parent_of, type_of = build_parent_map(gff_path)
    written = 0
    with open(gff_path) as fh, open(out_gff, "w") as out:
        out.write("##gff-version 3\n")
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 9:
                continue
            idm = re.search(r"ID=([^;]+)", c[8])
            fid = idm.group(1) if idm else None
            top = resolve_top_gene(fid, parent_of, type_of, gene_types) if fid else None
            if top in missing_gene_ids:
                out.write(line if line.endswith("\n") else line + "\n")
                written += 1
    return written


def run_gffcompare(existing_gff, liftoff_gff, outdir):
    if shutil.which("gffcompare") is None:
        print("gffcompare not found on PATH — skipping cross-check "
              "(module load gffcompare on the HPC). bedtools result is unaffected.",
              file=sys.stderr)
        return
    gc_dir = os.path.join(outdir, "gffcompare")
    os.makedirs(gc_dir, exist_ok=True)
    prefix = os.path.join(gc_dir, "cmp")
    subprocess.run(["gffcompare", "-r", existing_gff, "-o", prefix, liftoff_gff], check=True)
    print(f"gffcompare output written to {gc_dir} (see cmp.stats for class-code summary)")


def main():
    ap = argparse.ArgumentParser(
        description="Find existing-annotation genes with no overlap in a Liftoff GFF.")
    ap.add_argument("existing", help="Existing annotation GFF3 (e.g. FungiDB download)")
    ap.add_argument("liftoff", help="Liftoff output GFF3")
    ap.add_argument("--outdir", default="missing_genes_out")
    ap.add_argument("--type", default="gene",
                     help="Comma-separated top-level gene feature type(s) (default: gene)")
    ap.add_argument("--gffcompare", action="store_true",
                     help="Also run gffcompare as a structural cross-check")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    gene_types = set(args.type.split(","))

    existing_bed = os.path.join(args.outdir, "existing_genes.bed")
    liftoff_bed = os.path.join(args.outdir, "liftoff_genes.bed")
    n_existing = gff_to_bed(args.existing, gene_types, existing_bed)
    n_liftoff = gff_to_bed(args.liftoff, gene_types, liftoff_bed)
    print(f"Existing annotation genes: {n_existing}")
    print(f"Liftoff annotation genes:  {n_liftoff}")

    missing_bed = os.path.join(args.outdir, "missing_genes.bed")
    run_bedtools_intersect_v(existing_bed, liftoff_bed, missing_bed)
    missing_ids = missing_ids_from_bed(missing_bed)
    print(f"Existing genes with NO overlap in Liftoff: {len(missing_ids)}")

    out_gff = os.path.join(args.outdir, "missing_genes.gff3")
    written = extract_gene_records(args.existing, missing_ids, gene_types, out_gff)
    print(f"Wrote {len(missing_ids)} full gene records ({written} GFF lines) to {out_gff}")

    summary_path = os.path.join(args.outdir, "missing_genes_summary.txt")
    with open(summary_path, "w") as s:
        s.write(f"existing_genes\t{n_existing}\n")
        s.write(f"liftoff_genes\t{n_liftoff}\n")
        s.write(f"missing_no_overlap\t{len(missing_ids)}\n")
        pct = 100 * len(missing_ids) / n_existing if n_existing else 0
        s.write(f"pct_of_existing_missing\t{pct:.2f}\n")
    print(f"Wrote {summary_path}")

    if args.gffcompare:
        run_gffcompare(args.existing, args.liftoff, args.outdir)

    print("\nNext step: review missing_genes.gff3, then concatenate its records "
          "into the Liftoff GFF (sort by coordinate afterwards) to complete the annotation.")


if __name__ == "__main__":
    main()
