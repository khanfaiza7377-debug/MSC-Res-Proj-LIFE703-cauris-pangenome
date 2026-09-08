#!/usr/bin/env python3
"""
pangenome_classify.py — classify pan-gene clusters as core / soft-core / shell /
cloud from a gene presence-absence matrix (e.g. produced by GET_PANGENES).

Input: a TSV/CSV presence-absence matrix where
  * the first column is the gene/cluster identifier, and
  * each remaining column is a genome, with values that are either 0/1 or
    counts (any value > 0 is treated as "present").

Categories (Vernikos et al., 2015 convention; thresholds configurable):
  core       present in 100% of genomes
  soft-core  present in >= --soft (default 95%) but < 100%
  shell      present in > 1 genome but < --soft
  cloud      present in exactly 1 genome

Outputs
  pangenome_category_counts.tsv   counts and percentages per category
  pangenome_per_cluster.tsv       each cluster with its category and frequency
  pangenome_composition.png       (with --plot) bar chart of the categories

Usage
  python3 pangenome_classify.py presence_absence.tsv \
      [--soft 0.95] [--outdir results/pangenome] [--plot]

Author: Faiza Khan (201936907) — MSc Bioinformatics, University of Liverpool
"""

import argparse
import csv
import os


def read_matrix(path):
    """Read a presence-absence matrix; return (genomes, [(cluster, [present...])])."""
    delim = "\t" if path.endswith((".tsv", ".txt")) else ","
    rows = []
    with open(path, newline="") as fh:
        reader = csv.reader(fh, delimiter=delim)
        header = next(reader)
        genomes = header[1:]
        for r in reader:
            if not r:
                continue
            cluster = r[0]
            present = []
            for v in r[1:]:
                v = v.strip()
                try:
                    present.append(float(v) > 0)
                except ValueError:
                    present.append(v not in ("", "-", "NA", "0", "FALSE", "False"))
            rows.append((cluster, present))
    return genomes, rows


def classify(n_present, n_genomes, soft):
    if n_genomes == 0 or n_present == 0:
        return "absent"
    frac = n_present / n_genomes
    if frac >= 1.0:
        return "core"
    if frac >= soft:
        return "soft-core"
    if n_present == 1:
        return "cloud"
    return "shell"


def main():
    ap = argparse.ArgumentParser(description="Classify pan-gene clusters (core/shell/cloud).")
    ap.add_argument("matrix", help="Presence-absence matrix (TSV/CSV): clusters x genomes")
    ap.add_argument("--soft", type=float, default=0.95, help="Soft-core fraction (default 0.95)")
    ap.add_argument("--outdir", default="pangenome_out", help="Output directory")
    ap.add_argument("--plot", action="store_true", help="Write a composition bar chart (needs matplotlib)")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    genomes, rows = read_matrix(args.matrix)
    n_genomes = len(genomes)
    print(f"Genomes: {n_genomes} | clusters: {len(rows)}")

    order = ["core", "soft-core", "shell", "cloud"]
    counts = {c: 0 for c in order}
    per_cluster_path = os.path.join(args.outdir, "pangenome_per_cluster.tsv")
    with open(per_cluster_path, "w") as out:
        out.write("cluster\tn_present\tn_genomes\tfrequency\tcategory\n")
        for cluster, present in rows:
            n_present = sum(present)
            cat = classify(n_present, n_genomes, args.soft)
            counts[cat] = counts.get(cat, 0) + 1
            freq = round(n_present / n_genomes, 4) if n_genomes else 0
            out.write(f"{cluster}\t{n_present}\t{n_genomes}\t{freq}\t{cat}\n")

    total = sum(counts.values()) or 1
    counts_path = os.path.join(args.outdir, "pangenome_category_counts.tsv")
    with open(counts_path, "w") as out:
        out.write("category\tn_clusters\tpercent\n")
        print("\nPan-genome composition:")
        for c in order:
            pct = round(100 * counts[c] / total, 1)
            out.write(f"{c}\t{counts[c]}\t{pct}\n")
            print(f"  {c:10s}: {counts[c]} ({pct}%)")
        print(f"  {'TOTAL':10s}: {total}")
    print(f"\nWrote {counts_path} and {per_cluster_path}")

    if args.plot:
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            vals = [counts[c] for c in order]
            plt.figure(figsize=(7, 5))
            plt.bar(order, vals, color=["#1F4E79", "#5B9BD5", "#A9CCE3", "#D6EAF8"])
            for i, v in enumerate(vals):
                plt.text(i, v, str(v), ha="center", va="bottom")
            plt.ylabel("Number of gene clusters")
            plt.title(f"C. auris pan-genome composition ({n_genomes} genomes)")
            plt.tight_layout()
            out_png = os.path.join(args.outdir, "pangenome_composition.png")
            plt.savefig(out_png, dpi=300)
            print(f"Wrote {out_png}")
        except ImportError:
            print("matplotlib not available — skipping plot (TSVs still written).")


if __name__ == "__main__":
    main()
