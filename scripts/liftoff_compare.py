#!/usr/bin/env python3
"""
liftoff_compare.py — compare Liftoff lift-over efficiency between parameter sets.

Answers the question Eve posed: did the adjusted parameters affect Liftoff
efficiency? For each parameter-set directory it counts, per genome:
  * mapped genes   = number of 'gene' features in <base>_liftoff.gff
  * unmapped genes = number of IDs in <base>_unmapped.txt
  * efficiency (%) = mapped / (mapped + unmapped) * 100

Point it at the folders from Eve's Drive download, e.g.:
  python3 liftoff_compare.py faiza_liftoff_adriana_param faiza_liftoff_new_param \
      --outdir results/liftoff_compare --plot

It searches each directory recursively for *_liftoff.gff files and finds the
matching *_unmapped.txt next to each one. Standard library only; matplotlib is
needed only for --plot.

Author: Faiza Khan (201936907) — MSc Bioinformatics, University of Liverpool
"""

import argparse
import glob
import os


def count_genes(gff_path):
    """Count 'gene' features in a (Liftoff) GFF3 file."""
    n = 0
    with open(gff_path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.split("\t")
            if len(cols) >= 3 and cols[2] == "gene":
                n += 1
    return n


def count_unmapped(path):
    """Count non-empty, non-comment lines in an unmapped-genes file."""
    if not os.path.exists(path):
        return 0
    n = 0
    with open(path) as fh:
        for line in fh:
            if line.strip() and not line.startswith("#"):
                n += 1
    return n


def scan_dir(directory):
    """Return {genome_base: (mapped, unmapped)} for one parameter-set folder."""
    results = {}
    for gff in sorted(glob.glob(os.path.join(directory, "**", "*_liftoff.gff"), recursive=True)):
        base = os.path.basename(gff)[:-len("_liftoff.gff")]
        unmapped = os.path.join(os.path.dirname(gff), base + "_unmapped.txt")
        results[base] = (count_genes(gff), count_unmapped(unmapped))
    return results


def main():
    ap = argparse.ArgumentParser(description="Compare Liftoff efficiency across parameter sets.")
    ap.add_argument("dirs", nargs="+", help="Parameter-set directories (one per Liftoff run)")
    ap.add_argument("--outdir", default="liftoff_compare_out", help="Output directory")
    ap.add_argument("--plot", action="store_true", help="Write a grouped bar chart of % mapped (needs matplotlib)")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    per_set = {os.path.basename(os.path.normpath(d)): scan_dir(d) for d in args.dirs}
    genomes = sorted({g for d in per_set.values() for g in d})

    out_tsv = os.path.join(args.outdir, "liftoff_efficiency.tsv")
    with open(out_tsv, "w") as out:
        out.write("parameter_set\tgenome\tmapped_genes\tunmapped_genes\ttotal\tpct_mapped\n")
        print(f"{'param_set':28s} {'genome':22s} {'mapped':>7} {'unmap':>6} {'%map':>6}")
        for pset, data in per_set.items():
            for g in genomes:
                mapped, unmapped = data.get(g, (0, 0))
                total = mapped + unmapped
                pct = round(100 * mapped / total, 2) if total else 0.0
                out.write(f"{pset}\t{g}\t{mapped}\t{unmapped}\t{total}\t{pct}\n")
                print(f"{pset:28s} {g:22s} {mapped:7d} {unmapped:6d} {pct:6.1f}")
    print(f"\nWrote {out_tsv}")

    if args.plot:
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            sets = list(per_set.keys())
            x = range(len(genomes))
            width = 0.8 / max(len(sets), 1)
            plt.figure(figsize=(max(7, 1.2 * len(genomes)), 5))
            for i, pset in enumerate(sets):
                vals = []
                for g in genomes:
                    m, u = per_set[pset].get(g, (0, 0))
                    vals.append(100 * m / (m + u) if (m + u) else 0)
                plt.bar([xx + i * width for xx in x], vals, width=width, label=pset)
            plt.xticks([xx + width * (len(sets) - 1) / 2 for xx in x], genomes, rotation=30, ha="right")
            plt.ylabel("% genes mapped")
            plt.ylim(0, 100)
            plt.title("Liftoff efficiency by parameter set")
            plt.legend()
            plt.tight_layout()
            out_png = os.path.join(args.outdir, "liftoff_efficiency.png")
            plt.savefig(out_png, dpi=300)
            print(f"Wrote {out_png}")
        except ImportError:
            print("matplotlib not available — skipping plot (TSV still written).")


if __name__ == "__main__":
    main()
