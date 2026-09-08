#!/usr/bin/env python3
"""
ref_vs_lifted_scatter.py — sanity-check the Liftoff transfer.

For each gene, compares its length in the REFERENCE annotation with the length
of its lifted-over copy in each target genome. Liftoff preserves gene IDs, so
genes are matched by ID. If the transfer is faithful, points lie on the y=x
diagonal (same gene model → same length).

Outputs
  ref_vs_lifted_scatter.png     scatter panels (reference vs lifted gene length)
  ref_vs_lifted_summary.tsv     per-genome: matched genes, % identical length, correlation

Usage
  python3 ref_vs_lifted_scatter.py <reference.gff> <lifted1.gff> [lifted2.gff ...] \
      [--outdir results] [--type gene]

Author: Faiza Khan (201936907) — MSc Bioinformatics, University of Liverpool
"""

import argparse, os, re, math


def gene_lengths(gff, gene_types):
    """Return {gene_id: length_bp} for gene features in a GFF3 file."""
    lengths = {}
    with open(gff) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 9 or c[2] not in gene_types:
                continue
            m = re.search(r"ID=([^;]+)", c[8])
            if not m:
                continue
            lengths[m.group(1)] = int(c[4]) - int(c[3]) + 1
    return lengths


def pearson(xs, ys):
    n = len(xs)
    if n < 2:
        return float("nan")
    mx, my = sum(xs) / n, sum(ys) / n
    num = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    dx = math.sqrt(sum((x - mx) ** 2 for x in xs))
    dy = math.sqrt(sum((y - my) ** 2 for y in ys))
    return num / (dx * dy) if dx and dy else float("nan")


def main():
    ap = argparse.ArgumentParser(description="Reference vs lifted gene-length scatter.")
    ap.add_argument("reference", help="Reference GFF3 (the -g annotation used by Liftoff)")
    ap.add_argument("lifted", nargs="+", help="Lifted GFF3 file(s)")
    ap.add_argument("--outdir", default="ref_vs_lifted_out")
    ap.add_argument("--type", default="gene,protein_coding_gene",
                    help="Comma-separated gene feature types (default gene,protein_coding_gene)")
    ap.add_argument("--plot", action="store_true", default=True)
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    gtypes = set(args.type.split(","))
    ref = gene_lengths(args.reference, gtypes)
    print(f"Reference genes: {len(ref)}")

    per = []  # (label, xs, ys, identical, matched)
    for gff in args.lifted:
        label = os.path.basename(gff).replace("_liftoff.gff", "").replace(".gff", "")
        lift = gene_lengths(gff, gtypes)
        xs, ys, identical = [], [], 0
        for gid, llen in lift.items():
            if gid in ref:
                xs.append(ref[gid]); ys.append(llen)
                if ref[gid] == llen:
                    identical += 1
        per.append((label, xs, ys, identical, len(xs)))
        print(f"{label}: matched {len(xs)}, identical length {identical} "
              f"({100*identical/len(xs):.1f}%), r={pearson(xs,ys):.4f}")

    with open(os.path.join(args.outdir, "ref_vs_lifted_summary.tsv"), "w") as o:
        o.write("genome\tmatched_genes\tidentical_length\tpct_identical\tpearson_r\n")
        for label, xs, ys, ident, matched in per:
            pct = 100 * ident / matched if matched else 0
            o.write(f"{label}\t{matched}\t{ident}\t{pct:.2f}\t{pearson(xs,ys):.4f}\n")

    if args.plot:
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            n = len(per)
            cols = 3
            rows = (n + cols - 1) // cols
            fig, axes = plt.subplots(rows, cols, figsize=(4.2 * cols, 4.0 * rows))
            axes = axes.flatten() if n > 1 else [axes]
            for ax, (label, xs, ys, ident, matched) in zip(axes, per):
                ax.scatter(xs, ys, s=4, alpha=0.3, color="#1C7293")
                lim = max(max(xs, default=1), max(ys, default=1))
                ax.plot([0, lim], [0, lim], color="#B85042", lw=1)  # y=x
                ax.set_xlabel("Reference gene length (bp)")
                ax.set_ylabel("Lifted gene length (bp)")
                ax.set_title(f"{label}\n{matched} genes, {100*ident/matched:.1f}% identical")
            for ax in axes[n:]:
                ax.axis("off")
            plt.tight_layout()
            out = os.path.join(args.outdir, "ref_vs_lifted_scatter.png")
            plt.savefig(out, dpi=200)
            print("Wrote", out)
        except ImportError:
            print("matplotlib not available — TSV written, plot skipped.")


if __name__ == "__main__":
    main()
