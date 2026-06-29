#!/usr/bin/env python3
"""
gene_stats.py — general gene statistics for Candida auris genome annotations.

For one or more GFF3 annotation files, computes per genome:
  * gene count
  * gene lengths (bp)        -> mean / median / min / max / total
  * exon counts per gene     -> mean / median / max, single-exon gene count
  * intron counts per gene   -> mean / total   (introns = exons - 1 per gene)
  * number of mRNA transcripts

Outputs
  <label>_gene_stats.tsv   per-gene table (one row per gene), per input file
  genome_summary.tsv       one row per genome, for cross-genome comparison
  gene_length_boxplot.png  (optional, with --plot) box plot of gene lengths

Only the Python standard library is required for the statistics; matplotlib is
needed only for --plot. The script therefore runs on a laptop or on the HPC.

Usage
  python3 gene_stats.py annotation1.gff3 [annotation2.gff3 ...] \
      [--outdir stats_out] [--plot]

Author: Faiza Khan (201936907) — MSc Bioinformatics, University of Liverpool
"""

import argparse
import os
import statistics as st
from collections import defaultdict

GENE_TYPES = {"gene", "protein_coding_gene", "ncRNA_gene", "pseudogene"}
MRNA_TYPES = {"mrna", "transcript"}


def parse_attributes(attr_field):
    attrs = {}
    for item in attr_field.strip().split(";"):
        if not item or "=" not in item:
            continue
        key, val = item.split("=", 1)
        attrs[key.strip()] = val.strip()
    return attrs


def strip_prefix(value):
    if value is None:
        return None
    if ":" in value:
        value = value.split(":", 1)[1]
    return value


def parse_gff(gff_file):
    genes = {}
    mrnas_by_gene = defaultdict(list)
    exons_by_parent = defaultdict(int)
    cds_by_parent = defaultdict(int)
    with open(gff_file) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            seqid, _src, ftype, start, end, _score, strand, _phase, attr = fields
            ftype_l = ftype.lower()
            attrs = parse_attributes(attr)
            fid = strip_prefix(attrs.get("ID"))
            parent = strip_prefix(attrs.get("Parent"))
            if ftype in GENE_TYPES or ftype_l.endswith("_gene") or ftype_l == "gene":
                if fid:
                    genes[fid] = {"seqid": seqid, "start": int(start),
                                  "end": int(end), "strand": strand}
            elif ftype_l in MRNA_TYPES:
                if fid and parent:
                    mrnas_by_gene[parent].append(fid)
            elif ftype_l == "exon":
                if parent:
                    exons_by_parent[parent] += 1
            elif ftype == "CDS":
                if parent:
                    cds_by_parent[parent] += 1
    return genes, mrnas_by_gene, exons_by_parent, cds_by_parent


def gene_table(genes, mrnas_by_gene, exons_by_parent, cds_by_parent):
    rows = []
    for gene_id, g in genes.items():
        length = g["end"] - g["start"] + 1
        transcripts = mrnas_by_gene.get(gene_id, [])
        if transcripts:
            exon_counts = [exons_by_parent.get(t, 0) for t in transcripts]
            cds_counts = [cds_by_parent.get(t, 0) for t in transcripts]
        else:
            exon_counts = [exons_by_parent.get(gene_id, 0)]
            cds_counts = [cds_by_parent.get(gene_id, 0)]
        n_exons = max(exon_counts) if exon_counts else 0
        exon_source = "exon"
        if n_exons == 0 and max(cds_counts, default=0) > 0:
            n_exons = max(cds_counts)
            exon_source = "CDS"
        rows.append({"gene_id": gene_id, "seqid": g["seqid"], "strand": g["strand"],
                     "gene_length": length, "n_transcripts": len(transcripts),
                     "n_exons": n_exons, "n_introns": max(n_exons - 1, 0),
                     "exon_source": exon_source})
    rows.sort(key=lambda r: (r["seqid"], r["gene_id"]))
    return rows


def summarise(label, rows):
    lengths = [r["gene_length"] for r in rows]
    exons = [r["n_exons"] for r in rows]
    introns = [r["n_introns"] for r in rows]
    n_genes = len(rows)
    single_exon = sum(1 for e in exons if e <= 1)

    def safe(fn, data, default=0):
        return fn(data) if data else default

    return {"genome": label, "n_genes": n_genes,
            "n_transcripts": sum(r["n_transcripts"] for r in rows),
            "gene_len_mean": round(safe(st.mean, lengths), 1),
            "gene_len_median": safe(st.median, lengths),
            "gene_len_min": safe(min, lengths), "gene_len_max": safe(max, lengths),
            "gene_len_total": sum(lengths),
            "exons_per_gene_mean": round(safe(st.mean, exons), 2),
            "exons_per_gene_max": safe(max, exons),
            "introns_per_gene_mean": round(safe(st.mean, introns), 2),
            "introns_total": sum(introns), "single_exon_genes": single_exon,
            "single_exon_pct": round(100 * single_exon / n_genes, 1) if n_genes else 0}


def write_tsv(path, header, rows):
    with open(path, "w") as out:
        out.write("\t".join(header) + "\n")
        for r in rows:
            out.write("\t".join(str(r[h]) for h in header) + "\n")


def main():
    ap = argparse.ArgumentParser(description="Gene statistics from GFF3 annotation(s).")
    ap.add_argument("gff", nargs="+", help="One or more GFF3 annotation files")
    ap.add_argument("--outdir", default="gene_stats_out", help="Output directory")
    ap.add_argument("--plot", action="store_true", help="Write a gene-length box plot (needs matplotlib)")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    per_gene_header = ["gene_id", "seqid", "strand", "gene_length",
                       "n_transcripts", "n_exons", "n_introns", "exon_source"]
    summary_header = ["genome", "n_genes", "n_transcripts", "gene_len_mean",
                      "gene_len_median", "gene_len_min", "gene_len_max",
                      "gene_len_total", "exons_per_gene_mean", "exons_per_gene_max",
                      "introns_per_gene_mean", "introns_total", "single_exon_genes",
                      "single_exon_pct"]
    summaries = []
    lengths_by_genome = {}
    for gff in args.gff:
        label = os.path.splitext(os.path.basename(gff))[0]
        genes, mrnas, exons, cds = parse_gff(gff)
        rows = gene_table(genes, mrnas, exons, cds)
        write_tsv(os.path.join(args.outdir, f"{label}_gene_stats.tsv"), per_gene_header, rows)
        summ = summarise(label, rows)
        summaries.append(summ)
        lengths_by_genome[label] = [r["gene_length"] for r in rows]
        print(f"\n=== {label} ===")
        for k in summary_header[1:]:
            print(f"  {k:22s}: {summ[k]}")
    write_tsv(os.path.join(args.outdir, "genome_summary.tsv"), summary_header, summaries)
    print(f"\nWrote per-gene tables and genome_summary.tsv to {args.outdir}/")

    if args.plot:
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            labels = list(lengths_by_genome.keys())
            data = [lengths_by_genome[l] for l in labels]
            plt.figure(figsize=(max(6, 1.5 * len(labels)), 6))
            plt.boxplot(data, labels=labels, showfliers=False)
            plt.ylabel("Gene length / bp")
            plt.title("Gene length distribution per genome")
            plt.xticks(rotation=30, ha="right")
            plt.tight_layout()
            out_png = os.path.join(args.outdir, "gene_length_boxplot.png")
            plt.savefig(out_png, dpi=300)
            print(f"Wrote {out_png}")
        except ImportError:
            print("matplotlib not available — skipping plot (stats TSVs still written).")


if __name__ == "__main__":
    main()
