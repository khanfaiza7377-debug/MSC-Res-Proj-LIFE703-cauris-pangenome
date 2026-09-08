#!/usr/bin/env python3
"""
Phase 1 — Protein set comparison for Candida auris B8441.

Compares the annotated protein sets from three sources:
  - FungiDB v68
  - GenBank v2  (GCA_002759435.2)
  - GenBank v3  (GCA_002759435.3)

Produces shared/unique counts and a 3-way Venn diagram.

Fixes applied (flagged by Dr Evelina Basenko):
  1. Correct B8441 accessions — use GCA (not GCF).
  2. FungiDB sequences end with '*' (stop codon) while NCBI sequences do not;
     strip the trailing '*' so identical proteins aren't counted as unique.

Author: Faiza Khan (201936907) — MSc Bioinformatics, University of Liverpool
"""

import argparse
from Bio import SeqIO
from matplotlib_venn import venn3
import matplotlib.pyplot as plt


def load_seqs(fasta_path):
    """Return a set of unique protein sequences (trailing '*' removed)."""
    seqs = set()
    for record in SeqIO.parse(fasta_path, "fasta"):
        seqs.add(str(record.seq).rstrip("*"))
    return seqs


def main():
    ap = argparse.ArgumentParser(description="Compare C. auris B8441 protein sets.")
    ap.add_argument("--fungidb", required=True, help="FungiDB v68 proteins FASTA")
    ap.add_argument("--genbank_v2", required=True, help="GenBank v2 (GCA_002759435.2) FASTA")
    ap.add_argument("--genbank_v3", required=True, help="GenBank v3 (GCA_002759435.3) FASTA")
    ap.add_argument("--out", default="venn_diagram.png", help="Output PNG path")
    args = ap.parse_args()

    fungidb = load_seqs(args.fungidb)
    gb_v2 = load_seqs(args.genbank_v2)
    gb_v3 = load_seqs(args.genbank_v3)

    print(f"FungiDB v68 unique sequences:  {len(fungidb)}")
    print(f"GenBank v2 unique sequences:   {len(gb_v2)}")
    print(f"GenBank v3 unique sequences:   {len(gb_v3)}")
    print(f"Shared by all three:           {len(fungidb & gb_v2 & gb_v3)}")
    print(f"FungiDB only:                  {len(fungidb - gb_v2 - gb_v3)}")
    print(f"GenBank v2 only:               {len(gb_v2 - fungidb - gb_v3)}")
    print(f"GenBank v3 only:               {len(gb_v3 - fungidb - gb_v2)}")

    plt.figure(figsize=(8, 8))
    venn3([fungidb, gb_v2, gb_v3], set_labels=("FungiDB v68", "GenBank v2", "GenBank v3"))
    plt.title("C. auris B8441 protein set comparison")
    plt.savefig(args.out, dpi=300, bbox_inches="tight")
    print(f"Venn diagram saved to {args.out}")


if __name__ == "__main__":
    main()
