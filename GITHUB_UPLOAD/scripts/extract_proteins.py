#!/usr/bin/env python3
# Protein extraction and translation for Candida auris genomes.
# Uses NCBI Translation Table 12 (Alternative Yeast Nuclear / CTG clade):
#   CTG -> Ser (not Leu as in the standard code)
#
# Provided by Dr Evelina Basenko; used in Phase 2 to translate CDS features
# from a GFF3 + genome FASTA into a protein FASTA.

import sys
from collections import defaultdict


def parse_fasta(fasta_file):
    """Parse FASTA file, returning {seq_id: sequence} dict."""
    sequences = {}
    current_id = None
    current_seq = []
    with open(fasta_file) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_id:
                    sequences[current_id] = ''.join(current_seq)
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line.upper())
        if current_id:
            sequences[current_id] = ''.join(current_seq)
    return sequences


def reverse_complement(seq):
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N',
                  'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W',
                  'K': 'M', 'M': 'K', 'B': 'V', 'V': 'B',
                  'D': 'H', 'H': 'D'}
    return ''.join(complement.get(base, base) for base in reversed(seq))


def translate(dna_seq):
    """Translate a CDS using NCBI Translation Table 12 (CTG -> Ser)."""
    codon_table = {
        'TTT':'F','TTC':'F','TTA':'L','TTG':'L','TCT':'S','TCC':'S','TCA':'S','TCG':'S',
        'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','TGT':'C','TGC':'C','TGA':'*','TGG':'W',
        'CTT':'L','CTC':'L','CTA':'L','CTG':'S',
        'CCT':'P','CCC':'P','CCA':'P','CCG':'P',
        'CAT':'H','CAC':'H','CAA':'Q','CAG':'Q','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
        'ATT':'I','ATC':'I','ATA':'I','ATG':'M','ACT':'T','ACC':'T','ACA':'T','ACG':'T',
        'AAT':'N','AAC':'N','AAA':'K','AAG':'K','AGT':'S','AGC':'S','AGA':'R','AGG':'R',
        'GTT':'V','GTC':'V','GTA':'V','GTG':'V','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
        'GAT':'D','GAC':'D','GAA':'E','GAG':'E','GGT':'G','GGC':'G','GGA':'G','GGG':'G',
    }
    protein = []
    for i in range(0, len(dna_seq) - 2, 3):
        protein.append(codon_table.get(dna_seq[i:i+3].upper(), 'X'))
    return ''.join(protein)


def extract_proteins(gff_file, genome_file, output_file):
    print("Loading genome...", file=sys.stderr)
    genome = parse_fasta(genome_file)
    print(f"Loaded {len(genome)} sequences", file=sys.stderr)

    print("Parsing GFF...", file=sys.stderr)
    transcripts = defaultdict(list)
    with open(gff_file) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
            seqid, source, feature_type, start, end, score, strand, phase, attributes = fields
            if feature_type != 'CDS':
                continue
            parent = None
            for item in attributes.split(';'):
                if item.startswith('Parent='):
                    parent = item.split('=')[1]
                    if parent.startswith('transcript:'):
                        parent = parent.replace('transcript:', '') 
                    break
            if not parent:
                continue
            try:
                parsed_phase = int(phase)
            except ValueError:
                parsed_phase = 0
            transcripts[parent].append({'seqid': seqid, 'start': int(start),
                                        'end': int(end), 'strand': strand,
                                        'phase': parsed_phase})
    print(f"Found {len(transcripts)} transcripts", file=sys.stderr)

    print("Extracting and translating...", file=sys.stderr)
    count = 0
    warned_internal_stop = 0
    with open(output_file, 'w') as out:
        for transcript_id, cds_list in sorted(transcripts.items()):
            strand = cds_list[0]['strand']
            seqid = cds_list[0]['seqid']
            if seqid not in genome:
                print(f"  WARNING: {seqid} not found in genome (transcript {transcript_id})",
                      file=sys.stderr)
                continue
            cds_list.sort(key=lambda x: x['start'], reverse=(strand == '-'))
            first_phase = cds_list[0]['phase']
            cds_sequences = [genome[c['seqid']][c['start'] - 1:c['end']] for c in cds_list]
            full_cds = ''.join(cds_sequences)
            if strand == '-':
                full_cds = reverse_complement(full_cds)
            if first_phase > 0:
                full_cds = full_cds[first_phase:]
            remainder = len(full_cds) % 3
            if remainder > 0:
                full_cds = full_cds[:-remainder]
            if not full_cds:
                continue
            protein = translate(full_cds)
            if protein.endswith('*'):
                protein = protein[:-1]
            if '*' in protein:
                warned_internal_stop += 1
                print(f"  WARNING: internal stop codon in {transcript_id} "
                      f"(pos {protein.index('*') + 1})", file=sys.stderr)
            out.write(f">transcript:{transcript_id}\n")
            for i in range(0, len(protein), 60):
                out.write(protein[i:i+60] + '\n')
            count += 1
            if count % 5000 == 0:
                print(f"  Processed {count}...", file=sys.stderr)

    print(f"✓ Translated {count} proteins", file=sys.stderr)
    if warned_internal_stop:
        print(f"  ⚠ {warned_internal_stop} transcripts had internal stop codons",
              file=sys.stderr)


if __name__ == '__main__':
    if len(sys.argv) != 4:
        print("Usage: python3 extract_proteins.py <gff> <genome.fasta> <output.fasta>")
        sys.exit(1)
    extract_proteins(sys.argv[1], sys.argv[2], sys.argv[3])
<import> 
#the following command sequence is used to run the script with the provided arguments
# python3 extract_proteins.py <gff> <genome.fasta> <output.fasta>
#sequence fasta B86441 and gff file are used as input to extract and translate proteins, and the output is written to a specified fasta file.