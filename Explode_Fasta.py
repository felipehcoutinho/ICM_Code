#! /usr/bin/env python3

from collections import defaultdict
from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--input", help="File of genomic sequences", type =str)
parser.add_argument("--in_format", help="Input format of genomic sequences (Default: fasta)", type =str, default="fasta")
parser.add_argument("--out_format", help="Output format of genomic sequences (Default: fasta)", type =str, default="fasta")
args = parser.parse_args()

#2D dictionary to store all relevant information about sequences
seq_info = defaultdict(dict)

counter = 0
for seqobj in SeqIO.parse(args.input, args.in_format):
    counter += 1
    with open(f"{seqobj.id}.{args.out_format}", 'w') as OUT:
        SeqIO.write(seqobj, OUT, args.out_format)
    if ((counter % 1000) == 0):
        print(f"Processed {counter} Sequences")


