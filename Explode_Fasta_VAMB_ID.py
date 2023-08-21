#! /usr/bin/env python3

from collections import defaultdict
from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--assembly", help="File of genomic sequences", type =str)
parser.add_argument("--in_format", help="Input format of genomic sequences", type =str, default="fasta")
args = parser.parse_args()

#2D dictionary to store all relevant information about sequences
seq_info = defaultdict(dict)

for seqobj in SeqIO.parse(args.assembly, args.in_format):
    (sample_id,seq_id) = seqobj.id.split('C')
    with open(f"{sample_id}.fasta", 'a') as OUT:
        SeqIO.write(seqobj, OUT, "fasta")



