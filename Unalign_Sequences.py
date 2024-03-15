
from Bio import SeqIO
from Bio.Seq import Seq
import argparse
import re

parser = argparse.ArgumentParser()
parser.add_argument("--input", help="fasta_file", type=str)
parser.add_argument("--output", help="fasta_file", type=str)
args = parser.parse_args()

with open(args.output, 'w', newline='') as OUT:           
    for seqobj in SeqIO.parse(args.input, "fasta"):
        seqobj.seq = Seq(re.sub(r'-', '', str(seqobj.seq).upper()))
        SeqIO.write(seqobj, OUT, "fasta")
