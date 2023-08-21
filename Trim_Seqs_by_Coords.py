#! /usr/bin/env python3

from collections import defaultdict
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import pandas as pd
import argparse


parser = argparse.ArgumentParser()
parser.add_argument("--blacklist", help="txt file containing the list of sequences to be completely removed" , type=str)
parser.add_argument("--coords", help="tsv file containing the coordinates of the sequences to be trimmed" , type=str)
parser.add_argument("--genome_files", help="Fasta files containing the sequences to be trimmed or removed",  nargs="+", type=str)
parser.add_argument("--output", help="Output fasta file to write the passed sequences",  default="Trimmed_Sequences.fasta", type=str)

args = parser.parse_args()

blacklist = []
if (args.blacklist):
    with open(args.blacklist) as file:
        blacklist = file.readlines()
    blacklist = [seqid.strip() for seqid in blacklist] 

coords_df = []
if (args.coords):
    print(f'Reading info from {args.coords}')
    coords_df = pd.read_csv(args.coords, sep='\t',index_col=0,header=0)


seq_counter = 0
with open(args.output, 'w', newline='') as OUT:
    for file in args.genome_files:
        print(f"Processing {file}")
        for seqobj in SeqIO.parse(file,"fasta"):
            seq_counter += 1
            if (seqobj.id in blacklist):
                print(f"Removing {seqobj.id}")
            elif (seqobj.id in coords_df.index):
                print(f"Trimming {seqobj.id}")
                full_length = len(seqobj.seq)
                index_length = full_length - 1
                full_range = list(range(0,index_length))
                trim_coords = coords_df.loc[[seqobj.id]]
                for i,coord in trim_coords.iterrows():
                    #print(coord)
                    start_index = int(coord["nucleotide start"]) - 1
                    stop_index = int(coord["nucleotide stop"]) - 1
                    full_range = [posit for posit in full_range if posit not in list(range(start_index,stop_index))]
                #print(full_range)
                subseq = [seqobj.seq[posit] for posit in full_range]
                trimmed_seq = SeqRecord(seq=Seq("".join(subseq)), id=seqobj.id, name = seqobj.name, description=seqobj.description)
                print(trimmed_seq)
                SeqIO.write(trimmed_seq, OUT, "fasta")
            else:
                SeqIO.write(seqobj, OUT, "fasta")