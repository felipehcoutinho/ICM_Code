#! /usr/bin/env python3
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
#from ete2 import NCBITaxa
import pandas as pd
import argparse
import re
import subprocess
import json

parser = argparse.ArgumentParser()
parser.add_argument("--input", help="List of fasta files containing sequences to be concatenated. One sequence will be geenrated per file", type =str, nargs="+")
parser.add_argument("--format_input", help="Input files sequence format", default="fasta", type =str)
parser.add_argument("--output", help="The output fasta file with the concatenatd sequences. A single file containign all the concatenated sequences will be generated", default="Concatenated.fasta", type =str)
parser.add_argument("--format_output", help="Output file sequence format", default="fasta", type =str)
parser.add_argument("--info_output", help="The sequence info output table file to be generated ", default="Concatenated_Sequence_Info.tsv", type = str)
parser.add_argument("--character",help="Character to put between the concatenated sequences. Default: N", default="N", type = str)
parser.add_argument("--reps",help="Charcater to put between the concatenated sequences. Number of times the character between the concatenated sequences will eb repeate: Default: 100", default=100, type = int)
args = parser.parse_args()

def index_seqs(in_seq_files=[]):
    OUTSEQ = open(args.output,'w', newline='')
    seq_counter = 0
    for seq_file in in_seq_files:
        #Define sequence ID based on file name with extension removed
        out_seq_id_terms = re.sub('(.)+/','',seq_file).split(".")
        out_seq_id_terms.pop()
        # out_seq_id_terms.append(args.format_output)
        out_seq_id = ".".join(out_seq_id_terms)
        print (f"Concatenatinging sequences from {seq_file} and printing to {out_seq_id}")
        # seen_cds_ids = dict()
        seen_scaff_ids = dict()
        concat_dict = defaultdict(dict)
        concat_seq_info["Original_File"][out_seq_id] = seq_file
        concat_seq_info["Original_Sequences"][out_seq_id] = []
        concat_list = []

        for seqobj in SeqIO.parse(seq_file, args.format_input):
            scaffold_id = seqobj.id
            #Initialize the genome info for the specific genomic sequence if it does not exist yet
            if (scaffold_id in seen_scaff_ids):
                print(f"WARNING! Repeated ID {scaffold_id} in {seq_file}. Make sure sequences are not duplicated")
            else:
                seen_scaff_ids[scaffold_id] = True
            concat_seq_info["Original_Sequences"][out_seq_id].append(scaffold_id)
            concat_list.append(str(seqobj.seq))
        
        spacer_seq = args.character * args.reps
        concat_seq = spacer_seq.join(concat_list)
        concat_seq_info['Length'][out_seq_id] = len(concat_seq)
        seqobj = SeqRecord(id= out_seq_id, seq=Seq(concat_seq), description="")
        SeqIO.write(seqobj, OUTSEQ, args.format_output)


def print_results():
    print(f"Printing concats seq info to {args.info_output}")
    concat_info_df = pd.DataFrame.from_dict(concat_seq_info )
    concat_info_df.index.name = 'Concatenated_Sequence'
    concat_info_df.to_csv(args.info_output,sep="\t",na_rep="NA")


# cds_seq_info = defaultdict(dict)
concat_seq_info = defaultdict(dict)
# genome_seq_info = defaultdict(dict)
seq_file_info = defaultdict(dict)

index_seqs(in_seq_files=args.input)

print_results()