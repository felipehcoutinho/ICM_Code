#! /usr/bin/env python3
from Bio import SeqIO
import argparse
import re

parser = argparse.ArgumentParser()
parser.add_argument("--input_sequences", help="File containing the input sequences", type=str, nargs="+")
parser.add_argument("--protein", help="IDs are genomes but input file is of protein sequences", default=False, type=bool)
parser.add_argument("--matched_output_sequences", help="File to write matched output sequences", default='Matched_Sequences.fasta',type=str)
parser.add_argument("--unmatched_output_sequences", help="File to write unmatched output sequences", default='Unmatched_Sequences.fasta',type=str)
parser.add_argument("--list", help="Text file containing a list of IDs to bre retrieved from the input_sequences file", type=str)
parser.add_argument("--format_input", help="Format of the input sequences file", default='fasta', type=str)
parser.add_argument("--format_output", help="Format of the output sequences file", default='fasta', type=str)
parser.add_argument("--min_length", help="Minimum length of sequences to be include in the output file", default=0, type=int)
parser.add_argument("--max_length", help="Minimum length of sequences to be include in the output file", default=999999999, type=int)
parser.add_argument("--id_split_sep", help="Optional separator to split the Sequence IDs in the input fasta file", default="", type=str)
parser.add_argument("--id_split_pos", help="Optional index number used to match the split Sequence ID field to the IDs in the list file. Counts from 0", default=0, type=int)
parser.add_argument("--fetch_all", help="Flag to fetch all sequences regaless of them being in the list. Useful for simply filtering by length.", default=False, type=bool)
args = parser.parse_args()


def central():
    ids_list = ()
    if (args.list):
        ids_list = read_list(args.list)
    else:
        print("No file with list of sequence IDs provided. Will fetch all sequences within the length cutoffs.")
    fetch_seqs(args.input_sequences,args.format_input,args.matched_output_sequences,args.unmatched_output_sequences,args.format_output,ids_list,args.min_length,args.max_length,args.fetch_all)

def read_list(file):    
    ids = [line.rstrip() for line in open(file)]
    print('Obtained',len(ids),'Sequence identifiers from',file)
    if (len(set(ids)) != len(ids)):
        print('Warning! Repeated IDs.',len(set(ids)),'unique IDs')
    return set(ids)
    
def fetch_seqs(input_files_list,input_format,matched_output_file,unmatched_output_file,output_format,ids_list,min_length,max_length,fetch_all=False):
    if (len(ids_list) == 0):
        fetch_all = True
    total_seq_counter = 0
    passed_seq_counter = 0
    passed_seq_ids = []
    shouter = 1
    with open(matched_output_file, 'w', newline='') as OUT:
         with open(unmatched_output_file, 'w', newline='') as OUT2:
            for input_file in input_files_list:
                print(f"Processing {input_file}")
                for seqobj in SeqIO.parse(input_file, input_format):
                    total_seq_counter += 1
                    id = seqobj.id
                    if (args.protein):
                        id = re.sub('_(\d)+$','',id)
                    if (args.id_split_sep):
                        id = id.split(args.id_split_sep)[args.id_split_pos]
                    seq_length = len(seqobj.seq)
                    if (((id in ids_list) or (fetch_all == True)) and ((seq_length >= min_length) and (seq_length <= max_length))):
                        passed_seq_counter += 1
                        SeqIO.write(seqobj, OUT, output_format)
                        passed_seq_ids.append(id)
                    else:
                        SeqIO.write(seqobj, OUT2, output_format)
                    if (total_seq_counter == shouter):
                        print('Processed',total_seq_counter,'sequences. ',passed_seq_counter,' sequences passed.')
                        shouter += total_seq_counter
            print('Processed',total_seq_counter,'Sequences.',passed_seq_counter,'passed.')
    if ((args.protein == False) and (args.id_split_sep == "")):
        print("Checking if all IDs were retrieved...")
        for seqid in ids_list:
            if (seqid not in passed_seq_ids):
                print(f'{seqid} not retrieved!')
    if (passed_seq_counter != len(passed_seq_ids)):
        print('Warning! Passed sequence counter and number of IDs in the output file do not match!')

    
central()