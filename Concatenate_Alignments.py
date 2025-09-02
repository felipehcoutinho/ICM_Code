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
parser.add_argument("--aln_files", help="Files of aligned  sequences", type =str, nargs="+")
parser.add_argument("--in_format", help="Input Sequence file format", default="fasta", type =str)
parser.add_argument("--fasta_output", help="The output fasta file with the concatenatd alingments", default="Concatenated.fasta", type =str)
parser.add_argument("--info_cds_output", help="The CDS info output table file to be generated", default="CDS_Info.tsv", type =str)
parser.add_argument("--info_genome_output", help="The Genome info output table file to be generated ", default="Genome_Info.tsv", type =str)
parser.add_argument("--min_matches", help="Minimum number of sequences a concatenate must include to be printed to the output fasta file", default=0, type=int)
parser.add_argument("--clean_gvdb", help="Parses considering the sequences come from GVDB, i.e. and considers the genome identifer as everything before the first | or the first .fna or .fa or .fasta or .fas", default=False, type=bool)
args = parser.parse_args()

def index_seqs(in_seq_files=[]):
    # cds_to_file = defaultdict(dict)
    concat_dict = defaultdict(dict)
    cds_rm = re.compile("_(\d)+$")
    nw_rm = re.compile("\W")
    bar_rm = re.compile("\|.+$")
    ext_rm = re.compile("(\.fna.+)|(\.fa.+)")
    seq_counter = 0
    for seq_file in in_seq_files:
        print ("Indexing sequences from",seq_file)
        seen_cds_ids = dict()
        seen_scaff_ids = dict()
        for seqobj in SeqIO.parse(seq_file, "fasta"):
            cds_id = seqobj.id


            if seq_file not in seq_file_info['Alignment_Length']:
                seq_file_info['Alignment_Length'][seq_file] = len(seqobj.seq)
            if len(seqobj.seq) != seq_file_info['Alignment_Length'][seq_file]:
                exit(f"FATAL ERROR! Sequence length mismatch in {seq_file} for {cds_id}: {len(seqobj.seq)} vs {seq_file_info['Alignment_Length'][seq_file]}")

            scaffold_id = re.sub(cds_rm,"",cds_id)
            scaff_oid = scaffold_id

            if (args.clean_gvdb == True):
                # scaffold_id.split("|")[0]
                scaffold_id = re.sub(ext_rm,"",scaffold_id)
                scaffold_id = re.sub(bar_rm,"",scaffold_id)
                if (scaffold_id not in concat_seq_info['Original_Scaffold_IDs']):
                    concat_seq_info['Original_Scaffold_IDs'][scaffold_id] = set()
                concat_seq_info['Original_Scaffold_IDs'][scaffold_id].add(scaff_oid)

            #Initialize the genome info for the specific genomic sequence if it does not exist yet
            if (scaffold_id in seen_scaff_ids):
                print(f"WARNING! multiple matches to {scaffold_id} in {seq_file}. {cds_id} replacing last seen match")
            else:
                seen_scaff_ids[scaffold_id] = True

            #keep track of the number of sequences assigned to each scaffold in each file
            # if scaffold_id not in concat_seq_info[seq_file]:
            #     concat_seq_info[seq_file][scaffold_id] = 0 
            # concat_seq_info[seq_file][scaffold_id] += 1

            # keep track of the nCDS id assigned to each scaffold in each file
            if scaffold_id not in concat_seq_info[seq_file]:
                concat_seq_info[seq_file][scaffold_id] = "NA" 
            concat_seq_info[seq_file][scaffold_id] = cds_id
            
            #create a unique id combining file name and cds_id
            nid = seq_file + "|" + cds_id
            if (nid in cds_seq_info['Original_File']):
                print(f"WARNING! {cds_id} ids cuplicated in {seq_file}")
            else:
                cds_seq_info['Original_File'][nid] = set()
            
            cds_seq_info['Length'][nid] = len(seqobj.seq)
            cds_seq_info['Original_File'][nid].add(seq_file)
            
            concat_dict[scaffold_id][seq_file] = str(seqobj.seq)


    #print concatenated sequences
    print(f"Concatenating sequences and printing output to {args.fasta_output}")
    OUTSEQ = open(args.fasta_output,'w', newline='')
    for scaffold_id in concat_dict:
        concat_seq = ""
        matched_count = 0
        unmatched_count = 0
        concat_seq_info['Matched_Files'][scaffold_id] = set()
        for seq_file in in_seq_files:
            if (seq_file not in concat_dict[scaffold_id]):
                concat_seq += "-" * seq_file_info['Alignment_Length'][seq_file]
                unmatched_count += 1
            else:
                concat_seq += concat_dict[scaffold_id][seq_file]
                matched_count += 1
                concat_seq_info['Matched_Files'][scaffold_id].add(seq_file)
        
        concat_seq_info['Length'][scaffold_id] = len(concat_seq)
        concat_seq_info['Matched'][scaffold_id] = matched_count
        concat_seq_info['Unmatched'][scaffold_id] = unmatched_count
        if (matched_count >= args.min_matches):
            seqobj = SeqRecord(id= scaffold_id, seq=Seq(concat_seq))
            SeqIO.write(seqobj, OUTSEQ, "fasta")


def print_results(cds_output_df_file="CDS_Info.tsv",genome_output_df_file="Genome_Info.tsv"):
    print(f"Printing CDS info to {cds_output_df_file}")
    cds_info_df = pd.DataFrame.from_dict(cds_seq_info)
    cds_info_df.index.name = 'Sequence'
    cds_info_df.to_csv(cds_output_df_file,sep="\t",na_rep='NA')
    print(f"Printing concats seq info to Concat_Seq_Info.tsv")
    concat_info_df = pd.DataFrame.from_dict(concat_seq_info )
    concat_info_df.index.name = 'Concatenated_Sequence'
    concat_info_df.to_csv("Concat_Seq_Info.tsv",sep="\t",na_rep="NA")
    print(f"Printing file info to File_Info.tsv")
    file_info_df = pd.DataFrame.from_dict(seq_file_info)
    file_info_df.index.name = 'File'
    file_info_df.to_csv("File_Info.tsv",sep="\t",na_rep="NA")

cds_seq_info = defaultdict(dict)
concat_seq_info = defaultdict(dict)
genome_seq_info = defaultdict(dict)
seq_file_info = defaultdict(dict)

index_seqs(in_seq_files=args.aln_files)

print_results(cds_output_df_file=args.info_cds_output,genome_output_df_file=args.info_genome_output)