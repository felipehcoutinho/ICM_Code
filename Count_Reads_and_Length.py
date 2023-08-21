#! /usr/bin/env python3
#This script takes a list of directories as input, then iterates through set of fastq.gz files in each directory, counts the total number or reads, total base pairs, and calculates average read length of each fastq file

from collections import defaultdict
from Bio import SeqIO
from Bio.SeqUtils import GC
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd
import argparse
import gzip
import glob


parser = argparse.ArgumentParser()
parser.add_argument("--format", help="Input Sequence file format", default="fastq", type =str)
parser.add_argument("--directories", help="Directory containing the fastq files", type =str, required=True, nargs='+')
parser.add_argument("--extension", help="Extension of the fastq files", type =str, required=True)
parser.add_argument("--output", help="The output table file generated", default="Sample_Info.tsv", type =str)
args = parser.parse_args()


def count_seqs():
    sample_info_df = pd.DataFrame(data=None, columns=['Total_Reads','Total_Basepairs','Average_Read_Length'])
    for input_dir in args.directories:
        print(f"Processing {input_dir}")
        input_files = glob.glob(f'{input_dir}/*{args.extension}')
        for fq_file in input_files:
            print(f"\tProcessing {fq_file}")
            seq_counter = 0
            base_counter = 0
            with gzip.open(fq_file, "rt") as SEQHANDLE:
                for seqobj in SeqIO.parse(SEQHANDLE, args.format):
                    seq_counter += 1
                    base_counter += len(seqobj.seq)
                    if (seq_counter % 10000000 == 0):
                        print(f"\t\tProcessed {seq_counter} sequences")
            sample_info_df.loc[fq_file,'Total_Reads'] = seq_counter
            sample_info_df.loc[fq_file,'Total_Basepairs'] = base_counter
            sample_info_df.loc[fq_file,'Average_Read_Length'] = base_counter/seq_counter
    sample_info_df.to_csv(args.output,sep="\t",na_rep='NA')

count_seqs()