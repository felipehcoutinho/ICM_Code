#! /usr/bin/env python3
from collections import defaultdict
from Bio import SeqIO
from Bio import SearchIO
from Bio.Seq import Seq
import pandas as pd
import argparse
import subprocess
import re
import glob
import os

parser = argparse.ArgumentParser()
parser.add_argument("--m8_files", help="File of genomic sequences", type =str, nargs="+")
parser.add_argument("--max_evalue", help="Maximum e-value", default=0.00001, type=float)
parser.add_argument("--min_identity", help="Minimum identity (0-100 for blast and 0-1 for mmseqs)", default=30.0, type=float)
parser.add_argument("--min_bitscore", help="Minimum bitscore", default=50, type=int)
parser.add_argument("--min_alignment", help="Minimum alignment length", default=30, type=int)
parser.add_argument("--output", help="Output table file", default="Match_Counts.tsv", type=str)
args = parser.parse_args()

def parse_m8(m8_files):
    scores = defaultdict(lambda: defaultdict(int))
    #Iterate over the blastn output files
    for blast_m8_file in m8_files:
        print ('Parsing BLASTN output',blast_m8_file)
        seen_queries = defaultdict(dict)
        #Iterare over queries
        for qresult in SearchIO.parse(blast_m8_file, 'blast-tab'):
            #Iterate over hits in the query
            for hit in qresult.hits:
                #Iterate over HSPs in the hits
                for hsp in hit.hsps:
                    #Check if the HSP passes established blast cutoffs
                    is_valid = check_match_cutoff(hsp,args.max_evalue,args.min_bitscore,args.min_identity,args.min_alignment)
                    if (is_valid):
                        #Only go on if the query/hit pair has not been processed before
                        if (qresult.id not in seen_queries.keys()):
                            seen_queries[qresult.id] = 1
                            scores[blast_m8_file][hit.id] += 1
                            #initialize variables in the scores dictionary if they are not there already
    return(scores)
                            
def check_match_cutoff(hsp,max_evalue,min_bitscore,min_ident,min_ali):
    if ((hsp.bitscore < min_bitscore) or (hsp.evalue > max_evalue) or (hsp.ident_pct < min_ident) or (hsp.aln_span < min_ali)):
        return False
    else:
        return True

scores = parse_m8(args.m8_files)

info_dataframe = pd.DataFrame.from_dict(scores)
info_dataframe.index.name = "File"
info_dataframe.to_csv(args.output,sep="\t",na_rep=0)