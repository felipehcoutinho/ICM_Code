#! /usr/bin/env python3
from collections import defaultdict
# from Bio import SeqIO
# from Bio import SearchIO
# from Bio.Seq import Seq
# import pandas as pd
import argparse
# import subprocess
# import re
# import glob
# import os

parser = argparse.ArgumentParser()
parser.add_argument("--cds_list", help="List of CDS to include in output. If not provided all CDS are included if they pass the scoring criteria", type =str, default=None)
parser.add_argument("--input", help="Input m8 file", type =str, required=True)
parser.add_argument("--max_evalue", help="Maximum e-value", default=0.00001, type=float)
parser.add_argument("--min_identity", help="Minimum identity (0-100 for blast like output and 0-1 for mmseqs)", default=30.0, type=float)
parser.add_argument("--min_bitscore", help="Minimum bitscore", default=50, type=int)
parser.add_argument("--min_alignment", help="Minimum alignment length", default=30, type=int)
parser.add_argument("--output", help="Output table file", default="Filtered_m8.tsv", type=str)
args = parser.parse_args()

def parse_m8(m8_file):
    scores = defaultdict(lambda: defaultdict(int))
    #Iterate over the blastn output files
    use_list = False
    if (args.cds_list):
        use_list = True
    with open(args.output, 'w', newline='') as OUT:
        print ('Parsing m8 output',m8_file)
        seen_queries = defaultdict(dict)
        #Iterare over m8_file one line at a time using a filehandle
        with open(m8_file) as blast_m8:
            for line in blast_m8:
                #Split the line into fields
                fields = line.strip().split("\t")
                #Extract the query and subject IDs
                query_id = fields[0]
                subject_id = fields[1]
                #Extract the e-value, bitscore, identity, and alignment length
                evalue = float(fields[10])
                bitscore = float(fields[11])
                ident_pct = float(fields[2])
                aln_span = int(fields[3])
                #Check if the match passes the established cutoffs
                if ((bitscore >= args.min_bitscore) or (evalue <= args.max_evalue) or (ident_pct >= args.min_identity) or (aln_span >= args.min_alignment)):
                    if ((use_list == False) or (query_id in valid_cds)):
                       print (line.strip(), file=OUT) 
                         

valid_cds = dict()
if (args.cds_list):
    print ('CDS list provided. Only CDS in the list will be included in the output')
    with open(args.cds_list) as CDS_LIST:
        for line in CDS_LIST:
            cds = line.strip()
            valid_cds[cds] = 1
    print("Obtained",len(valid_cds),"unique CDS IDs from the list")
        
parse_m8(args.input)

