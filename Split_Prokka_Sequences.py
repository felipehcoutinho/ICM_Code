from collections import defaultdict
from Bio import SeqIO
import pandas as pd
import argparse
import re

parser = argparse.ArgumentParser()
parser.add_argument("--directories", help="Directories containing Prokka output files", nargs="+", type=str)

args = parser.parse_args()

def get_prefix(file):
    word_list = file.split("/")
    prefix_file = word_list[(len(word_list) - 2)]
    #print(f"Prefix is {prefix_file}")
    prefix_file = re.sub("Prokka_","",prefix_file)
    return prefix_file

for directory in args.directories:  
    print(f"Splitting sequences from {directory} according to Prokka annotation")
    directory_prefix = get_prefix(directory)
    info_table = directory+directory_prefix+".tsv"
    print(f'Reading info from {info_table}')
    info_data_frame = pd.read_csv(info_table, sep="\t",index_col="locus_tag",header=0)
    gene_file = directory+directory_prefix+".ffn"
    print(f"Processing {gene_file}")
    for seqobj in SeqIO.parse(gene_file, "fasta"):
        categ = info_data_frame.loc[seqobj.id]['ftype']
        product = info_data_frame.loc[seqobj.id]['product']
        product = re.sub("\s","_",product)
        split_file_name = f"Prokka_All_{categ}_{product}.fna"
        with open(split_file_name, 'a', newline='') as OUT:
            SeqIO.write(seqobj, OUT, "fasta")
            
    prot_file = directory+directory_prefix+".faa"
    print(f"Processing {prot_file}")
    for seqobj in SeqIO.parse(prot_file, "fasta"):
        categ = info_data_frame.loc[seqobj.id]['ftype']
        #product = info_data_frame.loc[seqobj.id]['product']
        #product = re.sub("\s","_",product)
        split_file_name = f"Prokka_All_{categ}.faa"
        with open(split_file_name, 'a', newline='') as OUT:
            SeqIO.write(seqobj, OUT, "fasta")
