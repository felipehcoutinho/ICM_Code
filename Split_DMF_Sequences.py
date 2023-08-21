from collections import defaultdict
from Bio import SeqIO
import pandas as pd
import argparse
import re

parser = argparse.ArgumentParser()
parser.add_argument("--min_score", help="Minimum score", default=0, type=float)
parser.add_argument("--assemblies_files", help="Fasta files containing scaffolds", nargs="+", type=str)
parser.add_argument("--info_table", help=".tsv file with scaffold scores from DeepMicrobeFinder", type=str)

args = parser.parse_args()

def get_prefix(file,file_format):
    prefix_file = re.sub(f".{file_format}","",file)
    prefix_file = re.sub("(.)+/","",prefix_file)
    return prefix_file

print(f'Reading info from {args.info_table}')
info_data_frame = pd.read_csv(args.info_table, sep="\t",index_col="Sequence Name",header=0)
info_data_frame['Highest_Score_Class'] = info_data_frame.idxmax(axis=1)
info_data_frame['Highest_Score'] = info_data_frame.max(axis=1)
print(info_data_frame.describe())

for assembly_file in args.assemblies_files:
    print(f"Splitting sequences from {assembly_file} according to DMF classification")
    assembly_file_prefix = get_prefix(assembly_file,"fasta")
    split_assemblies_files = []
    for seqobj in SeqIO.parse(assembly_file, "fasta"):
        taxon = info_data_frame.loc[seqobj.id]['Highest_Score_Class']
        high_score = info_data_frame.loc[seqobj.id]['Highest_Score']
        if (high_score >= args.min_score):
            split_file_name = f"Split_{assembly_file_prefix}_{taxon}.fasta"
            if (split_file_name not in split_assemblies_files):
                split_assemblies_files.append(split_file_name)
                #print("New split file:",split_file_name)
            with open(split_file_name, 'a', newline='') as OUT:
                SeqIO.write(seqobj, OUT, "fasta")
