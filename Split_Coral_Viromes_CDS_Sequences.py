from collections import defaultdict
from Bio import SeqIO
import pandas as pd
import argparse
import re

parser = argparse.ArgumentParser()
parser.add_argument("--assemblies_files", help="Fasta files containing scaffolds", nargs="+", type=str)

args = parser.parse_args()

def get_prefix(file,file_format):
    prefix_file = re.sub(f".{file_format}","",file)
    prefix_file = re.sub("(.)+/","",prefix_file)
    return prefix_file

for assembly_file in args.assemblies_files:
    print(f"Splitting sequences from {assembly_file} according to ID classification")
    assembly_file_prefix = get_prefix(assembly_file,"fasta")
    split_assemblies_files = []
    for seqobj in SeqIO.parse(assembly_file, "fasta"):
        taxon = seqobj.id.rsplit("_",1)[0]
        #print(f"{seqobj.id} {taxon}")
        split_file_name = f"{taxon}.faa"
        if (split_file_name not in split_assemblies_files):
            split_assemblies_files.append(split_file_name)
                #print("New split file:",split_file_name)
        with open(split_file_name, 'a', newline='') as OUT:
            SeqIO.write(seqobj, OUT, "fasta")
