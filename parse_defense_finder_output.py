#! /usr/bin/env python3

from collections import defaultdict
import pandas as pd
import argparse
import subprocess
import re
import glob
import os


parser = argparse.ArgumentParser()
parser.add_argument("--output_table", help="The output .tsv format table file generated", default="DenfenseFinder_Results.tsv", type =str)
parser.add_argument("--input_directory", help="Directory containing DefenseFinder output directories", type=str)
args = parser.parse_args()

list_dfinder_tables = glob.glob(f"{args.input_directory}/**/*defense_finder_systems.tsv", recursive=True)

genome_info = defaultdict(lambda: defaultdict(int))

for dfinder_tbl_file in list_dfinder_tables:
    print(f'Reading info from {dfinder_tbl_file}')
    dfinder_df = pd.read_csv(dfinder_tbl_file, sep="\t",index_col="sys_id",header=0)
    for sys_id,row in dfinder_df.iterrows():
        (sys_type,subtype,sys_beg,sys_end,protein_in_syst,genes_count,name_of_profiles_in_sys) = row
        sys_type = "Defense_System|"+sys_type
        genome_id = re.sub("_(\d)+$","",sys_beg)
        genome_info[genome_id][sys_type] += 1
        if ("List_of_Unique_Genes_in_Systems" not in genome_info[genome_id]):
            genome_info[genome_id]["List_of_Unique_Genes_in_Systems"] = set()
            genome_info[genome_id]["List_of_Defense_SubSystems_in_Genome"] = set()
        genome_info[genome_id][sys_type] += 1
        genome_info[genome_id]["List_of_Defense_SubSystems_in_Genome"].add(f"{sys_type}|{subtype}")
        for unique_gene in protein_in_syst.split(','):
            genome_info[genome_id]["List_of_Unique_Genes_in_Systems"].add(unique_gene)
            
print("Calculating metrics")
for genome_id in genome_info.keys():
    genome_info[genome_id]["Count_of_Unique_Genes_in_Systems"] = len(genome_info[genome_id]["List_of_Unique_Genes_in_Systems"])
    genome_info[genome_id]["Count_of_Unique_Defense_Systems"] = len(genome_info[genome_id].keys()) - 2
    genome_info[genome_id]["List_of_Unique_Genes_in_Systems"] = ",".join(genome_info[genome_id]["List_of_Unique_Genes_in_Systems"])
    genome_info[genome_id]["List_of_Defense_SubSystems_in_Genome"] = ",".join(genome_info[genome_id]["List_of_Defense_SubSystems_in_Genome"])

print(f"Printing output to {args.output_table}")
genome_info_df = pd.DataFrame.from_dict(genome_info)
genome_info_df.index.name = 'MAG'
genome_info_df_transp = genome_info_df.T
genome_info_df_transp.index.name = 'MAG'
genome_info_df_transp.sort_index(axis=1, inplace=True)
genome_info_df_transp.index.name = 'MAG'
genome_info_df_transp.to_csv(args.output_table,sep="\t",na_rep=0)