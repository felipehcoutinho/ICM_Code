#! /usr/bin/env python3

from collections import defaultdict
import pandas as pd
import argparse
import re
import numpy

parser = argparse.ArgumentParser()
parser.add_argument("--seq_info_table", type=str)
parser.add_argument("--metab_path", type=str)
parser.add_argument("--output", default='Seq_Info+AMG_Metab_Pathway_Counts.tsv', type=str)
args = parser.parse_args()


path_df = pd.read_csv(args.metab_path, sep='\t',index_col='KEGG Entry',header=0)


seq_df = pd.read_csv(args.seq_info_table, sep='\t',index_col='Scaffold',header=0,na_values="DUMMY")

seq_df['Metabolisms'] = numpy.empty((len(seq_df), 0)).tolist()
seq_df['Pathways'] = numpy.empty((len(seq_df), 0)).tolist()

print(seq_df.index)

print(seq_df.columns)

quote_match_obj = re.compile("\'")
bracket_match_obj = re.compile("(\[)|(\])")

for scaffold,row in seq_df.iterrows():
    if (numpy.isnan(row['AMG_Count']) == False):
        amg_string = row['AMG_List']
        #print(scaffold,amg_string)
        amg_string = re.sub(quote_match_obj,"",amg_string)
        amg_string = re.sub(bracket_match_obj,"",amg_string)
        #print(scaffold,amg_string)
        amg_list = amg_string.split(' ')
        print(scaffold,amg_list)
        for amgko in amg_list:
            for entry,rowb in path_df.iterrows():
                entry_ko_string = rowb['Present AMG KOs']
                entry_ko_list = entry_ko_string.split(',')
                if (amgko in entry_ko_list):
                    print('Matched',amgko,'to',entry)
                    seq_df.loc[scaffold,'Metabolisms'].append(rowb['Metabolism'])
                    seq_df.loc[scaffold,'Pathways'].append(rowb['Pathway'])


seq_df.to_csv(args.output,sep="\t",na_rep='NA')