#! /usr/bin/env python3
import pandas as pd
import argparse
import re

parser = argparse.ArgumentParser()
parser.add_argument("--dictionary", help="List of matrixes to merge", type=str, required=True)
parser.add_argument("--index_dictionary", help="Column to use as the Index of the dictionary table", type=str, required=True)
#parser.add_argument("--nid_dictionary", help="Column to retrieve the New ID from the dictionary table", type=str, required=True)
parser.add_argument("--info", help="Seq Info table", type=str , required=True)
parser.add_argument("--index_info", help="Column to use as the Index of the info table", type=str, required=True)
parser.add_argument("--output", help="Output file", default="Renamed_IDs_Info.tsv", type=str)

args = parser.parse_args()


print(f'Reading from {args.dictionary}')
dict_df = pd.read_csv(args.dictionary, sep="\t",index_col=args.index_dictionary,header=0,low_memory=False)

oid2nid = dict()
print(f'Building dictionary')
for seq_id,row in dict_df.iterrows():
    desc = dict_df.loc[seq_id,"Description"]
    oid = dict_df.loc[seq_id,"Original_ID"]
    frag_match_obj = re.search("_fragment_(\d)+$",desc)
    if (frag_match_obj != None):
        oid = oid+frag_match_obj.group(0)
        #print(seq_id,oid,desc)
    oid2nid[oid] = seq_id
    

print(f'Reading from {args.info}')
info_df = pd.read_csv(args.info, sep="\t",index_col=args.index_info,header=0,low_memory=False)

print(f'Changing IDs')
info_df.rename(index=oid2nid,inplace=True)

print(f'Printing to {args.output}')
info_df.to_csv(args.output,sep="\t",na_rep='NA',index=True)