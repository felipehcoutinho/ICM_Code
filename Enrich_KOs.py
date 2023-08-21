#! /usr/bin/env python3
import json
import pandas as pd
import argparse
from collections import defaultdict

parser = argparse.ArgumentParser()
parser.add_argument("--ko", help="Table containing KOs", type =str)
parser.add_argument("--json", help=".json file from KEGG", type =str)
args = parser.parse_args()

kegg_info =  defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(dict)))))

json_fh = open(args.json,)
kegg_data = json.load(json_fh)

valid_categs = ("09100 Metabolism","09120 Genetic Information Processing","09130 Environmental Information Processing","09140 Cellular Processe")

for category in kegg_data["children"]:
    print("Category",category['name'])
    categ_name = category['name']
    if (categ_name in valid_categs):
        for hier1 in category['children']:
            #print("\tLevel 1",hier1['name'])
            h1_name = hier1['name']
            if (hier1['children']):
                for hier2 in hier1['children']:
                    h2_name = hier2['name']
                    #print("\t\tLevel 2",hier2['name'])
                    if (hier2['children']):
                        for hier3 in hier2['children']:
                            #print("\t\t\tLevel 3",hier3['name'])
                            h3_name = hier3['name']
                            (ko,desc) = hier3['name'].split(' ',1)
                            #print("\t\t\tKO",ko,"Description",desc)
                            kegg_info[ko][categ_name][h1_name][h2_name][h3_name]["Match"] = True
                            kegg_info[ko][categ_name][h1_name][h2_name][h3_name]["Description"] = desc
  

json_fh.close()

ko_df = pd.read_csv(args.ko, sep="\t",index_col=0,header=0)

with open("KO_Enriched.tsv", 'w', newline='') as OUT:
    #OUT.write('Query_Scaffold\tSubject_Scaffold\tAAI\tMatched_CDS\tPerc_Matched_CDS\n')
    for i,row in ko_df.iterrows():
        ko = row['Query']
        #print(ko)
        if (ko in kegg_info.keys()):
            for categ in kegg_info[ko].keys():
                for h1_name in kegg_info[ko][categ].keys():
                    for h2_name in kegg_info[ko][categ][h1_name].keys():
                        for h3_name in kegg_info[ko][categ][h1_name][h2_name].keys():
                            desc = kegg_info[ko][categ][h1_name][h2_name][h3_name]["Description"]
                            OUT.write("\t".join(row.astype("string")))
                            OUT.write(f"\t{categ}\t{h1_name}\t{h2_name}\t{h3_name}\t{desc}\n")