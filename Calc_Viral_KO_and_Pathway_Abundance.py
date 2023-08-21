#! /usr/bin/env python3

from collections import defaultdict
import pandas as pd
import argparse
import re

parser = argparse.ArgumentParser()
parser.add_argument("--abundance", help="Genomic Sequence Abundance matrix", type=str)
parser.add_argument("--ko_to_pathway", help="ko_to_pathway table generated by VIBRANT",  type=str)
parser.add_argument("--amg", help="Table with the list of AMGs in each scaffold. Each genomic sequence takes a single line of the table with a list of KOs in the AMG_List column. Format example: ['K09474', 'K21480', 'K00643', 'K09882']", type=str)
parser.add_argument("--vibrant_amg", help="AMG Individuals table generated by VIBRANT", type=str)
args = parser.parse_args()

abund_df = pd.read_csv(args.abundance, sep='\t',index_col='Sequence',header=0)

kopath_df = pd.read_csv(args.ko_to_pathway, sep='\t',index_col='KEGG Entry',header=0)

ko_path = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))

for map_id,row in kopath_df.iterrows():
    ko_list = kopath_df.loc[map_id,'Present AMG KOs'].split(',')
    path = kopath_df.loc[map_id,'Pathway']
    metab = kopath_df.loc[map_id,'Metabolism']
    #print(metab,ko_list)
    for ko in ko_list:
        ko_path[ko]['Metabolisms'][metab] = 1
        ko_path[ko]['Pathways'][path] = 1


ko_abund = defaultdict(lambda: defaultdict(int))
met_abund = defaultdict(lambda: defaultdict(int))
path_abund = defaultdict(lambda: defaultdict(int))

if (args.amg):
    info_df = pd.read_csv(args.amg, sep='\t',index_col='Scaffold',header=0)
    for scaffold,row in abund_df.iterrows():
        amgs_string = str(info_df.loc[scaffold,'AMG_List'])
        if (amgs_string != 'nan'):
            amgs_string = re.sub("'","",amgs_string)
            amgs_string = re.sub("\[","",amgs_string)
            amgs_string = re.sub("\]","",amgs_string)
            amgs_string = re.sub("\s","",amgs_string)
            amgs_list = amgs_string.split(',')
            #print(scaffold,amgs_list,len(amgs_list))
            for ko in amgs_list:
                metab_list = ko_path[ko]['Metabolisms'].keys()
                metab_count = len(metab_list)
                path_list = ko_path[ko]['Pathways'].keys()
                path_count = len(path_list)
                #print('Processing KO',ko,'Metabolisms:',metab_count,'Pathways:',path_count)
                for sample in abund_df.columns:
                    ko_abund[sample][ko] += abund_df.loc[scaffold,sample]
                    for metab in metab_list:
                        met_abund[sample][metab] += (abund_df.loc[scaffold,sample] / metab_count)
                    for path in path_list:
                        path_abund[sample][path] += (abund_df.loc[scaffold,sample] / path_count)
elif(args.vibrant_amg):
    info_df = pd.read_csv(args.vibrant_amg, sep='\t',index_col='protein',header=0)
    for protein,row in info_df.iterrows():
        ko = row['AMG KO']
        scaffold = row['scaffold']
        metab_list = ko_path[ko]['Metabolisms'].keys()
        metab_count = len(metab_list)
        path_list = ko_path[ko]['Pathways'].keys()
        path_count = len(path_list)
        for sample in abund_df.columns:
            ko_abund[sample][ko] += abund_df.loc[scaffold,sample]
            for metab in metab_list:
                met_abund[sample][metab] += (abund_df.loc[scaffold,sample] / metab_count)
            for path in path_list:
                path_abund[sample][path] += (abund_df.loc[scaffold,sample] / path_count)
else:
    print("Must provide --amg or --vibrant_amg!")

ko_abund_df = pd.DataFrame.from_dict(ko_abund)
ko_abund_df.index.name = 'KO'
ko_abund_df.to_csv('Viral_KO_Abundance.tsv',sep="\t",na_rep=0)

met_abund_df = pd.DataFrame.from_dict(met_abund)
met_abund_df.index.name = 'Metabolism'
met_abund_df.to_csv('Viral_Metabolism_Abundance.tsv',sep="\t",na_rep=0)

path_abund_df = pd.DataFrame.from_dict(path_abund)
path_abund_df.index.name = 'Pathway'
path_abund_df.to_csv('Viral_Pathway_Abundance.tsv',sep="\t",na_rep=0)