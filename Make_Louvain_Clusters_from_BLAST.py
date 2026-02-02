#! /usr/bin/env python3
from collections import defaultdict
from Bio import SearchIO
import pandas as pd
import argparse
import re
import networkx as nx

parser = argparse.ArgumentParser()
parser.add_argument("--blast_input", help="M8 format file of blast hits", type =str)
parser.add_argument("--output", default="Louvain_Communities.tsv", help="File to write community assignments", type =str)
args = parser.parse_args()

#graph = nx.read_gml(args.blast_input)

#G = nx.petersen_graph()

graph = nx.read_weighted_edgelist(args.blast_input, comments='#', delimiter="\t")

lv_comms = nx.community.louvain_communities(graph, seed=666)

comms_df = pd.DataFrame(lv_comms) #, columns=['Communities']

print(f"Printing Community Info to  {args.output}")
#comms_df = pd.DataFrame.from_dict(cds_seq_info)
#comms_df.index.name = 'Sequence'
comms_df.to_csv(args.output,sep="\t",na_rep='NA')

