from collections import defaultdict
import pandas as pd
import argparse
import re

parser = argparse.ArgumentParser()
parser.add_argument("--assembly_stats", help=".tsv Tables containing Bin/MAG assembly data", nargs="+", type=str)
parser.add_argument("--checkm_stats", help=".tsv Tables containing Bin/MAG CheckM data", nargs="+", type=str)
parser.add_argument("--gtdbtk_stats", help=".tsv Tables containing Bin/MAG GTDBtk data", nargs="+", type=str)
parser.add_argument("--mag_info_table", help=".tsv file to write MAG info", default="MAG_Info.tsv", type=str)

args = parser.parse_args()

bin_info = defaultdict(dict)

all_gtdbtk_dfs = []
    
for gtdbtk_file in args.gtdbtk_stats:
    print(f"Reading {gtdbtk_file}")
    gtdbtk_df = pd.read_csv(gtdbtk_file, sep="\t",index_col="user_genome",header=0)
    all_gtdbtk_dfs.append(gtdbtk_df)

info_df = pd.concat(all_gtdbtk_dfs,axis=0)
info_df.index.name = 'Genome'

#print(info_df.describe())

rank_dict = dict(d = "Domain" ,  p = "Phylum" , c = "Class" , o = "Order" , f = "Family" , g= "Genus" , s = "Species")

for rank in rank_dict:
    info_df[rank_dict[rank]] = "NA"

for genome,row in info_df.iterrows():
    #print(genome,row)
    full_tax = row["classification"].split(";")
    for level in full_tax:
        (rank,taxon) = level.split("__")
        #print(rank,taxon)
        if (taxon != ""):
            info_df[rank_dict[rank]][genome] = taxon



info_df.to_csv(args.mag_info_table,sep="\t",na_rep='NA')