#! /usr/bin/env python3

from collections import defaultdict
import pandas as pd
import glob

length_info = dict()
map_info = dict()
raw_abund_info = defaultdict(dict)
rpkm_abund_info = defaultdict(dict)

counts_files = glob.glob('*Counts.tsv')

for file in counts_files:
    print(f"Processing {file}")
    sample = file.split("x")[0]
    with open(file, 'r') as IN:
        for line in IN.readlines():
            vals = line.split('\t')
            if (vals[0] != "*"):
                (scaffold,length,read_count,dummy) = vals
                length_info[scaffold] = int(length)
                raw_abund_info[sample][scaffold] = int(read_count)
                map_info[sample] = (int(read_count) + map_info.get(sample,0))
                

            
for sample in map_info:
    print(f"Calculating RPKM for {sample}")  
    for scaffold in length_info:
        #print("Scaffold",scaffold,"Raw abundance",raw_abund_info[sample][scaffold],"Total mapped reads",map_info[sample])
        rpkm_abund_info[sample][scaffold] = ((((raw_abund_info[sample][scaffold] / length_info[scaffold]) * 1000) / map_info[sample]) * 1000000)
        
rpkm_matrix = pd.DataFrame.from_dict(rpkm_abund_info)
rpkm_matrix.index.name = 'Sequence'
rpkm_matrix.to_csv("RPKM.tsv",sep="\t",na_rep='NA')