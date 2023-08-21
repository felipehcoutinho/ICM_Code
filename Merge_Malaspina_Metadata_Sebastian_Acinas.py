#! /usr/bin/env python3

from collections import defaultdict
import pandas as pd
import argparse
import re

acinas_df = pd.read_csv('/mnt/lustre/scratch/fcoutinho/Malaspina_Virus/Metadata/Malaspina_Table_S1+Lara+Sebastian.tsv', sep='\t',index_col='Sample',header=0)

seb_df = pd.read_csv('/mnt/lustre/scratch/fcoutinho/Malaspina_Virus/Metadata/envTable_withKmeansRatioLR.tsv', sep='\t',index_col='MPcode',header=0)



acinas_df['cluster'] = 'NA'

match_count = 0

for uid,row in seb_df.iterrows():
	seb_station = row['station']
	for buid,brow in acinas_df.iterrows():
		acinas_station = brow['station']
		if (seb_station == acinas_station):
			match_count += 1
			acinas_df.loc[buid,'cluster'] = seb_df.loc[uid,'cluster']
			

print('Matched',match_count)				
acinas_df.to_csv('Malaspina_Table_S1+Lara+Sebastian2.tsv',sep="\t",na_rep='NA')