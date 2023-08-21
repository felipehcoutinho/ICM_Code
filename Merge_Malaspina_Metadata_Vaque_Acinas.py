#! /usr/bin/env python3

from collections import defaultdict
import pandas as pd
import argparse
import re

acinas_df = pd.read_csv('/mnt/lustre/scratch/fcoutinho/Malaspina_Virus/Metadata/Malaspina_Table_S1+Lara+Sebastian2.tsv', sep='\t',index_col='Sample',header=0)

vaque_df = pd.read_csv('/mnt/lustre/scratch/fcoutinho/Malaspina_Virus/Metadata/Vaque_Malaspina_Metadata.tsv', sep='\t',index_col=None,header=0)


vals_to_try = [-15,-14,-11,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,0]
match_count = 0
for column in vaque_df.columns:
	if column not in acinas_df.columns:
		acinas_df[column] = 'NA'
		for v,row in vaque_df.iterrows():
			v_depth = row['DEPTH']
			v_station = row['STATION']
			for val in vals_to_try:
				match_depth = val + v_depth
				for buid,brow in acinas_df.iterrows():
					acinas_station = brow['station']
					acinas_depth = brow['Depth'] * -1
					if ((acinas_station == v_station) and (acinas_depth == match_depth)):
						print('Match! Station',acinas_station,'Depth',match_depth,'Variable',column,'Value',row[column],'Malaspina Sample',buid)
						acinas_df.loc[buid,column] = row[column]

			

print('Matched',match_count)				
acinas_df.to_csv('Malaspina_Table_S1+Lara+Sebastian2+Vaque.tsv',sep="\t",na_rep='NA')