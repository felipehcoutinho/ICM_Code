#! /usr/bin/env python3

import argparse
from collections import defaultdict
import pandas as pd
import glob
import re

parser = argparse.ArgumentParser()
parser.add_argument("--input_dir", help="Directory containing VIBRANT output directories", type=str)
parser.add_argument("--output", help="Output table .tsv file", type=str, default="VIBRANT_Merged_Info.tsv")
args = parser.parse_args()


#This is the main function that calls all the other functions according to the user specified parameters
def central():
	full_dataframe = pd.DataFrame({'A' : []})
	subdirs = glob.glob(f'{args.input_dir}/VIBRANT_results_*')
	for subdir in subdirs:
		#print('Processing',subdir)
		prefix = subdir
		prefix = re.sub(f"{args.input_dir}","",prefix)
		prefix = re.sub(f"VIBRANT_results_","",prefix)
		amg_file = 'VIBRANT_AMG_individuals_'+prefix+'.tsv'
		quality_file = 'VIBRANT_genome_quality_'+prefix+'.tsv'
		#print('\tAMG:',amg_file)
		#print('\tQuality:',quality_file)
		parse_vibrant(f'{subdir}/{quality_file}',f'{subdir}/{amg_file}')
	
	info_dataframe = pd.DataFrame.from_dict(seq_info)
	info_dataframe.index.name = 'Sequence'
	info_dataframe.to_csv(args.output,sep="\t",na_rep='NA')

def index_info(table_file,index_col_name):
	print(f'Reading info from {table_file}')
	info_data_frame = pd.read_csv(table_file,sep='\t',index_col=index_col_name,header=0)
	return info_data_frame

def collect_qual_info(row):
	#print(row)
	scaffold = row.name
	frag_match = re.search(compiled_obj,scaffold)
	clean_name = scaffold.split(' ')[0]
	#print(scaffold,clean_name)
	#print('Match:',frag_match)
	if (frag_match):
		pass
		#clean_name = clean_name+frag_match[0]
		#print('Matched a fragment!',i,frag_match,'Updated Clean Name:',clean_name)

	seq_info['type'][clean_name] = row['type']
	seq_info['Quality'][clean_name] = row['Quality']
	seq_info['Is_Virus'][clean_name] = True
	
def collect_amg_info(row):
	#print(row)
	protein = row.name
	scaffold = row['scaffold']
	frag_match = re.search(compiled_obj,scaffold)
	clean_name = scaffold.split(' ')[0]
	#print(scaffold,clean_name)
	#print('Match:',frag_match)
	if (frag_match):
		pass
		#clean_name = clean_name+frag_match[0]
		#print('Matched a fragment!',i,frag_match,'Updated Clean Name:',clean_name)

	if (clean_name not in seq_info['AMG_Count']):
		seq_info['AMG_Count'][clean_name] = 0
		seq_info['AMG_List'][clean_name] = set()
	
	seq_info['AMG_Count'][clean_name]+= 1
	seq_info['AMG_List'][clean_name].add(row['AMG KO'])

def parse_vibrant(vibrant_out_quality_file,vibrant_out_amg_file):
	if (vibrant_out_quality_file != 'NA'):
		vibrant_info_data_frame = index_info(vibrant_out_quality_file,'scaffold')
		vibrant_info_data_frame = vibrant_info_data_frame[vibrant_info_data_frame.Quality != 'complete circular']
		print('Processing',vibrant_out_quality_file)
		vibrant_info_data_frame.apply(collect_qual_info,axis='columns')
	
	if (vibrant_out_amg_file != 'NA'):
		vibrant_info_data_frame = index_info(vibrant_out_amg_file,'protein')
		print('Processing',vibrant_out_amg_file)
		vibrant_info_data_frame.apply(collect_amg_info,axis='columns')
	return 0

seq_info = defaultdict(dict)
compiled_obj = re.compile('_fragment_(\d)+$')
central()
