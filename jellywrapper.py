#! /usr/bin/env python3

from collections import defaultdict
import pandas as pd
import argparse
import subprocess
import re
import glob

parser = argparse.ArgumentParser()
parser.add_argument("--metagenomes_dir", help="Directories fastq files to be used for abundance calculation",  nargs="+", type=str)
parser.add_argument("--metagenomes_extension", help="Extension of the fastq files in metagenomes_dir", default="fastq.gz", type=str)
parser.add_argument("--kmer_length", help="Kmer Length", default="fastq.gz", type=str)
parser.add_argument("--threads", help="Threads to use", default=1, type=int)
parser.add_argument("--parse_only", help="Flag to skip running any programs and only parse their output", default=False, type=bool)
args = parser.parse_args()


def central():
	merged_kmer_df = pd.DataFrame({'Dummy' : []})
	samples_index = index_samples(args.metagenomes_dir,args.metagenomes_extension)
	for sample in samples_index.keys():
		if (samples_index[sample]['Is_Paired'] == True):
			r1_file = samples_index[sample]['R1']
			r2_file = samples_index[sample]['R2']
			print(f'Processing paired reads of {sample}',r1_file,r2_file)
			if (args.parse_only == False):
				command = f'jellyfish count --mer-len {args.kmer_length} --size 1M --threads {args.threads} --Files 2 --canonical -o {sample}.jf <(zcat {r1_file}) <(zcat {r2_file})'
				subprocess.call(command, executable='/bin/bash', shell=True)
				command = f'jellyfish dump -c --tab {sample}.jf > Kmer_Count_{sample}.tsv'
				subprocess.call(command,  executable='/bin/bash', shell=True)
				command = f'rm -f {sample}.jf'
				#subprocess.call(command,  executable='/bin/bash', shell=True)
			kmer_df = pd.read_csv(f'Kmer_Count_{sample}.tsv', sep='\t',index_col=0,header=None)
			kmer_df.columns = [sample]
			frames = [merged_kmer_df,kmer_df]
			merged_kmer_df = pd.concat(frames,axis=1)
		else:
			#print(f'{sample} is unpaired and will be skipped')
			pass

	merged_kmer_df.drop(['Dummy'], axis=1, inplace=True)
	merged_kmer_df.to_csv('Kmer_Counts.tsv',sep="\t",na_rep=0,index=True,index_label='Kmer')

def index_samples(metagenomes_dir,metagenomes_extension):
	samples_index = defaultdict(dict)
	for dir in metagenomes_dir:
		files = glob.glob(f'{dir}/*{metagenomes_extension}')
		suffix = '_split_t_paired'
		#print(f'File suffix inferred from {files[0]}: ',suffix)
		
		for file in files:
			id = file
			if (re.search(f'1{suffix}(\.)+{metagenomes_extension}$',id)):
				pair = id
				id = re.sub('(.)+/','',id)
				id = re.sub(f'1{suffix}','',id)
				id = re.sub(f'(\.)+{metagenomes_extension}','',id)
				pair = re.sub(f'1{suffix}',f'2{suffix}',pair)
				samples_index[id]['R1'] =  file
				if (pair in files):
					samples_index[id]['Is_Paired'] = True
				else:
					samples_index[id]['Is_Paired'] = False
				print('R1',id,file,pair,'Is Paired = ',samples_index[id]['Is_Paired'])
			elif (re.search(f'2{suffix}(\.)+{metagenomes_extension}$',id)):
				pair = id
				id = re.sub('(.)+/','',id)
				id = re.sub(f'2{suffix}','',id)
				id = re.sub(f'(\.)+{metagenomes_extension}','',id)
				pair = re.sub(f'2{suffix}',f'1{suffix}',pair)
				samples_index[id]['R2'] =  file
				if (pair in files):
					samples_index[id]['Is_Paired'] = True
				else:
					samples_index[id]['Is_Paired'] = False
				print('R2',id,file,pair,'Is Paired = ',samples_index[id]['Is_Paired'])
			else:
				id = re.sub(f'(\.)+{metagenomes_extension}','',id)
				#print('Unpaired',id,file)
				samples_index[id]['R1'] =  file
				samples_index[id]['Is_Paired'] = False
				
	return(samples_index)

central()