#! /usr/bin/env python3
import pandas as pd
import argparse


parser = argparse.ArgumentParser()
parser.add_argument("--matrixes", help="List of matrixes to merge", nargs="+", required=True)
parser.add_argument("--output", help="Output file", default="Merged_Matrixes.tsv", type=str)
parser.add_argument("--sep", help="Field Separator in the input matrixes", default="\t", type=str)
parser.add_argument("--index", help="Index variable", type=str)
parser.add_argument("--axis", help="Axis to perform the merge. 0 for rows 1 for colums", type=int, default=1)
args = parser.parse_args()

merged_data_frame = pd.DataFrame({'Dummy' : []})

for matrix_file in args.matrixes:
	print(f'Reading from {matrix_file}')
	if (args.index):
		data_frame = pd.read_csv(matrix_file, sep=args.sep,index_col=args.index,header=0,low_memory=False)
	else:
		data_frame = pd.read_csv(matrix_file, sep=args.sep,index_col=0,header=0,low_memory=False)
	
	frames = [merged_data_frame,data_frame]
	merged_data_frame = pd.concat(frames,axis=args.axis)

merged_data_frame.drop(['Dummy'], axis=1, inplace=True)
merged_data_frame = merged_data_frame.loc[:,~merged_data_frame.columns.duplicated()]

merged_data_frame.index.name = args.index
merged_data_frame.to_csv(args.output,sep="\t",na_rep='NA',index=True)

#if (args.index):
#	merged_data_frame.to_csv(args.output,sep="\t",na_rep='NA',index=True)
#else:
#	merged_data_frame.to_csv(args.output,sep="\t",na_rep='NA',index=False)