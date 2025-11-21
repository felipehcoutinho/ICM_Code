from Bio import Entrez
from Bio import SeqIO
import argparse
import pandas as pd
from collections import defaultdict

parser = argparse.ArgumentParser()
parser.add_argument("--list", help="The file with the list of sequence identifiers to be downloaded")
parser.add_argument("--output", help="The file to which the records will be written to", default="Sequences")
parser.add_argument("--email", help="The email to be used for Entrez queries")
parser.add_argument("--database", help="The database from which to retrieve the queries")
parser.add_argument("--format", help="The file format in which the queries should be retrieved", default="gb")
args = parser.parse_args()

#nothing is in this script seems even close to working

def read_list_file(input_file):
	print ('Reading file',input_file)
	with open(input_file) as input:
		return input.read().split('\n')


#Entrez.email = args.email
Entrez.email = "fhernandes@icm.csic.es"
acc_nums = read_list_file(args.list)
seq_info = defaultdict(dict)


handle = Entrez.efetch(db="genome", id="ASM236995", rettype="gb", retmode="text")
print(handle.readline().strip())
# acc_nums = ["GCF_002369955.1"]
# print(f"Fetching records and printing output sequences to {args.output}")
# with open(args.output, 'w', newline='') as OUT:
# 	try:
# 		handle = Entrez.efetch(db="nucleotide", id="186972394", rettype="gb", retmode="text")
# 		for SeqRec in SeqIO.read(handle,"fasta"):
# 			print (f"ID: {SeqRec.id}")
# 			SeqIO.write(SeqRec,OUT,"fasta")
# 	except:
# 		print (f"Error fetching 186972394")

# print(f"Fetching records and printing output sequences to {args.output}")
# with open(args.output, 'w', newline='') as OUT:
# 	for acc_num in acc_nums:
# 		try:
# 			handle = Entrez.efetch(db=args.database, id=acc_nums, rettype=args.format, retmode="text")
# 			#handle = Entrez.efetch(db="nucleotide", id="186972394", rettype="gb", retmode="text")
# 			for SeqRec in SeqIO.read(handle,args.format):
# 				print (f"ID: {SeqRec.id}")
# 				# print (f"Description: {SeqRec.description}")
# 				#print (f"Annotations: {SeqRec.annotations}")
# 				annot_dict = SeqRec.annotations
# 				SeqIO.write(SeqRec,OUT,args.format)
# 				for key in annot_dict:
# 					seq_info[key][SeqRec.id] = annot_dict[key]
# 		except:
# 			print (f"Error fetching {acc_num}")
# 		continue
		

print(f"Printing record info to Fetched_Seq_Info.tsv")
info_dataframe = pd.DataFrame.from_dict(seq_info)
info_dataframe.index.name = 'Sequence'
info_dataframe.to_csv("Fetched_Seq_Info.tsv",sep="\t",na_rep='NA')