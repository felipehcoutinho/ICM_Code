from Bio import Entrez
from Bio import SeqIO
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--list", help="The file with the list of sequence identifiers to be downloaded")
parser.add_argument("--output", help="The file to which the records will be written to", default="Sequences")
parser.add_argument("--email", help="The email to be used for Entrez queries")
parser.add_argument("--database", help="The database from which to retrieve the queries")
parser.add_argument("--format", help="The file format in which the queries should be retrieved", default="gb")
args = parser.parse_args()


def read_list_file(input_file):
	print ('Reading file',input_file)
	with open(input_file) as input:
		return input.read().split('\n')


Entrez.email = args.email
acc_nums = read_list_file(args.list)

with open(args.output, 'w', newline='') as OUT:
	for acc_num in acc_nums:
		handle = Entrez.efetch(db=args.database, id=acc_num, rettype=args.format, retmode="text")
		SeqRec = SeqIO.read(handle,args.format)
		SeqIO.write(SeqRec,OUT,args.format)



