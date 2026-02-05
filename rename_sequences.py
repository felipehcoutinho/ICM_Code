from Bio import SeqIO
import re
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--input", help="File of genomic sequences", type =str, required=True)
parser.add_argument("--output", help="File of genomic sequences", type =str, required=True)
parser.add_argument("--string_rename", help="String to use when renaming genomic sequences", default='_Seq_', type=str)
parser.add_argument("--min_length", help="Minimum length of sequences to be include in the output file", default=0, type=int)
parser.add_argument("--max_length", help="Minimum length of sequences to be include in the output file", default=999999999, type=int)
args = parser.parse_args()


# prefix_file = re.sub('(.)+/','',args.input)
# out_seq_file = "Renamed_"+prefix_file

out_seq_file = args.output

OUT = open(out_seq_file,'w', newline='')
seq_counter = 0
for seqobj in SeqIO.parse(args.input, "fasta"):
    seq_length = len(seqobj.seq)
    seq_counter += 1
    new_id = args.string_rename+"_Seq_"+str(seq_counter)
    seqobj.id = new_id
    if ((seq_length >= args.min_length) and (seq_length <= args.max_length)):
        SeqIO.write(seqobj, OUT, "fasta")
OUT.close()
    