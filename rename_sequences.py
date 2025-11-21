from Bio import SeqIO
import re
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--input", help="File of genomic sequences", type =str, required=True)
#parser.add_argument("--input", help="File of genomic sequences", type =str, required=True)
parser.add_argument("--string_rename", help="String to use when renaming genomic sequences", default='_Seq_', type=str)
args = parser.parse_args()


prefix_file = re.sub('(.)+/','',args.input)
out_seq_file = "Renamed_"+prefix_file

OUT = open(out_seq_file,'w', newline='')
seq_counter = 0
for seqobj in SeqIO.parse(args.input, "fasta"):
    seq_counter += 1
    new_id = args.string_rename+"_Seq_"+str(seq_counter)
    seqobj.id = new_id
    SeqIO.write(seqobj, OUT, "fasta")
OUT.close()
    