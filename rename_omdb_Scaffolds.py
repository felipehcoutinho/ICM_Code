from Bio import SeqIO
import re
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--input", help="File of genomic sequences", type =str, required=True)
args = parser.parse_args()


prefix_file = re.sub('(.)+/','',args.input)
out_seq_file = "Renamed_"+prefix_file

OUT = open(out_seq_file,'w', newline='')

for seqobj in SeqIO.parse(args.input, "fasta"):
    new_id = re.sub(r"(.+rep=)|(;.+$)","",seqobj.description)
    #print(f"OID {seqobj.id} NID {new_id}")
    seqobj.description = "oid="+str(seqobj.id)+";"+str(re.sub(r"^.+\s","",seqobj.description))
    seqobj.id = new_id
    SeqIO.write(seqobj, OUT, "fasta")
OUT.close()
    