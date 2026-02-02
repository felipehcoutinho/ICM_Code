from Bio import SeqIO
import re
import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--input", help="File of genomic sequences", type =str, required=True)
parser.add_argument("--dictionary", help="two column tsv file with no header with instructions on how to rename contigs. first column with original name and second column with new names", type =str, required=True)
parser.add_argument("--string_rename", help="String to use when renaming genomic sequences", default='_Seq_', type=str)
args = parser.parse_args()


if args.dictionary:
    names_df = pd.read_csv(args.dictionary, sep="\t",index_col=0,header=0)


prefix_file = re.sub('(.)+/','',args.input)
out_seq_file = "Renamed_"+prefix_file

OUT = open(out_seq_file,'w', newline='')
seq_counter = 0
for seqobj in SeqIO.parse(args.input, "fasta"):
    seq_counter += 1
    new_id = None 
    if (args.dictionary):
        #using seqobj.id as the idnex, get the value of the first column from names_df and assing it to new_id
        try:
            new_id = names_df.loc[seqobj.id].values[0]
        except KeyError:
            new_id = args.string_rename+"_Seq_"+str(seq_counter)
    else:
        new_id = args.string_rename+"_Seq_"+str(seq_counter)
    seqobj.id = new_id
    SeqIO.write(seqobj, OUT, "fasta")
OUT.close()
    