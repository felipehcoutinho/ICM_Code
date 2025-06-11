#This script scans the sequences of protein encoding genes for start and stop codons
#! /usr/bin/env python3
#Virathon: Genomic Analysis of Viruses of Archaea and Bacteria
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--gene_files", help="File of CDS DNA sequences", type =str, nargs="+")
args = parser.parse_args()

valid_start_codons = ["ATG", "GTG", "TTG"]
valid_stop_codons = ["TAA", "TAG", "TGA"]

# for i in range(1,15-3,3):
#     print(i)
# exit()

for seq_file in args.gene_files:
    #print ("Indexing sequences from",seq_file)
    #Iterate over sequences in the file. Collect basic Info
    for seqobj in SeqIO.parse(seq_file,"fasta"):
        seq_length = len(seqobj.seq)
        seq_id = seqobj.id
        #Check first and last codons
        if (seqobj.seq[0:3] not in valid_start_codons):
            print("Warning: Sequence\t",seq_id,"\tdoes not start with a valid start codon:",seqobj.seq[0:3])
        if (seqobj.seq[-3:] not in valid_stop_codons):
            print("Warning: Sequence\t",seq_id,"\tdoes not end with a valid stop codon:",seqobj.seq[-3:])
        #Check each codon in the middle for internal stop codons
        codon_counter = 0
        for i in range(1,seq_length-3,3):
            codon = seqobj.seq[i-1:i+2]
            codon_counter += 1
            if (codon in valid_stop_codons):
                print("Warning: Sequence\t",seq_id,"\t of length\t",seq_length,"\tcontains an internal stop codon\t:",codon,"codon number:\t",codon_counter,"\tat positions",i,"-",i+2)
