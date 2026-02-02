#! /usr/bin/env python3
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
#from ete2 import NCBITaxa
import pandas as pd
import argparse
import re


parser = argparse.ArgumentParser()
parser.add_argument("--aln_files", help="Files of aligned  sequences", type =str, nargs="+")
parser.add_argument("--in_format", help="Input Sequence file format", default="fasta", type =str)
parser.add_argument("--fasta_output", help="The output fasta file with the concatenatd alingments", default="Concatenated.fasta", type =str)
parser.add_argument("--info_genome_output", help="The Genome info output table file to be generated ", default="Nucleotide_Info.tsv", type =str)
parser.add_argument("--min_prevalence", help="Minimum number of sequences a concatenate must include to be printed to the output fasta file", default=1, type=int)
args = parser.parse_args()

def index_seqs(in_seq_files=[]):
    for seq_file in in_seq_files:
        ext_rm = re.compile("\..+$")
        out_tbl_file = re.sub(ext_rm,"",seq_file)
        pref_rm = re.compile("^.*/")
        out_tbl_file = re.sub(pref_rm,"",out_tbl_file)
        out_id =  "Consensus_"+out_tbl_file
        out_tbl_file = out_tbl_file + "_Base_Counts_per_Position.tsv"
        out_seq_file = out_id+".fasta"
        OUTSEQ = open(out_seq_file,'w', newline='')
        print ("Analysing",seq_file)
        consen_dict = defaultdict(dict)
        # seen_cds_ids = dict()
        # seen_scaff_ids = dict()
        seq_counter = 0
        for seqobj in SeqIO.parse(seq_file, "fasta"):
            # cds_id = seqobj.id
            # print ("\t Seq:",seqobj.id)
            seq_length = len(seqobj.seq)
            for i in range(1,seq_length,1):
                nuc = seqobj.seq[i]
                # print ("\t Position:",i)
                # if (nuc not in consen_dict[i]):
                #     consen_dict[i][nuc] = 0
                # consen_dict[i][nuc] += 1
                if (i not in consen_dict[nuc]):
                    consen_dict[nuc][i] = 0
                consen_dict[nuc][i] += 1
    
        print(f"\tPrinting consensus info to {out_tbl_file}")
        concat_info_df = pd.DataFrame.from_dict(consen_dict)
        #concat_info_df = concat_info_df.transpose()
        concat_info_df.index.name = 'Position'
        concat_info_df = concat_info_df.drop(columns='-')
        # print(concat_info_df.columns)
        concat_info_df["Consensus"] = concat_info_df.idxmax(axis="columns")
        concat_info_df.to_csv(out_tbl_file,sep="\t",na_rep=0)
        concatamer= "".join(concat_info_df["Consensus"])
        print("\tWriting consensus to",out_seq_file,"Length of concatatmer: ",len(concatamer))
        seqobj = SeqRecord(id = out_id, seq=Seq(concatamer), description ="")
        SeqIO.write(seqobj, OUTSEQ, "fasta")

index_seqs(in_seq_files=args.aln_files)

