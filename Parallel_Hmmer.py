import multiprocessing as mp
from joblib import Parallel, delayed
import os
import subprocess
#import networkx
#import obonet
#import csv
#import re
import argparse
#from Bio.Seq import Seq
from Bio import SeqIO
#from Bio import SearchIO
#from collections import defaultdict

#Get input parameters from user
parser = argparse.ArgumentParser()
parser.add_argument("--threads", help="Number of threads to use per chunk", type =int, default=4)
parser.add_argument("--chunks", help="Number of chunks to split the query file into", type =int, default = 5)
parser.add_argument("--database", help="Hmmer database file", type =str)
parser.add_argument("--query", help="Query fasta file", type =str)
parser.add_argument("--program", help="Program to run: hmmscan (Default) or hmmsearch", type =str, default="hmmscan")

args = parser.parse_args()

#Split query will count the number of sequences in the query files, then proceed to split it into smaller chunks according to the specified number of threads
def split_query():
    #Count the number of seqs in the query file
    command = f'grep -c ">" {args.query}'
    query_seq_count = int(subprocess.getoutput(command))
    #Calculate the number of sequences in each chunk (split files)
    chunk_size = round(query_seq_count / args.chunks) + 1
    print('Splitting',args.query,'into',args.chunks,'chunks of',chunk_size,'sequences.')
    seq_counter = 0
    split_counter = 1
    split_name = 'Split_'+str(split_counter)+'.faa'
    split_names_list.append(split_name)
    split_handle = open(split_name,"w")
    #Iterate over all sequences in the query file. Keep track of the number of sequences added to each split file in seq_counter
    for seqobj in SeqIO.parse(args.query,format="fasta"):
        SeqIO.write(seqobj,split_handle, "fasta")
        seq_counter += 1
        #If we reached the number of sequences in the split file, close the handle and open a new split file handle and reset the seq_counter. Keep track of split file names in split_names_list
        if (seq_counter == chunk_size):
            split_counter += 1
            split_name = 'Split_'+str(split_counter)+'.faa'
            split_names_list.append(split_name)
            split_handle.close()
            split_handle = open(split_name,"w")
            seq_counter = 0
    split_handle.close()
        
#para_loop calls the call_hmmer function using multiprocessing for each split file
def para_loop():
    if __name__ == "__main__":
         results = Parallel(n_jobs=args.chunks)(delayed(call_hmmer)(i) for i in split_names_list)

#call_hmmer calls hmmscan or hmmsearch using a single cpu for the specified in_file (which will be a split fasta) and the args.tdabase file
def call_hmmer(in_file):
    command = f'{args.program} -o {in_file}xDB --noali --cpu {args.threads} {args.database} {in_file}'
    subprocess.call(command, shell=True)
    return 1

split_names_list = list()
split_query()
para_loop()
