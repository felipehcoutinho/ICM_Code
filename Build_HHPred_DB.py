#! /usr/bin/python3
import os
import pandas as pd
import subprocess
import re
import argparse
from Bio.Seq import Seq
from Bio import SeqIO
from Bio import SearchIO
from Bio.SeqRecord import SeqRecord
from collections import defaultdict

parser = argparse.ArgumentParser()
parser.add_argument("--protein", help="The fasta file of proteins to be used for building the Database", type =str)
parser.add_argument("--prefix", help="Prefix of the output files", default="DB", type=str)
parser.add_argument("--min_cluster_size", help="The minimum number of proteins in a MMseqs cluster to be included in the DB", default=3, type=int)
parser.add_argument("--mmseqs_c", help="Value of -c for MMseqs cluster", default=0.7, type=float)
parser.add_argument("--mmseqs_m", help="Value of --min-seq-id for MMseqs cluster", default=0.35, type=float)
parser.add_argument("--threads", help="The number of threads to be used per data chunk", default=1, type=int)
parser.add_argument("--onlyparse", help="Skip analysis and only parse the output files", default=False, type=bool)

args = parser.parse_args()

def central():
    protein_clusters_file = cluster_proteins(args.protein)
    (cluster_info,protein_info) = split_clusters(protein_clusters_file)
    cluster_info = align_and_build(cluster_info)
    og_info_file = args.prefix + '_OG_Info.tsv'
    print_results(cluster_info,og_info_file,'NA','OG')
    protein_info_file = args.prefix + '_Protein_Info.tsv'
    print_results(protein_info,protein_info_file,'NA','Protein')
    

def split_clusters(protein_clusters_file):
    print(f"Reading {protein_clusters_file}")
    protein_og_df = pd.read_csv(protein_clusters_file, sep="\t",index_col=1,header=None)
    protein_og_df.index.name = 'Protein'
    protein_og_df.rename(columns={0: "OG"},inplace=True)
    cluster_info = defaultdict(dict)
    protein_info = defaultdict(dict)
    print("Calculating cluster protein counts")
    for protein,row in protein_og_df.iterrows():
        protein_cluster = row["OG"]
        protein_info['OG'][protein] = protein_cluster
        if protein_cluster not in cluster_info['Members']:
            cluster_info['Members'][protein_cluster] = 0
        cluster_info['Members'][protein_cluster] += 1
    
    print('Spliting proteins into cluster files')
    #Create a directory to store unaligned fasta file of clusters
    unaligned_dir = args.prefix +'_Unaligned_OG'
    if args.onlyparse == False:
        create_dir(unaligned_dir)
    #Iterate over each sequence of the proteins file
    for seqobj in SeqIO.parse(args.protein, 'fasta'):
        prot_cluster = protein_info['OG'][seqobj.id]
        #Print the sequence to its cluster file if the cluster meets the criteria for minimum number of members
        if (cluster_info['Members'][prot_cluster] >= args.min_cluster_size):
            cluster_info['Unaligned_Protein_File'][prot_cluster] = args.prefix +'_Unaligned_OG'+'/'+prot_cluster+'.faa'
            if (args.onlyparse == False):
                fasta_cluster_name =  cluster_info['Unaligned_Protein_File'][prot_cluster]
                fasta_handle = open(fasta_cluster_name, "a+")
                SeqIO.write(seqobj,fasta_handle, "fasta")
                fasta_handle.close()

    return(cluster_info,protein_info)

    
def align_and_build(cluster_info):   
    print('Aligning Clusters and building relevant files')
    aligned_og_dir = args.prefix +'_Aligned_OG'
    a3m_dir =  args.prefix +'_a3m_Aligned_OG'
    hhm_dir =  args.prefix +'_HHM_Aligned_OG'
    
    if args.onlyparse == False:
        #Create a directory to store aligned fasta file of clusters
        create_dir(aligned_og_dir)
        create_dir(a3m_dir)
        create_dir(hhm_dir)
        #Iterate over all clusters
    for prot_cluster in cluster_info['Members'].keys():
        #Skip clusters that do not meet the criteria for minimum size
        if (cluster_info['Members'][prot_cluster] >= args.min_cluster_size):
            #Align cluster with muscle
            unaligned_file = cluster_info['Unaligned_Protein_File'][prot_cluster]
            aligned_file =  args.prefix +'_Aligned_OG'+'/'+prot_cluster+'.faa'
            cluster_info['Aligned_Protein_File'][prot_cluster] = aligned_file
            command = f'muscle -in {unaligned_file} -out {aligned_file} -quiet'
            if args.onlyparse == False:
                subprocess.call(command, shell=True)
            a3m_file = f'{a3m_dir}/Aligned_{prot_cluster}.a3m'
            cluster_info['Aligned_a3m_File'][prot_cluster] = a3m_file
            command = f'reformat.pl fas a3m {aligned_file} {a3m_file}'
            if args.onlyparse == False:
                subprocess.call(command, shell=True)
            #Convert a3m file to HHM
            hhm_file = f'{hhm_dir}/Aligned_{prot_cluster}.hhm'
            #cluster_info['Aligned_HHM_File'][prot_cluster] = hhm_file
            command = f'hhmake -i {a3m_file} -o {hhm_file}'
            if args.onlyparse == False:
                #subprocess.call(command, shell=True)
    
    
    if args.onlyparse == False:
        print("Building Indexes")
        command = f'ffindex_build -s {args.prefix}_Aligned.ff{data,index} {args.prefix}_Aligned_OG/'
        #subprocess.call(command, shell=True)
        command = f'ffindex_build -s {args.prefix}_a3m.ff{data,index} {args.prefix}_a3m_Aligned_OG/'
        #subprocess.call(command, shell=True)
    return(cluster_info)


def cluster_proteins(protein_file):
    if args.onlyparse == False:
        #Builds and Runs the MMseqs clustering command
        print('Clustering proteins from',protein_file)
        command = f'mmseqs easy-linclust {args.protein} {args.prefix} mmseqs_tmp --threads {args.threads} -c {args.mmseqs_c} --min-seq-id {args.mmseqs_m} --cov-mode 0'
        subprocess.call(command, shell=True)

    protein_clusters_file = args.prefix + '_cluster.tsv'
    return(protein_clusters_file)


def get_prefix(file,format):
    prefix_file = re.sub(f'.{format}','',file)
    prefix_file = re.sub('(.)+/','',prefix_file)
    return prefix_file


def print_results(info_dict,out_file,missing_val,index_name):
    print(f'Printing results to {out_file}')
    info_df = pd.DataFrame.from_dict(info_dict)
    if (index_name):
        info_df.index.name = index_name
    info_df.to_csv(out_file,sep="\t",na_rep=missing_val)
    return True


def create_dir(directory):
    try:
        os.mkdir(directory)
    except OSError:
        raise Exception(f'Output directory {directory} already exists!')

central()