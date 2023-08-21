#! /usr/bin/python3
#Viral words version 0.3
import os
import multiprocessing as mp
from joblib import Parallel, delayed
import numpy as np
import pandas as pd
import subprocess
import csv
import re
import argparse
import gensim
from Bio.Seq import Seq
from Bio import SeqIO
from Bio import SearchIO
from Bio.SeqRecord import SeqRecord
from sklearn.cluster import DBSCAN
from sklearn.metrics.pairwise import pairwise_distances
from collections import defaultdict

parser = argparse.ArgumentParser()
parser.add_argument("--prefix", help="Prefix of the output files", default="VWV", type=str)
parser.add_argument("--output_directory", help="Directory to write output files", default="VWV_Output", type=str)
parser.add_argument("--genome", help="The fasta file of genomes to be used for training", type =str)
parser.add_argument("--protein", help="The OPTIONAL fasta file of proteins to be used for training. Protein IDs must follow the PRODIGAL naming scheme. Will skip gene calling if provided", type =str)
parser.add_argument("--map", help="The OPTIONAL map file of proteins to be used for training. Will skip gene calling and clustering if provided", type =str)
parser.add_argument("--use_score_map", help="Use a Genome Map based on HMM searches instead of the one derived from just protein clustering (slower)", default=False, type=bool)
parser.add_argument("--taxonomy", help="Optional .tsv file with the taxonomic classification of the genomes", type =str)
parser.add_argument("--min_cluster_size", help="The minimum number of proteins in a MMseqs cluster", default=3, type=int)
parser.add_argument("--mmseqs_c", help="Value of -c for MMseqs cluster", default=0.7, type=float)
parser.add_argument("--mmseqs_s", help="Value of -s for MMseqs cluster", default=7.5, type=float)
parser.add_argument("--mmseqs_m", help="Value of --min-seq-id for MMseqs cluster", default=0.35, type=float)
parser.add_argument("--window_size", help="The window size to use when training the Word2Vec algorithm", type=int, default = 5)
parser.add_argument("--epsilon", help="The epsilon value to be used by the DBSCAN algorithm", type=float, default = 0.025)
parser.add_argument("--min_samples", help="The min_samples value be used by the DBSCAN algorithm", default=2, type=int)
parser.add_argument("--metric", help="The distance metric to be used by the DBSCAN algorithm", default='cosine', type=str)
parser.add_argument("--correct", help="Wether to correct distances with the Coutinho approach", default=False, type=bool)
parser.add_argument("--alpha", help="Value of alpha to be use by the Coutinho correction approach", default=1.5, type=float)
parser.add_argument("--hhsuite_query_db", help="Directory containing the pre-computed hhsuite database to use for annotating the OGs", type=str)
parser.add_argument("--fasta_query_db", help="Fasta file containing the protein sequences to query against the OGs for annotation", type=str)
parser.add_argument("--annotation_query_db", help="Optional .tsv dictionary of query DB entry annotations", type=str)
parser.add_argument("--evalue", help="The maximum e-value of protein searches to consider a match as valid", default = 0.001, type=float)
parser.add_argument("--score", help="The minimum score of protein searches to consider a match as valid", default=50, type=int)
parser.add_argument("--optimize", help="Flag to run the optimization module", default=False, type=bool)
parser.add_argument("--onlyparse", help="Skip analysis and only parse the output files", default=False, type=bool)
parser.add_argument("--chunks", help="Number of chunks to split the query files when performing protein search", type =int, default = 1)
parser.add_argument("--threads", help="The number of threads to be used per data chunk", default=1, type=int)


args = parser.parse_args()

def central():
    if args.onlyparse == False:
        create_dir(args.output_directory)
    if ((args.genome) and (not args.protein)):
        args.protein = call_phanotate(args.genome)
    if (args.protein):
        protein_clusters_file = cluster_proteins(args.protein)
        (cluster_info,protein_info,genome_info,genome_og_info) = collect_cluster_info(protein_clusters_file)
        (genome_info,cluster_info,protein_info,rev_alias_dict,map_file) = print_map(genome_info,cluster_info,protein_info,False)
    if (args.hhsuite_query_db or args.fasta_query_db or args.use_score_map):
       (cluster_info,protein_info) = split_clusters(cluster_info,protein_info)
       cluster_info = align_and_build(cluster_info)
    if (args.use_score_map):    
        (cluster_info,protein_info) = align_fasta_to_hmm(cluster_info,protein_info)
    if (args.hhsuite_query_db):
        (cluster_info,protein_info) = align_hhm_to_hhmdb(cluster_info,protein_info,args.hhsuite_query_db)
    if (args.use_score_map == True):
        (genome_info,cluster_info,protein_info,rev_alias_dict,map_file) = print_map(genome_info,cluster_info,protein_info,True)
    cluster_info = make_modules(cluster_info,map_file,rev_alias_dict)
    
    #Print OG info and Protein info to output files
    og_info_file = args.output_directory+'/'+ args.prefix + '_OG_Info.tsv'
    print_results(cluster_info,og_info_file,'NA','OG')
    protein_info_file = args.output_directory+'/'+ args.prefix + '_Protein_Info.tsv'
    print_results(protein_info,protein_info_file,'NA','Protein')

def align_hhm_to_hhmdb(cluster_info,protein_info,hhsuite_query_db):
    print(f'Querying OGs against {hhsuite_query_db}')
    db_prefix = get_prefix(hhsuite_query_db,"No_Extension")
    output_directory = args.output_directory+'/'+ args.prefix +f'_HHSearch_Alignments_{db_prefix}'
    if args.onlyparse == False:
        print('Creating',output_directory)
        #create_dir(output_directory)
    valid_prot_clusters = []
    for prot_cluster in cluster_info['Members'].keys():
        #Skip clusters that do not meet the criteria for minimum size
        if (cluster_info['Members'][prot_cluster] >= args.min_cluster_size):
            valid_prot_clusters.append(cluster_info['Aligned_HHM_File'][prot_cluster])

    result_files = Parallel(n_jobs=args.chunks)(delayed(call_hhsearch)(query_file=prot_cluster) for prot_cluster in valid_prot_clusters)
    
    query_dict_df = []
    if (args.annotation_query_db):
        print('Reading annotation info from',args.annotation_query_db)
        query_dict_df =  pd.read_csv(args.annotation_query_db,sep='\t',na_values='',header=0,index_col=0)
    file_counter = 0    
    print(f'Parsing HHSearch files')
    for outfile in result_files:
        with open(outfile, 'r') as input:
            for line in input.readlines():
                vals = line.split('\t')
                vals[11] = float(vals[11].rstrip())
                protein_id = vals[0]
                protein_cluster = protein_info['OG'][protein_id]
                if protein_cluster not in cluster_info['Best_Hit_DB_Score'].keys():
                    cluster_info['Best_Hit_DB_Domain'][protein_cluster] = vals[1]
                    cluster_info['Best_Hit_DB_Percent_Matched'][protein_cluster] = float(vals[2]) * 100
                    cluster_info['Best_Hit_DB_Score'][protein_cluster] = vals[11]
                    cluster_info['Best_Hit_DB_evalue'][protein_cluster] = vals[10]
                    if (args.annotation_query_db):
                        for category in query_dict_df.columns:
                            cluster_info['Best_Hit_DB_'+category][protein_cluster] = query_dict_df.loc[vals[1],category]
                if vals[11] > cluster_info['Best_Hit_DB_Score'][protein_cluster]:
                    cluster_info['Best_Hit_DB_Domain'][protein_cluster] = vals[1]
                    cluster_info['Best_Hit_DB_Percent_Matched'][protein_cluster] = float(vals[2]) * 100
                    cluster_info['Best_Hit_DB_Score'][protein_cluster] = vals[11]
                    cluster_info['Best_Hit_DB_evalue'][protein_cluster] = vals[10]
                    if (args.annotation_query_db):
                        for category in query_dict_df.columns:
                            cluster_info['Best_Hit_'+category][protein_cluster] = query_dict_df.loc[vals[1],category]
        file_counter += 1
        if (file_counter % 100 == 0):
            print('Processed HHPred files:',file_counter)

    return(cluster_info,protein_info)
                    
def call_phanotate(genomes_file):
    file_prefix = get_prefix(genomes_file,"fasta")
    if args.onlyparse == False:
        command = f"python3 phanotate.py -o {file_prefix}.faa -f fasta {genomes_file}"
        print('Running',command)
        subprocess.call(command, shell=True)
        
    return(f"{file_prefix}.faa")

def call_hhsearch(query_file):
    db_prefix = get_prefix(args.hhsuite_query_db,"No_Extension")
    in_prefix = get_prefix(query_file,"hhm")
    hhr_outfile = args.output_directory+'/'+ args.prefix +f'_HHSearch_Alignments_{db_prefix}/{in_prefix}.xDB.hhr'
    m8_outfile = args.output_directory+'/'+ args.prefix +f'_HHSearch_Alignments_{db_prefix}/{in_prefix}.xDB.m8'
    if args.onlyparse == False:
        command = f'hhsearch -cpu {args.threads} -i {query_file} -d {args.hhsuite_query_db} -o {hhr_outfile} -blasttab {m8_outfile}'
        print('Running',command)
        subprocess.call(command, shell=True)
    return(m8_outfile)

def make_modules(cluster_info,map_file,rev_alias_dict):
    genome_tokens = list (read_genome_structure(map_file))

    print('Setting up Word2Vec model')
    model = gensim.models.Word2Vec (genome_tokens, vector_size=100, window=args.window_size, min_count=args.min_cluster_size, workers=args.threads, sg=1)

    print('Training Word2Vec model')
    model.train(genome_tokens,total_examples=len(genome_tokens),epochs=10)

    #words = sorted(model.wv.index_to_key)
    #scaled_data = [model.wv[w] for w in words]
    #cluster_distrib_distance = pairwise_distances(scaled_data, metric=args.metric)
    print('Calculating distribution distance among clusters')
    cluster_distrib_distance = pairwise_distances(model.wv.get_normed_vectors(), metric=args.metric)
    print('Dimensions of cluster distribution distance data frame:',cluster_distrib_distance.shape)
    
    transformed_cluster_distance = cluster_distrib_distance
    
    if args.correct == True:
        print('Reading',genus_range_output)
        genus_tax_range = pd.read_csv(genus_range_output,sep='\t')
        nrow = len(genus_tax_range.iloc[:,0])
        for row in range(0,nrow):
            protein_cluster = str(genus_tax_range['Variable'][row])
            alias = cluster_info['Alias'][protein_cluster]
            word = words[row]
            if (word != alias):
                print('Fatal error for row:',row,alias,word,genus_range_output,'and',map_file,'are not in the same order')
        del genus_tax_range['Variable']
        print('Calculating taxonomic distance among clusters')
        cluster_tax_distance = pairwise_distances(genus_tax_range, metric=args.metric)
        print('Dimensions of cluster taxonomic distance data frame:',cluster_tax_distance.shape)
        transformed_cluster_distance = ((cluster_distrib_distance) / (args.alpha * (cluster_tax_distance + 1)))
        print('Dimensions of transformed cluster distance data frame:',transformed_cluster_distance.shape)
    
    print('Printing OG word vector distance matrix')
    output_wv_dist_file = args.output_directory+'/'+ args.prefix +'_OG_WV_Distances.tsv'
    #pd.DataFrame(transformed_cluster_distance).to_csv(output_wv_dist_file,sep="\t",na_rep='NA',header=[rev_alias_dict[w] for w in model.wv.index_to_key])
    dbscan = DBSCAN(eps=args.epsilon,n_jobs=args.threads,metric='precomputed',min_samples=args.min_samples)
    #Perform DBSCAN clustering on the scaled data
    dbscan.fit(transformed_cluster_distance)
    
    dbscan_outfile = args.output_directory+'/'+ args.prefix+'_DBSCAN_Clustering.tsv'
    out_array = dbscan.labels_
    np.savetxt(dbscan_outfile,out_array,delimiter='\t')
    print('Median Module value is',np.median(dbscan.labels_))

    print('Collecting Module info from',dbscan_outfile)
    words = model.wv.index_to_key
    with open(dbscan_outfile, 'r') as input:
        counter = 0
        for line in input.readlines():
            prot_cluster = rev_alias_dict[str(words[counter])]
            module_num = line.rstrip()
            #print(prot_cluster,type(prot_cluster),module_num,type(module_num),type(cluster_info))
            cluster_info['Module'][prot_cluster] = module_num
            counter += 1

    return(cluster_info)
   

def split_clusters(cluster_info,protein_info):
    print('Spliting proteins into cluster files')
    #Create a directory to store unaligned fasta file of clusters
    unaligned_dir = args.output_directory+'/'+ args.prefix +'_Unaligned_OG'
    if args.onlyparse == False:
        create_dir(unaligned_dir)
    #Iterate over each sequence of the proteins file
    for seqobj in SeqIO.parse(args.protein, 'fasta'):
        #Collect description and legth information and store it in protein_info
        protein_info['Description'][seqobj.id] = seqobj.description
        protein_info['Length'][seqobj.id] = len(seqobj.seq)
        prot_cluster = protein_info['OG'][seqobj.id]
        #Print the sequence to its cluster file if the cluster meets the criteria for minimum number of members
        if (cluster_info['Members'][prot_cluster] >= args.min_cluster_size):
            cluster_info['Unaligned_Protein_File'][prot_cluster] = args.output_directory+'/'+ args.prefix +'_Unaligned_OG'+'/'+prot_cluster+'.faa'
            if (args.onlyparse == False):
                fasta_cluster_name =  cluster_info['Unaligned_Protein_File'][prot_cluster]
                fasta_handle = open(fasta_cluster_name, "a+")
                SeqIO.write(seqobj,fasta_handle, "fasta")
                fasta_handle.close()

    return(cluster_info,protein_info)

def align_fasta_to_hmm(cluster_info,protein_info):
    genome_hmm_scores = defaultdict(dict)
    #split the files according to the specified chunk number
    split_files = split_fasta(args.protein)
    #Align proteins against the generated hmm and parse the result
    hmmer_dbs = [f'{args.output_directory}/All_HMM_{args.prefix}.hmm']
    if __name__ == "__main__":
         result_files = Parallel(n_jobs=args.chunks)(delayed(call_hmmsearch)(query_file=file,db_file=hmmer_db) for file in split_files for hmmer_db in hmmer_dbs)
    #Parse the output files
    for file in result_files:
        print(f'Parsing {file}')
        #Iterate over each query
        for qresult in SearchIO.parse(file, 'hmmer3-text'):
            #Iterate over each Hit
            for hit in qresult.hits:
                #Iterate over each HSP
                for hsp in hit.hsps:
                    protein = hit.id
                    genome = hit.id
                    genome = re.sub('_(\d)+$','',genome)
                    is_valid = check_search_cutoffs(hsp)
                    if is_valid:
                        if (protein not in protein_info['Best_Hit_ID']):
                            protein_info['Best_Hit_ID'][protein] = qresult.id
                            protein_info['Best_Hit_Score'][protein] = hsp.bitscore
                        elif (hsp.bitscore > protein_info['Best_Hit_Score'][protein]):
                            protein_info['Best_Hit_ID'][protein] = qresult.id
                            protein_info['Best_Hit_Score'][protein] = hsp.bitscore
                        ori_score = genome_hmm_scores[qresult.id].get(genome,0)
                        if (ori_score < hsp.bitscore):
                            genome_hmm_scores[qresult.id][genome] = hsp.bitscore
                    elif (genome not in genome_hmm_scores[qresult.id]):
                        genome_hmm_scores[qresult.id][genome] = 0
    
    genome_HMM_output = args.output_directory+'/'+ args.prefix + '_GenomesxOG_Scores.tsv'
    print_results(genome_hmm_scores,genome_HMM_output,'0','Genome')
    return(cluster_info,protein_info)

def call_hmmsearch(query_file='NA',db_file='NA'):
    query_prefix = get_prefix(query_file,'faa')
    out_file =  args.output_directory+'/'+ args.prefix +'_'+query_prefix+'xAll_HMM'
    if args.onlyparse == False:
        print(f'Aligning proteins from {query_file} to {db_file}')
        command = f'hmmsearch -o {out_file} --noali --cpu {args.threads} {db_file} {query_file}'
        subprocess.call(command, shell=True)
    return(out_file)


def split_fasta(protein_file):
    #Count the number of seqs in the query file
    split_names_list = []
    command = f'grep -c ">" {protein_file}'
    query_seq_count = int(subprocess.getoutput(command))
    #Calculate the number of sequences in each chunk (split files)
    chunk_size = round(query_seq_count / args.chunks)
    print('Splitting',protein_file,'into',args.chunks,'chunks of',chunk_size,'sequences.')
    seq_counter = 0
    split_counter = 1
    split_name = args.output_directory+'/'+ args.prefix +'_Split_'+str(split_counter)+'.faa'
    split_names_list.append(split_name)
    split_handle = open(split_name,"w")
    #Iterate over all sequences in the query file. Keep track of the number of sequences added to each split file in seq_counter
    for seqobj in SeqIO.parse(protein_file,format="fasta"):
        SeqIO.write(seqobj,split_handle, "fasta")
        seq_counter += 1
        #If we reached the number of sequences in the split file, close the handle and open a new split file handle and reset the seq_counter. Keep track of split file names in split_names_list
        if (seq_counter == chunk_size):
            split_counter += 1
            split_name = args.output_directory+'/'+ args.prefix +'_Split_'+str(split_counter)+'.faa'
            split_names_list.append(split_name)
            split_handle.close()
            split_handle = open(split_name,"w")
            seq_counter = 0
    split_handle.close()
    
    return(split_names_list)
    
def align_and_build(cluster_info):   
    print('Aligning Clusters and building relevant files')
    aligned_og_dir = args.output_directory+'/'+ args.prefix +'_Aligned_OG'
    hmm_og_dir = args.output_directory+'/'+ args.prefix +'_HMM_OG'
    a3m_dir =  args.output_directory+'/'+ args.prefix +'_a3m_Aligned_OG'
    hhm_dir =  args.output_directory+'/'+ args.prefix +'_HHM_Aligned_OG'
    
    if args.onlyparse == False:
        #Create a directory to store aligned fasta file of clusters
        create_dir(aligned_og_dir)
        if (args.fasta_query_db or args.use_score_map):
            #Create a directory to store the HMMs
            create_dir(hmm_og_dir)
        if (args.hhsuite_query_db):
            #Creat direcotires to store files to run hhpred and the output of the runs
            db_prefix = get_prefix(args.hhsuite_query_db,"No_Extension")
            hhpred_results_dir =  args.output_directory+'/'+ args.prefix +f'_HHSearch_Alignments_{db_prefix}'
            create_dir(a3m_dir)
            create_dir(hhm_dir)
            create_dir(hhpred_results_dir)
        #Iterate over all clusters
    for prot_cluster in cluster_info['Members'].keys():
        #Skip clusters that do not meet the criteria for minimum size
        if (cluster_info['Members'][prot_cluster] >= args.min_cluster_size):
            #Align cluster with muscle
            unaligned_file = cluster_info['Unaligned_Protein_File'][prot_cluster]
            aligned_file =  args.output_directory+'/'+ args.prefix +'_Aligned_OG'+'/'+prot_cluster+'.faa'
            cluster_info['Aligned_Protein_File'][prot_cluster] = aligned_file
            command = f'muscle -in {unaligned_file} -out {aligned_file} -quiet'
            if args.onlyparse == False:
                subprocess.call(command, shell=True)
            #Build HMM from alignment
            if (args.fasta_query_db or args.use_score_map):
                hmmer_file =  args.output_directory+'/'+ args.prefix +'_HMM_OG'+'/'+prot_cluster+'.hmm'
                cluster_info['HMM_File'][prot_cluster] = hmmer_file
                command = f'hmmbuild -n {prot_cluster} {hmmer_file} {aligned_file}'
                if args.onlyparse == False:
                    subprocess.call(command, shell=True)
            if (args.hhsuite_query_db):
                a3m_file = f'{a3m_dir}/Aligned_{prot_cluster}.a3m'
                cluster_info['Aligned_a3m_File'][prot_cluster] = a3m_file
                command = f'reformat.pl fas a3m {aligned_file} {a3m_file}'
                if args.onlyparse == False:
                    subprocess.call(command, shell=True)
                 #Convert a3m file to HHM
                hhm_file = f'{hhm_dir}/Aligned_{prot_cluster}.hhm'
                cluster_info['Aligned_HHM_File'][prot_cluster] = hhm_file
                command = f'hhmake -i {a3m_file} -o {hhm_file}'
                if args.onlyparse == False:
                    subprocess.call(command, shell=True)


    if args.onlyparse == False:
        #Merge all HMMs into a single file
        command = f'cat {hmm_og_dir}/*.hmm > {args.output_directory}/All_HMM_{args.prefix}.hmm'
        subprocess.call(command, shell=True)
        #Run HMMpress on the merged hmm file
        print('Building HMM database')
        command = f'hmmpress {args.output_directory}/All_HMM_{args.prefix}.hmm'
        subprocess.call(command, shell=True)

    return(cluster_info)


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

def cluster_proteins(protein_file):
    if args.onlyparse == False:
        #Builds and Runs the MMseqs clustering command
        print('Clustering proteins from',protein_file)
        tmp_dir = args.output_directory+'/tmp'+ args.prefix
        command = f'mmseqs easy-cluster {args.protein} {args.prefix} {tmp_dir} --threads {args.threads} -c {args.mmseqs_c} -s {args.mmseqs_s} --min-seq-id {args.mmseqs_m} --cov-mode 0'
        subprocess.call(command, shell=True)
        #MMSeqs will write to CWD by default so the files need to be transferred to the specified output Dir 
        if (args.output_directory != '/.'):
            command = f'mv {args.prefix}_cluster.tsv {args.prefix}_all_seqs.fasta {args.prefix}_rep_seq.fasta {args.output_directory}/'
            subprocess.call(command, shell=True)
    protein_clusters_file = args.output_directory+'/'+ args.prefix + '_cluster.tsv'
    return(protein_clusters_file)
    
def collect_cluster_info(protein_clusters_file):
    #Parses the _clusters.tsv output file from MMseqs and assigns the proteins to their respecive clusters in the  protein_info dictionary    print('Collecting cluster info from',protein_clusters_file)
    print('Collecting Cluster Info from',protein_clusters_file)
    cluster_info = defaultdict(dict)
    protein_info = defaultdict(dict)
    genome_info = defaultdict(dict)
    genome_og_info = defaultdict(lambda: defaultdict(int))
    
    compiled_obj = re.compile('_(\d)+$')
    with open(protein_clusters_file, 'r') as input:
        for line in input.readlines():
            vals = line.split('\t')
            protein_cluster = 'Cluster_'+vals[0]
            protein = vals[1].rstrip()
            protein_info['OG'][protein] = protein_cluster
            genome = protein
            genome = re.sub(compiled_obj,'',genome)
            protein_info['Genome'][protein] = genome
            genome_og_info[protein_cluster][genome] += 1
            if genome not in genome_info['Total_Proteins']:
                genome_info['Total_Proteins'][genome] = 0
            genome_info['Total_Proteins'][genome] += 1
            if protein_cluster not in cluster_info['Members']:
                cluster_info['Members'][protein_cluster] = 0
            cluster_info['Members'][protein_cluster] += 1
    
    #Build the DF of Genome x OG counts and print to file
    genome_og_info_df = pd.DataFrame.from_dict(genome_og_info)
    genome_og_info_df.index.name = 'Genome'
    genome_og_info_file = args.output_directory+'/'+args.prefix + "_GenomexOG_Counts.tsv"
    print(f"Printing info to {genome_og_info_file}")
    #genome_og_info_df.to_csv(genome_og_info_file,sep="\t",na_rep=0)
    
    return(cluster_info,protein_info,genome_info,genome_og_info)

def get_alias(oid):
    #Convert numbers in a string to letters according to the dictionary below
    num_dict = {'0':'A','1':'B','2':'C','3':'D','4':'E','5':'F','6':'G','7':'H','8':'I','9':'J','.':'_'}
    protein_cluster_alias = list(oid)
    for char_posit in range(0,len(protein_cluster_alias)):
        char = protein_cluster_alias[char_posit]
        if char in num_dict.keys():
            protein_cluster_alias[char_posit] = num_dict[char]
    protein_cluster_alias = ''.join(protein_cluster_alias)
    
    return(protein_cluster_alias)
    
def print_map(genome_info,cluster_info,protein_info,score_map):
    map_file = args.output_directory+'/'+args.prefix + "_OG_Map.txt"
    if (score_map == True):
        map_file = args.output_directory+'/'+args.prefix + "_OG_Score_Map.txt"
    #Word vectos does not allow numbers so the OG ids are assigned an alias and we keep track of it in re_valias_dict and cluster_info
    rev_alias_dict = {}
    print('Printing genome map to',map_file)
    with open(map_file, 'w', newline='') as csvfile:
        tablewriter = csv.writer(csvfile, delimiter='\t')
        for genome in genome_info['Total_Proteins'].keys():
            ordered_hmms = []
            for prot_posit in range(1,(int(genome_info['Total_Proteins'][genome]) + 1)):
                protein = genome+'_'+str(prot_posit)
                protein_cluster = protein_info['OG'][protein]
                if ((score_map == True) and (protein in protein_info['Best_Hit_ID'])):
                    protein_cluster = protein_info['Best_Hit_ID'][protein]
                protein_cluster_alias = get_alias(protein_cluster)
                cluster_info['Alias'][protein_cluster] = protein_cluster_alias
                rev_alias_dict[protein_cluster_alias] = protein_cluster 
                ordered_hmms.append(protein_cluster_alias)
            tablewriter.writerow(ordered_hmms)
    
    return(genome_info,cluster_info,protein_info,rev_alias_dict,map_file)


def check_search_cutoffs(hsp):
    if ((hsp.bitscore < args.score) or (hsp.evalue > args.evalue)):
        return False
    else:
        return True


def read_genome_structure(map_file):
    """This method reads the Genome Map file which is plain text format"""
    print("Reading {0}".format(map_file))
    with open (map_file, 'r') as f:
        for i, line in enumerate (f): 
            if (i%100==0):
                print("Read {0} entries".format (i))
            #yield gensim.utils.simple_preprocess(line)
            yield list(gensim.utils.tokenize(line,deacc=False,lowercase=False))


def get_prefix(file,format):
    prefix_file = re.sub(f'.{format}','',file)
    prefix_file = re.sub('(.)+/','',prefix_file)
    return prefix_file


central()