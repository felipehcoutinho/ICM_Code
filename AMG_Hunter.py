#! /usr/bin/env python3

from collections import defaultdict
from Bio import SeqIO
from Bio import SearchIO
#from ete2 import NCBITaxa
import pandas as pd
import argparse
import re
import subprocess
import json

parser = argparse.ArgumentParser()
parser.add_argument("--cds", help="Input fasta file of CDS protein sequences", type =str)
parser.add_argument("--ref", help="Optional Fasta file of reference protein sequences", type =str)
parser.add_argument("--min_seq_len", help="The minimum sequence length to be included in the phylogeny", type =int, default=0)
parser.add_argument("--max_seq_len", help="The maximum sequence length to be included in the phylogeny", type =int, default=999999)
parser.add_argument("--pfam_db", help="The pfam Database file formated for Hmmer", type =str, default="/mnt/netapp2/Store_uni/COMPARTIDO/pfamdb/DB/Pfam-A.hmm")
parser.add_argument("--pfam_hits_file", help="The output of the hmmer search against Pfam (Will use this instead of searching if provided)", type =str, default="NA")
parser.add_argument("--pfam_min_score", help="The minimum score to consider a pfam hit", type =int, default=50)
parser.add_argument("--pfam_max_evalue", help="The maximum e-value to consider a pfam hit", type =float, default=0.00001)
parser.add_argument("--info_pfam_output", help="The Pfam hits info table file to be generated", type =str, default="Pfam_Info.tsv")
parser.add_argument("--skip_pfam", help="Flag to skip Pfam search and parsing", type = bool, default=False)
parser.add_argument("--kegg_db", help="The KEGG Database (KOfam) file formated for Hmmer", type =str, default="/mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/Databases/KOfam/All_KOfam.hmm")
parser.add_argument("--kegg_hits_file", help="The output of the hmmer search against KEGG (Will use this instead of searching if provided)", type =str, default="NA")
parser.add_argument("--kegg_min_score", help="The minimum score to consider a KEGG hit", type =int, default=50)
parser.add_argument("--kegg_max_evalue", help="The maximum e-value to consider a KEGG hit", type =float, default=0.00001)
parser.add_argument("--info_kegg_output", help="The KEGG hits info table file to be generated", type =str, default="KEGG_Info.tsv")
parser.add_argument("--kegg_json", help=".json file from KEGG", type =str, default="/mnt/smart/users/fcoutinho/Databases/KEGG/ko00001.json")
parser.add_argument("--skip_kegg", help="Flag to skip KEGG search and parsing", type = bool, default=False)
parser.add_argument("--uniref_db", help="The UniRef100 Database formated for DIAMOND", type =str, default="NA")
parser.add_argument("--uniref_hits_file", help="The Uniref m8 generated by BLASTp (or similar)", type =str, default="NA")
parser.add_argument("--enriched_uniref_hits_file", help="he Uniref m8 generated by BLASTp (or similar) + Enrich_m8", type =str)
parser.add_argument("--uniref_min_score", help="The minimum score to consider a Uniref hit", type =int, default=50)
parser.add_argument("--uniref_max_evalue", help="The maximum e-value to consider a UniRef hit", type =float, default=0.00001)
parser.add_argument("--enrich_uniref", help="Flag to enrich uniref annotations with taxonomic data before parsing", type = bool, default=False)
parser.add_argument("--skip_uniref", help="Flag to skip uniref search and parsing", type = bool, default=False)
parser.add_argument("--info_cds_output", help="The CDS info output table file to be generated", default="CDS_Info.tsv", type =str)
parser.add_argument("--info_genome_output", help="The Genome info output table file to be generated ", default="Genome_Info.tsv", type =str)
parser.add_argument("--genome_abundance", help="Optional Genomic Sequence Abundance matrix used to calculate KO, Pathway and Module abundances", type=str)
parser.add_argument("--gene_abundance", help="Optional Gene Sequence Abundance matrix used to calculate KO, Pathway and Module abundances", type=str)
parser.add_argument("--parse_only", help="Flag to skip running any programs and only parse their output", default=False, type=bool)
parser.add_argument("--threads", help="Number of threads to use during analysis", default=1, type=int)
parser.add_argument("--annotate", help="Flag to run the annotation module", default=False, type=bool)
parser.add_argument("--kegg_phylogeny", help="Flag to run the phylogeny module starting from KEGG hits", default=False, type=bool)
parser.add_argument("--kegg_ko", help="KO to be used as reference for the KEGG phylogeny", type=str)
parser.add_argument("--ko_hmm_dir", help="Directory where KEGG KO hmm models are located", default="/mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/Databases/KOfam/profiles/", type=str)
parser.add_argument("--derep_id", help="Identity threshold to dereplicate sequences when building KEGG phylogeny", default=1, type=float)

args = parser.parse_args()

def central():
    if (args.annotate == True):
        index_seqs(cds_file=args.cds)
        if (args.skip_uniref == False):
            if (args.uniref_hits_file != "NA"):
                uniref_hits_file = args.uniref_hits_file
            else:
                uniref_hits_file = call_diamond(cds_file=args.cds,db_file=args.uniref_db,program="blastp")
            print(f"Using Uniprot hits file: {uniref_hits_file}")
            if (args.enrich_uniref == True):
                enriched_uniref_hits_file = enrich_m8(m8_file=uniref_hits_file)
                parse_uniref_enriched(uniref_hits_file=enriched_uniref_hits_file,min_score=args.uniref_min_score,max_evalue=args.uniref_max_evalue)
            else:
                parse_uniref_raw(uniref_hits_file=uniref_hits_file,min_score=args.uniref_min_score,max_evalue=args.uniref_max_evalue)
        if (args.skip_pfam == False):
            if (args.pfam_hits_file != "NA"):
                pfam_hits_file = args.pfam_hits_file
            else:
                pfam_hits_file = call_hmmer(cds_file=args.cds,db_file=args.pfam_db,program="hmmsearch")
            print(f"Using pfam hits file: {pfam_hits_file}")
            (pfam_genome_hmm_scores,pfam_pairwise_scores) = parse_hmmer_output(hmmer_out_file=pfam_hits_file,min_score=args.pfam_min_score,max_evalue=args.pfam_max_evalue,multi_hits=True,prefix="Pfam",output_df_file="Pfam_Info.tsv")
        if (args.skip_kegg == False):
            if (args.kegg_hits_file != "NA"):
                kegg_hits_file = args.kegg_hits_file
            else:
                kegg_hits_file = call_hmmer(cds_file=args.cds,db_file=args.kegg_db,program="hmmsearch")
            print(f"Using KEGG hits file: {kegg_hits_file}")
            (kegg_genome_hmm_scores,kegg_pairwise_scores) = parse_hmmer_output(hmmer_out_file=kegg_hits_file,min_score=args.kegg_min_score,max_evalue=args.kegg_max_evalue,multi_hits=False,prefix="KEGG",output_df_file="KEGG_Info.tsv")
        print_results(cds_output_df_file=args.info_cds_output,genome_output_df_file=args.info_genome_output)
    if (args.genome_abundance):
        calc_kegg_abundance_from_genome(abundance_file=args.genome_abundance,genome_mode=True)
    if (args.gene_abundance):
        calc_kegg_abundance_from_gene(abundance_file=args.gene_abundance,genome_mode=False)    
    if (args.kegg_phylogeny):
        build_phylogeny_from_kegg_hits(cds_file=args.cds,ref_file=args.ref,ko=args.kegg_ko,ko_hmm_dir=args.ko_hmm_dir)

def parse_hmmer_output(hmmer_out_file="NA",max_evalue=0.00001,min_score=50,multi_hits=False,prefix="NA",output_df_file="Info.tsv",print_df=True):
    genome_hmm_scores = defaultdict(dict)
    pairwise_scores = defaultdict(dict)
    hsp_count = 0
    #Parse the output
    print(f'Parsing {hmmer_out_file}')
    #Iterate over each query
    for qresult in SearchIO.parse(hmmer_out_file, 'hmmer3-text'):
        #Iterate over each Hit
        for hit in qresult.hits:
            #Iterate over each HSP
            for hsp in hit.hsps:
                genome = hit.id
                genome = re.sub('_(\d)+$','',genome)
                #print(qresult.id,genome,hit.id)
                is_valid = check_hmmer_match_cutoff(hsp,max_evalue,min_score)
                if is_valid:
                    hsp_count += 1
                    pairwise_scores['Genome'][hsp_count] = genome
                    pairwise_scores['Query'][hsp_count] = qresult.id
                    pairwise_scores['Subject'][hsp_count] = hit.id
                    pairwise_scores['Score'][hsp_count] = hsp.bitscore
                    pairwise_scores['e-value'][hsp_count] = hsp.evalue
                    pairwise_scores['Subject_Description'][hsp_count] = hit.description
                    #pairwise_scores['Query_Description'][hsp_count] = qresult.description #Is always empty foe KOfam DB
                    #Add kegg hierarchy info the the pairwise scores dictionary everytime a hit is found, regardless of being the best one
                    ko = qresult.id
                    pairwise_scores['KEGG_Function'][hsp_count] = ""
                    pairwise_scores['KEGG_Modules'][hsp_count] = set()
                    pairwise_scores['KEGG_Pathways'][hsp_count] = set()
                    for categ in kegg_info[ko].keys():
                        for h1_name in kegg_info[ko][categ].keys():
                            for h2_name in kegg_info[ko][categ][h1_name].keys():
                                for h3_name in kegg_info[ko][categ][h1_name][h2_name].keys():
                                    pairwise_scores['KEGG_Function'][hsp_count] = kegg_info[ko][categ][h1_name][h2_name][h3_name]["Description"]
                                    pairwise_scores['KEGG_Modules'][hsp_count].add(h1_name)
                                    pairwise_scores['KEGG_Pathways'][hsp_count].add(h2_name)
                    #Keep track of best hit data for each query
                    genome_ori_score = genome_hmm_scores[qresult.id].get(genome,0)
                    if (genome_ori_score < hsp.bitscore):
                        genome_hmm_scores[qresult.id][genome] = hsp.bitscore
                    if ((multi_hits == True) and (prefix == "Pfam")):
                        if (hit.id not in cds_seq_info[prefix+'_Subjects']):
                            cds_seq_info[prefix+'_Subjects'][hit.id] = set()
                        cds_seq_info[prefix+'_Subjects'][hit.id].add(qresult.id)
                    elif ((multi_hits == False) and (prefix == "KEGG")):
                        #Add best hit data to cds info dictionary
                        if (hsp.bitscore > cds_seq_info[prefix+'_Best_Subject_Score'][hit.id]):
                            cds_seq_info[prefix+'_Best_Subject'][hit.id] = qresult.id
                            cds_seq_info[prefix+'_Best_Subject_Score'][hit.id] = hsp.bitscore                                            
    #Add the information of best hits only to genome_seq_info 
    print(f"Indexing {prefix} hits information")
    #Iterate over each query cds sequence for which a Best_Subject was found
    for cds in cds_seq_info['Length']:
        #Get the genome id
        scaffold_id = re.sub('_(\d)+$','',cds)
        #Get the best subject ID, and associted KEGG info for the cds Best Hit KO
        if (prefix == "KEGG"):
            ko = cds_seq_info['KEGG_Best_Subject'][cds]
            #When dealing with KEGG hits, get the KEGG info for the best hit of the cds if it was assigned to a KO. Note that KEGG_Best_Subject_Function holds a single function associated with the best hit of the CDS, while KEGG_Matched_Functions holds the set of the functons associated with all the KOS matched to the scaffold.
            if (ko != "NA"):
                for categ in kegg_info[ko].keys():
                    for h1_name in kegg_info[ko][categ].keys():
                        for h2_name in kegg_info[ko][categ][h1_name].keys():
                            for h3_name in kegg_info[ko][categ][h1_name][h2_name].keys():
                                #Add the KEGG hierarchy info to the cds info dictionary
                                desc = kegg_info[ko][categ][h1_name][h2_name][h3_name]["Description"]
                                cds_seq_info['KEGG_Best_Subject_Function'][cds] = desc
                                cds_seq_info['KEGG_Best_Subject_Pathways'][cds].add(h2_name)
                                cds_seq_info['KEGG_Best_Subject_Modules'][cds].add(h1_name)
                                #Add the KEGG hierarchy info to the genome info dictionary
                                genome_seq_info['KEGG_Matched_Subject_IDs'][scaffold_id].add(cds_seq_info['KEGG_Best_Subject'][cds])
                                genome_seq_info['KEGG_Matched_Functions'][scaffold_id].add(cds_seq_info['KEGG_Best_Subject_Function'][cds])
                                genome_seq_info['KEGG_Matched_Pathways'][scaffold_id].update(cds_seq_info['KEGG_Best_Subject_Pathways'][cds])
                                genome_seq_info['KEGG_Matched_Modules'][scaffold_id].update(cds_seq_info['KEGG_Best_Subject_Modules'][cds])                         
        #Get the best subject ID for the cds Pfam domain hits
        elif ((prefix == "Pfam") and (cds in cds_seq_info['Pfam_Subjects'])):
            #print(f"Indexing Pfam hits for scaffold: {scaffold_id}.")
            #print(f"CDS: {cds}.")
            #print(f"Pfam_Subjects: {cds_seq_info['Pfam_Subjects'][cds]}")
            genome_seq_info['Pfam_Matched_Subjects'][scaffold_id].update(cds_seq_info['Pfam_Subjects'][cds])                        
    if (print_df == True):
        print(f"Printing info to {output_df_file}")
        info_df = pd.DataFrame.from_dict(pairwise_scores)
        info_df.index.name = 'Hit_#'
        info_df.to_csv(output_df_file,sep="\t",na_rep='NA') 
    return(genome_hmm_scores,pairwise_scores)
             

def index_seqs(cds_file="NA"):
    seen_cds_ids = dict()
    seen_scaff_ids = dict()
    print ("Indexing CDS sequences from",cds_file)
    compiled_obj = re.compile("_(\d)+$")
    seq_counter = 0
    for seqobj in SeqIO.parse(cds_file, "fasta"):
        cds_id = seqobj.id
        #Do not allow duplicated IDs
        if (cds_id in seen_cds_ids):
            raise Exception(f'Duplicated ID: {seqobj.id} in {cds_file}')
        seen_cds_ids[cds_id] = True
        cds_seq_info['Length'][cds_id] = len(seqobj.seq)
        cds_seq_info['KEGG_Best_Subject'][cds_id] = "NA"
        cds_seq_info['KEGG_Best_Subject_Score'][cds_id] = 0
        cds_seq_info['KEGG_Best_Subject_Function'][cds_id] = "NA"
        cds_seq_info['KEGG_Best_Subject_Pathways'][cds_id] = set()
        cds_seq_info['KEGG_Best_Subject_Modules'][cds_id] = set()
        cds_seq_info['Pfam_Subjects'][cds_id] = set()
        
        scaffold_id = re.sub(compiled_obj,"",cds_id)
        #Initialize the genome info for the specific genomic sequence if it does not exist yet
        if (scaffold_id not in seen_scaff_ids):
            seen_scaff_ids[scaffold_id] = True
            genome_seq_info['CDS_Count'][scaffold_id] = 0
            genome_seq_info['KEGG_Matched_Subject_IDs'][scaffold_id] = set()
            genome_seq_info['KEGG_Matched_Functions'][scaffold_id] = set()
            genome_seq_info['KEGG_Matched_Pathways'][scaffold_id] = set()
            genome_seq_info['KEGG_Matched_Modules'][scaffold_id] = set()
            genome_seq_info['Pfam_Matched_Subjects'][scaffold_id] = set()
        #Increment cds count of the scaffold
        genome_seq_info['CDS_Count'][scaffold_id] += 1
        seq_counter += 1
        if (seq_counter % 100000 == 0):
            print(f"Processed {seq_counter} sequences")

def calc_kegg_abundance_from_gene(abundance_file=None,genome_mode=False):
    gene_abund_df = pd.read_csv(abundance_file, sep='\t',index_col='Sequence',header=0)

    ko_abund = defaultdict(lambda: defaultdict(int))
    met_abund = defaultdict(lambda: defaultdict(int))
    path_abund = defaultdict(lambda: defaultdict(int))
    
    #Iterate over all gene ids in the abundance file
    for gene_id in gene_abund_df.index:
        #Get the KEGG KO id for the gene
        ko = cds_seq_info['KEGG_Best_Subject'][gene_id]
        #Iterate over all samples in the abundance file
        for sample in gene_abund_df.columns:
            #Add the gene abundance to the KO abundance
            ko_abund[sample][ko] += gene_abund_df.loc[gene_id,sample]
            #Get the set of KEGG pathways for the gene
            path_set = set(cds_seq_info['KEGG_Best_Subject_Pathways'][gene_id])
            #Get the number of pathways for the gene
            path_count = len(path_set)
            #Iterate over all pathways for the gene
            for path in path_set:
                #Add the gene abundance to the pathway abundance
                path_abund[sample][path] += (gene_abund_df.loc[gene_id,sample] / path_count)
            #Get the set of KEGG modules for the gene
            met_set = set(cds_seq_info['KEGG_Best_Subject_Modules'][gene_id])
            #Get the number of modules for the gene
            met_count = len(met_set)
            #Iterate over all modules for the gene
            for met in met_set:
                #Add the gene abundance to the module abundance
                met_abund[sample][met] += (gene_abund_df.loc[gene_id,sample] / met_count)

    ko_abund_df = pd.DataFrame.from_dict(ko_abund)
    ko_abund_df.index.name = 'KO'
    ko_abund_df.to_csv('KO_Abundance_from_Gene.tsv',sep="\t",na_rep=0)
    
    met_abund_df = pd.DataFrame.from_dict(met_abund)
    met_abund_df.index.name = 'Metabolism'
    met_abund_df.to_csv('Metabolism_Abundance_from_Gene.tsv',sep="\t",na_rep=0)

    path_abund_df = pd.DataFrame.from_dict(path_abund)
    path_abund_df.index.name = 'Pathway'
    path_abund_df.to_csv('Pathway_Abundance_from_Gene.tsv',sep="\t",na_rep=0)

def calc_kegg_abundance_from_genome(abundance_file=None,genome_mode=True):
    abund_df = pd.read_csv(abundance_file, sep='\t',index_col='Sequence',header=0)

    ko_abund = defaultdict(lambda: defaultdict(int))
    met_abund = defaultdict(lambda: defaultdict(int))
    path_abund = defaultdict(lambda: defaultdict(int))
    
    print("Calculating functional abundances")
    for scaffold,row in abund_df.iterrows():
        for cds_num in range(1,(genome_seq_info['CDS_Count'][scaffold] + 1)):
            cds_id = scaffold+"_"+str(cds_num)
            if (cds_id in cds_seq_info['KEGG_Best_Subject']):
                ko = cds_seq_info['KEGG_Best_Subject'][cds_id]
                for sample in abund_df.columns:
                    ko_abund[sample][ko] += abund_df.loc[scaffold,sample]
                    path_set = set(cds_seq_info['KEGG_Best_Subject_Pathways'][cds_id])
                    path_count = len(path_set)
                    for path in path_set:
                        path_abund[sample][path] += (abund_df.loc[scaffold,sample] / path_count)
                    met_set = set(cds_seq_info['KEGG_Best_Subject_Modules'][cds_id])
                    met_count = len(met_set)
                    for met in met_set:
                        met_abund[sample][met] += (abund_df.loc[scaffold,sample] / met_count)
                
    ko_abund_df = pd.DataFrame.from_dict(ko_abund)
    ko_abund_df.index.name = 'KO'
    ko_abund_df.to_csv('KO_Abundance.tsv',sep="\t",na_rep=0)
    
    met_abund_df = pd.DataFrame.from_dict(met_abund)
    met_abund_df.index.name = 'Metabolism'
    met_abund_df.to_csv('Metabolism_Abundance.tsv',sep="\t",na_rep=0)

    path_abund_df = pd.DataFrame.from_dict(path_abund)
    path_abund_df.index.name = 'Pathway'
    path_abund_df.to_csv('Pathway_Abundance.tsv',sep="\t",na_rep=0)

def build_phylogeny_from_kegg_hits(cds_file=None,ref_file=None,ko=None,ko_hmm_dir=None):
    print("Building phylogeny from KEGG hits")
    if(args.parse_only == False):
        print(f"Making copy of {ko}.hmm ")
        command = f"cp {ko_hmm_dir}/{ko}.hmm ."
        subprocess.call(command, shell=True)
        print(f"Pressing {ko}.hmm ")
        command = f"hmmpress {ko}.hmm"
        subprocess.call(command, shell=True)

    list_protein_files = []
    if (ref_file):
        list_protein_files = [cds_file,ref_file]
    else:
        list_protein_files = [cds_file]
    
    list_nr_derep_protein_files = []
    for protein_file in list_protein_files:
        protein_vs_ko_outfile = call_hmmer(cds_file=protein_file,db_file=f"{ko}.hmm",program="hmmsearch")
        (protein_ko_genome_hmm_scores,protein_ko_pairwise_scores) = parse_hmmer_output(hmmer_out_file=protein_vs_ko_outfile,min_score=args.kegg_min_score,max_evalue=args.kegg_max_evalue,multi_hits=False,prefix="KEGG",print_df=False)
        #matched_cds_ids = list(protein_ko_pairwise_scores['Subject'].values())
        matched_cds_ids = {protein:True for protein in protein_ko_pairwise_scores['Subject'].values()}
        #print(f"Fetching: {matched_cds_ids}")
        sub_protein_file = f"Matched_{ko}_"+get_prefix(protein_file,"(fasta)|(faa)|(fa)")+".fasta"
        subset_seqs(input_file=protein_file,output_file=sub_protein_file,ids_dict=matched_cds_ids,min_length=args.min_seq_len,max_length=args.max_seq_len)
        nr_sub_protein_file = "NR_"+get_prefix(sub_protein_file,"fasta")+".fasta"
        derep_seqs(input_file=sub_protein_file,output_file=nr_sub_protein_file,cluster_id=args.derep_id)
        list_nr_derep_protein_files.append(nr_sub_protein_file)
    
    merged_nr_matched_file = f"Merged_NR_Matched_{ko}.faa"
    command = f"cat "+" ".join(list_nr_derep_protein_files)+f" > {merged_nr_matched_file}"
    subprocess.call(command, shell=True)
    
    aln_merged_nr_matched_file = "Aligned_"+merged_nr_matched_file
    align_seqs(input_file=merged_nr_matched_file,output_file=aln_merged_nr_matched_file)
    
    tree_file = aln_merged_nr_matched_file+".newick"
    build_tree(input_file=aln_merged_nr_matched_file,output_file=tree_file)

def build_tree(input_file=None,output_file=None):
    command = f"FastTreeMP -nosupport -out {output_file} {input_file}"
    subprocess.call(command, shell=True)  
    
def align_seqs(input_file=None,output_file=None):
    command = f"muscle -in {input_file} -out {output_file} -maxhours 24 -maxiters 3"
    subprocess.call(command, shell=True)    

def derep_seqs(input_file=None,output_file=None,cluster_id=1):
    command = f"cd-hit -i {input_file} -o {output_file} -c {cluster_id} -n 5 -M 100000 -T {args.threads}"
    subprocess.call(command, shell=True)


def enrich_m8(m8_file=None):
    m8_file_prefix = get_prefix(m8_file,'blastp')
    outfile = f"Enriched_{m8_file_prefix}.blastp"
    command = f"perl /mnt/lustre/bio/users/fcoutinho/Scripts/Enrich_m8_Table.pl --names /mnt/lustre/repos/bio/databases/public/NCBI/taxonomy/taxonomy_2021-04-06/names.dmp --nodes /mnt/lustre/repos/bio/databases/public/NCBI/taxonomy/taxonomy_2021-04-06/nodes.dmp --acc2taxid_dir /mnt/lustre/repos/bio/databases/public/UniProtKB/UniProt_Release_2021-04-07/ --description /mnt/lustre/bio/users/fcoutinho/Databases/UniRef100/UniRef100_Release_2021_04_07_Id_to_Desc.txt --input {m8_file} --output {outfile}"
    if (args.parse_only == False):
        pass
        subprocess.call(command, shell=True)
    return outfile



def call_diamond (cds_file="NA",db_file=args.uniref_db,program="blastp"):
    cds_file_prefix = get_prefix(cds_file,'(faa)|(fasta)|(fa)')
    db_file_prefix = get_prefix(db_file,'dmnd')
    outfile = 'NA'
    if (program == 'blastp'):
        outfile = cds_file_prefix+'x'+db_file_prefix+'.blastp'
        if (args.parse_only == False):
            print(f'Querying {cds_file} against {db_file}')
            command = f'diamond blastp --index-chunks 1 --matrix BLOSUM45 --threads {args.threads} --more-sensitive --db {db_file} --outfmt 6 --query {cds_file} --max-target-seqs 10 --evalue 0.001 --out {outfile}'
            subprocess.call(command, shell=True)    
    else:
        print('Not a valid DIAMOND program!')
        
    return outfile


def call_hmmer (cds_file="NA",db_file="NA",program="hmmscan"):
    cds_file_prefix = get_prefix(cds_file,'(faa)|(fasta)|(fa)')
    db_file_prefix = get_prefix(db_file,'hmm')
    outfile = 'NA'
    if (program == 'hmmscan'):
        outfile = cds_file_prefix+'x'+db_file_prefix+'.hmmscan'
        if (args.parse_only == False):
            print(f'Querying {db_file} against {cds_file}')
            command = f'hmmscan -o {outfile} --noali --cpu {args.threads} {db_file} {cds_file}'
            subprocess.call(command, shell=True)
    elif (program == 'hmmsearch'):
        outfile = cds_file_prefix+'x'+db_file_prefix+'.hmmsearch'
        if (args.parse_only == False):
            print(f'Querying {cds_file} against {db_file}')
            command = f'hmmsearch -o {outfile} --noali --cpu {args.threads} {db_file} {cds_file}'
            subprocess.call(command, shell=True)    
    else:
        print('Not a valid Hmmer program!')
        
    return outfile

def parse_uniref_raw(uniref_hits_file="NA",min_score=50,max_evalue=0.00001):
    print(f"Reading {uniref_hits_file}")
    with open(uniref_hits_file) as UNIREF_FH:
        compiled_obj = re.compile('(Uncharacterized)|(hypothetical)')
        for line in UNIREF_FH:
            line = line.rstrip("\n")
            [cds,subject,perc_id,ali_len,mismatches,gap_open,qstart,qend,sstarts,send,evalue,bitscore] = line.split("\t")
            bitscore = float(bitscore)
            evalue = float(evalue)
            if ((bitscore >= min_score) and (evalue <= max_evalue)):
                if (cds not in cds_seq_info['UniRef_Best_Hit_ID']):
                    cds_seq_info['UniRef_Best_Hit_ID'][cds] = subject
                    cds_seq_info['UniRef_Best_Hit_Score'][cds] = bitscore  
                if (bitscore > cds_seq_info['UniRef_Best_Hit_Score'][cds]):
                    cds_seq_info['UniRef_Best_Hit_ID'][cds] = subject
                    cds_seq_info['UniRef_Best_Hit_Score'][cds] = bitscore
                if (not re.search(compiled_obj,subject)):
                    if (cds not in cds_seq_info['UniRef_Best_NH_Hit_ID']):
                        cds_seq_info['UniRef_Best_Non_Hypothetical_Hit_ID'][cds] = subject
                        cds_seq_info['UniRef_Best_Non_Hypothetical_Hit_Score'][cds] = bitscore
                    if (bitscore > cds_seq_info['UniRef_Best_Non_Hypothetical_Hit_Score'][cds]):
                        cds_seq_info['UniRef_Best_Non_Hypothetical_Hit_ID'][cds] = subject
                        cds_seq_info['UniRef_Best_Non_Hypothetical_Hit_Score'][cds] = bitscore


def parse_uniref_enriched(uniref_hits_file="NA",min_score=50,max_evalue=0.00001):
    print(f"Reading {uniref_hits_file}")
    with open(uniref_hits_file) as UNIREF_FH:
        compiled_obj = re.compile('(Uncharacterized)|(hypothetical)')
        for line in UNIREF_FH:
            line = line.rstrip("\n")
            [cds,subject,perc_id,ali_len,mismatches,gap_open,qstart,qend,sstarts,send,evalue,bitscore,sdesc,staxid,sdomain,sphylum,sclass,sorder,sfamily,sgenus,sspecies] = line.split("\t")
            bitscore = float(bitscore)
            evalue = float(evalue)
            if ((bitscore >= min_score) and (evalue <= max_evalue)):
                if (cds not in cds_seq_info['UniRef_Best_Hit_ID']):
                    cds_seq_info['UniRef_Best_Hit_ID'][cds] = subject
                    cds_seq_info['UniRef_Best_Hit_Description'][cds] = sdesc
                    cds_seq_info['UniRef_Best_Hit_Score'][cds] = bitscore
                    cds_seq_info['UniRef_Best_Hit_Domain'][cds] = sdomain
                    cds_seq_info['UniRef_Best_Hit_Phylum'][cds] = sphylum
                    cds_seq_info['UniRef_Best_Hit_Class'][cds] = sclass
                    cds_seq_info['UniRef_Best_Hit_Order'][cds] = sorder
                    cds_seq_info['UniRef_Best_Hit_Family'][cds] = sfamily
                    cds_seq_info['UniRef_Best_Hit_Genus'][cds] = sgenus
                    cds_seq_info['UniRef_Best_Hit_Species'][cds] = sspecies
                if (bitscore > cds_seq_info['UniRef_Best_Hit_Score'][cds]):
                    cds_seq_info['UniRef_Best_Hit_ID'][cds] = subject
                    cds_seq_info['UniRef_Best_Hit_Description'][cds] = sdesc
                    cds_seq_info['UniRef_Best_Hit_Score'][cds] = bitscore
                    cds_seq_info['UniRef_Best_Hit_Domain'][cds] = sdomain
                    cds_seq_info['UniRef_Best_Hit_Phylum'][cds] = sphylum
                    cds_seq_info['UniRef_Best_Hit_Class'][cds] = sclass
                    cds_seq_info['UniRef_Best_Hit_Order'][cds] = sorder
                    cds_seq_info['UniRef_Best_Hit_Family'][cds] = sfamily
                    cds_seq_info['UniRef_Best_Hit_Genus'][cds] = sgenus
                    cds_seq_info['UniRef_Best_Hit_Species'][cds] = sspecies
                if (not re.search(compiled_obj,sdesc)):
                    if (cds not in cds_seq_info['UniRef_Best_NH_Hit_ID']):
                        cds_seq_info['UniRef_Best_NH_Hit_ID'][cds] = subject
                        cds_seq_info['UniRef_Best_NH_Hit_Description'][cds] = sdesc
                        cds_seq_info['UniRef_Best_NH_Hit_Score'][cds] = bitscore
                        cds_seq_info['UniRef_Best_NH_Hit_Domain'][cds] = sdomain
                        cds_seq_info['UniRef_Best_NH_Hit_Phylum'][cds] = sphylum
                        cds_seq_info['UniRef_Best_NH_Hit_Class'][cds] = sclass
                        cds_seq_info['UniRef_Best_NH_Hit_Order'][cds] = sorder
                        cds_seq_info['UniRef_Best_NH_Hit_Family'][cds] = sfamily
                        cds_seq_info['UniRef_Best_NH_Hit_Genus'][cds] = sgenus
                        cds_seq_info['UniRef_Best_NH_Hit_Species'][cds] = sspecies
                    if (bitscore > cds_seq_info['UniRef_Best_NH_Hit_Score'][cds]):
                        cds_seq_info['UniRef_Best_NH_Hit_ID'][cds] = subject
                        cds_seq_info['UniRef_Best_NH_Hit_Description'][cds] = sdesc
                        cds_seq_info['UniRef_Best_NH_Hit_Score'][cds] = bitscore
                        cds_seq_info['UniRef_Best_NH_Hit_Domain'][cds] = sdomain
                        cds_seq_info['UniRef_Best_NH_Hit_Phylum'][cds] = sphylum
                        cds_seq_info['UniRef_Best_NH_Hit_Class'][cds] = sclass
                        cds_seq_info['UniRef_Best_NH_Hit_Order'][cds] = sorder
                        cds_seq_info['UniRef_Best_NH_Hit_Family'][cds] = sfamily
                        cds_seq_info['UniRef_Best_NH_Hit_Genus'][cds] = sgenus
                        cds_seq_info['UniRef_Best_NH_Hit_Species'][cds] = sspecies


def parse_uniref_old(uniref_hits_file="NA",min_score=50,max_evalue=0.00001):
    print(f"Reading {uniref_hits_file}")
    uniref_df = pd.read_csv(uniref_hits_file, sep="\t",index_col=False,header=None)
    print(f"Processing {uniref_hits_file}")
    compiled_obj = re.compile('(Uncharacterized)|(hypothetical)')
    for i,row in uniref_df.iterrows():
        [cds,subject,perc_id,ali_len,mismatches,gap_open,qstart,qend,sstarts,send,evalue,bitscore,sdesc,staxid,sdomain,sphylum,sclass,sorder,sfamily,sgenus,sspecies] = row
        if ((bitscore >= min_score) and (evalue <= max_evalue)):
            if (cds not in cds_seq_info['UniRef_Best_Hit_ID']):
                cds_seq_info['UniRef_Best_Hit_ID'][cds] = subject
                cds_seq_info['UniRef_Best_Hit_Description'][cds] = sdesc
                cds_seq_info['UniRef_Best_Hit_Score'][cds] = bitscore
                cds_seq_info['UniRef_Best_Hit_Domain'][cds] = sdomain
                cds_seq_info['UniRef_Best_Hit_Phylum'][cds] = sphylum
                cds_seq_info['UniRef_Best_Hit_Class'][cds] = sclass
                cds_seq_info['UniRef_Best_Hit_Order'][cds] = sorder
                cds_seq_info['UniRef_Best_Hit_Family'][cds] = sfamily
                cds_seq_info['UniRef_Best_Hit_Genus'][cds] = sgenus
                cds_seq_info['UniRef_Best_Hit_Species'][cds] = sspecies
            if (bitscore > cds_seq_info['UniRef_Best_Hit_Score'][cds]):
                cds_seq_info['UniRef_Best_Hit_ID'][cds] = subject
                cds_seq_info['UniRef_Best_Hit_Description'][cds] = sdesc
                cds_seq_info['UniRef_Best_Hit_Score'][cds] = bitscore
                cds_seq_info['UniRef_Best_Hit_Domain'][cds] = sdomain
                cds_seq_info['UniRef_Best_Hit_Phylum'][cds] = sphylum
                cds_seq_info['UniRef_Best_Hit_Class'][cds] = sclass
                cds_seq_info['UniRef_Best_Hit_Order'][cds] = sorder
                cds_seq_info['UniRef_Best_Hit_Family'][cds] = sfamily
                cds_seq_info['UniRef_Best_Hit_Genus'][cds] = sgenus
                cds_seq_info['UniRef_Best_Hit_Species'][cds] = sspecies
            if (not re.search(compiled_obj,sdesc)):
                if (cds not in cds_seq_info['UniRef_Best_NH_Hit_ID']):
                    cds_seq_info['UniRef_Best_NH_Hit_ID'][cds] = subject
                    cds_seq_info['UniRef_Best_NH_Hit_Description'][cds] = sdesc
                    cds_seq_info['UniRef_Best_NH_Hit_Score'][cds] = bitscore
                    cds_seq_info['UniRef_Best_NH_Hit_Domain'][cds] = sdomain
                    cds_seq_info['UniRef_Best_NH_Hit_Phylum'][cds] = sphylum
                    cds_seq_info['UniRef_Best_NH_Hit_Class'][cds] = sclass
                    cds_seq_info['UniRef_Best_NH_Hit_Order'][cds] = sorder
                    cds_seq_info['UniRef_Best_NH_Hit_Family'][cds] = sfamily
                    cds_seq_info['UniRef_Best_NH_Hit_Genus'][cds] = sgenus
                    cds_seq_info['UniRef_Best_NH_Hit_Species'][cds] = sspecies
                if (bitscore > cds_seq_info['UniRef_Best_NH_Hit_Score'][cds]):
                    cds_seq_info['UniRef_Best_NH_Hit_ID'][cds] = subject
                    cds_seq_info['UniRef_Best_NH_Hit_Description'][cds] = sdesc
                    cds_seq_info['UniRef_Best_NH_Hit_Score'][cds] = bitscore
                    cds_seq_info['UniRef_Best_NH_Hit_Domain'][cds] = sdomain
                    cds_seq_info['UniRef_Best_NH_Hit_Phylum'][cds] = sphylum
                    cds_seq_info['UniRef_Best_NH_Hit_Class'][cds] = sclass
                    cds_seq_info['UniRef_Best_NH_Hit_Order'][cds] = sorder
                    cds_seq_info['UniRef_Best_NH_Hit_Family'][cds] = sfamily
                    cds_seq_info['UniRef_Best_NH_Hit_Genus'][cds] = sgenus
                    cds_seq_info['UniRef_Best_NH_Hit_Species'][cds] = sspecies
            
def check_hmmer_match_cutoff(hsp,max_evalue,min_bitscore):
    if ((hsp.bitscore < min_bitscore) or (hsp.evalue > max_evalue)):
        return False
    else:
        return True


def print_results(cds_output_df_file="CDS_Info.tsv",genome_output_df_file="Genome_Info.tsv"):
    print(f"Printing CDS info to {cds_output_df_file}")
    cds_info_df = pd.DataFrame.from_dict(cds_seq_info)
    cds_info_df.index.name = 'Sequence'
    cds_info_df.to_csv(cds_output_df_file,sep="\t",na_rep='NA')
    print(f"Printing Genome info to {genome_output_df_file}")
    genome_info_df = pd.DataFrame.from_dict(genome_seq_info)
    genome_info_df.index.name = 'Sequence'
    genome_info_df.to_csv(genome_output_df_file,sep="\t",na_rep="NA")

    
def get_prefix(file,format):
    prefix_file = re.sub(f'.{format}','',file)
    prefix_file = re.sub('(.)+/','',prefix_file)
    return prefix_file

def load_kegg_info(json_file=None):
    print("Loading KEGG Info")
    kegg_info =  defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(dict)))))
    
    json_fh = open(json_file,)
    kegg_data = json.load(json_fh)

    valid_categs = ("09100 Metabolism","09120 Genetic Information Processing","09130 Environmental Information Processing","09140 Cellular Processes")

    comp_obj_1 = re.compile("^(\d)+\s")
    comp_obj_2 = re.compile("\s\[(.)+\]$")

    for category in kegg_data["children"]:
        categ_name = category['name']
        if (categ_name in valid_categs):
            print("Indexing info for category",category['name'])
            for hier1 in category['children']:
                #print("\tLevel 1",hier1['name'])
                h1_name = hier1['name']
                if (hier1['children']):
                    for hier2 in hier1['children']:
                        h2_name = hier2['name']
                        #print("\t\tLevel 2",hier2['name'])
                        if (hier2['children']):
                            for hier3 in hier2['children']:
                                #print("\t\t\tLevel 3",hier3['name'])
                                h3_name = hier3['name']
                                (ko,desc) = hier3['name'].split(' ',1)
                                #print("\t\t\tKO",ko,"Description",desc)
                                h1_name = re.sub(comp_obj_1,"",h1_name)
                                h2_name = re.sub(comp_obj_1,"",h2_name)
                                h2_name = re.sub(comp_obj_2,"",h2_name)
                                kegg_info[ko][categ_name][h1_name][h2_name][h3_name]["Match"] = True
                                #use regular expression to remove leading spaces from desc                                                    
                                desc = re.sub("^(\s)","",desc)
                                kegg_info[ko][categ_name][h1_name][h2_name][h3_name]["Description"] = desc
      

    json_fh.close()
    return(kegg_info)

def subset_seqs(input_file=None,input_format="fasta",output_file=None,output_format="fasta",ids_dict=None,min_length=0,max_length=999999,protein=False):
    total_seq_counter = 0
    passed_seq_counter = 0
    shouter = 1
    compiled_obj = re.compile('_(\d)+$')
    with open(output_file, 'w', newline='') as OUT:
        for seqobj in SeqIO.parse(input_file, input_format):
            total_seq_counter += 1
            if (total_seq_counter == shouter):
                print('Processed',total_seq_counter,'sequences')
                shouter += total_seq_counter
            seqid = seqobj.id
            if (protein):
                seqid = re.sub(compiled_obj,'',seqid)
            if (seqid in ids_dict):
                seq_length = len(seqobj.seq)
                if ((seq_length >= min_length) and (seq_length <= max_length)):
                    passed_seq_counter += 1
                    SeqIO.write(seqobj, OUT, output_format)
    print('Processed',total_seq_counter,'Sequences.',passed_seq_counter,'passed.')
    
    
#2D dictionaries to store all relevant information about Genomic and CDS sequences
cds_seq_info = defaultdict(dict)
genome_seq_info = defaultdict(dict)
#Multidimensional dictionary that stores KEGG annotation data
kegg_info = load_kegg_info(json_file=args.kegg_json)

central()
