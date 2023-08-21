#! /usr/bin/env python3
from Bio import SearchIO
from collections import defaultdict
import pandas as pd
import argparse
import re
import json

parser = argparse.ArgumentParser()
parser.add_argument("--abundance", help="Genomic Sequence Abundance matrix", type=str)
parser.add_argument("--kegg_json", help=".json file from KEGG", type =str, default="/mnt/lustre/bio/users/fcoutinho/KOfam/ko00001.json")
parser.add_argument("--kegg_hits_file", help="The KEGG Enriched pairwise table generated by Virathon + Enrich_KO", type =str)
parser.add_argument("--kegg_min_score", help="The minimum score to consider a KEGG hit", type =int, default=50)
parser.add_argument("--kegg_max_evalue", help="The maximum e-value to consider a KEGG hit", type =float, default=0.00001)
args = parser.parse_args()

def load_kegg_info(json_file=None):
    print("Loading KEGG Info")
    kegg_info =  defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(dict)))))
    
    json_fh = open(json_file,)
    kegg_data = json.load(json_fh)

    valid_categs = ("09100 Metabolism","09120 Genetic Information Processing","09130 Environmental Information Processing","09140 Cellular Processes")

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
                                kegg_info[ko][categ_name][h1_name][h2_name][h3_name]["Match"] = True
                                kegg_info[ko][categ_name][h1_name][h2_name][h3_name]["Description"] = desc
      

    json_fh.close()
    return(kegg_info)

def parse_hmmer_output(hmmer_out_file,max_evalue=0.00001,min_score=50,multi_hits=False,prefix="NA"):
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
                    pairwise_scores['Query_Description'][hsp_count] = qresult.description
                    genome_ori_score = genome_hmm_scores[qresult.id].get(genome,0)
                    if (genome_ori_score < hsp.bitscore):
                        genome_hmm_scores[qresult.id][genome] = hsp.bitscore
                    if (multi_hits == True):
                        if (hit.id not in cds_seq_info[prefix+'_Subjects']):
                            cds_seq_info[prefix+'_Subjects'][hit.id] = []
                        cds_seq_info[prefix+'_Subjects'][hit.id].append(qresult.id)
                    else:
                        if (hit.id not in cds_seq_info[prefix+'_Best_Subject']):
                            cds_seq_info[prefix+'_Best_Subject'][hit.id] = "NA"
                            cds_seq_info[prefix+'_Best_Score'][hit.id] = 0
                            cds_seq_info[prefix+'_Best_Subject_Function'][hit.id] = "NA"
                            cds_seq_info[prefix+'_Best_Subject_Pathways'][hit.id] = []
                        elif(hsp.bitscore > cds_seq_info[prefix+'_Best_Score'][hit.id]):
                            cds_seq_info[prefix+'_Best_Subject'][hit.id] = qresult.id
                            cds_seq_info[prefix+'_Best_Score'][hit.id] = hsp.bitscore
                            ko = qresult.id
                            for categ in kegg_info[ko].keys():
                                for h1_name in kegg_info[ko][categ].keys():
                                    for h2_name in kegg_info[ko][categ][h1_name].keys():
                                        for h3_name in kegg_info[ko][categ][h1_name][h2_name].keys():
                                            desc = kegg_info[ko][categ][h1_name][h2_name][h3_name]["Description"]
                                            cds_seq_info[prefix+'_Best_Subject_Function'][hit.id] = desc
                                            cds_seq_info[prefix+'_Best_Subject_Pathways'][hit.id].append(h2_name)
                            
    return(genome_hmm_scores,pairwise_scores)

def check_hmmer_match_cutoff(hsp,max_evalue,min_bitscore):
    if ((hsp.bitscore < min_bitscore) or (hsp.evalue > max_evalue)):
        return False
    else:
        return True

genome_seq_info = defaultdict(dict)
cds_seq_info = defaultdict(dict)

kegg_info = load_kegg_info(json_file=args.kegg_json)
abund_df = pd.read_csv(args.abundance, sep='\t',index_col='Sequence',header=0)

(kegg_genome_hmm_scores,kegg_pairwise_scores) = parse_hmmer_output(hmmer_out_file=args.kegg_hits_file,min_score=args.kegg_min_score,max_evalue=args.kegg_max_evalue,multi_hits=False,prefix="KEGG")

print("Assigning best hit KOs to genomic sequences")
for cds in cds_seq_info.keys():
    genome = re.sub("_(\d)+$","",cds)
    best_ko = cds_seq_info['KEGG_Best_Subject'][cds]
    if genome not in genome_seq_info:
        genome_seq_info[genome]['KOs'] = []
    genome_seq_info[genome]['KOs'].append(best_ko)

ko_abund = defaultdict(lambda: defaultdict(int))
met_abund = defaultdict(lambda: defaultdict(int))
path_abund = defaultdict(lambda: defaultdict(int))

print("Calculating functional abundances")
for scaffold,row in abund_df.iterrows():
    ko_list = genome_seq_info[genome]['KOs']
    for sample in abund_df.columns:
        ko_abund[sample][ko] += abund_df.loc[scaffold,sample]
        path_count = len(cds_seq_info[prefix+'_Best_Subject_Pathways'][cds])
        

ko_abund_df = pd.DataFrame.from_dict(ko_abund)
ko_abund_df.index.name = 'KO'
ko_abund_df.to_csv('KO_Abundance.tsv',sep="\t",na_rep=0)

met_abund_df = pd.DataFrame.from_dict(met_abund)
met_abund_df.index.name = 'Metabolism'
met_abund_df.to_csv('Metabolism_Abundance.tsv',sep="\t",na_rep=0)

path_abund_df = pd.DataFrame.from_dict(path_abund)
path_abund_df.index.name = 'Pathway'
path_abund_df.to_csv('Pathway_Abundance.tsv',sep="\t",na_rep=0)