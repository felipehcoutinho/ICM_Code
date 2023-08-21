import numpy as np
import pandas as pd
import subprocess
import csv
import re
import argparse
import gensim
import multiprocessing as mp
from joblib import Parallel, delayed
from sklearn.cluster import DBSCAN
from sklearn.cluster import KMeans
from sklearn.metrics.pairwise import pairwise_distances
from sklearn.metrics import adjusted_rand_score
from collections import defaultdict
from mango import scheduler, Tuner
from scipy.stats import uniform

parser = argparse.ArgumentParser()
parser.add_argument("--map_file", help="The genome map file", type =str)
parser.add_argument("--min_cluster_size", help="The minimum number of proteins in a MMseqs cluster", default=3, type=int)
parser.add_argument("--tax_distance_file", help="The optional OG x OG taxonomic distance .tsv file to be used when performing correction for taxonomic distance", type =str)
parser.add_argument("--classification_file", help="The optional OG x Classification .tsv file to be used when performing ARI measures", type =str)
parser.add_argument("--classification_variable", help="The column of the classification file to be used to measure ARI of superclusters", type =str)
parser.add_argument("--unknown_value", help="The string on unknown values ", type =str, default = 'NA')
parser.add_argument("--threads", help="The number of threads to be used", default=1, type=int)
parser.add_argument("--dist_metric", help="Distance metric to be used", default='cosine', type=str)
parser.add_argument("--dbscan_algorithm", help="DBSCAN algorithm to be used", default='auto', type=str)

args = parser.parse_args()

def read_genome_structure():
    """This method reads the input file which is plain text format"""
    #print("Reading {0}".format(args.map_file))
    with open (args.map_file, 'r') as f:
        for i, line in enumerate (f): 
            if (i%100==0):
                pass #print("Read {0} entries".format (i))
            #yield gensim.utils.simple_preprocess(line)
            yield list(gensim.utils.tokenize(line,deacc=False,lowercase=False))

def fit_model(eps_i,min_sample_i,alpha,dist_metric,window_i,transformed_cluster_distance,sub_cluster_class_info,scaled_data):
    epsilon = eps_i
    iter_uid = 'Epsilon_'+str(epsilon)+'_Min_Samples_'+str(min_sample_i)+'_Alpha_'+str(alpha)+'_Correction_'+'NA'+'_Window_Size_'+str(window_i)+'_Metric_'+str(dist_metric)
    dbscan = DBSCAN(eps=epsilon,n_jobs=1,metric='precomputed',min_samples=min_sample_i,algorithm=args.dbscan_algorithm)
    dbscan.fit(transformed_cluster_distance)
    
    #kmeans_kwargs = {"init": "random","n_init": 500,"max_iter": 300,"random_state": 42,}
    #kmeans = KMeans(n_clusters=25, **kmeans_kwargs)
    #Perform k-means clustering on the scaled data
    #kmeans.fit(scaled_data)
    #cluster_labels = kmeans.labels_

    cluster_labels = dbscan.labels_
    
    ari = 'NA'
    if args.classification_file:
        valid_posits = sub_cluster_class_info[args.classification_variable].notnull()
        #valid_posits = (sub_cluster_class_info[args.classification_variable] != NaN)
        #print(valid_posits)
        #print('Valids:',valids.value_counts())
        #na_count_cci = np.count_nonzero(np.isnan(cluster_class_info[valids]))
        #na_count_dsl = np.count_nonzero(np.isnan(dbscan.labels_[valids]))
        #print('NAs CCI:',na_count_cci,'DSL:',na_count_dsl)
        sliced_sub_cluster_class_info = sub_cluster_class_info[args.classification_variable]
        cci_valids = sliced_sub_cluster_class_info[valid_posits]
        dsl_valids = cluster_labels[valid_posits]
        #print('CCI valids:',len(cci_valids),'DSL valids:',len(dsl_valids))
        #print(sliced_sub_cluster_class_info,cci_valids)
        ari = adjusted_rand_score(cci_valids,dsl_valids)
    
    return (iter_uid, epsilon, len(set(cluster_labels)), np.count_nonzero(cluster_labels == -1), np.median(cluster_labels),ari)

def get_alias(protein_cluster):
    num_dict = {'0':'A','1':'B','2':'C','3':'D','4':'E','5':'F','6':'G','7':'H','8':'I','9':'J','.':'_'}
    protein_cluster_alias = list(protein_cluster)
    for char_posit in range(0,len(protein_cluster_alias)):
        char = protein_cluster_alias[char_posit]
        if char in num_dict.keys():
            protein_cluster_alias[char_posit] = num_dict[char]
    alias = ''.join(protein_cluster_alias)
    rev_alias[alias] = protein_cluster
    return alias
    
def get_class_info():
    print('Reading',args.classification_file)
    info = pd.read_csv(args.classification_file,sep='\t',na_values='',index_col=0,header=0)
    return info

@scheduler.parallel(n_jobs=args.threads)
def make_super_clusters(window_i=5,alpha_i=0,eps_i=1):
    #(window_i,alpha_i,correct,min_sample_i,eps_i) = params
    #(window_i,alpha_i,correct,eps_i) = params
    genome_tokens = list (read_genome_structure())
    
    #print('Setting up Word2Vec model')
    model = gensim.models.Word2Vec (genome_tokens, vector_size=100, window=window_i, min_count=3, workers=1, sg=1)
    #print('Training Word2Vec model with window size',window_i)
    model.train(genome_tokens,total_examples=len(genome_tokens),epochs=10)
        

    alpha = alpha_i
    #print('Using alpha',alpha_i)

    #dist_metric = 'manhattan'
    #print('Calculating distances using',dist_metric)
    words = model.wv.index_to_key
    scaled_data = model.wv.get_normed_vectors()
    oids = [rev_alias[w] for w in words]
    sub_cluster_class_info =  cluster_class_info.loc[oids,]
    sub_cluster_class_info = sub_cluster_class_info.reindex(oids)
    #print('Original shape',cluster_class_info.shape,'Subset shape',sub_cluster_class_info.shape)
    #scaled_data = [model.vw[w] for w in words]
    #print('Calculating distribution distance among clusters')
    cluster_distrib_distance = pairwise_distances(scaled_data, metric=args.dist_metric)
    #print('Dimensions of cluster distribution distance data frame:',cluster_distrib_distance.shape)
    transformed_cluster_distance = cluster_distrib_distance
    if alpha > 0:
        sub_cluster_tax_distance = cluster_tax_distance.loc[oids,oids]
        #print('Dimensions of cluster taxonomic distance data frame:',cluster_tax_distance.shape)
        transformed_cluster_distance = ((cluster_distrib_distance) / (alpha * (sub_cluster_tax_distance + 1)))
        #print('Dimensions of transformed cluster distance data frame:',transformed_cluster_distance.shape)
    
    #print('NaN check:',transformed_cluster_distance.isnull().values.any())    
    #print('Training DBSCAN models')
    results = fit_model(eps_i,2,alpha,args.dist_metric,window_i,transformed_cluster_distance,sub_cluster_class_info,scaled_data)
    iter_uid = results[0]
    clustering_info['Epsilon'][iter_uid] = results[1]
    clustering_info['Cluster_Count'][iter_uid] = results[2]
    clustering_info['Outlier_Count'][iter_uid] = results[3]
    clustering_info['Median_Cluster_Num'][iter_uid] = results[4]
    #clustering_info['Min_Samples'][iter_uid] = 2
    #clustering_info['Distance_Correction'][iter_uid] = correct
    clustering_info['Distance_Metric'][iter_uid] = args.dist_metric
    clustering_info['Alpha'][iter_uid] = alpha
    clustering_info['Window_Size'][iter_uid] = window_i
    clustering_info['ARI'][iter_uid] = results[5]
    #print('Alpha',alpha,'Window_Size',window_i,'Epsilon',eps_i,'ARI',results[-1],'Cluster_Count',results[2],'Outlier_Count',results[3],'Median_Cluster_Num',results[4])
    
    with open('Mango_VW_Optimization_Parameters.tsv', 'a+', newline='') as csvfile:
        tablewriter = csv.writer(csvfile, delimiter='\t')
        line = [results[0],results[1],results[2],results[3],results[4],results[5],alpha,window_i,args.dist_metric,args.classification_variable,args.dbscan_algorithm]
        tablewriter.writerow(line)
    
    return(results[-1] * -1)


cluster_class_info = []
cluster_tax_distance = []
clustering_info = defaultdict(dict)
rev_alias = {}

if (args.tax_distance_file):
    cluster_tax_distance =  pd.read_csv(args.tax_distance_file,sep='\t',na_values='',index_col=0,header=0)
    cluster_tax_distance = cluster_tax_distance.fillna(0)
    print('Shape of Cluster Tax Distance',cluster_tax_distance.shape)

#cluster_class_info = defaultdict(dict)
cluster_class_info = get_class_info()
for i in cluster_class_info.index.values.tolist():
    alias = get_alias(i)

#Order: window_i,alpha_i,correct,eps_i
param_space = dict(window_i=range(3,51),eps_i=uniform(0,1),alpha_i=uniform(-1,0))
conf_dict = dict(num_iteration=1000, initial_random=5, domain_size=1000)
tuner = Tuner(param_space, make_super_clusters, conf_dict) 
results = tuner.minimize()
print(f'Optimal value of parameters: {results["best_params"]} and objective: {results["best_objective"]}')



