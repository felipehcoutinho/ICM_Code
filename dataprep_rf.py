from iphop.modules import utility
from iphop.modules import dataprep
import sys
import os
import logging
logger = logging.getLogger(__name__)
import math
import subprocess as sp
import pandas as pd
pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
import numpy as np
import csv
from itertools import combinations
import click
from Bio import Phylo
csv.field_size_limit(sys.maxsize) ## Need this for GTDB file which has some very very long fields

step = 7.5

def aggregate_rf(args):
	check = check_already_done(args)
	if check:
		logger.info(f"[{step}/Skip] We already found all the expected files, we skip...")
	else:
		logger.info(f"[{step}] Aggregating all results and formatting for RF...")
		## Load all info from parsed file (should fit easily in memory)
		## Load all observations from label
		df_labels = load_current_obs(args)
		## Load info blast from blast_parsed
		df_blast = load_signal_blast(args)
		## Load info crispr from crispr_parsed
		df_crispr = load_signal_crispr(args)
		## Fill in the matrix
		compute_matrices(df_blast,df_crispr,df_labels,args)

def check_already_done(args):
	args["matrix_labels"] = os.path.join(args["tmp"], "matrixlabels.csv")
	tag = 0
	if not os.path.exists(args["matrix_labels"]): # If we don't even have a label file, we stop there
		tag = 1
	for tool in args["list_tools_rf"]:
		if args["list_tools_rf"][tool] == 1: # we check if we wanted this tool (that should always be 1, but just to check)
			# logger.debug(f"We want to use {tool}")
			code = "matrix_"+tool+"_rf"
			args[code] = os.path.join(args["tmp"], "matrix"+tool+"_rf.csv")
			if os.path.exists(args[code]): ## we have the righ matrix
				args["list_tools_rf"][tool] = 2 # If so, we check it as done
				# logger.debug(f"We want to use {tool} --> {args[code]}")
			else:
				# logger.debug(f"We don't have the matrix for {tool} ==> {args[code]}")
				tag = 1 # If we did want this, that means we don't have all the matrix we wanted, so we need to load this
	if tag == 0:
		return True # If we reach here, that means we did not find an issue, so we think everything is here and correct
	else:
		## so we don't have all the matrices, we need to reload everything
		for tool in args["list_tools_rf"]:
			if args["list_tools_rf"][tool] == 2:
				args["list_tools_rf"][tool] = 1
		return False


def load_signal_blast(args):
	if not os.path.exists(args["blastparsed"]):
			logger.critical(f"Pblm - trying to use blast, but the file {args['blastparsed']} does not exist ? ...")
			sys.exit("")
	## New idea, we will store a dict of (consistent) pandas dataframe
	logger.debug(f"Loading blast data -- {args['blastparsed']}")
	df_blast = pd.read_csv(args["blastparsed"],delimiter=',',quotechar='"',dtype={'Virus': "string",'Host': "string",'Host contig': "string"})
	# filter and transform
	df_blast.drop(df_blast[df_blast['N match'] < args["thresholds"]["min_blast_hit"]].index, inplace=True) ## dropping all rows with blast hit < 500bp
	df_blast.drop(df_blast[df_blast['Id %'] < args["thresholds"]["min_blast_pcent"]].index, inplace=True) ## dropping all rows with percentaged id under the threshold
	# after filtering, for the RF, we actually re-aggregated the separate hits, then transform (the transform is only to get the same sorting as the original dataprep)
	### We actually should not group by here, these will be summed later
	# df_blast = df_blast.groupby(['Virus','Host','V Cover','H Cover','Match length','Id %','Hit number (for duplications)'], as_index=False).agg({'N match':'sum'})
	df_blast['N_match_transformed'] = df_blast['N match'].apply(lambda x: round(x/500) - 1)
	# tmp = df_blast.loc[df_blast["Virus"] == "MT889366"]
	return df_blast

def load_signal_crispr(args):
	if not os.path.exists(args["crisprparsed"]):
			print(f"Pblm - trying to use crispr, but the file {args['crisprparsed']} does not exist ? ...")
			sys.exit("")
	df_crispr = pd.read_csv(args["crisprparsed"],delimiter=',',quotechar='"',dtype={'Virus': "string",'Host': "string",'Spacer': "string"})
	df_crispr.drop(df_crispr[df_crispr['Array size'] < args["thresholds"]["min_n_spacer"]].index, inplace=True) ## dropping all rows with arrays including less than 2 spacers
	df_crispr.drop(df_crispr[df_crispr['Spacer length'] < args["thresholds"]["min_len_crispr"]].index, inplace=True) ## dropping all crispr spacers shorter than 10bp
	df_crispr.drop(df_crispr[df_crispr['Pcent Id'] < args["thresholds"]["min_crispr_pcent"]].index, inplace=True) ## dropping all rows with hits < 75%
	df_crispr.drop(df_crispr[df_crispr['N mismatch'] > args["thresholds"]["max_mis_crispr"]].index, inplace=True) ## dropping all rows with > 4 mismatches (too noisy)
	return df_crispr


def load_current_obs(args):
	if not os.path.exists(args["matrix_labels"]):
		logger.critical(f"Pblm - trying to load label info, but the file {args['matrix_labels']} does not exist ? ...")
		sys.exit("")
	df_labels = pd.read_csv(args["matrix_labels"],delimiter=',',quotechar='"',dtype={'Observation_n': "int",'Virus': "string",'Repr host': "string",'Repr host taxonomy':"string"})
	df_labels = df_labels.rename(columns={'Repr host':'Host'})
	return df_labels

def compute_matrices(df_blast,df_crispr,df_labels,args):
	logger.debug("Loading trees")
	tab_trees = dataprep.load_trees(args)
	logger.debug("Load taxo and repr")
	check_host = df_blast.groupby(['Host']).Host.unique()
	host_info = {}
	dataprep.load_taxo_repr(args,check_host,host_info)
	# logger.debug(host_info)
	### Now blast
	df_blast["Repr"] = df_blast["Host"].apply(lambda x: dataprep.get_repr(x,host_info))
	df_blast = df_blast[~df_blast.Repr.isin(["Unknown"])]
	## We take the top 50 hits for each virus, otherwise we have hits in df that are not in the labels file
	tmp = df_blast.groupby(['Virus','Repr','Host']).sum().reset_index() ## For each Repr, we take the sum of all hits
	tmp = tmp.groupby(['Virus','Repr','N_match_transformed']).first().reset_index()
	tmp = tmp.groupby(['Virus']).apply(lambda x: x.nlargest(n=50,columns='N_match_transformed',keep='all')).reset_index(drop=True)
	# logger.debug(tmp)
	df_blast = df_blast[df_blast[['Virus','Host']].apply(tuple,axis=1).isin(tmp[['Virus','Host']].apply(tuple,axis=1))]
	## Prep a dummy row to append if not enough rows
	dummy = {'Virus': 'DummyVirus', 'Host': 'DummyHost', 'V Cover': None, 'H Cover': None, 'Match length': None, 'Id %': None, 'N match': None}
	## Prep out file
	code = "matrix_blast_rf"
	rfb = args[code]
	with open(rfb, 'w', newline='') as frf:
		outwriter = csv.writer(frf, delimiter=',',quotechar='"', quoting=csv.QUOTE_MINIMAL, lineterminator="\n")
		outwriter.writerow(["Observation_n","D1","D2","D3","D4","D5","L1","L2","L3","L4","L5","Id1","Id2","Id3","Id4","Id5"])
		for virus, all_labels in df_labels.groupby(['Virus']):
			store_dist = {}
			logger.debug(f"Processing data for virus {virus}")
			# logger.debug(f"This is linked to {all_labels}")
			# logger.debug(f"Loading all distances")
			list_hosts = all_labels.groupby(['Host']).Host.unique()
			# with click.progressbar(combinations(list_hosts,2),label = print("Listing all relevant distances")) as bar:
			# logger.debug(list_hosts)
			dataprep.load_pair_dist(list_hosts,store_dist,tab_trees,args)
			# logger.debug(store_dist)
			## Now compute the 5 top blast for this virus
			selected_blast = df_blast.loc[df_blast["Virus"] == virus]
			# full_blast = df_blast.loc[df_blast["Virus"] == virus]
			if selected_blast.empty:
				logger.debug(f"No features at all, we can't do anything")
			else:
				# logger.debug(f"And the corresponding Blast info is {selected_blast}")
				selected_blast = selected_blast.nlargest(5,['N match'], keep='all').reset_index(drop=True) ## Select 5 best based on match total # (i.e. combination of % id and Match length)
				while len(selected_blast) < 5:
					# logger.debug("We need to add one more empty rows at the end")
					selected_blast = selected_blast.append(dummy,ignore_index=True)
					# logger.debug(f"selected blast is now {selected_blast}")
					# logger.debug(len(selected_blast))
				# logger.debug(f"Which is then filtered to {selected_blast}")
				for obs, obs_info in all_labels.groupby(['Observation_n']):
					host_pivot = obs_info.iloc[0]["Host"]
					# logger.debug(f"Looking at observation {obs} -- {host_pivot}")
					## Calculate host distance
					selected_blast['Dist'] = selected_blast['Repr'].apply(lambda x: update_dist(host_pivot,x,store_dist))
					selected_blast = selected_blast.sort_values(by = ["Dist","N match","Id %"], ascending = [False,False,False])
					vec_dist = selected_blast['Dist'].transpose()
					vec_len = selected_blast['Match length'].transpose()
					vec_id = selected_blast['Id %'].transpose()
					outwriter.writerow([obs,vec_dist[0],vec_dist[1],vec_dist[2],vec_dist[3],vec_dist[4],vec_len[0],vec_len[1],vec_len[2],vec_len[3],vec_len[4],vec_id[0],vec_id[1],vec_id[2],vec_id[3],vec_id[4]])
	## Now CRISPR
	logger.debug("Load taxo and repr")
	check_host = df_crispr.groupby(['Host']).Host.unique()
	host_info = {}
	dataprep.load_taxo_repr(args,check_host,host_info)
	df_crispr["Repr"] = df_crispr["Host"].apply(lambda x: dataprep.get_repr(x,host_info))
	df_crispr = df_crispr[~df_crispr.Repr.isin(["Unknown"])]
	## We take the top 30 hits for each virus, otherwise we have hits in df that are not in the labels file
	tmp = df_crispr.groupby(['Virus','Repr','Host']).min().reset_index() ## For each Repr, we take the sum of all hits
	tmp['N mismatch'] = tmp['N mismatch'].apply(lambda x: 4 - x )
	tmp = tmp.groupby(['Virus','Repr','Spacer length','N mismatch']).first().reset_index()
	tmp = tmp.groupby(['Virus']).apply(lambda x: x.nlargest(30,['N mismatch'], keep='all')).reset_index(drop=True)
	df_crispr = df_crispr[df_crispr[['Virus','Host']].apply(tuple,axis=1).isin(tmp[['Virus','Host']].apply(tuple,axis=1))]
	## Prep a dummy row to append if not enough rows
	dummy = {'Virus': 'DummyVirus', 'Host': 'DummyHost', 'Spacer': 'DummySpacer', 'Array size': 2, 'Spacer complexity score': None, 'Spacer length': None, 'Pcent Id': None, 'N mismatch': None, 'Repr': None, 'Dist': None}
	## Prep out file
	code = "matrix_crispr_rf"
	rfc = args[code]
	with open(rfc, 'w', newline='') as frf:
		outwriter = csv.writer(frf, delimiter=',',quotechar='"', quoting=csv.QUOTE_MINIMAL, lineterminator="\n")
		outwriter.writerow(["Observation_n","D1","D2","D3","D4","D5","C1","C2","C3","C4","C5","L1","L2","L3","L4","L5","Nm1","Nm2","Nm3","Nm4","Nm5"])
		for virus, all_labels in df_labels.groupby(['Virus']):
			store_dist = {}
			logger.debug(f"Processing data for virus {virus}")
			# logger.debug(f"This is linked to {all_labels}")
			## Now compute the 5 top CRISPR for this virus
			selected_crispr = df_crispr.loc[df_crispr["Virus"] == virus]
			full_crispr = df_crispr.loc[df_crispr["Virus"] == virus]
			if selected_crispr.empty:
				logger.debug(f"No CRISPR at all, we can't do anything")
			else:
				# logger.debug(f"And the corresponding CRISPR info is {selected_crispr}")
				# print(f"Loading all distances")
				list_hosts = all_labels.groupby(['Host']).Host.unique()
				# with click.progressbar(combinations(list_hosts,2),label = print("Listing all relevant distances")) as bar:
				dataprep.load_pair_dist(list_hosts,store_dist,tab_trees,args)
				##
				selected_crispr = selected_crispr.nsmallest(5,['N mismatch'], keep='all').reset_index(drop=True) ## Select 5 best based on mismatch
				while len(selected_crispr) < 5:
					selected_crispr = selected_crispr.append(dummy,ignore_index=True)
				# print(f"Which is then filtered to {selected_crispr}")
				for obs, obs_info in all_labels.groupby(['Observation_n']):
					# print(f"Looking at observation {obs}")
					host_pivot = obs_info.iloc[0]["Host"]
					# logger.debug(f"Looking at observation {obs} with host pivot {host_pivot}")
					selected_crispr['Dist'] = selected_crispr['Repr'].apply(lambda x: update_dist(host_pivot,x,store_dist))
					# print(f"We have now updated selected crispr as:  {selected_crispr}")
					selected_crispr.sort_values(by = ["N mismatch","Dist","Spacer complexity score","Spacer length"], ascending = [True,False,True,True], inplace = True, ignore_index = True)
					# logger.debug(f"We have now sorted selected crispr as:  {selected_crispr}")
					vec_dist = selected_crispr['Dist'].transpose()
					vec_comp = selected_crispr['Spacer complexity score'].transpose()
					vec_len = selected_crispr['Spacer length'].transpose()
					vec_mismatch = selected_crispr['N mismatch'].transpose()
					outwriter.writerow([obs,vec_dist[0],vec_dist[1],vec_dist[2],vec_dist[3],vec_dist[4],vec_comp[0],vec_comp[1],vec_comp[2],vec_comp[3],vec_comp[4],vec_len[0],vec_len[1],vec_len[2],vec_len[3],vec_len[4],vec_mismatch[0],vec_mismatch[1],vec_mismatch[2],vec_mismatch[3],vec_mismatch[4]])

def update_dist(host_1,host_2,store_dist): ## special for RF, need to deal with NaN
	dist = 0
	code = str(host_1)+"_vs_"+str(host_2)
	if host_2 is None:
		dist = None
	elif host_2 != host_2: ## Shoul detect the NaN
		dist = None
	elif host_1 == host_2:
		dist = 0
	else:
		dist = store_dist[code]
	return dist
