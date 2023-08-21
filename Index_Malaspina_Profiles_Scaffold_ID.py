from collections import defaultdict
from Bio import SeqIO
import pandas as pd
import argparse
import re
import glob

parser = argparse.ArgumentParser()
#Index and filtering options
parser.add_argument("--assemblies_dir", help="Directory containing fasta files of assembled contigs", type=str)
parser.add_argument("--bins_dir", help="Directory containing fasta files of binned genomes", type=str)
parser.add_argument("--format", help="Assembly files format", default="fasta", type=str)
parser.add_argument("--table_file", help=".tsv file with intermediate scaffold info", default="Scaffolds_Info.tsv", type=str)
parser.add_argument("--scaffold_info_table", help=".tsv file to write scaffold info", default="Scaffolds_Info.tsv", type=str)
parser.add_argument("--bin_info_table", help=".tsv file to write Bin info", default="Bin_Info.tsv", type=str)
args = parser.parse_args()

def central():
    seq_info = defaultdict(dict)
    #Index assembled seq info
    (seq_info,bin_info) = index_assemblies()
    #Convert the 2d dictionary info_dict into a pandas dataframe and print it to output_dataframe_file in .tsv format
    print(f"Printing scaffold info to {args.scaffold_info_table}")
    info_dataframe = pd.DataFrame.from_dict(seq_info)
    info_dataframe.index.name = 'Scaffold'
    info_dataframe.to_csv(args.scaffold_info_table,sep="\t",na_rep='NA')
    
    print(f"Printing Bin info to {args.bin_info_table}")
    info_dataframe = pd.DataFrame.from_dict(bin_info)
    info_dataframe.index.name = 'Bin'
    info_dataframe.to_csv(args.bin_info_table,sep="\t",na_rep='NA')

def index_assemblies():
    print(f'Reading info from {args.table_file}')
    info_data_frame = pd.read_csv(args.table_file, sep="\t",index_col="Scaffold",header=0)
    info_data_frame['Description'] = info_data_frame['Description'].str.replace(r'^S(\d)+C', '', regex=True)
    info_data_frame = info_data_frame.reset_index().set_index('Description')
    #print(info_data_frame)
    seq_info = defaultdict(dict)
    bin_info = defaultdict(dict)
    print ('Indexing scaffolds')
    match_string = re.compile("MP(\d)+(bis)*")
    assemblies_list = glob.glob(f'{args.assemblies_dir}/*fasta')
    asb_count = len(assemblies_list)
    bins_list = glob.glob(f'{args.bins_dir}/*fa')
    bins_count = len(bins_list)
    print(f"Procesing {asb_count} assemblines and {bins_count} bins")
    for assembly_file in assemblies_list:
        print(f'Processing assembly {assembly_file}')
        seq_counter = 0
        #Iterate over genomic sequences in the file. Collect basic Info
        for seqobj in SeqIO.parse(assembly_file, "fasta"):
            seq_counter += 1
            match_obj = re.search(match_string,seqobj.id)
            seq_info['Sample'][seqobj.id] = match_obj.group()
            seq_info['Description'][seqobj.id] = seqobj.description
            seq_info['Length'][seqobj.id] = len(seqobj.seq)
            seq_info['Original_File'][seqobj.id] = assembly_file
            if (seqobj.id in info_data_frame.index):
                seq_info['ID_2'][seqobj.id] = info_data_frame.loc[seqobj.id]['Original_ID']
                seq_info['ID_3'][seqobj.id] = info_data_frame.loc[seqobj.id]['Scaffold']
                
    match_string = re.compile("^S(\d)+C")
    for bin_file in bins_list:
        print(f'Processing bin {bin_file}')
        bin_id = get_prefix(bin_file,"fa")
        bin_info['File'][bin_id] = bin_file
        #Iterate over genomic sequences in the file. Collect basic Info
        for seqobj in SeqIO.parse(bin_file, "fasta"):
            oid = seqobj.id
            nid = re.sub(match_string,"",oid)
            if (nid not in seq_info['Bin_file']):
                seq_info['Bin_file'][nid] = []
                seq_info['Bin'][nid] = []
            seq_info['Bin_file'][nid].append(bin_file)
            seq_info['Bin'][nid].append(bin_id)
            #bin_info['Sample'][bin_id] = seq_info['Sample'][nid]

    return (seq_info,bin_info)

def get_prefix(file,format):
    prefix_file = re.sub(f'.{format}','',file)
    prefix_file = re.sub('(.)+/','',prefix_file)
    return prefix_file
    
central()