#! /usr/bin/env python3
import pandas as pd
from collections import defaultdict
from Bio import SeqIO
from Bio import SearchIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import subprocess
import argparse
import re
import glob


parser = argparse.ArgumentParser()
parser.add_argument("--mags_extension", help="Extension of the Mags", default="fasta", type=str)
parser.add_argument("--mags_dir", help="Directory containing MAGs", type=str)
parser.add_argument("--viral_genomes", help="Multifasta file containing all viral genomes to match", type=str)
parser.add_argument("--threads", help="The number of threads to be used", default=1, type=int)
parser.add_argument("--parse_only", help="Flag to only parse results", default=False, type=bool)
parser.add_argument("--sample_specific", help="Flag to only allow matches between viruses and MAGs from the same sample", default=False, type=bool)
parser.add_argument("--scaffold_specific", help="Flag to only allow matches between viruses and MAGs from the same scaffold", default=False, type=bool)
parser.add_argument("--virus_scaffold_info_file", help=".tsv File linking viral scaffolds to unique sample IDs", type=str)
parser.add_argument("--host_scaffold_info_file", help=".tsv File linking host scaffolds to unique sample IDs", type=str)
parser.add_argument("--remove_virus", help="Flag to produce MAGs with viral sequences removed", default=False, type=bool)
args = parser.parse_args()

mag_files = glob.glob(f"{args.mags_dir}/*{args.mags_extension}")

mag_info = defaultdict(lambda: defaultdict(int))
seq_info = defaultdict(dict)
coord_info = defaultdict(dict)
seen_pairs = defaultdict(dict)
hsinfo = "NA"
vsinfo = "NA"

if ((args.host_scaffold_info_file) and (args.virus_scaffold_info_file)):
    print(f"Reading in info from: {args.host_scaffold_info_file}")
    hsinfo = pd.read_csv(args.host_scaffold_info_file, sep="\t",index_col="Sequence",header=0)
    print(f"Reading in info from: {args.virus_scaffold_info_file}")
    vsinfo = pd.read_csv(args.virus_scaffold_info_file, sep="\t",index_col="Sequence",header=0)


def check_match_cutoff(hsp,max_evalue,min_bitscore,min_ident,min_ali):
    if ((hsp.bitscore < min_bitscore) or (hsp.evalue > max_evalue) or (hsp.ident_pct < min_ident) or (hsp.aln_span < min_ali)):
        return False
    else:
        return True

def get_prefix(file,format):
    prefix_file = re.sub(f'.{format}','',file)
    prefix_file = re.sub('(.)+/','',prefix_file)
    return prefix_file

with open('All_MAG_Seqs.fasta', 'w', newline='') as OUT:
    for mag in mag_files:
        #print('Processing',mag)
        for seqobj in SeqIO.parse(mag, 'fasta'):
            mag_info['Sequence_Count'][mag] += 1
            mag_info['Bp_Count'][mag] += len(seqobj.seq)
            mag_info['Virus_List'][mag] = set()
            mag_info['Virus_Regions'][mag] = []
            seq_info["Source"][seqobj.id] = "Host"
            seq_info["MAG"][seqobj.id] = mag
            seq_info["Length"][seqobj.id] = len(seqobj.seq)
            if (args.parse_only == False):
                SeqIO.write(seqobj, OUT, "fasta")
            
for seqobj in SeqIO.parse(args.viral_genomes, 'fasta'):
    seq_info["Source"][seqobj.id] = "Virus"
    seq_info["Length"][seqobj.id] = len(seqobj.seq)


#Build BLAST DB of MAG genomes
command = f'makeblastdb -in All_MAG_Seqs.fasta -dbtype nucl -title DB_All_MAG_Seqs -out DB_All_MAG_Seqs'
if (args.parse_only == False):
    subprocess.call(command, shell=True)

command = f"blastn -db DB_All_MAG_Seqs -query {args.viral_genomes} -out Viral_GenomesxAll_MAG_Seqs.blastn -outfmt 6 -evalue 0.001 -perc_identity 90 -max_target_seqs 100 -num_threads {args.threads}"
if (args.parse_only == False):
    subprocess.call(command, shell=True)

print ('Parsing BLASTN output')
sample_re_obj = re.compile("MP(\d)+")

valid_count = 0
#Iterare over queries
for qresult in SearchIO.parse("Viral_GenomesxAll_MAG_Seqs.blastn", 'blast-tab'):
    #Iterate over hits in the query
    for hit in qresult.hits:
        #Iterate over HSPs in the hits
        for hsp in hit.hsps:
            #Check if the HSP passes established blast cutoffs
            viral_genome = qresult.id
            mag_seq = hit.id
            is_valid = check_match_cutoff(hsp,0.00001,500,100,seq_info["Length"][viral_genome])
            #prophage_frac = seq_info["Length"][viral_genome] / seq_info["Length"][mag_seq]
            #if ((is_valid) and  (prophage_frac <= 0.5)):
            #print(seen_pairs[mag].keys())
            if args.sample_specific:
                #match_obj = re.search(sample_re_obj,mag_seq)
                #mag_sample = match_obj.group()
                #match_obj = re.search(sample_re_obj,viral_genome)
                #vir_sample = match_obj.group()
                mag_sample = hsinfo["Sample"][mag_seq]
                vir_sample = vsinfo["Sample"][viral_genome]
                print("Checking match for:",mag_seq,mag_sample,viral_genome,vir_sample)
                if (mag_sample != vir_sample):
                    is_valid = False
            if args.scaffold_specific:
                clean_mag_scaff_name =  re.sub("S(\d)+C","",mag_seq)
                clean_mag_scaff_name =  re.sub("_Scaffold","",clean_mag_scaff_name)
                match_obj = re.search("MP(\d)+_(\d)+",viral_genome)
                clean_vir_scaff_name =  match_obj.group()
                #print(mag_seq,clean_mag_scaff_name,viral_genome,clean_vir_scaff_name)
                if (clean_mag_scaff_name != clean_vir_scaff_name):
                    is_valid = False
                
            if ((is_valid == True) and (viral_genome not in seen_pairs[mag].keys())):
                mag = seq_info["MAG"][hit.id]
                mag_info['Prophage_Count'][mag] += 1
                mag_info['Prophage_Base_Pairs'][mag] += seq_info["Length"][viral_genome]
                mag_info['Virus_List'][mag].add(viral_genome)
                mag_info['Virus_Regions'][mag].append(f"{mag_seq}|{hsp.hit_start}|{hsp.hit_end}")
                seen_pairs[mag][viral_genome] = True
                print("Valid",viral_genome,seq_info["Length"][viral_genome],mag,mag_seq,hsp.hit_start,hsp.hit_end)
                valid_count += 1
                coord_info["Host_Genome"][valid_count] = mag
                coord_info["Host_Sequence"][valid_count] = mag_seq
                coord_info["Viral_Genome"][valid_count] = viral_genome
                coord_info["Start_Host_Sequence"][valid_count] = int(hsp.hit_start) + 1
                coord_info["End_Host_Sequence"][valid_count] = hsp.hit_end
                full_sequence = False
                if (hsp.aln_span >= seq_info["Length"][mag_seq]):
                    full_sequence = True
                coord_info["Full_Viral_Sequence"][valid_count] = full_sequence
                

print("Printing results to MAG_Info.tsv")
mag_info_df = pd.DataFrame.from_dict(mag_info)
mag_info_df.index.name = 'MAG'
mag_info_df.to_csv("MAG_Info.tsv",sep="\t",na_rep=0)

print("Printing results to Coord_Info.tsv")
coord_info_df = pd.DataFrame.from_dict(coord_info)
coord_info_df.index.name = 'Genome'
coord_info_df.to_csv("Coord_Info.tsv",sep="\t",na_rep=0)

if (args.remove_virus):
    print("Removing Viral Sequences from MAGs")
    for mag in mag_files:
        print('Processing',mag)
        mag_nid = "No_Virus_" + get_prefix(mag,"DUMMY")
        with open(mag_nid, 'w', newline='') as OUT:
            for seqobj in SeqIO.parse(mag, 'fasta'):
                original_seq = seqobj.seq
                full_length = len(seqobj.seq)
                index_length = full_length - 1
                sub_coord = coord_info_df[coord_info_df["Host_Sequence"] == seqobj.id]
                if (len(sub_coord) >= 1):
                    for i,coord in sub_coord.iterrows():
                        print(coord)
                        start_index = int(coord["Start_Host_Sequence"]) - 1
                        stop_index = int(coord["End_Host_Sequence"]) - 1
                        full_range = list(range(0,index_length))
                        trim_range = list(range(start_index,stop_index))
                        blacklist = {posit : True for posit in trim_range}
                        trimmed_seq_nucs = []
                        for nuc_posit in full_range:
                            if nuc_posit in blacklist:
                                trimmed_seq_nucs.append("X")
                            else:
                                trimmed_seq_nucs.append(original_seq[nuc_posit])
                        original_seq = "".join(trimmed_seq_nucs)
                    split_seqs = re.split("(X)+",original_seq)
                    frag_count = 0 
                    for split in split_seqs:
                        frag_count += 1
                        if (len(split) > 1):
                            split_seq = SeqRecord(seq=Seq(split), id=seqobj.id+f"_Split_{frag_count}", name = seqobj.name, description=seqobj.description)
                            #seqobj.seq = Seq(original_seq)
                            SeqIO.write(split_seq, OUT, "fasta")
                else:
                    SeqIO.write(seqobj, OUT, "fasta")