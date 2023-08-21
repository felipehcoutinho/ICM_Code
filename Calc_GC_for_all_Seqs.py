#This script iterates through multiple fasta files and calculates the GC content for the entire file, considering all sequences together
#It does NOT calculate the GC content for each sequence in the file
#Import modules
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import GC
import glob

#Get options for directory containing the input files, file extensions, and output file name
import argparse
parser = argparse.ArgumentParser(description='This script iterates through multiple fasta files and calculates the GC content for the entire file, considering all sequences together. It does NOT calculate the GC content for each sequence in the file.')
parser.add_argument('-i', '--input_dir', help='Directory containing the input files', required=True)
parser.add_argument('-e', '--file_ext', help='File extension of the input files', required=True)
parser.add_argument('-o', '--output_file', help='Name of the output file', required=True)
args = parser.parse_args()

#Glob the files in the directory according to the specified extension

files = glob.glob(args.input_dir + "/*" + args.file_ext)

def calc_all_seqs_gc(input_file):
    """This function receives as input of a fasta file, iterates through all sequences in the file, and returns the GC content for the entire file, considering all sequences together"""
    #Create an empty list to store the sequences
    ind_seqs = []
    #Iterate through the sequences in the file and add them to the list
    for seq_record in SeqIO.parse(input_file, "fasta"):
        ind_seqs.append(str(seq_record.seq))
    #Concatenate all sequences in the list into one sequence
    concat_seq = ''.join(ind_seqs)
    #Calculate the GC content for the entire file
    gc_content = round(GC(concat_seq),2)
    #print(f"GC for {input_file} is {gc_content}")
    return gc_content

#Apply the calc_all_seqs_gc function to all files in the directory using the apply function, and assign the output to a dictionary in which the file names (minus the full path) are the keys and the GC content is the value
gc_dict = dict(zip([file.split("/")[-1] for file in files], [calc_all_seqs_gc(file) for file in files]))

#Convert gc_dict into a pandas data frame in which the file names are the index and the GC content is the column
import pandas as pd
gc_df = pd.DataFrame.from_dict(gc_dict, orient='index')

#Write the data frame to a tsv file with the specified output file name
gc_df.to_csv(args.output_file, sep='\t', header=False)

#Print a message to the user to let them know the script is done
print("Done!")

#Exit the script
exit()



