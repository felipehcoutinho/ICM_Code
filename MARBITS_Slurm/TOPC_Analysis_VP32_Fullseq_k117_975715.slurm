#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=00-06:00:00
#SBATCH --job-name=VP32_topc
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=100G
#SBATCH --output=/mnt/smart/users/fcoutinho/Job_Logs/jobLog_%J.out
#SBATCH --error=/mnt/smart/users/fcoutinho/Job_Logs/jobLog_%J.err
#SBATCH --mail-type=All
#SBATCH --mail-user=fhernandes@icm.csic.es

##############################################
#Move to the relvant directory
# cd /mnt/smart/users/fcoutinho/TARA_Polar/VP32/

#Load the prodigla module
# module load prodigal/2.6.3

# #Run prodigal
# prodigal -q -p meta -a VP32_original_Scaffold_k117_975715.faa -d VP32_original_Scaffold_k117_975715.fna -f gff -i VP32_original_Scaffold_k117_975715.fasta -o VP32_original_Scaffold_k117_975715.gff

# module purge

###Query proteins against Pfam
#Set up database
cd /mnt/smart/users/fcoutinho/Databases/
mkdir PfamA
cd PfamA

wget https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz

gunzip Pfam-A.hmm.gz

module load hmmer

hmmpress Pfam-A.hmm

#run hmmsearch
cd /mnt/smart/users/fcoutinho/TARA_Polar/VP32/

hmmsearch -o VP32xPfamA --tblout VP32xPfamA.tsv -T 50 --noali --cpu 23 /mnt/smart/users/fcoutinho/Databases/PfamA/Pfam-A.hmm VP32_original_Scaffold_k117_975715.faa

###Create IMGvr db to find close relatives
module load mmseqs2/13-45111

cd /mnt/smart/users/fcoutinho/Databases/IMGVR/

mmseqs createdb /mnt/smart/users/fcoutinho/Databases/IMGVR/IMGVR_all_proteins-high_confidence.faa.gz IMGvrDB
mmseqs createindex IMGvrDB tmp

rm -fr tmp/

###Query Proteins Against IMGvr db to find close relatives
cd /mnt/smart/users/fcoutinho/TARA_Polar/VP32/

mmseqs easy-search VP32_original_Scaffold_k117_975715.faa /mnt/smart/users/fcoutinho/Databases/IMGVR/IMGvrDB VP32_original_Scaffold_k117_975715xIMGvr.tsv tmp --threads 23 --max-seqs 1000 --min-seq-id 0.3 --min-aln-len 30
rm -fr tmp/

#Calculate pairwise protein scores
module load python
python3 /mnt/smart/users/fcoutinho/Repos/virathon/Virathon.py --pairwise_protein_scores True --genome VP32_original_Scaffold_k117_975715.fasta --cds VP32_original_Scaffold_k117_975715.faa --pps_hits_table VP32_original_Scaffold_k117_975715xIMGvr.tsv

#Select close relatives with min 50% AAI and 50% shared genes
awk -F "\t" '{ if (($3 >= 50) && ($5 >= 50)) { print } }'  VP32_original_Scaffold_k117_975715xIMGvr.tsv.Pairwise_Protein_Scores.tsv | cut -f 2 > List_of_IMGvr_Close_Relatives_AAI_50_Shared_Genes_50_Scaffold_k117_975715.txt 

###Query proteins against uniprot using mmseqs2
module load mmseqs2/13-45111

mmseqs createdb /mnt/smart/users/fcoutinho/Databases/Uniprot_02-10-2024/uniprot_sprot.fasta uniprot_sprotDB
mmseqs createindex uniprot_sprotDB tmp
mmseqs easy-search VP32_original_Scaffold_k117_975715.faa uniprot_sprotDB VP32_original_Scaffold_k117_975715xUniprot.tsv tmp --threads 24 --max-seqs 10 --min-seq-id 0.3 --min-aln-len 30 --format-mode 0 --format-output query,theader,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits

module purge

#Separate target description from ID
sed -i 's/ /'$'\t''/'  VP32_original_Scaffold_k117_975715xUniprot.tsv