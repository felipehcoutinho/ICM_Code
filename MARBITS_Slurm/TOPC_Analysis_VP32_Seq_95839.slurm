#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=00-12:00:00
#SBATCH --job-name=iphop_topc_vir
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=50G
#SBATCH --output=/mnt/smart/users/fcoutinho/Job_Logs/jobLog_%J.out
#SBATCH --error=/mnt/smart/users/fcoutinho/Job_Logs/jobLog_%J.err
#SBATCH --mail-type=All
#SBATCH --mail-user=fhernandes@icm.csic.es

##############################################
#Load the prodigla module
module load prodigal/2.6.3

#Move to the relvant directory
cd /mnt/smart/users/fcoutinho/TARA_Polar/VP32/

#Run prodigal
prodigal -q -p meta -a TOPC_Viral_Seq_95839.faa -d TOPC_Viral_Seq_95839.fna -f gff -i TOPC_Viral_Seq_95839.fasta -o TOPC_Viral_Seq_95839.gff

module purge

#Query proteins against uniprot using mmseqs2
module load mmseqs2/13-45111

mmseqs createdb /mnt/smart/users/fcoutinho/Databases/Uniprot_02-10-2024/uniprot_sprot.fasta uniprot_sprotDB
mmseqs createindex uniprot_sprotDB tmp
mmseqs easy-search TOPC_Viral_Seq_95839.faa uniprot_sprotDB TOPC_Viral_Seq_95839xUniprot.tsv tmp --threads 24 --max-seqs 10 --min-seq-id 0.3 --min-aln-len 30 --format-mode 0 --format-output query,theader,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits

#Separate target description from ID
sed -i 's/ /'$'\t''/'  TOPC_Viral_Seq_95839xUniprot.tsv

#Query Proteins Against IMGvr dt to find close relatives
mmseqs easy-search TOPC_Viral_Seq_95839.faa /mnt/smart/users/fcoutinho/Databases/IMGVR/IMGVR_all_proteins-high_confidence.faa.gz TOPC_Viral_Seq_95839xIMGvr.tsv tmp --threads 24 --max-seqs 1000 --min-seq-id 0.3 --min-aln-len 30


#Calculate pairwise protein scores
module load python
python3 /mnt/smart/users/fcoutinho/Repos/virathon/Virathon.py --pairwise_protein_scores True --genome TOPC_Viral_Seq_95839.fasta --cds TOPC_Viral_Seq_95839.faa --pps_hits_table TOPC_Viral_Seq_95839xIMGvr.tsv

#elect close relatives with min 50% AAI and 50% shared genes

###Select MQ and HQ Bins
awk -F "\t" '{ if (($3 >= 50) && ($5 >= 50)) { print } }' TOPC_Viral_Seq_95839xIMGvr.tsv | cut -f 2 > OctoMicro_MQ_and_HQ_Bins_R2_List.txt 
