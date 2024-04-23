#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=05-00:00:00
#SBATCH --job-name=Antarctic_lagoon_ASV_Phylogeny
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=50G
#SBATCH --output=/mnt/smart/users/fcoutinho/Job_Logs/jobLog_%J.out
#SBATCH --error=/mnt/smart/users/fcoutinho/Job_Logs/jobLog_%J.err
#SBATCH --mail-type=All
#SBATCH --mail-user=fhernandes@icm.csic.es

##############################################
cd /mnt/smart/users/fcoutinho/Antarctic_Lagoons/DADA_Output/

# module load muscle

# muscle -super5 All_ASV_Sequences.fasta -output Aligned_All_ASV_Sequences.fasta -threads 48

module load iqtree/2.0.6

#Using model finder
iqtree2 -s Aligned_All_ASV_Sequences.fasta --mem 45G --threads-max 47 -m MFP

mv Aligned_All_ASV_Sequences.fasta.treefile Aligned_All_ASV_MF.nwk
