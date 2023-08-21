#!/bin/bash
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=03-00:00:00
#SBATCH --job-name=VWV_DB
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=50G
#SBATCH --output=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.out
#SBATCH --error=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=felipehcoutinho@gmail.com

##############################################
module load mmseqs2

####efam
#tar -zxvf 10_aligned_fasta.tar.gz

####SEED
mmseqs easy-linclust /mnt/lustre/bio/users/fcoutinho/Databases/SEED/subsystems_proteins.all.unique.faa.gz SEED_Clustered tmp --threads 48 -c 0.7 --min-seq-id 0.35 --cov-mode 0
