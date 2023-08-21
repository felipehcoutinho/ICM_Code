#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=00-02:00:00
#SBATCH --job-name=Build_MMSeqs_DB_UniRef100
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=400G
#SBATCH --output=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.out
#SBATCH --error=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=felipehcoutinho@gmail.com

##############################################
module load mmseqs2/15c7762
mmseqs createdb /mnt/lustre/repos/bio/databases/public/UniProtKB/UniProt_Release_2021-04-07/uniref100.fasta /mnt/lustre/scratch/fcoutinho/UniProt100/DB_UniRef100
mmseqs createindex /mnt/lustre/scratch/fcoutinho/UniProt100/DB_UniRef100 tmp --split 0
