#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=07-00:00:00
#SBATCH --job-name=Download_Tara_Viromes
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --output=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.out
#SBATCH --error=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=felipehcoutinho@gmail.com

##############################################
module load perl
perl /mnt/lustre/bio/users/fcoutinho/Scripts/download_SRA.pl --list /mnt/lustre/scratch/fcoutinho/StG/List_TARA_Oceans_All_Viromes.txt