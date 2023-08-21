#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=1-00:00:00
#SBATCH --job-name=PHIST_Profiles_Malaspina
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=50G
#SBATCH --output=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.out
#SBATCH --error=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=felipehcoutinho@gmail.com

##############################################
module load python

python3 /mnt/lustre/bio/users/fcoutinho/PHIST/phist.py -t 48 /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Viruses/Split_Profiles_Malaspina_Viruses_CheckV_Trimmed_Genomes/ /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Viruses/MAGs_Without_Viruses/ ./PHIST_Output/PHIST_Kmers.csv ./PHIST_Output/PHIST_Predictions.csv
