#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=03-00:00:00
#SBATCH --job-name=QueroBins_Malaspina
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=510G
#SBATCH --output=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.out
#SBATCH --error=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=felipehcoutinho@gmail.com

##############################################
module load python/3.8.5
module load spades

python3 /mnt/lustre/bio/users/fcoutinho/Scripts/QueroBins.py --coassembly_table /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies/Malaspina_Profiles_Coassembly_Dictionary.tsv --threads 48