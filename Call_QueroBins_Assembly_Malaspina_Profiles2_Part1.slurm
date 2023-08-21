#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=20-00:00:00
#SBATCH --job-name=QueroBins_Malaspina2_Part1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=750G
#SBATCH --output=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.out
#SBATCH --error=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=felipehcoutinho@gmail.com

##############################################
module load python/3.8.5
module load spades

python3 /mnt/lustre/bio/users/fcoutinho/Scripts/QueroBins.py --coassembly_table /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Redo_Profiles_Malaspina_Assembly_Dictionary_Eps_0.0125_Groups_Non_Zero_Part_1.tsv --threads 48 --assemble True --max_memory 750