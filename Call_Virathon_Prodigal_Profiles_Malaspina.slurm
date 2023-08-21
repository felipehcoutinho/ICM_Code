#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=02-00:00:00
#SBATCH --job-name=Virathon_Prodigal_Profiles_Malaspina
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=150G
#SBATCH --output=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.out
#SBATCH --error=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=felipehcoutinho@gmail.com

##############################################
module load python/3.8.5
module load prodigal/2.6.3 

#python3 /mnt/lustre/bio/users/fcoutinho/Scripts/virathon_v0.1.py --genome_files /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Viruses/Profiles_Malaspina_Viruses_CheckV_Trimmed_Genomes.fasta --call_prodigal True

module load diamond/2.0.7
module load hmmer
module load perl

python3 /mnt/lustre/bio/users/fcoutinho/Scripts/AMG_Hunter.py --annotate True --threads 24 --cds /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Viruses/Profiles_Malaspina_Viruses_CheckV_Trimmed_Genomes.faa