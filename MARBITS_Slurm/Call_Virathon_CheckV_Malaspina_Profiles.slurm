#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=00-48:00:00
#SBATCH --job-name=Virathon_CheckV_Malaspina_Profiles_Prok_Viruses
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=20G
#SBATCH --output=/mnt/lustre/scratch/fcoutinho/Job_Logs/%J.Log.out
#SBATCH --error=/mnt/lustre/scratch/fcoutinho/Job_Logs/%J.Log.err
#SBATCH --mail-type=All
#SBATCH --mail-user=fhernandes@icm.csic.es

##############################################
module load python/3.8.5 checkv/0.7.0

export CHECKVDB="/mnt/lustre/repos/bio/databases/public/checkv/checkv-db-v1.0/"

cd /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Prok_Viruses/

python3 /mnt/lustre/bio/users/fcoutinho/Repos/virathon/Virathon.py --genome_files /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Prok_Viruses/Profiles_Malaspina_Viruses_CheckV_Trimmed_Genomes.fasta --call_checkv True --threads 47