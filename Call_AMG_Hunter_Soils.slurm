#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=05-00:00:00
#SBATCH --job-name=AMG_Hunter_Soil
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=150G
#SBATCH --output=/mnt/lustre/scratch/fcoutinho/Job_Logs/%J.Log.out
#SBATCH --error=/mnt/lustre/scratch/fcoutinho/Job_Logs/%J.Log.err
#SBATCH --mail-type=All
#SBATCH --mail-user=fhernandes@icm.csic.es

##############################################
module load diamond/2.0.7
module load python/3.8.5
module load hmmer
module load perl/5.28
module load prodigal

cd /mnt/lustre/scratch/fcoutinho/Soil_Elena/
python3 /mnt/lustre/bio/users/fcoutinho/Repos/virathon/Virathon.py --threads 23 --genome_files /mnt/lustre/scratch/fcoutinho/Soil_Elena/Post_CheckV_Trimmed_True_Virus_Soil_Genomes.fasta --call_prodigal True

cd /mnt/lustre/scratch/fcoutinho/Soil_Elena/Annotation/
python3 /mnt/lustre/bio/users/fcoutinho/Scripts/AMG_Hunter.py --annotate True --threads 23 --cds /mnt/lustre/scratch/fcoutinho/Soil_Elena/Post_CheckV_Trimmed_True_Virus_Soil_Genomes.faa

#python3 /mnt/lustre/bio/users/fcoutinho/Scripts/AMG_Hunter.py --annotate True --threads 23 --cds /mnt/lustre/scratch/fcoutinho/Soil_Elena/True_Virus_Soil_PEGs.faa --uniref_hits_file /mnt/lustre/scratch/fcoutinho/Soil_Elena/Annotation/True_Virus_Soil_PEGsxuniref100.blastp

