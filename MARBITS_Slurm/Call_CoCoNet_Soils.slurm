#!/bin/bash
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=3-00:00:00
#SBATCH --job-name=CoCoNet_Soil
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=100G
#SBATCH --output=/mnt/lustre/scratch/fcoutinho/Job_Logs/%J.Log.out
#SBATCH --error=/mnt/lustre/scratch/fcoutinho/Job_Logs/%J.Log.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=fhernandes@icm.csic.es

##############################################
module load coconet/1.1.0-064cb8c
source activate
export PYTHONUSERBASE=intentionally-disabled
coconet run --fasta /mnt/lustre/scratch/fcoutinho/Soil_Elena/Post_CheckV_Trimmed_True_Virus_Soil_Genomes.fasta --output /mnt/lustre/scratch/fcoutinho/Soil_Elena/CoCoNet_Soil/ --threads 23 --bam /mnt/lustre/scratch/fcoutinho/Soil_Elena/Bowtie_Abundance/*bam