#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=vir
#SBATCH --time=03-00:00:00
#SBATCH --job-name=BlanesVir_Assembly
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=250G
#SBATCH --output=/mnt/smart/scratch/vir/felipe/Job_Logs/jobLog_%J.out
#SBATCH --error=/mnt/smart/scratch/vir/felipe/Job_Logs/jobLog_%J.err
#SBATCH --mail-type=All
#SBATCH --mail-user=fhernandes@icm.csic.es

##############################################
cd /mnt/smart/scratch/vir/felipe/Workshops/Assemblies/

module load spades/4.2.0

spades.py -1 /mnt/smart/scratch/vir/felipe/Workshops/Raw_Reads/Blanes_Viromes/fastq/ABR_S3_R1_001.fastq.gz -2 /mnt/smart/scratch/vir/felipe/Workshops/Raw_Reads/Blanes_Viromes/fastq/ABR_S3_R2_001.fastq.gz -o Assembly_ABR_S3 --threads 24 --memory 250 --meta