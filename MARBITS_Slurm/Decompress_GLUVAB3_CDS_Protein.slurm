#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=00-04:00:00
#SBATCH --job-name=Decompress_GLUVAB3_CDS
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=5G
#SBATCH --output=/mnt/lustre/scratch/fcoutinho/Malaspina_Virus/jobLog_%J.out
#SBATCH --error=/mnt/lustre/scratch/fcoutinho/Malaspina_Virus/jobLog_%J.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=felipehcoutinho@gmail.com

##############################################
tar -zxvf GLUVAB3_CDS_Protein.tgz