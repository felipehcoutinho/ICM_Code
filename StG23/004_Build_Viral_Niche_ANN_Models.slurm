#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=00-03:00:00
#SBATCH --job-name=VirNicheStG23
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=20G
#SBATCH --output=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.out
#SBATCH --error=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=fhernandes@icm.csic.es

##############################################
module load R

Rscript /mnt/lustre/bio/users/fcoutinho/StG_23/Scripts/004_StG_23_ANN_Niche_Modelling.R /mnt/lustre/scratch/fcoutinho/StG_23/Viruses/Abundance/RPKM_Abundance_NA.tsv /mnt/lustre/scratch/fcoutinho/MaLME/Marine_Viral_Communities_Sample_Metadata.tsv 100 1 5 Temperature,Salinity,Oxygen,Chlorophyll_A,Iron_5m,Ammonium_5m 50 Viruses

Rscript /mnt/lustre/bio/users/fcoutinho/StG_23/Scripts/004_Build_gVOD_ANN_Models.R Viral_Abundance Depth,Prokaryote_Abundance,Temperature,Salinity 1000 0.7 10 FALSE
#Rscript /mnt/lustre/bio/users/fcoutinho/StG_23/Scripts/004_Build_gVOD_ANN_Models.R Prokaryote_Abundance Depth,Temperature,Salinity 1000 0.7 10 FALSE
