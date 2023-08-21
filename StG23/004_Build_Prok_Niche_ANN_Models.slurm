#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=00-06:00:00
#SBATCH --job-name=ProkNicheStG23
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=50G
#SBATCH --output=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.out
#SBATCH --error=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=fhernandes@icm.csic.es

##############################################
module load R

cd /mnt/lustre/scratch/fcoutinho/StG_23/Prokaryotes/Niche_Models/

Rscript /mnt/lustre/bio/users/fcoutinho/StG_23/Scripts/004_StG_23_ANN_Niche_Modelling.R /mnt/lustre/scratch/fcoutinho/StG_23/Prokaryotes/Abundance/RPKM_Abundance_OceanDNA_All_Species_Rep_MAGs_by_Species_Cluster.tsv /mnt/lustre/scratch/fcoutinho/StG_23/Metadata/TARA_Salazar_19_MGs_Metadata.tsv 5 100 1 0.5 Temperature,Oxygen,ChlorophyllA,Salinity,Iron.5m,Ammonium.5m 50 Prokaryotes

cd /mnt/lustre/scratch/fcoutinho/StG_23/Viruses/Niche_Models/

Rscript /mnt/lustre/bio/users/fcoutinho/StG_23/Scripts/004_StG_23_ANN_Niche_Modelling.R /mnt/lustre/scratch/fcoutinho/StG_23/Viruses/Abundance/RPKM_Abundance_NA.tsv /mnt/lustre/scratch/fcoutinho/MaLME/Marine_Viral_Communities_Sample_Metadata.tsv 5 100 1 0.5 Temperature,Salinity,Oxygen,Chlorophyll_A,Iron_5m,Ammonium_5m 50 Viruses



