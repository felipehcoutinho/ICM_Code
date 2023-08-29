#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --time=00-03:00:00
#SBATCH --job-name=VirRERFStG23
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=61
#SBATCH --mem=20G
#SBATCH --output=/mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/Job_Logs/%J.out
#SBATCH --error=/mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/Job_Logs/%J.err
#SBATCH --mail-type=All
#SBATCH --mail-user=fhernandes@icm.csic.es

##############################################
module load cesga/system miniconda3/22.11.1-1
conda activate /mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/Conda_Envs/basic

Rscript /mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/repos/ICM_Code/StG23/004_Build_RF_Reverse_Ecology_Models.R /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/StG23/Metadata/Marine_Viral_Communities_Sample_Metadata_Updated_with_Guidi.tsv Temperature,Ammonium_5m,Nitrate_5m,Nitrite_5m,Phosphate,Silicate,Iron_5m,Oxygen,Salinity,Chlorophyll_A,POC,PIC,Mean.Flux.at.150m,NPP.8d.VGPM..mgC.m2.day. /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/StG23/Viruses/Abundance/StG23_Viruses_Min_Prev_50_Raw_Abundance.tsv 1000 0.5 5000 500 1 60

