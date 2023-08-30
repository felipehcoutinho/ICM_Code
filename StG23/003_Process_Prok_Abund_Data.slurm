#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --time=00-02:30:00
#SBATCH --job-name=FilterProkAbdStG23
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=50G
#SBATCH --output=/mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/Job_Logs/%J.out
#SBATCH --error=/mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/Job_Logs/%J.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=fhernandes@icm.csic.es

##############################################
module load cesga/system miniconda3/22.11.1-1
conda activate /mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/Conda_Envs/basic
cd /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/StG23/Prokaryotes/Abundance/

#Rscript /mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/repos/ICM_Code/StG23/003_Process_Prok_Abund_Data.R /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/StG23/Prokaryotes/Abundance/RPKM_Abundance_OceanDNA_All_Species_Rep_MAGs_Scaffolds.tsv /mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/Databases/OceanDNA/Sub_OceanDNA_MAGs_Info.tsv MAG_ID 50 /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/StG23/Prokaryotes/Abundance/RPKM_Abundance_OceanDNA_All_Species_Rep_MAGs_by_MAG_ID.tsv

Rscript /mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/repos/ICM_Code/StG23/003_Process_Prok_Abund_Data.R /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/StG23/Prokaryotes/Abundance/Percentage_Abundance_OceanDNA_All_Species_Rep_MAGs_Scaffolds.tsv /mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/Databases/OceanDNA/Sub_OceanDNA_MAGs_Info.tsv MAG_ID 50 /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/StG23/Prokaryotes/Abundance/Percentage_Abundance_OceanDNA_All_Species_Rep_MAGs_by_MAG_ID.tsv

Rscript /mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/repos/ICM_Code/StG23/003_Process_Prok_Abund_Data.R /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/StG23/Prokaryotes/Abundance/Raw_Abundance_OceanDNA_All_Species_Rep_MAGs_Scaffolds.tsv /mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/Databases/OceanDNA/Sub_OceanDNA_MAGs_Info.tsv MAG_ID 50 /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/StG23/Prokaryotes/Abundance/Raw_Abundance_OceanDNA_All_Species_Rep_MAGs_by_MAG_ID.tsv
