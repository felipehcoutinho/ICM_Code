#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --time=00-03:00:00
#SBATCH --job-name=Prok_RERF_StG23
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

cd /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/StG23/Prokaryotes/RE_RF_Models

Rscript /mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/repos/ICM_Code/StG23/004_Build_RF_Reverse_Ecology_Models.R /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/StG23/Metadata/TARA_MGs_Paired_With_Salazar_Data_Info+Metadata_Updated_with_Guidi.tsv Temperature,Oxygen,ChlorophyllA,Salinity,PO4,NO3,HCO3,Si,Ammonium.5m,NO2NO3,Iron.5m,Mean.Flux.at.150m,NPP.8d.VGPM..mgC.m2.day. /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/StG23/Prokaryotes/Abundance/RPKM_Abundance_OceanDNA_All_Species_Rep_MAGs_by_MAG_ID.tsv 5 0.7 1000 500 1 47 3 Prokaryotes_RPKM_All_Seqs

Rscript /mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/repos/ICM_Code/StG23/007_Importance_Analysis_RF_Reverse_Ecology_Models.R /mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/Databases/OceanDNA/Sub_OceanDNA_MAGs_Info_Species_Rep_only.tsv /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/StG23/Prokaryotes/RE_RF_Models/Prokaryotes_RPKM_All_Seqs_Reverse_Ecology_Best_RF_Importance.tsv Phylum 0.01 1 Prokaryotes_RPKM_All_Seqs

Rscript /mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/repos/ICM_Code/StG23/004_Build_RF_Reverse_Ecology_Models.R /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/StG23/Metadata/TARA_MGs_Paired_With_Salazar_Data_Info+Metadata_Updated_with_Guidi.tsv Temperature,Oxygen,ChlorophyllA,Salinity,PO4,NO3,HCO3,Si,Ammonium.5m,NO2NO3,Iron.5m,Mean.Flux.at.150m,NPP.8d.VGPM..mgC.m2.day. /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/StG23/Prokaryotes/Abundance/Percentage_Abundance_OceanDNA_All_Species_Rep_MAGs_by_MAG_ID.tsv 5 0.7 1000 500 1 47 3 Prokaryotes_Percentage_All_Seqs

Rscript /mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/repos/ICM_Code/StG23/007_Importance_Analysis_RF_Reverse_Ecology_Models.R /mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/Databases/OceanDNA/Sub_OceanDNA_MAGs_Info_Species_Rep_only.tsv /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/StG23/Prokaryotes/RE_RF_Models/Prokaryotes_Percentage_All_Seqs_Reverse_Ecology_Best_RF_Importance.tsv Phylum 0.01 1 Prokaryotes_Percentage_All_Seqs

#Rscript /mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/repos/ICM_Code/StG23/004_Build_RF_Reverse_Ecology_Models.R /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/StG23/Metadata/TARA_MGs_Paired_With_Salazar_Data_Info+Metadata_Updated_with_Guidi.tsv Temperature,Oxygen,ChlorophyllA,Salinity,PO4,NO3,HCO3,Si,Ammonium.5m,NO2NO3,Iron.5m,Mean.Flux.at.150m,NPP.8d.VGPM..mgC.m2.day. /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/StG23/Prokaryotes/Abundance/Raw_Abundance_OceanDNA_All_Species_Rep_MAGs_by_MAG_ID.tsv 5 0.7 1000 500 1 47 3 Prokaryotes_Raw_All_Seqs

#Rscript /mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/repos/ICM_Code/StG23/004_Build_RF_Reverse_Ecology_Models.R /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/StG23/Metadata/TARA_MGs_Paired_With_Salazar_Data_Info+Metadata_Updated_with_Guidi.tsv Temperature,Mean.Flux.at.150m,NPP.8d.VGPM..mgC.m2.day. /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/StG23/Prokaryotes/Abundance/Raw_Abundance_OceanDNA_All_Species_Rep_MAGs_by_MAG_ID.tsv 5 0.7 5000 500 1 60 3

#Rscript /mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/repos/ICM_Code/StG23/007_Importance_Analysis_RF_Reverse_Ecology_Models.R /mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/Databases/OceanDNA/Sub_OceanDNA_MAGs_Info_Species_Rep_only.tsv /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/StG23/Prokaryotes/RE_RF_Models/Reverse_Ecology_Best_RF_Importance.tsv Phylum 0.01 1 Prokaryotes