#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --time=00-01:30:00
#SBATCH --job-name=VirRERFStG23
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=50G
#SBATCH --output=/mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/Job_Logs/%J.out
#SBATCH --error=/mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/Job_Logs/%J.err
#SBATCH --mail-type=All
#SBATCH --mail-user=fhernandes@icm.csic.es

##############################################
module load cesga/system miniconda3/22.11.1-1
conda activate /mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/Conda_Envs/basic

cd /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/StG23/Viruses/RE_RF_Models/

#Fit RF models for selected variables, using perc abundances, and all virals equences regardles sof prevalence across samples
    Rscript /mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/repos/ICM_Code/StG23/004_Build_RF_Reverse_Ecology_Models.R /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/StG23/Metadata/Marine_Viral_Communities_Sample_Metadata_Updated_with_Guidi.tsv Chlorophyll_A,Mean.Flux.at.150m,NPP.8d.VGPM..mgC.m2.day. /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/StG23/Viruses/Abundance/Percentage_Abundance_NA.tsv 5 0.7 1000 500 1 47 1 NPP+Flux+Chloro_Viruses_Percentage_All_Seqs

    Rscript /mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/repos/ICM_Code/StG23/007_Importance_Analysis_RF_Reverse_Ecology_Models.R /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/StG23/Viruses/Info_Sequences/Full_Info_IMGVR_Scaffolds_Comp_75-100_Max_Conta_0_vOTU_Representatives.tsv /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/StG23/Viruses/RE_RF_Models/Test_Viruses_Percentage_All_Seqs_Reverse_Ecology_Best_RF_Importance.tsv Host_Phylum_PHIST 0.01 1 NPP+Flux+Chloro_Viruses_Percentage_All_Seqs

#Rscript /mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/repos/ICM_Code/StG23/004_Build_RF_Reverse_Ecology_Models.R /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/StG23/Metadata/Marine_Viral_Communities_Sample_Metadata_Updated_with_Guidi.tsv Temperature,Ammonium_5m,Nitrate_5m,Nitrite_5m,Phosphate,Silicate,Iron_5m,Oxygen,Salinity,Chlorophyll_A,POC,PIC,Mean.Flux.at.150m,NPP.8d.VGPM..mgC.m2.day. /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/StG23/Viruses/Abundance/StG23_Viruses_Min_Prev_50_Raw_Abundance.tsv 5 0.7 1000 500 1 47 10 Viruses_Raw_Min_Prev_50

#Rscript /mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/repos/ICM_Code/StG23/007_Importance_Analysis_RF_Reverse_Ecology_Models.R /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/StG23/Viruses/Info_Sequences/Full_Info_IMGVR_Scaffolds_Comp_75-100_Max_Conta_0_vOTU_Representatives.tsv /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/StG23/Viruses/RE_RF_Models/Viruses_Raw_Reverse_Ecology_Best_RF_Importance.tsv Host_Phylum_PHIST 0.01 1 Viruses_Raw_Min_Prev_50

#Temperature,Ammonium_5m,Nitrate_5m,Nitrite_5m,Phosphate,Silicate,Iron_5m,Oxygen,Salinity,Chlorophyll_A,POC,PIC,Mean.Flux.at.150m,NPP.8d.VGPM..mgC.m2.day.
#Ammonium_5m,Nitrate_5m,Nitrite_5m,Iron_5m,Chlorophyll_A,Mean.Flux.at.150m,NPP.8d.VGPM..mgC.m2.day.

#Rscript /mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/repos/ICM_Code/StG23/004_Build_RF_Reverse_Ecology_Models.R /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/StG23/Metadata/Marine_Viral_Communities_Sample_Metadata_Updated_with_Guidi.tsv Mean.Flux.at.150m /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/StG23/Viruses/Abundance/RPKM_Abundance_NA.tsv 5 0.7 1000 500 1 47 1 Test_Viruses_RPKM_All_Seqs

#Rscript /mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/repos/ICM_Code/StG23/004_Build_RF_Reverse_Ecology_Models.R /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/StG23/Metadata/Marine_Viral_Communities_Sample_Metadata_Updated_with_Guidi.tsv Ammonium_5m,Nitrate_5m,Nitrite_5m,Iron_5m,Chlorophyll_A,Mean.Flux.at.150m,NPP.8d.VGPM..mgC.m2.day. /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/StG23/Viruses/Abundance/RPKM_Abundance_NA.tsv 5 0.7 1000 500 1 47 3 Viruses_RPKM_All_Seqs

#Rscript /mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/repos/ICM_Code/StG23/007_Importance_Analysis_RF_Reverse_Ecology_Models.R /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/StG23/Viruses/Info_Sequences/Full_Info_IMGVR_Scaffolds_Comp_75-100_Max_Conta_0_vOTU_Representatives.tsv /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/StG23/Viruses/RE_RF_Models/Viruses_RPKM_All_Seqs_Reverse_Ecology_Best_RF_Importance.tsv Host_Phylum_PHIST 0.01 1 Viruses_RPKM_All_Seqs

#Rscript /mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/repos/ICM_Code/StG23/004_Build_RF_Reverse_Ecology_Models.R /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/StG23/Metadata/Marine_Viral_Communities_Sample_Metadata_Updated_with_Guidi.tsv Mean.Flux.at.150m /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/StG23/Viruses/Abundance/Percentage_Abundance_NA.tsv 5 0.7 1000 500 1 47 1 Test_Viruses_Percentage_All_Seqs

#Rscript /mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/repos/ICM_Code/StG23/007_Importance_Analysis_RF_Reverse_Ecology_Models.R /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/StG23/Viruses/Info_Sequences/Full_Info_IMGVR_Scaffolds_Comp_75-100_Max_Conta_0_vOTU_Representatives.tsv /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/StG23/Viruses/RE_RF_Models/Test_Viruses_Percentage_All_Seqs_Reverse_Ecology_Best_RF_Importance.tsv Host_Phylum_PHIST 0.01 1 Test_Viruses_Percentage_All_Seqs

#It seems changing from RPKM to Percentage does not change the relativ eimportance of Nanoarchaeota infecting viruses in the RF models for POC flux.
#But RMSE and cor metrics are ALWAYS identical for Raw, Perc and RPKM inout tables. Something is wrong.

#Rscript /mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/repos/ICM_Code/StG23/004_Build_RF_Reverse_Ecology_Models.R /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/StG23/Metadata/Marine_Viral_Communities_Sample_Metadata_Updated_with_Guidi.tsv Ammonium_5m,Nitrate_5m,Nitrite_5m,Iron_5m,Chlorophyll_A,Mean.Flux.at.150m,NPP.8d.VGPM..mgC.m2.day. /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/StG23/Viruses/Abundance/Percentage_Abundance_NA.tsv 5 0.7 1000 500 1 47 3 Viruses_Percentage_All_Seqs

#Rscript /mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/repos/ICM_Code/StG23/007_Importance_Analysis_RF_Reverse_Ecology_Models.R /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/StG23/Viruses/Info_Sequences/Full_Info_IMGVR_Scaffolds_Comp_75-100_Max_Conta_0_vOTU_Representatives.tsv /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/StG23/Viruses/RE_RF_Models/Viruses_Percentage_All_Seqs_Reverse_Ecology_Best_RF_Importance.tsv Host_Phylum_PHIST 0.01 1 Viruses_Percentage_All_Seqs

#Rscript /mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/repos/ICM_Code/StG23/004_Build_RF_Reverse_Ecology_Models.R /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/StG23/Metadata/Marine_Viral_Communities_Sample_Metadata_Updated_with_Guidi.tsv Ammonium_5m,Nitrate_5m,Nitrite_5m,Iron_5m,Chlorophyll_A,Mean.Flux.at.150m,NPP.8d.VGPM..mgC.m2.day. /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/StG23/Viruses/Abundance/Raw_Abundance_NA.tsv 5 0.7 1000 500 1 47 3 Viruses_Raw_All_Seqs

#Rscript /mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/repos/ICM_Code/StG23/007_Importance_Analysis_RF_Reverse_Ecology_Models.R /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/StG23/Viruses/Info_Sequences/Full_Info_IMGVR_Scaffolds_Comp_75-100_Max_Conta_0_vOTU_Representatives.tsv /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/StG23/Viruses/RE_RF_Models/Viruses_Raw_All_Seqs_Reverse_Ecology_Best_RF_Importance.tsv Host_Phylum_PHIST 0.01 1 Viruses_Raw_All_Seqs