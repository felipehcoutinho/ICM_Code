#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO

#SBATCH --time=00-00:30:00
#SBATCH --job-name=ANN_Niche_Prok_RPKM
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --output=/mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/Job_Logs/%J.out
#SBATCH --error=/mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/Job_Logs/%J.err
#SBATCH --mail-type=All
#SBATCH --mail-user=fhernandes@icm.csic.es

##############################################
module load cesga/system miniconda3/22.11.1-1
conda activate /mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/Conda_Envs/basic

cd /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/ICCOM/

resp_file="/mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/StG23/Prokaryotes/Abundance/Transposed_RPKM_Abundance_OceanDNA_All_Species_Rep_MAGs_by_MAG_ID.tsv" #"/mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/ICCOM/Test_100.tsv" # #"/mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/StG23/Prokaryotes/Abundance/RPKM_Abundance_OceanDNA_All_Species_Rep_MAGs_by_MAG_ID.tsv" #"/mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/StG23/Viruses/Abundance/Head_10K_RPKM_Abundance_NA.tsv"  # # RPKM_Abundance_NA.tsv
pred_file="/mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/StG23/Metadata/TARA_MGs_Paired_With_Salazar_Data_Info+Metadata_Updated_with_Guidi.tsv" #"/mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/StG23/Metadata/Marine_Viral_Communities_Sample_Metadata_Updated_with_Guidi.tsv"
prefix="TestProkRPKM"
pred_vars="Temperature Oxygen ChlorophyllA Salinity Iron.5m Ammonium.5m"
resp_vars="OceanDNA-a1001  OceanDNA-a1038  OceanDNA-a1061  OceanDNA-a1068  OceanDNA-a1070  OceanDNA-a1075  OceanDNA-a1079  OceanDNA-a1082  OceanDNA-a1128  OceanDNA-a113  OceanDNA-a1143"

time python3 /mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/repos/ICM_Code/EcoModels.py \
--response_table $resp_file \
--predictor_table $pred_file \
--predictor_vars $pred_vars \
--prefix $prefix \
--index Sample_UID

#--response_vars $resp_vars \