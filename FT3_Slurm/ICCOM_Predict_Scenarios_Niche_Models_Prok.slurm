#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO

#SBATCH --time=00-00:30:00
#SBATCH --job-name=ANN_Niche_Prok_Pred_Scenarios
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=4G
#SBATCH --output=/mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/Job_Logs/%J.out
#SBATCH --error=/mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/Job_Logs/%J.err
#SBATCH --mail-type=All
#SBATCH --mail-user=fhernandes@icm.csic.es

##############################################
module load cesga/system miniconda3/22.11.1-1
conda activate /mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/Conda_Envs/ml

cd /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/ICCOMM/Prokaryotes/Niche_Models/

resp_file="/mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/StG23/Prokaryotes/Abundance/Transposed_RPKM_Abundance_OceanDNA_All_Species_Rep_MAGs_by_MAG_ID.tsv" #"/mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/StG23/Prokaryotes/Abundance/Transposed_RPKM_Abundance_OceanDNA_All_Species_Rep_MAGs_by_MAG_ID.tsv" #"/mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/StG23/Prokaryotes/Abundance/Top100_Transposed_RPKM_Abundance_OceanDNA_All_Species_Rep_MAGs_by_MAG_ID.tsv" #"/mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/StG23/Prokaryotes/Abundance/Transposed_RPKM_Abundance_OceanDNA_All_Species_Rep_MAGs_by_MAG_ID.tsv"
pred_file="/mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/StG23/Metadata/TARA_MGs_Paired_With_Salazar_Data_Info+Metadata_Updated_with_Guidi.tsv"
test_file="/mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/StG23/Metadata/TARA_MGs_Paired_With_Salazar_Data_Info+Metadata_Updated_with_Guidi_Scenarios.tsv"
pred_vars="Temperature Oxygen ChlorophyllA Salinity Iron.5m Ammonium.5m"
resp_vars="OceanDNA-b25897 OceanDNA-b45131 OceanDNA-b25736 OceanDNA-b41655 OceanDNA-b25900 OceanDNA-b25738 OceanDNA-b25795 OceanDNA-b34523 OceanDNA-b16309 OceanDNA-b28264"
#resp_vars="OceanDNA-a1795 OceanDNA-a1744 OceanDNA-a1301"
prd_idx="Group"
rsp_idx="Sample_UID"
prefix="ProkZ"
threads=10

#Perform predictions on pre-computed models
python3 /mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/repos/ICM_Code/EcoModels.py \
--build_models ANN \
--load_model /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/ICCOMM/Prokaryotes/Niche_Models/ProkZ_ANN_MOR_Model.pkl \
--test_predictor_table $test_file \
--skip_predictor_importance True \
--response_table $resp_file \
--predictor_table $pred_file \
--predictor_vars $pred_vars \
--prefix $prefix \
--min_resp_prev 50 \
--z_transform True \
--debug True \
--threads $threads