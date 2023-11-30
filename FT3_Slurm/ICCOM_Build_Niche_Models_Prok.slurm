#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO

#SBATCH --time=00-02:30:00
#SBATCH --job-name=ANN_Niche_Prok_RPKM
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=4G
#SBATCH --output=/mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/Job_Logs/%J.out
#SBATCH --error=/mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/Job_Logs/%J.err
#SBATCH --mail-type=All
#SBATCH --mail-user=fhernandes@icm.csic.es

##############################################
module load cesga/system miniconda3/22.11.1-1
conda activate /mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/Conda_Envs/basic

cd /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/ICCOMM/

resp_file="/mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/StG23/Prokaryotes/Abundance/Transposed_RPKM_Abundance_OceanDNA_All_Species_Rep_MAGs_by_MAG_ID.tsv" #"/mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/StG23/Prokaryotes/Abundance/Transposed_RPKM_Abundance_OceanDNA_All_Species_Rep_MAGs_by_MAG_ID.tsv" #"/mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/StG23/Prokaryotes/Abundance/Top100_Transposed_RPKM_Abundance_OceanDNA_All_Species_Rep_MAGs_by_MAG_ID.tsv" #"/mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/StG23/Prokaryotes/Abundance/Transposed_RPKM_Abundance_OceanDNA_All_Species_Rep_MAGs_by_MAG_ID.tsv"
pred_file="/mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/StG23/Metadata/TARA_MGs_Paired_With_Salazar_Data_Info+Metadata_Updated_with_Guidi.tsv"
test_file="/mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/StG23/Metadata/TARA_MGs_Paired_With_Salazar_Data_Info+Metadata_Updated_with_Guidi_Scenario_Warming+3.tsv"
prefix="ANNProkZ"
pred_vars="Temperature Oxygen ChlorophyllA Salinity Iron.5m Ammonium.5m"
#resp_vars="OceanDNA-a1001  OceanDNA-a1038  OceanDNA-a1061  OceanDNA-a1068  OceanDNA-a1070  OceanDNA-a1075  OceanDNA-a1079  OceanDNA-a1082  OceanDNA-a1128  OceanDNA-a113  OceanDNA-a1143 OceanDNA-a1795 OceanDNA-a1744 OceanDNA-a1301"
#resp_vars="OceanDNA-a1795 OceanDNA-a1744 OceanDNA-a1301"
prd_idx="Group"
rsp_idx="Sample_UID"
threads=23

python3 /mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/repos/ICM_Code/EcoModels.py \
--response_table $resp_file \
--predictor_table $pred_file \
--predictor_vars $pred_vars \
--prefix $prefix \
--build_models ANN \
--min_resp_prev 50 \
--z_transform True \
--test_predictor_table $test_file \
--threads $threads

#--test_predictor_table "/mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/StG23/Metadata/TARA_MGs_Paired_With_Salazar_Data_Info+Metadata_Updated_with_Guidi_Scenario_Warming+3.tsv" \
#--response_vars $resp_vars \
#--z_transform True \
#--run_hpt True \