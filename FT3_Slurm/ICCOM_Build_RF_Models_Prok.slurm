#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO

#SBATCH --time=00-04:00:00
#SBATCH --job-name=ProkRF
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=15G
#SBATCH --output=/mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/Job_Logs/%J.out
#SBATCH --error=/mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/Job_Logs/%J.err
#SBATCH --mail-type=All
#SBATCH --mail-user=fhernandes@icm.csic.es

##############################################
module load cesga/system miniconda3/22.11.1-1
conda activate /mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/Conda_Envs/ml

cd /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/ICCOMM/Prokaryotes/RF/

pred_file="/mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/StG23/Prokaryotes/Abundance/Transposed_RPKM_Abundance_OceanDNA_All_Species_Rep_MAGs_by_MAG_ID.tsv" #"/mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/StG23/Prokaryotes/Abundance/Top100_Transposed_RPKM_Abundance_OceanDNA_All_Species_Rep_MAGs_by_MAG_ID.tsv" 
resp_file="/mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/StG23/Metadata/TARA_MGs_Paired_With_Salazar_Data_Info+Metadata_Updated_with_Guidi.tsv"

resp_vars="Temperature Oxygen ChlorophyllA Salinity Iron.5m Ammonium.5m Mean.Flux.at.150m NPP.8d.VGPM..mgC.m2.day." #resp_vars="Mean.Flux.at.150m NPP.8d.VGPM..mgC.m2.day."
prd_idx="Sample_UID"
rsp_idx="Group"
prefix="ProkZ"
threads=47

python3 /mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/repos/ICM_Code/EcoModels.py \
--response_table $resp_file \
--predictor_table $pred_file \
--response_vars $resp_vars \
--response_index $rsp_idx \
--predictor_index $prd_idx \
--build_models RF \
--threads $threads \
--z_transform True \
--prefix $prefix


prefix="ProkZTop100"
threads=6
pred_vars="OceanDNA-b16309 OceanDNA-b25904 OceanDNA-b25897 OceanDNA-b45014 OceanDNA-b45131 OceanDNA-b16286 OceanDNA-a449 OceanDNA-b16691 OceanDNA-b43155 OceanDNA-b25900 OceanDNA-b25738 OceanDNA-b16363 OceanDNA-b16317 OceanDNA-b16405 OceanDNA-b25301 OceanDNA-b28717 OceanDNA-b25781 OceanDNA-b22023 OceanDNA-a113 OceanDNA-b43160 OceanDNA-b8887 OceanDNA-b16578 OceanDNA-b25736 OceanDNA-a434 OceanDNA-b26056 OceanDNA-b41711 OceanDNA-b2444 OceanDNA-b44852 OceanDNA-b25417 OceanDNA-a633 OceanDNA-b18515 OceanDNA-b25819 OceanDNA-a117 OceanDNA-b41655 OceanDNA-b24487 OceanDNA-b20104 OceanDNA-b42077 OceanDNA-a93 OceanDNA-b5053 OceanDNA-b28264 OceanDNA-b34523 OceanDNA-b42461 OceanDNA-b26068 OceanDNA-b43072 OceanDNA-b25795 OceanDNA-a455 OceanDNA-b19757 OceanDNA-b44718 OceanDNA-b27101 OceanDNA-b16552 OceanDNA-b25306 OceanDNA-b2181 OceanDNA-b38823 OceanDNA-b41651 OceanDNA-b2483 OceanDNA-b25793 OceanDNA-b2171 OceanDNA-b25218 OceanDNA-b37046 OceanDNA-b28311 OceanDNA-b9944 OceanDNA-b42068 OceanDNA-b24350 OceanDNA-b8720 OceanDNA-b44823 OceanDNA-b12280 OceanDNA-b1154 OceanDNA-b26561 OceanDNA-b39627 OceanDNA-b43201 OceanDNA-b42483 OceanDNA-b24481 OceanDNA-b44719 OceanDNA-b43110 OceanDNA-a198 OceanDNA-a640 OceanDNA-b45032 OceanDNA-b44721 OceanDNA-b19229 OceanDNA-b9950 OceanDNA-b32485 OceanDNA-b25912 OceanDNA-b25899 OceanDNA-a190 OceanDNA-b7699 OceanDNA-b39704 OceanDNA-b42527 OceanDNA-b25141 OceanDNA-b9322 OceanDNA-a98 OceanDNA-b28721 OceanDNA-b28720 OceanDNA-b42286 OceanDNA-b2136 OceanDNA-b36751 OceanDNA-b42706 OceanDNA-b8277 OceanDNA-b24489 OceanDNA-b39211 OceanDNA-b45059"

python3 /mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/repos/ICM_Code/EcoModels.py \
--predictor_vars $pred_vars \
--response_table $resp_file \
--predictor_table $pred_file \
--response_vars $resp_vars \
--response_index $rsp_idx \
--predictor_index $prd_idx \
--build_models RF \
--threads $threads \
--z_transform True \
--prefix $prefix

#--skip_predictor_importance True \
#--test_predictor_table $pred_file \
#--predictor_vars $pred_vars \
#--z_transform False \
#--run_hpt False \
#--hpt_k_cv 3 \
#--min_pred_prev 10 \