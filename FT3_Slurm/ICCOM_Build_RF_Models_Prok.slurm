#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO

#SBATCH --time=00-03:00:00
#SBATCH --job-name=ProkRF
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=5G
#SBATCH --output=/mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/Job_Logs/%J.out
#SBATCH --error=/mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/Job_Logs/%J.err
#SBATCH --mail-type=All
#SBATCH --mail-user=fhernandes@icm.csic.es

##############################################
module load cesga/system miniconda3/22.11.1-1
conda activate /mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/Conda_Envs/basic

cd /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/ICCOMM/

pred_file="/mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/StG23/Prokaryotes/Abundance/Transposed_RPKM_Abundance_OceanDNA_All_Species_Rep_MAGs_by_MAG_ID.tsv" 
resp_file="/mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/StG23/Metadata/TARA_MGs_Paired_With_Salazar_Data_Info+Metadata_Updated_with_Guidi.tsv"
prefix="Prok"
resp_vars="Temperature Oxygen ChlorophyllA Salinity Iron.5m Ammonium.5m Mean.Flux.at.150m NPP.8d.VGPM..mgC.m2.day." #resp_vars="Mean.Flux.at.150m NPP.8d.VGPM..mgC.m2.day."
prd_idx="Sample_UID"
rsp_idx="Group"
threads=23

python3 /mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/repos/ICM_Code/EcoModels.py \
--response_table $resp_file \
--predictor_table $pred_file \
--response_vars $resp_vars \
--response_index $rsp_idx \
--predictor_index $prd_idx \
--build_models RF \
--threads $threads \
--run_hpt True \
--z_transform True \
--prefix $prefix

#--z_transform False \
#--run_hpt False \
#--hpt_k_cv 3 \]
#--min_pred_prev 10 \