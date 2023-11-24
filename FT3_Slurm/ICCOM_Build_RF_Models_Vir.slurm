#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO

#SBATCH --time=00-00:30:00
#SBATCH --job-name=RF_Niche_Vir_RPKM
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=10G
#SBATCH --output=/mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/Job_Logs/%J.out
#SBATCH --error=/mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/Job_Logs/%J.err
#SBATCH --mail-type=All
#SBATCH --mail-user=fhernandes@icm.csic.es

##############################################
module load cesga/system miniconda3/22.11.1-1
conda activate /mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/Conda_Envs/basic

cd /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/ICCOMM/

pred_file="/mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/StG23/Viruses/Abundance/RPKM_Abundance_NA.tsv" #"/mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/StG23/Viruses/Abundance/Head_10K_RPKM_Abundance_NA.tsv" #
resp_file="/mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/StG23/Metadata/Marine_Viral_Communities_Sample_Metadata_Updated_with_Guidi.tsv"
prefix="TestVirRPKM"
resp_vars="Temperature Salinity Oxygen Chlorophyll_A Iron_5m Ammonium_5m Mean.Flux.at.150m NPP.8d.VGPM..mgC.m2.day."
pred_vars=""
prd_idx="Sequence"
rsp_idx="MG_UID"
threads=47

python3 /mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/repos/ICM_Code/EcoModels.py \
--response_table $resp_file \
--predictor_table $pred_file \
--response_vars $resp_vars \
--response_index $rsp_idx \
--predictor_index $prd_idx \
--transpose_predictors True \
--build_rf True \
--run_hpt False \
--threads $threads \
--min_pred_prev 10 \
--hpt_k_cv 0 \
--prefix $prefix

#--z_transform False \