#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO

#SBATCH --time=00-06:00:00
#SBATCH --job-name=ANN_Niche_Vir_RPKM
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=10G
#SBATCH --output=/mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/Job_Logs/%J.out
#SBATCH --error=/mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/Job_Logs/%J.err
#SBATCH --mail-type=All
#SBATCH --mail-user=fhernandes@icm.csic.es

##############################################
module load cesga/system miniconda3/22.11.1-1
conda activate /mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/Conda_Envs/ml

cd /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/ICCOMM/Viruses/Niche_Models/

resp_file="/mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/StG23/Viruses/Abundance/RPKM_Abundance_NA.tsv" #"/mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/StG23/Viruses/Abundance/Head_10K_RPKM_Abundance_NA.tsv" #
pred_file="/mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/StG23/Metadata/Marine_Viral_Communities_Sample_Metadata_Updated_with_Guidi.tsv"
prefix="VirZ"
pred_vars="Temperature Salinity Oxygen Chlorophyll_A Iron_5m Ammonium_5m"
resp_vars=""
rsp_idx="Sequence"
prd_idx="MG_UID"
threads=47

python3 /mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/repos/ICM_Code/EcoModels.py \
--predictor_table $pred_file \
--predictor_index $prd_idx \
--predictor_vars $pred_vars \
--response_table $resp_file \
--response_index $rsp_idx \
--transpose_response True \
--build_models ANN \
--min_resp_prev 50 \
--z_transform True \
--threads $threads \
--prefix $prefix

#--response_vars $resp_vars \
