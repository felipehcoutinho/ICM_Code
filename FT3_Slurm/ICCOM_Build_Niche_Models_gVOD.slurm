#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO

#SBATCH --time=00-00:10:00
#SBATCH --job-name=ANN_Niche_Prok+Vir_gVOD
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --output=/mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/Job_Logs/%J.out
#SBATCH --error=/mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/Job_Logs/%J.err
#SBATCH --mail-type=All
#SBATCH --mail-user=fhernandes@icm.csic.es

##############################################
module load cesga/system miniconda3/22.11.1-1
conda activate /mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/Conda_Envs/ml

cd /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/ICCOMM/gVOD/


resp_pred_file="/mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/Databases/gVOD/gVOD_Sub_Data_Renamed_Vars.tsv"
test_file="/mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/Databases/gVOD/gVOD_Sub_Data_Renamed_Vars_Scenarios.tsv"
prefix="gVOD_V+P"
pred_vars="Depth_m Temperature Salinity" # Chlorophyll_a_g_per_l  PO4_M
resp_vars="Viruses_per_ml Prokaryotes_per_ml" #Prokaryotes_per_ml   #Synechococcus_per_ml Prochlorococcus_per_ml
rspprd_idx="No"
threads=10

python3 /mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/repos/ICM_Code/EcoModels.py \
--resp_pred_table $resp_pred_file \
--response_vars $resp_vars \
--predictor_vars $pred_vars \
--resp_pred_index $rspprd_idx \
--test_predictor_table $test_file \
--z_transform True \
--threads $threads \
--build_models ANN \
--prefix $prefix
