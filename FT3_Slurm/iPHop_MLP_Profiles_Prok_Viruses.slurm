#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO

#SBATCH --time=00-04:00:00
#SBATCH --job-name=Test_iphop
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=36
#SBATCH --mem=100G
#SBATCH --output=/mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/Job_Logs/%J.out
#SBATCH --error=/mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/Job_Logs/%J.err
#SBATCH --mail-type=All
#SBATCH --mail-user=fhernandes@icm.csic.es

##############################################
module load cesga/system miniconda3/22.11.1-1

conda activate /mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/Conda_Envs/iphop_env

cd /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/Databases/

mkdir iphop_test_results

iphop predict --num_threads 35 --fa_file test_input_phages.fna --db_dir iphop_db/Test_db_rw/ --out_dir iphop_test_results/test_input_phages_iphop