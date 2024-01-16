#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO

#SBATCH --time=04-00:00:00
#SBATCH --job-name=iphop_Blanes
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=170G
#SBATCH --output=/mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/Job_Logs/%J.out
#SBATCH --error=/mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/Job_Logs/%J.err
#SBATCH --mail-type=All
#SBATCH --mail-user=fhernandes@icm.csic.es

##############################################
module load cesga/system miniconda3/22.11.1-1

conda activate /mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/Conda_Envs/iphop_env

cd /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/Blanes_Vir_Nestor/

mkdir iPHoP_Results

iphop predict --num_threads 47 --fa_file /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/Blanes_Vir_Nestor/Vir_MAGs_Sequences/dsDNAphage_Blanes_virus.fasta --db_dir /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/Databases/iPHoP/Aug_2023_pub_rw/ --out_dir /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/Blanes_Vir_Nestor/iPHoP_Results/