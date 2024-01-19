#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO

#SBATCH --time=00-06:00:00
#SBATCH --job-name=Antarctic_Lagoon_DADA2
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=25G
#SBATCH --output=/mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/Job_Logs/%J.out
#SBATCH --error=/mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/Job_Logs/%J.err
#SBATCH --mail-type=All
#SBATCH --mail-user=fhernandes@icm.csic.es

##############################################
module load bioconductor

cd /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/Antarctic_Lagoon/DADA/

Rscript /mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/repos/ICM_Code/Antarctic_Lagoons/DADA2_Antarctic_Lagoons.R