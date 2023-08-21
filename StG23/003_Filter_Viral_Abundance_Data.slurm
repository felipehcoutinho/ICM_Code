#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=00-00:10:00
#SBATCH --job-name=FilterVirAbdStG23
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=5G
#SBATCH --output=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.out
#SBATCH --error=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=fhernandes@icm.csic.es

##############################################
module load R

Rscript /mnt/lustre/bio/users/fcoutinho/StG_23/Scripts/003_Process_Prok_Abund_Data.R /mnt/lustre/scratch/fcoutinho/StG_23/Prokaryotes/Abundance/RPKM_Abundance_OceanDNA_All_Species_Rep_MAGs_Scaffolds.tsv /mnt/lustre/bio/users/fcoutinho/OceanDNA_MAGs/Sub_OceanDNA_MAGs_Info.tsv 100 /mnt/lustre/scratch/fcoutinho/StG_23/Prokaryotes/Abundance/RPKM_Abundance_OceanDNA_All_Species_Rep_MAGs_by_Species_Cluster.tsv
