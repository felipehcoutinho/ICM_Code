#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=30-00:00:00
#SBATCH --job-name=Metabolic_Keystones
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=250G
#SBATCH --output=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%A_%a.out
#SBATCH --error=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=fhernandes@icm.csic.es
#SBATCH --array=7,8,9,69-99%5

##############################################
module load metabolic
source activate

d=/mnt/lustre/scratch/fcoutinho/Keystones/Sub_CDS_Files_${SLURM_ARRAY_TASK_ID}/

METABOLIC-G.pl -t 47 -in $d -o keystones_metabolic_faa_outputs_Sub_CDS_Files_${SLURM_ARRAY_TASK_ID}

