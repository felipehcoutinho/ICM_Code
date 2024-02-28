#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=10-00:00:00
#SBATCH --job-name=Metabolic_Blanes
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=250G
#SBATCH --output=/mnt/smart/users/fcoutinho/Job_Logs/jobLog_%A_%a.out
#SBATCH --error=/mnt/smart/users/fcoutinho/Job_Logs/jobLog_%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=fhernandes@icm.csic.es
#SBATCH --array=1-24%8

##############################################
module load metabolic
source activate

cd /mnt/smart/users/fcoutinho/Blanes_Vir_Nestor/Metabolic_Output/

d=/mnt/smart/users/fcoutinho/Blanes_Vir_Nestor/Metabolic_Output/Batch_${SLURM_ARRAY_TASK_ID}/

METABOLIC-G.pl -t 12 -kofam-db small -m-cutoff 0.5 -in $d -o Metabolic_Outputs_Batch_${SLURM_ARRAY_TASK_ID}

