#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=02-00:00:00
#SBATCH --job-name=Metabolic_Coral_Reefs
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=250G
#SBATCH --output=/mnt/smart/users/fcoutinho/Job_Logs/jobLog_%J.out
#SBATCH --error=/mnt/smart/users/fcoutinho/Job_Logs/jobLog_%J.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=fhernandes@icm.csic.es

##############################################

module load metabolic
source activate

cd /mnt/smart/users/fcoutinho/Coral_Reef_Viruses/

#METABOLIC-G.pl -t 23 -m-cutoff 0.5 -in /mnt/smart/users/fcoutinho/Coral_Reef_Viruses/viruses/proteins_per_bin/ -kofam-db small -o Coral_Reefs_vMAGs_METABOLIC_Results

METABOLIC-G.pl -t 47 -m-cutoff 0.5 -in /mnt/smart/users/fcoutinho/Coral_Reef_Viruses/bMAGs/Proteins/ -kofam-db small -o Coral_Reefs_bMAGs_METABOLIC_Results


