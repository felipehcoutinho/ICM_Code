#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=00-24:00:00
#SBATCH --job-name=Virathon_CheckV_Coral_Reef_Viruses
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=20G
#SBATCH --output=/mnt/smart/users/fcoutinho/Job_Logs/%J.Log.out
#SBATCH --error=/mnt/smart/users/fcoutinho/Job_Logs/%J.Log.err
#SBATCH --mail-type=All
#SBATCH --mail-user=fhernandes@icm.csic.es

##############################################
module load checkv/0.7.0 python/3.8.5

export CHECKVDB="/mnt/lustre/repos/bio/databases/public/checkv/checkv-db-v1.0/"

cd /mnt/smart/users/fcoutinho/Coral_Reef_Viruses/viruses/

python3 /mnt/smart/users/fcoutinho/virathon/Virathon.py --genome_files /mnt/smart/users/fcoutinho/Coral_Reef_Viruses/viruses/Coral_Reef_Viruses_All_vMAG_Sequences.fasta --call_checkv True --info_output CheckV_Info_All_Coral_Reef_Virus_Sequences.tsv --threads 47 