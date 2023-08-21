#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=25-00:00:00
#SBATCH --job-name=Asmbl_MLP_Profiles_Large_Fraction
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=750G
#SBATCH --output=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.out
#SBATCH --error=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=fhernandes@icm.csic.es

##############################################
module load python/3.8.5
module load spades/3.13.0

python3 /mnt/lustre/bio/users/fcoutinho/Scripts/Virathon.py --max_ram 730 --threads 48 --assemble True --raw_read_table /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Metadata/Malaspina_Profiles_PA_Metagenome_Info.tsv
