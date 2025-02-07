#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=6-00:00:00
#SBATCH --job-name=Asmbl_TOPC
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=250G
#SBATCH --output=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.out
#SBATCH --error=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=fhernandes@icm.csic.es

##############################################
module load python/3.8.5
module load spades/3.13.0

cd /mnt/smart/users/fcoutinho/TARA_Polar/VP32/Reassembly/

python3 /mnt/smart/shared/mlCourse/fcoutinho_data/Virathon.py --max_ram 240 --threads 23 --assemble True --raw_read_table /mnt/smart/users/fcoutinho/TARA_Polar/VP32/Reassembly/Assembly_Instructions.tsv 
