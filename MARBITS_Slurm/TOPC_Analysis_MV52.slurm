#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=01-12:00:00
#SBATCH --job-name=mv52_pipeline
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=250G
#SBATCH --output=/mnt/smart/users/fcoutinho/Job_Logs/jobLog_%J.out
#SBATCH --error=/mnt/smart/users/fcoutinho/Job_Logs/jobLog_%J.err
#SBATCH --mail-type=All
#SBATCH --mail-user=fhernandes@icm.csic.es

##############################################
#iPhop
module load iphop/iphop_env

conda activate /mnt/smart/apps/condaEnvs/iphop_env

cd /mnt/smart/users/fcoutinho/TARA_Polar/MV52/

mkdir iPhop_Results

iphop predict --num_threads 23 --fa_file /mnt/smart/users/fcoutinho/TARA_Polar/MV52/MV52_Genomic.fasta --min_score 75 --db_dir /mnt/smart/users/fcoutinho/Databases/iPHoP_db/Aug_2023_pub_rw --out_dir /mnt/smart/users/fcoutinho/TARA_Polar/MV52/iPhop_Results/