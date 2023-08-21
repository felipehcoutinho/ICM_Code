#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=03-00:00:00
#SBATCH --job-name=iPhop_Soil
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=150G
#SBATCH --output=/mnt/lustre/scratch/fcoutinho/Job_Logs/%J.Log.out
#SBATCH --error=/mnt/lustre/scratch/fcoutinho/Job_Logs/%J.Log.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=fhernandes@icm.csic.es

##############################################
module load iphop
conda activate iphop_env

iphop predict --fa_file /mnt/lustre/scratch/fcoutinho/Soil_Elena/Post_CheckV_Trimmed_True_Virus_Soil_Genomes.fasta --out_dir /mnt/lustre/scratch/fcoutinho/Soil_Elena/iPhop_Output/ --min_score 70 --num_threads 23 --db_dir /mnt/lustre/bio/users/fcoutinho/Databases/iPhop/Sept_2021_pub_rw/