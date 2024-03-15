#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=30-00:00:00
#SBATCH --job-name=iphop_mlp_profiles
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=750G
#SBATCH --output=/mnt/smart/users/fcoutinho/Job_Logs/jobLog_%J.out
#SBATCH --error=/mnt/smart/users/fcoutinho/Job_Logs/jobLog_%J.err
#SBATCH --mail-type=All
#SBATCH --mail-user=fhernandes@icm.csic.es

##############################################
module load iphop/iphop_env

conda activate /mnt/smart/apps/condaEnvs/iphop_env

cd /mnt/smart/users/fcoutinho/MLP_Prok_Virus/

mkdir iPhop_Results

iphop predict --num_threads 47 --fa_file /mnt/smart/users/fcoutinho/MLP_Prok_Virus/Sequences/Profiles_Malaspina_Viruses_CheckV_Trimmed_Genomes.fasta --db_dir /mnt/lustre/bio/users/fcoutinho/Databases/iPhop/Sept_2021_pub_rw/ --out_dir /mnt/smart/users/fcoutinho/MLP_Prok_Virus/iPhop_Results/