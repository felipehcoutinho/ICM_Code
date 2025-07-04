#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=12-00:00:00
#SBATCH --job-name=iphop_Blanes_Viromes
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=36
#SBATCH --mem=450G
#SBATCH --output=/mnt/smart/scratch/vir/felipe/Job_Logs/jobLog_%J.out
#SBATCH --error=/mnt/smart/scratch/vir/felipe/Job_Logs/jobLog_%J.err
#SBATCH --mail-type=All
#SBATCH --mail-user=fhernandes@icm.csic.es

##############################################
#iPhop
module load iphop/iphop_env

conda activate /mnt/smart/apps/condaEnvs/iphop_env

cd /mnt/smart/scratch/vir/felipe/Viromes_Blanes/

mkdir iPhop_Results

iphop predict --num_threads 35 --fa_file VP_Reps_Batch_1_2_Trimmed_Viral_Genomes.fasta --min_score 75 --db_dir /mnt/smart/scratch/vir/felipe/Databases/iPHoP_db/Aug_2023_pub_rw --out_dir iPhop_Results/