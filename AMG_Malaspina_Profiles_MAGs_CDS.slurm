#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=01-00:00:00
#SBATCH --job-name=Malaspina_Profiles_CDS_AMG_Hunter
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=50G
#SBATCH --output=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%A_%a.out
#SBATCH --error=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%A_%a.err
#SBATCH --mail-type=All
#SBATCH --mail-user=felipehcoutinho@gmail.com
#SBATCH --array=1-20%5

##############################################
module load diamond/2.0.7
module load python/3.8.5
module load hmmer
module load perl

file=$(ls *.fasta | sed -n ${SLURM_ARRAY_TASK_ID}p)
#python3 /mnt/lustre/bio/users/fcoutinho/Scripts/AMG_Hunter.py --annotate True --threads 12 --cds $file

python3 /mnt/lustre/bio/users/fcoutinho/Scripts/AMG_Hunter.py --annotate True --parse_only True --threads 12 --cds $file --info_cds_output CDS_Info_$file.tsv --info_genome_output Genome_Info_$file.tsv