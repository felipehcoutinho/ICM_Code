#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=01-00:00:00
#SBATCH --job-name=AMG_Abundance_Profiles_Malaspina_Virus
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=50G
#SBATCH --output=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.out
#SBATCH --error=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=felipehcoutinho@gmail.com

##############################################
module load python/3.8.5

python3 /mnt/lustre/bio/users/fcoutinho/Scripts/AMG_Hunter.py --annotate True --parse_only True --abundance /mnt/lustre/scratch/elopez/5_bowtie_results/After_checkV_output/RPKM.tsv --cds /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Viruses/Profiles_Malaspina_Viruses_CheckV_Trimmed_Genomes.faa