#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=00-03:00:00
#SBATCH --job-name=PPS_Malaspina_GLUVAB3_Virathon
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=300G
#SBATCH --output=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.out
#SBATCH --error=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.err
#SBATCH --mail-type=END
#SBATCH --mail-user=felipehcoutinho@gmail.com

##############################################
module load python/3.8.5
module load mmseqs2/15c7762

python3 /mnt/lustre/bio/users/fcoutinho/Scripts/virathon_v0.1.py --pairwise True --pps_subject_db /mnt/lustre/repos/bio/databases/local/GLUVAB3/MMSeqs2_CDS_DB/DB_GLUVAB3_CDS --cds /mnt/lustre/scratch/fcoutinho/Malaspina_Virus/Malaspina_CoAssembly_Viral_Scaffolds_Renamed.faa --threads 24 --genome_files /mnt/lustre/scratch/fcoutinho/Malaspina_Virus/Malaspina_CoAssembly_Viral_Scaffolds_Renamed.fasta
