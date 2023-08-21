#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=03-00:00:00
#SBATCH --job-name=Metabolic_Profiles_Malaspina
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=250G
#SBATCH --output=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.out
#SBATCH --error=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=felipehcoutinho@gmail.com

##############################################
module load metabolic
source activate

METABOLIC-G.pl -t 24 -m-cutoff 0.5 -in /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/CDS_by_MAG/ -kofam-db small -o Profiles_Malaspina_MAGs_METABOLIC_Results


