#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=01-00:00:00
#SBATCH --job-name=Salvage_Bins_Malaspina_Profiles
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=50G
#SBATCH --output=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.out
#SBATCH --error=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=felipehcoutinho@gmail.com

##############################################
module load python

python3 /mnt/lustre/bio/users/fcoutinho/Scripts/QueroBins2.py --refine_bins True --refinement_tax_level class --call_checkm True --parse_only True --assemblies /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Filtered_Renamed_S*fasta