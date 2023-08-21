#!/bin/bash
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=5-00:00:00
#SBATCH --job-name=VAMB_Profiles_Malaspina
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=750G
#SBATCH --output=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.out
#SBATCH --error=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=felipehcoutinho@gmail.com

##############################################
module load vamb
source activate

vamb --outdir /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/VAMB_Binning/VAMB_Bins_Round_1/ --fasta /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/VAMB_Binning/Profiles_Malaspina_Scaffolds_Filtered_2500_bp.fasta --jgi /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/VAMB_Binning/Corrected_Profiles_Malaspina_JGI_Depth.tsv -o C --minfasta 200000 -p 8