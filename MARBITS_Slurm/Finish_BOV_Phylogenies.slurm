#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=02-06:00:00
#SBATCH --job-name=Finish_BOV_Phylogenies
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=50G
#SBATCH --output=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.out
#SBATCH --error=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.err
#SBATCH --mail-type=All
#SBATCH --mail-user=felipehcoutinho@gmail.com

##############################################
#module load python
#module load hmmer
module load muscle
module load fasttree
#module load cd-hit


muscle -in Merged_NR_Matched_K14652.faa -out Aligned_Merged_NR_Matched_K14652.faa -maxhours 24 -maxiters 3

FastTreeMP -nosupport -out Tree_Merged_NR_Matched_K14652.newick Aligned_Merged_NR_Matched_K14652.faa




