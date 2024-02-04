#!/bin/bash
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=5-12:00:00
#SBATCH --job-name=VirMatcher_TOPC
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=250G
#SBATCH --output=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.out
#SBATCH --error=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=felipehcoutinho@gmail.com

##############################################
eval "$(conda shell.bash hook)"

conda activate /mnt/lustre/bio/users/fcoutinho/envs/VirMatcher
export PATH="/home/fcoutinho/WIsH:$PATH"

VirMatcher --preparer --gtdbtk-out-dir /mnt/lustre/scratch/fcoutinho/TARA_Polar/TOPC_GTDBtk_Output/ --gtdbtk-in-dir /mnt/lustre/scratch/fcoutinho/TARA_Polar/08_hq-bins-pooled-with-renamed-contigs/ -v /mnt/lustre/scratch/fcoutinho/TARA_Polar/Renamed_Viruses_Genomes_Filtered_All_Tara_Polar_Coassemblies_Contigs.fasta --threads 48 -o VirMatcher_Prep_Output --python-aggregator

#VirMatcher --help