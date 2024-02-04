#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=3-00:00:00
#SBATCH --job-name=VHMNet_Rumen
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=500G
#SBATCH --output=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.out
#SBATCH --error=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=felipehcoutinho@gmail.com

##############################################
module load python/3.8.5
module load virhostmatcher-net/1.0

#python3 /mnt/lustre/bio/users/fcoutinho/Scripts/virathon_v0.1.py --genome_files /mnt/lustre/scratch/fcoutinho/Rumen/Assembled/Rumen_Viruses_Renamed.fasta --cds /mnt/lustre/scratch/fcoutinho/Rumen/Assembled/Rumen_Viruses_Renamed.faa --call_vhm True --vhmnet_mode_short True --threads 48

VirHostMatcher-Net.py -q Split_Genomes_Rumen_Viruses_Renamed/ -o  VHMNet_Output_Rumen_Viruses_Renamed/ -t 24 -i tmp -n 10 --short-contig