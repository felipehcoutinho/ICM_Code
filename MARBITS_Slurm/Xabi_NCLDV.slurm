#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=3-00:00:00
#SBATCH --job-name=NCLDV_Phylogenies
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=250G
#SBATCH --output=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.out
#SBATCH --error=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.err
#SBATCH --mail-type=All
#SBATCH --mail-user=felipehcoutinho@gmail.com
#SBATCH --constraint=supermicro

##############################################
#module load python/3.8.5
#module load hmmer
#module load muscle
#module load fasttree
#module load mmseqs2/15c7762
#module load diamond/2.0.7
#odule load perl

#python3 /mnt/lustre/bio/users/fcoutinho/Scripts/virathon_v0.1.py --genome /mnt/lustre/scratch/xlopezal/00_NCLDV/NCLDV-2022-04-19/all_NCLDV.fasta --cds /mnt/lustre/scratch/xlopezal/00_NCLDV/NCLDV-2022-04-19/all_NCLDV.faa --og_phylogeny True --min_cluster 5 --threads 24

#python3 /mnt/lustre/bio/users/fcoutinho/Scripts/AMG_Hunter.py --annotate True --threads 24 --cds /mnt/lustre/scratch/xlopezal/00_NCLDV/NCLDV-2022-04-19/all_NCLDV.faa

#eval "$(conda shell.bash hook)"
#conda activate /mnt/lustre/bio/users/fcoutinho/envs/coconet
#coconet run --fasta /mnt/lustre/scratch/xlopezal/00_NCLDV/NCLDV-2022-04-19/all_NCLDV.fasta --output /mnt/lustre/scratch/fcoutinho/xlopez/NCLDV_Test/CoCoNet_NCLDV --threads 24

module load coconet
source activate

coconet run --fasta /mnt/lustre/scratch/xlopezal/00_NCLDV/NCLDV-2022-04-19/all_NCLDV.fasta --output /mnt/lustre/scratch/fcoutinho/xlopez/NCLDV_Test/CoCoNet_NCLDV --threads 24

