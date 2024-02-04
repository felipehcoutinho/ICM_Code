#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=01-00:00:00
#SBATCH --job-name=Vpop_Rumen_Virathon
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=25G
#SBATCH --output=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.out
#SBATCH --error=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=felipehcoutinho@gmail.com

##############################################
module load python/3.8.5
module load blast/2.7.1

#python3 /mnt/lustre/bio/users/fcoutinho/Scripts/virathon_v0.1.py --genome_files /mnt/lustre/scratch/fcoutinho/Rumen/Assembled/Rumen_Viruses_Renamed.fasta --gene /mnt/lustre/scratch/fcoutinho/Rumen/Assembled/Rumen_Viruses_Renamed.fna --make_pops True --threads 48

python3 /mnt/lustre/bio/users/fcoutinho/Scripts/virathon_v0.1.py --genome_files /mnt/lustre/scratch/fcoutinho/Rumen/Assembled/Rumen_Viruses_Renamed.fasta --gene /mnt/lustre/scratch/fcoutinho/Rumen/Assembled/Rumen_Viruses_Renamed.fna --make_pops True --threads 48 --parse_only True --info_output Rumen_Viruses_Vpop_Info.tsv