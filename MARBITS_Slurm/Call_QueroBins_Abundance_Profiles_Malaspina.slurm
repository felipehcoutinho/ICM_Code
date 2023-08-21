#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=00-03:00:00
#SBATCH --job-name=Abundance_Profiles_Malaspina_QueroBins_Test
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=36
#SBATCH --mem=10G
#SBATCH --output=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.out
#SBATCH --error=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=felipehcoutinho@gmail.com

##############################################
module load python/3.8.5
module load bowtie2
module load samtools

python3 /mnt/lustre/bio/users/fcoutinho/Scripts/QueroBins2.py --threads 36 --assemblies /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies/Group_5.fasta --metagenomes_extension fastq.gz --metagenomes_dir /mnt/lustre/repos/bio/projects/malaspina/release1/station019/Epi/0.22-3/metaDNA/MP0311/rawDatas/cleanReads/
