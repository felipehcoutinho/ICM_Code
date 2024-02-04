#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=01-12:00:00
#SBATCH --job-name=Malaspina_Profiles_Compress
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=5G
#SBATCH --output=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.out
#SBATCH --error=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=fhernades@icm.csic.es

##############################################
for i in Assembly_MP* ; do tar -zcvf Contigs_and_Scaffolds_$i.tgz $i/scaffolds.fasta $i/contigs.fasta ; mv Contigs_and_Scaffolds_$i.tgz /mnt/lustre/bio/users/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Scaffolds/Large_Fraction/ ; done ;


