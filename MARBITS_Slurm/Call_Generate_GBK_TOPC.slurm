#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=01-00:00:00
#SBATCH --job-name=GBK_TOPC
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=250G
#SBATCH --output=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.out
#SBATCH --error=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=felipehcoutinho@gmail.com

##############################################
module load perl

perl /mnt/lustre/bio/users/fcoutinho/Scripts/Generate_GBK3.pl --skip_unk --explode --annot /mnt/lustre/scratch/fcoutinho/TARA_Polar/Enriched_Renamed_Viruses_Genomes_Filtered_All_Tara_Polar_Coassemblies_CDSxUniRef100.blastp --protein /mnt/lustre/scratch/fcoutinho/TARA_Polar/Renamed_Viruses_Genomes_Filtered_All_Tara_Polar_Coassemblies_Contigs.faa --nucleotide /mnt/lustre/scratch/fcoutinho/TARA_Polar/Renamed_Viruses_Genomes_Filtered_All_Tara_Polar_Coassemblies_Contigs.fasta --gff /mnt/lustre/scratch/fcoutinho/TARA_Polar/Renamed_Viruses_Genomes_Filtered_All_Tara_Polar_Coassemblies_Contigs.gff