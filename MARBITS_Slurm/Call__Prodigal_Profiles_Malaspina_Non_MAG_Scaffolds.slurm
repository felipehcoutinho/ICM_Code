#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=03-00:00:00
#SBATCH --job-name=Prodigal_Profiles_Malaspina_Non_MAG_Scaffolds
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=150G
#SBATCH --output=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.out
#SBATCH --error=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=felipehcoutinho@gmail.com

##############################################
module load prodigal/2.6.3 

prodigal -q -p meta -a Profiles_Malaspina_Non_MAG_Scaffolds_CDS.faa -d Profiles_Malaspina_Non_MAG_Scaffolds_CDS.fna -f gff -i /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/Profiles_Malaspina_Non_MAG_Scaffolds.fasta -o Profiles_Malaspina_Non_MAG_Scaffolds_CDS.gff