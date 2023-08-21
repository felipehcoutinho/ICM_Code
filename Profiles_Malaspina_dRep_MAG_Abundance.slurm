#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=10-00:00:00
#SBATCH --job-name=Profiles_Malaspina_dRep_MAG_Abundance
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=44
#SBATCH --mem=150G
#SBATCH --output=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.out
#SBATCH --error=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.err
#SBATCH --mail-type=All
#SBATCH --mail-user=felipehcoutinho@gmail.com

##############################################
module load python/3.8.5
module load bowtie2/2.3.4.1
module load samtools

python3 /mnt/lustre/bio/users/fcoutinho/Scripts/virathon_v0.1.py --threads 48 --abundance_table True --abundance_rpkm True --raw_read_table /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Profiles_Malaspina_Metagenome_Info.tsv --genome /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/Dereplicated_MAGs_Original_and_Refined/dereplicated_genomes/*fa --merge True --index True

