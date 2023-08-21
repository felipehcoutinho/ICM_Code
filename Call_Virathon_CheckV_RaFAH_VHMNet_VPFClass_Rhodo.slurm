#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=1-12:00:00
#SBATCH --job-name=Rhodo_Virathon
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=250G
#SBATCH --output=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.out
#SBATCH --error=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=felipehcoutinho@gmail.com

##############################################
module load python/3.8.5
#module load checkv/0.7.0
#module load vpf-class
#module load prodigal
#module load blast/2.7.1
module load diamond/2.0.7
module load hmmer/3.3
#module load virhostmatcher-net/1.0
module load perl
#module load R


#export CHECKVDB="/mnt/lustre/repos/bio/databases/public/checkv/checkv-db-v1.0/"
#export PATH="/mnt/lustre/bio/users/fcoutinho/Scripts/:$PATH"

#python3 /mnt/lustre/bio/users/fcoutinho/Scripts/virathon_v0.1.py --genome_files /mnt/lustre/scratch/fcoutinho/Rhodo/Rhodo_Virus_Genomes.fasta --call_checkv True --call_vpf_class True --call_rafah True --threads 24

#python3 /mnt/lustre/bio/users/fcoutinho/Scripts/virathon_v0.1.py --genome_files /mnt/lustre/scratch/fcoutinho/Rhodo/Rhodo_Virus_Genomes.fasta --cds /mnt/lustre/scratch/fcoutinho/Rhodo/Rhodo_Virus_Genomes.faa --gene /mnt/lustre/scratch/fcoutinho/Rhodo/Rhodo_Virus_Genomes.fna --threads 24 --call_vhm True --vhmnet_mode_short True --make_pops True

python3 /mnt/lustre/bio/users/fcoutinho/Scripts/AMG_Hunter.py --annotate True --threads 48 --cds /mnt/lustre/scratch/fcoutinho/Rhodo/Rhodo_Virus_Genomes.faa --info_cds_output Rhodo_Virus_CDS_Info.tsv --info_genome_output Rhodo_Virus_Genome_Info.tsv

#python3 /mnt/lustre/bio/users/fcoutinho/Scripts/virathon_v0.1.py --genome_files /mnt/lustre/scratch/fcoutinho/Rhodo/Rhodo_Virus_Genomes.fasta --call_checkv True --call_vpf_class True --call_rafah True --cds /mnt/lustre/scratch/fcoutinho/Rhodo/Rhodo_Virus_Genomes.faa --gene /mnt/lustre/scratch/fcoutinho/Rhodo/Rhodo_Virus_Genomes.fna --call_vhm True --vhmnet_mode_short True --make_pops True --parse_only True