#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=01-12:00:00
#SBATCH --job-name=Metabolic_Antarctic_Lagoons
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=250G
#SBATCH --output=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.out
#SBATCH --error=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=felipehcoutinho@gmail.com

##############################################
cd /mnt/smart/users/fcoutinho/Antarctic_Lagoons/Binning/AL_MAGs_R1

module load prodigal

for genome in *fa; do prodigal -q -p single -a "CDS_${genome%.fa}.faa" -d "Genes_${genome%.fa}.fna" -f gff -i $genome -o "${genome%.fa}.gff" -t 4 ; done

# mkdir Genomes
# mkdir Gene_Calls
# mv *fa Genomes/
# mv CDS* *gff Genes* Gene_Calls/

module load metabolic
source activate

cd /mnt/smart/users/fcoutinho/Antarctic_Lagoons/Binning/AL_MAGs_R1

###Not working because metabolic files ar ein lustre
METABOLIC-G.pl -t 12 -m-cutoff 0.5 -in /mnt/smart/users/fcoutinho/Antarctic_Lagoons/Binning/AL_MAGs_R1/Proteins_for_Metabolic/ -kofam-db small -o Antarctic_Lagoons_MAGs_METABOLIC_Results
