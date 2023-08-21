#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=1-00:00:00
#SBATCH --job-name=GTDBtk_TOPC_Virathon
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=250G
#SBATCH --output=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.out
#SBATCH --error=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=felipehcoutinho@gmail.com

##############################################
module load gtdbtk/1.5
source activate

gtdbtk classify_wf --extension fasta --genome_dir /mnt/lustre/bio/shared/topc/tara.polar.mags.curated.fasta.sequences/output/08_hq-bins-pooled-with-renamed-contigs/ --out_dir /mnt/lustre/scratch/fcoutinho/TARA_Polar/TOPC_GTDBtk_Output --pplacer_cpus 12 --cpus 48 --prefix TOPC_GTDBtk
