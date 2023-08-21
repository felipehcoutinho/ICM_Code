#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=00-02:00:00
#SBATCH --job-name=MMSeqs_PETHxGLUVAB3
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=200G
#SBATCH --output=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.out
#SBATCH --error=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=felipehcoutinho@gmail.com

##############################################
module load mmseqs2/15c7762
mmseqs easy-search /mnt/lustre/scratch/fcoutinho/Plastisphere/PETH_IDESA.faa /mnt/lustre/repos/bio/databases/local/GLUVAB3/MMSeqs2_CDS_DB/DB_GLUVAB3_CDS /mnt/lustre/scratch/fcoutinho/Plastisphere/PETH_IDESAxGLUVAB3_CDS.blastp tmp --threads 24 --max-seqs 10000
