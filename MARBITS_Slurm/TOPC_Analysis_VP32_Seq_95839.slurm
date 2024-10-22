#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=00-12:00:00
#SBATCH --job-name=iphop_topc_vir
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=50G
#SBATCH --output=/mnt/smart/users/fcoutinho/Job_Logs/jobLog_%J.out
#SBATCH --error=/mnt/smart/users/fcoutinho/Job_Logs/jobLog_%J.err
#SBATCH --mail-type=All
#SBATCH --mail-user=fhernandes@icm.csic.es

##############################################
#Load the prodigla module
module load prodigal/2.6.3

#Move to the relvant directory
cd /mnt/smart/users/fcoutinho/TARA_Polar/VP32/

#Run prodigal
prodigal -q -p meta -a TOPC_Viral_Seq_95839.faa -d TOPC_Viral_Seq_95839.fna -f gff -i TOPC_Viral_Seq_95839.fasta -o TOPC_Viral_Seq_95839.gff