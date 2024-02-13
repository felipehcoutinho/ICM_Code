#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=05-00:00:00
#SBATCH --job-name=VirSorter_Blanes_MAGs
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=150G
#SBATCH --output=/mnt/smart/users/fcoutinho/Job_Logs/jobLog_%J.out
#SBATCH --error=/mnt/smart/users/fcoutinho/Job_Logs/jobLog_%J.err
#SBATCH --mail-type=All
#SBATCH --mail-user=fhernandes@icm.csic.es

##############################################
module load virsorter/2.2.3
conda activate vs2

cd /mnt/smart/users/fcoutinho/Blanes_Vir_Nestor/Sequences/Prok_MAGs/

#virsorter run --seqfile /mnt/smart/users/fcoutinho/Blanes_Vir_Nestor/Sequences/Prok_MAGs/bbmo_concatenated.fa --jobs 47 --prep-for-dramv --rm-tmpdir --min-score 0.5 all

virsorter run --seqfile /mnt/smart/users/fcoutinho/Blanes_Vir_Nestor/Sequences/Prok_MAGs/bbmo_concatenated.fa --jobs 47 --prep-for-dramv --rm-tmpdir --min-score 0.5 --include-groups dsDNAphage,NCLDV,RNA,ssDNA,lavidaviridae all