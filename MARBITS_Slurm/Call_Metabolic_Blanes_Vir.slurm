#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=02-00:00:00
#SBATCH --job-name=Metabolic_Blanes_Vir
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=250G
#SBATCH --output=/mnt/smart/users/fcoutinho/Job_Logs/jobLog_%J.out
#SBATCH --error=/mnt/smart/users/fcoutinho/Job_Logs/jobLog_%J.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=fhernandes@icm.csic.es

##############################################
# cd /mnt/smart/users/fcoutinho/Blanes_Vir_Nestor/Sequences/Viruses/
# mkdir Split_CDS
# cd Split_CDS

# module load python

# python3 /mnt/smart/users/fcoutinho/ICM_Code/Get_Seqs.py --fetch_all True --group True --id_split_sep DUMMY --protein True --id_split_pos 0 --input /mnt/smart/users/fcoutinho/Blanes_Vir_Nestor/Sequences/Viruses/dsDNAphage_Blanes_virus.faa

cd /mnt/smart/users/fcoutinho/Blanes_Vir_Nestor/
module load metabolic
source activate

METABOLIC-G.pl -t 47 -m-cutoff 0.5 -in /mnt/smart/users/fcoutinho/Blanes_Vir_Nestor/Sequences/Viruses/Split_CDS/ -kofam-db small -o Blanes_Viruses_METABOLIC_Results


