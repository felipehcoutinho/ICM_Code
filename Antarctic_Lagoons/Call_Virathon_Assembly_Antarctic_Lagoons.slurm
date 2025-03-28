#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=60-00:00:00
#SBATCH --job-name=Asmbl_Antartic_Lagoons
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=750G
#SBATCH --output=/mnt/smart/users/fcoutinho/Job_Logs/jobLog_%J.out
#SBATCH --error=/mnt/smart/users/fcoutinho/Job_Logs/jobLog_%J.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=fhernandes@icm.csic.es

##############################################
module load python/3.8.5
# module load spades/4.0.0
module load megahit

cd /mnt/smart/users/fcoutinho/Antarctic_Lagoons/WGS_Assemblies/

# python3 /mnt/smart/users/fcoutinho/Repos/virathon/Virathon.py --max_ram 730 --threads 48 --assemble True --raw_read_table /mnt/smart/users/fcoutinho/Antarctic_Lagoons/Metadata/WGS_Samples.tsv

python3 /mnt/smart/users/fcoutinho/Repos/virathon/Virathon.py --max_ram 730 --threads 47 --assembler megahit --raw_read_table /mnt/smart/users/fcoutinho/Antarctic_Lagoons/Metadata/WGS_Samples.tsv

perl /mnt/smart/users/fcoutinho/Repos/ICM_Code/assembly_stats_table.pl */*fa  > MG_Assembly_Stats.tsv