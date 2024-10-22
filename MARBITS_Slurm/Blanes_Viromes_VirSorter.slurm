#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=00-04:00:00
#SBATCH --job-name=vs2_blanes
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=150G
#SBATCH --output=/mnt/smart/users/fcoutinho/Job_Logs/jobLog_%J.out
#SBATCH --error=/mnt/smart/users/fcoutinho/Job_Logs/jobLog_%J.err
#SBATCH --mail-type=All
#SBATCH --mail-user=fhernandes@icm.csic.es

##############################################
#module load cesga/2020
#module load virsorter2/2.2.4-python-3.9.9

module load virsorter
conda activate vs2
#mkdir /mnt/smart/users/fcoutinho/OctoMicro/WGS_Assemblies
cd /mnt/smart/users/fcoutinho/Viromes_Blanes/

#virsorter config --set HMMSEARCH_THREADS=23
#rm -fr VS2_Output/
virsorter run --seqfile /mnt/smart/users/fcoutinho/Viromes_Blanes/Blanes_Viromes_All_Scaffolds_Min_Len_5kbp.fasta --working-dir VS2_Output --tmpdir VS2_tmp --exclude-lt2gene --jobs 47 --include-groups dsDNAphage,NCLDV,ssDNA,lavidaviridae

#--use-conda-off
