#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=01-12:00:00
#SBATCH --job-name=AMG_Hunter_Blanes
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=10G
#SBATCH --output=/mnt/smart/users/fcoutinho/Job_Logs/jobLog_%J.out
#SBATCH --error=/mnt/smart/users/fcoutinho/Job_Logs/jobLog_%J.err
#SBATCH --mail-type=All
#SBATCH --mail-user=fhernandes@icm.csic.es

##############################################
module load python
module load hmmer
cd /mnt/smart/users/fcoutinho/Blanes_Vir_Nestor/AMG_Hunter_Annotation/

# python3 /mnt/smart/users/fcoutinho/Repos/ICM_Code/AMG_Hunter.py --threads 24 --cds /mnt/smart/users/fcoutinho/Blanes_Vir_Nestor/Sequences/Viruses/VP_Reps_dsDNAphage_Blanes_virus.faa --info_cds_output CDS_Info_VP_Reps_dsDNAphage_Blanes_Virus.tsv --info_genome_output Genome_Info_VP_Reps_dsDNAphage_Blanes_Virus.tsv --kegg_annotate True

python3 /mnt/smart/users/fcoutinho/Repos/ICM_Code/AMG_Hunter.py --threads 24 --cds /mnt/smart/users/fcoutinho/Blanes_Vir_Nestor/Sequences/Viruses/VP_Reps_dsDNAphage_Blanes_virus.faa --info_cds_output CDS_Info_VP_Reps_dsDNAphage_Blanes_Virus.tsv --info_genome_output Genome_Info_VP_Reps_dsDNAphage_Blanes_Virus.tsv --kegg_annotate True --parse_only True --genome_abundance /mnt/smart/users/fcoutinho/Blanes_Vir_Nestor/AMG_Hunter_Annotation/VP_Abundances_with_Rep_Ids.tsv