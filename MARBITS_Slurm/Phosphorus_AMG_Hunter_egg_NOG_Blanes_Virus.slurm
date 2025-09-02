#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=07-12:00:00
#SBATCH --job-name=AMG_Hunter_Blanes
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=30G
#SBATCH --output=/mnt/smart/scratch/vir/felipe/Job_Logs/jobLog_%J.out
#SBATCH --error=/mnt/smart/scratch/vir/felipe/Job_Logs/jobLog_%J.err
#SBATCH --mail-type=All
#SBATCH --mail-user=fhernandes@icm.csic.es

##############################################
module load python
module load hmmer

###Setup db
# cd /mnt/smart/scratch/vir/felipe/Databases/eggNOG/
# wget http://eggnog5.embl.de/download/eggnog_5.0/per_tax_level/2/2_hmms.tar.gz
# tar -zxvf 2_hmms.tar.gz
# ls 2/ > hmmlist

# while read i; do cat 2/$i >> eggNOG_Bacteria.hmm ; done < hmmlist

# hmmpress eggNOG_Bacteria.hmm

cd /mnt/smart/scratch/vir/felipe/Blanes_Vir_Nestor/AMG_Hunter_Annotation/MSebastian_Phosphorus/

python3 /mnt/smart/scratch/vir/felipe/Repos/ICM_Code/AMG_Hunter.py --threads 24 --cds /mnt/smart/scratch/vir/felipe/Blanes_Vir_Nestor/Sequences/Viruses/dsDNAphage_Blanes_virus.faa --info_cds_output CDS_Info_dsDNAphage_Blanes_VirusxeggNOG_Bacteria.tsv --info_genome_output Genome_Info_nfo_dsDNAphage_Blanes_VirusxeggNOG_Bacteria --kegg_annotate True --kegg_db /mnt/smart/scratch/vir/felipe/Databases/eggNOG/eggNOG_Bacteria.hmm --kegg_min_score 100

grep -v "Role" /mnt/smart/scratch/vir/felipe/Blanes_Vir_Nestor/AMG_Hunter_Annotation/MSebastian_Phosphorus/COGs_Pmetabolism.txt | cut -f 1 | cut -d " " -f 1 > List_P_COGs.txt

head -n 1 CDS_Info_dsDNAphage_Blanes_VirusxeggNOG_Bacteria.tsv > P_Best_Hits_CDS_Info_dsDNAphage_Blanes_VirusxeggNOG_Bacteria.tsv

grep -w -f List_P_COGs.txt CDS_Info_dsDNAphage_Blanes_VirusxeggNOG_Bacteria.tsv >> P_Best_Hits_CDS_Info_dsDNAphage_Blanes_VirusxeggNOG_Bacteria.tsv
