#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=00-03:00:00
#SBATCH --job-name=DBSetup_Trembl+SwissProt
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=150G
#SBATCH --output=/mnt/smart/scratch/vir/felipe/Job_Logs/jobLog_%J.out
#SBATCH --error=/mnt/smart/scratch/vir/felipe/Job_Logs/jobLog_%J.err
#SBATCH --mail-type=All
#SBATCH --mail-user=fhernandes@icm.csic.es

##############################################
###Trembl set + Swiss Prot
cd /mnt/smart/scratch/vir/felipe/Databases/Trembl+SwissProt/

#Get swissprot files (release 2025-06-18)
# wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/uniprot_sprot_viruses.xml.gz
# wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/uniprot_sprot_archaea.xml.gz
# wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/uniprot_sprot_bacteria.xml.gz

# #Get trembl files (release 2025-06-18)
# wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/uniprot_trembl_archaea.xml.gz

# wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/uniprot_trembl_viruses.xml.gz

# wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/uniprot_trembl_bacteria.xml.gz

# gunzip *gz

# module load python

#Extract sequences and get the relevant info from uniprot xml files. because we are using the tax specific files, no need to filter out eukaryotes nd unclassified
# python3  /mnt/smart/scratch/vir/felipe/Repos/ICM_Code/Parse_Uniprot_XML.py --xml *xml --fasta_output SwissProt+Trembl_Arch+Bac+Vir_Sequences.fasta --info_output SwissProt+Trembl_Arch+Bac+Vir_Info.tsv

# rm -f *xml

module load mmseqs2
mmseqs createdb SwissProt+Trembl_Arch+Bac+Vir_Sequences.fasta DB_SwissProt+Trembl_Arch+Bac+Vir_Sequences
mmseqs createindex DB_SwissProt+Trembl_Arch+Bac+Vir_Sequences tmp



