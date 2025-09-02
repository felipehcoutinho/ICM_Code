#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=07-12:00:00
#SBATCH --job-name=Setup_Phylo_MV52
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=6G
#SBATCH --output=/mnt/smart/scratch/vir/felipe/Job_Logs/jobLog_%J.out
#SBATCH --error=/mnt/smart/scratch/vir/felipe/Job_Logs/jobLog_%J.err
#SBATCH --mail-type=All
#SBATCH --mail-user=fhernandes@icm.csic.es

##############################################


# cd /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/TARA_Polar/MV52

# module load cesga/system miniconda3/22.11.1-1
# conda activate /mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/Conda_Envs/basic

# prodigal -a MV52.faa -d MV52.fna -o MV52.gff -i MV52.fasta

# #Remove Uniref cluster represented solely by Eukaryotic sequences
# cd /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/Databases/Uniref50_03-06-2025/
# python3 Remove_Euk_Only_Cluster_Reps.py
# #output file: /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/Databases/Uniref50_03-06-2025/NonEuk_Clusters_Uniref50.fasta

# #Cluster IMG vr protein sequences by 50% ID
# sbatch /mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/repos/ICM_Code/FT3_Slurm/MV52_Phylo_Setup.slurm

# ###In marbits from now on:
# cd /mnt/smart/users/fcoutinho/Databases/MV52_PhyloDB

# #Merge files
# # cd /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/TARA_Polar/MV52

# # cat /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/TARA_Polar/MV52/MV52.faa /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/Databases/Uniref50_03-06-2025/NonEuk_Clusters_Uniref50.fasta /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/Databases/IMG_VR/IMGVR_HQ_DBclu_rep.fasta > PhyloDB.faa

# # module load cesga/system miniconda3/22.11.1-1
# # conda activate /mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/Conda_Envs/basic

# ##Copy this file to your folder: /mnt/smart/users/fcoutinho/Databases/MV52_PhyloDB/PhyloDB.faa
# #Add sequences with this command:
# cat ../../TARA_Polar/MV52/CDS_MV52.faa >> PhyloDB.faa

# #rUN amg hunter use the version in /mnt/smart/users/fcoutinho/Repos/ICM_Code/AMG_Hunter.py:
# cd /mnt/smart/users/fcoutinho/Databases/MV52_PhyloDB/
# module load python
# module load hmmer
# module load muscle
# module load fasttree
# #module load cd-hit

# #heme oxygenase
# #/mnt/smart/users/fcoutinho/Databases/MV52_PhyloDB/PhyloDB.faa
# python3 /mnt/smart/users/fcoutinho/Repos/ICM_Code/AMG_Hunter.py --threads 12 --kegg_phylogeny True --kegg_ko K21480 --cds PhyloDB.faa --derep_id 0 


# ###making the swiss prot subset 
# cd /mnt/smart/scratch/vir/felipe/Databases/MV52_PhyloDB/Swiss_Prot_Set

# #Get uniprot files (release 2025-06-18)
# wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/uniprot_sprot_viruses.xml.gz
# wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/uniprot_sprot_archaea.xml.gz
# wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/uniprot_sprot_bacteria.xml.gz

# gunzip *gz

# module load python

# #Extract sequences and get the relevant info from uniprot xml files. because we are using the tax specific files, no nee dto filter out eukaryotes nd unclassified
# python3  /mnt/smart/scratch/vir/felipe/Repos/ICM_Code/Parse_Uniprot_XML.py --xml *.xml --fasta_output Swiss_Prot_Arch+Bac+Vir_Sequences.fasta --info_output Swiss_Prot_Arch+Bac+Vir_Info.tsv

# rm -f *xml

# #This time, no clustering of IMG VR sequences 
# python3 /mnt/smart/scratch/vir/felipe/Repos/ICM_Code/Get_Seqs.py --fetch_all True --input Swiss_Prot_Arch+Bac+Vir_Sequences.fasta /mnt/smart/scratch/vir/felipe/Databases/IMGVR/IMGVR_Scaffolds_Comp_75-100_Max_Conta_0_vOTU_Representatives.faa /mnt/smart/scratch/vir/felipe/Databases/MV52_PhyloDB/Secs_ad_K*.faa /mnt/smart/scratch/vir/felipe/Databases/MV52_PhyloDB/MV52.faa --matched_out SwissProt_PhyloDB.faa

# module load R/4.2.2
# Rscript /mnt/smart/scratch/vir/felipe/Databases/MV52_PhyloDB/Make_UniRef_NonEuk_info.R

# mv Minimal_tree_Deco_Info.tsv Swiss_Prot_Set_Minimal_tree_Deco_Info.tsv

# sbatch /mnt/smart/scratch/vir/felipe/Repos/ICM_Code/MARBITS_Slurm/TOPC_Analysis_MV52.slurm


###Trembl set + Swiss Prot
cd /mnt/smart/scratch/vir/felipe/Databases/MV52_PhyloDB/Trembl_Set/

wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/uniprot_trembl_archaea.xml.gz
wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/uniprot_trembl_bacteria.xml.gz
# wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/uniprot_trembl_viruses.xml.gz

gunzip *gz

module load python

#Extract sequences and get the relevant info from uniprot xml files. because we are using the tax specific files, no need to filter out eukaryotes nd unclassified
python3  /mnt/smart/scratch/vir/felipe/Repos/ICM_Code/Parse_Uniprot_XML.py --xml /mnt/smart/scratch/vir/felipe/Databases/MV52_PhyloDB/Swiss_Prot_Set/*.xml /mnt/smart/scratch/vir/felipe/Databases/MV52_PhyloDB/Trembl_Set/*xml --fasta_output SwissProt+Trembl_Arch+Bac+Vir_Sequences.fasta --info_output SwissProt+Trembl_Arch+Bac+Vir_Info.tsv



