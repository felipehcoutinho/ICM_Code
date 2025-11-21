#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=00-03:00:00
#SBATCH --job-name=MV52
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=160G
#SBATCH --output=/mnt/smart/scratch/vir/felipe/Job_Logs/jobLog_%J.out
#SBATCH --error=/mnt/smart/scratch/vir/felipe/Job_Logs/jobLog_%J.err
#SBATCH --mail-type=All
#SBATCH --mail-user=fhernandes@icm.csic.es

##############################################
###Gene calling
# cd /mnt/smart/scratch/vir/felipe/TARA_Polar/MV52/Sequences/

# module load prodigal

# #NOTE: Very important to use prodigal in single mode for this virus, as the meta mode produces different results (possibly with truncated AMGs)
prodigal -q -p single -a MV52.faa -d MV52.fna -f gff -i MV52_Genomic.fasta -o MV52.gff

# /mnt/smart/scratch/vir/felipe/Software/prodigal-gv -q -p single -a MV52.faa -d MV52.fna -f gff -i ../MV52_Genomic.fasta -o MV52.gff


# ###Pharokka Annotation
# module load perl

# perl /mnt/smart/scratch/vir/felipe/Repos/ICM_Code/Generate_GBK3.pl --nucleotide /mnt/smart/scratch/vir/felipe/TARA_Polar/MV52/Sequences/MV52_Genomic.fasta --protein /mnt/smart/scratch/vir/felipe/TARA_Polar/MV52/Sequences/MV52.faa --gff /mnt/smart/scratch/vir/felipe/TARA_Polar/MV52/Sequences/MV52.gff

# mv Output.gbk MV52.gbk

# cd /mnt/smart/scratch/vir/felipe/TARA_Polar/MV52/Annotation

# module purge
# conda activate /mnt/smart/scratch/vir/felipe/envs/pharokka

# pharokka.py --genbank -i /mnt/smart/scratch/vir/felipe/TARA_Polar/MV52/Sequences/MV52.gbk -o Pharokka_MV52 -d /mnt/smart/scratch/vir/felipe/Databases/Pharokka/ --prefix Pharokka_MV52 --threads 4 --skip_mash --force

# conda deactivate

# ###Phold Annotation
# conda activate pholdENV

# phold run -i /mnt/smart/scratch/vir/felipe/TARA_Polar/MV52/Annotation/Pharokka_MV52/Pharokka_MV52.gbk -o MV52_phold -t 24 --database /mnt/smart/scratch/vir/felipe/Databases/Phold_DB/ --cpu --keep_tmp_files --force

###AMG hunter annotation
cd /mnt/smart/scratch/vir/felipe/TARA_Polar/MV52/Annotation
module load python
module load hmmer
module load mmseqs2

# python3 /mnt/smart/scratch/vir/felipe/Repos/ICM_Code/AMG_Hunter.py --threads 24 --cds /mnt/smart/scratch/vir/felipe/TARA_Polar/MV52/Sequences/MV52.faa --info_cds_output CDS_Info_MV52.tsv --info_genome_output Genome_Info_MV52.tsv --pfam_annotate True --uniref_annotate True --uniref_db /mnt/smart/scratch/vir/felipe/Databases/Trembl/TremblDB --kegg_annotate True

# python3 /mnt/smart/scratch/vir/felipe/Repos/ICM_Code/AMG_Hunter.py --threads 24 --cds /mnt/smart/scratch/vir/felipe/TARA_Polar/MV52/Sequences/MV52.faa --info_cds_output CDS_Info_MV52.tsv --info_genome_output Genome_Info_MV52.tsv --kegg_annotate True

python3 /mnt/smart/scratch/vir/felipe/Repos/ICM_Code/AMG_Hunter.py --threads 24 --cds /mnt/smart/scratch/vir/felipe/TARA_Polar/MV52/Sequences/MV52.faa --info_cds_output CDS_Info_MV52.tsv --info_genome_output Genome_Info_MV52.tsv --pfam_annotate True --uniref_annotate True --uniref_db /mnt/smart/scratch/vir/felipe/Databases/Trembl/TremblDB --kegg_annotate True --parse_only True

###iPhop
# module load iphop/iphop_env

# conda activate /mnt/smart/apps/condaEnvs/iphop_env

# cd /mnt/smart/scratch/vir/felipe/TARA_Polar/MV52/

# mkdir iPhop_Results

# iphop predict --num_threads 23 --fa_file /mnt/smart/scratch/vir/felipe/TARA_Polar/MV52/MV52_Genomic.fasta --min_score 75 --db_dir /mnt/smart/scratch/vir/felipe/Databases/iPHoP_db/Aug_2023_pub_rw --out_dir /mnt/smart/scratch/vir/felipe/TARA_Polar/MV52/iPhop_Results/

#Phylogenies
# cd /mnt/smart/scratch/vir/felipe/TARA_Polar/MV52/Phylo_SwissProt/Bitscore_50/

# module load python
# module load hmmer
# # module load muscle
# # module load fasttree
# #module load cd-hit

# #PhyloDb full set (used min score 50)
# #PhyloDB=/mnt/smart/scratch/vir/felipe/Databases/MV52_PhyloDB/PhyloDB.faa
# #PhyloDb swissprot subset (used min score 75)
# PhyloDB=/mnt/smart/scratch/vir/felipe/Databases/MV52_PhyloDB/Swiss_Prot_Set/SwissProt_PhyloDB.faa

# #Annotate PhyloDB with AMG_Hunter
# cd /mnt/smart/scratch/vir/felipe/Databases/MV52_PhyloDB/Swiss_Prot_Set/

# python3 /mnt/smart/scratch/vir/felipe/Repos/ICM_Code/AMG_Hunter.py --threads 48 --cds $PhyloDB --kegg_annotate True --pfam_annotate True

# #Aminolevulinate synthase K00643 Domains: PF00155 Aminotran_1_2 Assembly_MV_52_125
# python3 /mnt/smart/scratch/vir/felipe/Repos/ICM_Code/AMG_Hunter.py --threads 24 --kegg_phylogeny True --kegg_ko K00643 --cds $PhyloDB --derep_id 0 --kegg_min_score 50
#Because this gene has an aminotransferase domain shared with other proteins off different function (e.g. P26505 8-amino-7-oxononanoate synthase) it is necessary to better select sequences (not only they hit with this KO but htis needs to be the best hit and contain the domains)

#cluster with a viral clade among the sequences from viruses infecting alpha and gamma proteobacteria, with some cyanophages mixed in. Related viruses are marine and non-marine saline alkaline

# #Hydroxymethilbilane synthase K01749 Assembly_MV_52_212
# python3 /mnt/smart/scratch/vir/felipe/Repos/ICM_Code/AMG_Hunter.py --threads 24 --kegg_phylogeny True --kegg_ko K01749 --cds $PhyloDB --derep_id 0 --kegg_min_score 50
#clusters with a viral clade with a few other viruses of proteobacteria that is sister to a clade made mostly or archeal sequences. clusters  close to Rickettsiales_2026788_PBIG01000028.1_4_K00228 Very few other viruses also encode this gene.  Related viruses are marine

# #Uroporphirinogen decarboxylase K01599 Assembly_MV_52_214
# python3 /mnt/smart/scratch/vir/felipe/Repos/ICM_Code/AMG_Hunter.py --threads 24 --kegg_phylogeny True --kegg_ko K01599 --cds $PhyloDB --derep_id 0 --kegg_min_score 50
#Clusters as a viral clade sister to alphaproteobacteria. very few viruses have this gene. Related viruses are marine

# #Coproporphyrinogen III-oxidase / K00228 / Assembly_MV_52_213
# python3 /mnt/smart/scratch/vir/felipe/Repos/ICM_Code/AMG_Hunter.py --threads 24 --kegg_phylogeny True --kegg_ko K00228 --cds $PhyloDB --derep_id 0 --kegg_min_score 50 
#Clusters as a viral clade sister to alphaproteobacteria and close to Rickettsiales_2026788_PBIG01000028.1_4_K00228. very few viruses have this gene. Related viruses are marine 

# #heme oxygenase K21480 : Assembly_MV_52_160
# python3 /mnt/smart/scratch/vir/felipe/Repos/ICM_Code/AMG_Hunter.py --threads 24 --kegg_phylogeny True --kegg_ko K21480 --cds $PhyloDB --derep_id 0 --kegg_min_score 50
#all viral sequences cluster together and separated from cellular organisms, except from the cyanophages which cluster with their hosts. Many viruses have the gene. Related viruses are marine, freshwater, and non-marine saline alkaline

# #PFO K05371: Assembly_MV_52_158
# python3 /mnt/smart/scratch/vir/felipe/Repos/ICM_Code/AMG_Hunter.py --threads 24 --kegg_phylogeny True --kegg_ko K05371 --cds $PhyloDB --derep_id 0 --kegg_min_score 50
#all viral sequences cluster together and separated from cellular organisms, except from the cyanophages which cluster with their hosts. Only cyanobacteria have this gene, but viruses predicted to infect alpha, gamma, cyano, and bacteroidetes carry it. Related viruses are marine and non-marine saline alkaline




#Plot trees
# conda activate /mnt/smart/scratch/vir/felipe/envs/basic/
# Rscript /mnt/smart/scratch/vir/felipe/Repos/ICM_Code/Decorate_Trees.R Aligned_Merged_NR_Matched_K05371.faa.newick
# Rscript /mnt/smart/scratch/vir/felipe/Repos/ICM_Code/Decorate_Trees.R Aligned_Merged_NR_Matched_K00643.faa.newick
# Rscript /mnt/smart/scratch/vir/felipe/Repos/ICM_Code/Decorate_Trees.R Aligned_Merged_NR_Matched_K21480.faa.newick
# Rscript /mnt/smart/scratch/vir/felipe/Repos/ICM_Code/Decorate_Trees.R Aligned_Merged_NR_Matched_K01749.faa.newick
# Rscript /mnt/smart/scratch/vir/felipe/Repos/ICM_Code/Decorate_Trees.R Aligned_Merged_NR_Matched_K01599.faa.newick
# Rscript /mnt/smart/scratch/vir/felipe/Repos/ICM_Code/Decorate_Trees.R Aligned_Merged_NR_Matched_K00228.faa.newick
