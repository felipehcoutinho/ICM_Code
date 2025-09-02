#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=01-01:00:00
#SBATCH --job-name=Phylo_gvSAGs
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=50G
#SBATCH --output=/mnt/smart/scratch/vir/felipe/Job_Logs/jobLog_%J.out
#SBATCH --error=/mnt/smart/scratch/vir/felipe/Job_Logs/jobLog_%J.err
#SBATCH --mail-type=All
#SBATCH --mail-user=fhernandes@icm.csic.es

##############################################
module load python
module load hmmer
module load muscle
module load fasttree


###Setting up the database
#cd /mnt/smart/scratch/vir/felipe/gvSAGs/Phylogenies/
# #make the fasta file with all protein sequences
# python3 /mnt/smart/scratch/vir/felipe/Repos/ICM_Code/Get_Seqs.py --fetch_all True --input /mnt/smart/scratch/vir/felipe/gvSAGs/gvSAGs.faa /mnt/smart/scratch/vir/felipe/Databases/GVDB/GVDB_genomes/GVDB_prot/*faa --matched_out gvSAG_PhyloDB.faa

# #Remove characters from Seq IDs that will mess up analysis
# sed 's/:/_/g' /mnt/smart/scratch/vir/felipe/gvSAGs/Phylogenies/gvSAG_PhyloDB.faa > /mnt/smart/scratch/vir/felipe/gvSAGs/Phylogenies/Clean_Names_gvSAG_PhyloDB.faa

###To properly look for best hist need to do hmmpress combining # /mnt/smart/scratch/vir/felipe/Databases/GVDB/hmm/NCLDV_markers.hmm with at least one of (check which and if they overlap with markers):
# /mnt/smart/scratch/vir/felipe/Databases/GVDB/hmm/vogdb.hmm
# /mnt/smart/scratch/vir/felipe/Databases/GVDB/hmm/gvog.hmm

###Use amg hunter to identify markers and build individual alignments
#Markers in /mnt/smart/scratch/vir/felipe/Databases/GVDB/hmm/NCLDV_markers.hmm: A32 D5 mcp mRNAc PolB RNAPL RNAPS RNR SFII VLTF3
#Markers sued by Xabi: D5 RNAPL RNAPS RNR SFII VLTF3 mRNAc mcp polB
cd /mnt/smart/scratch/vir/felipe/gvSAGs/Concat_Phy/Round_3/

###REMEMBER to delete previous hmm_phylo_seqs.faa files before running, otherwise the sequences from the last run will be appended to the files
rm -f *hmm_phylo_seqs.faa ; python3 /mnt/smart/scratch/vir/felipe/Repos/ICM_Code/AMG_Hunter.py --threads 24 --hmm_phylogeny True --db_hmm_phylogeny /mnt/smart/scratch/vir/felipe/Databases/GVDB/hmm/NCLDV_markers.hmm --cds /mnt/smart/scratch/vir/felipe/gvSAGs/Phylogenies/Clean_Names_gvSAG_PhyloDB.faa --phylo_min_score 100 --derep_id 0 --hmm_ids A32 D5 mcp mRNAc PolB RNAPL RNAPS RNR SFII VLTF3 --phylo_single_cds True --clean_gvdb True

###Concatenate the alignments
python3 /mnt/smart/scratch/vir/felipe/Repos/ICM_Code/Concatenate_Alignments.py --aln_files *hmm_phylo_seqs_aln.faa --fasta_output gvSAGs_Concatenated_Alingment_min_2_matches.faa --clean_gvdb True --min_matches 2

###Build tree
FastTreeMP -nosupport -fastest gvSAGs_Concatenated_Alingment_min_2_matches.faa > gvSAGs_Concatenated_Alingment_min_2_matches.newick

###Decorate the tree with info from GVDB and print to pdf:
# /mnt/smart/scratch/vir/felipe/gvSAGs/Concat_Phy/Round_3/gvSAGs_Concatenated_Alingment_min_2_matches_Decorated_Tree_Rectangular.pdf
# /mnt/smart/scratch/vir/felipe/gvSAGs/Concat_Phy/Round_3/gvSAGs_Concatenated_Alingment_min_2_matches_Unrooted_Decorated_Tree.pdf

module load miniconda
source /home/apps/miniconda/4.5.x/etc/profile.d/conda.sh
conda activate /mnt/smart/scratch/vir/felipe/envs/basic/

Rscript /mnt/smart/scratch/vir/felipe/gvSAGs/Concat_Phy/gvSAGs_Decorate_Trees.R gvSAGs_Concatenated_Alingment_min_2_matches.newick /mnt/smart/scratch/vir/felipe/Databases/GVDB/GVDB_genomes/GVDB_genome_descriptions.tsv


###Ignore commands below
# ###Concatenate alignments
# ###NOTE: when cocatenating aligments, some GVDB genomes are represented by multiple contigs, e.g. ERX552255.28.dc.fa|contig_2127 
# # As it is the script currently considers ERX552255.28.dc.fa|contig_2127 ERX552255.28.dc.fa|contig_4798 as separate genomes and they show separately in the tree
# # GCA_000837185.1_ViralProj14097|AF250284.1_239 should be assigned to GCA.000837185.1.ViralProj14097 in the genome info table
# #No issue with the SAGs sequences as their naming is: AH2982_RM2_B7_NODE_3_length_27938_cov_591.035885_1 and there is no binning genoems ATM
# # GVMAG-S-1091232-186
# # GCA.003233895.1.ASM323389v1 > GCA_003233895.1_ASM323389v1|MG011689.1_51
# # medusavirus > medusavirus.fna
# cd /mnt/smart/scratch/vir/felipe/gvSAGs/Concat_Phy/
# module load python

# python3 /mnt/smart/scratch/vir/felipe/Repos/ICM_Code/Concatenate_Alignments.py --aln_files /mnt/smart/scratch/emm2/xlopezal/01_SAGs/05_marker_filogenies/*_phylogeny/*_modified.faa --fasta_output gvSAGs_Concatenated_Alingment_min_2_matches.faa --clean_gvdb True --min_matches 2 > concat_min2.log


# ###Testing for single matches vs markers only
# module load python
# module load hmmer
# module load muscle
# module load fasttree

# cd /mnt/smart/scratch/vir/felipe/gvSAGs/Concat_Phy/Test_set/

# python3 /mnt/smart/scratch/vir/felipe/Repos/ICM_Code/Get_Seqs.py --list test_seq_ids.txt --input /mnt/smart/scratch/vir/felipe/gvSAGs/Phylogenies/Clean_Names_gvSAG_PhyloDB.faa --protein True --matched_out Test_Seqs.faa

# ###REMEMBER to delete previous hmm_phylo_seqs.faa files before running, otherwise the sequences from the last run will be appended to the files
# rm -f *hmm_phylo_seqs.faa ; python3 /mnt/smart/scratch/vir/felipe/Repos/ICM_Code/AMG_Hunter.py --threads 24 --hmm_phylogeny True --db_hmm_phylogeny /mnt/smart/scratch/vir/felipe/Databases/GVDB/hmm/NCLDV_markers.hmm --cds Test_Seqs.faa --phylo_min_score 100 --derep_id 0 --hmm_ids mcp mRNAc PolB --phylo_single_cds True --clean_gvdb True

# python3 /mnt/smart/scratch/vir/felipe/Repos/ICM_Code/Concatenate_Alignments.py --aln_files *_hmm_phylo_seqs_aln.faa --fasta_output Test_gvSAGs_Concatenated_Alingment_min_2_matches.faa --clean_gvdb True --min_matches 2

# FastTreeMP -nosupport -fastest Test_gvSAGs_Concatenated_Alingment_min_2_matches.faa > Test_gvSAGs_Concatenated_Alingment_min_2_matches.newick

# conda activate basic

# Rscript /mnt/smart/scratch/vir/felipe/gvSAGs/Concat_Phy/gvSAGs_Decorate_Trees.R Test_gvSAGs_Concatenated_Alingment_min_2_matches.newick /mnt/smart/scratch/vir/felipe/Databases/GVDB/GVDB_genomes/GVDB_genome_descriptions.tsv

# #Build tree from concatenated alignment
# module load fasttree

# FastTreeMP -nosupport -fastest gvSAGs_Concatenated_Alingment_min_2_matches.faa > gvSAGs_Concatenated_Alingment_min_2_matches.newick

# #Decorate the tree with info from GVDB
# module load miniconda
# source /home/apps/miniconda/4.5.x/etc/profile.d/conda.sh
# conda activate basic

# Rscript /mnt/smart/scratch/vir/felipe/gvSAGs/Concat_Phy/gvSAGs_Decorate_Trees.R gvSAGs_Concatenated_Alingment_min_2_matches.newick /mnt/smart/scratch/vir/felipe/Databases/GVDB/GVDB_genomes/GVDB_genome_descriptions.tsv

