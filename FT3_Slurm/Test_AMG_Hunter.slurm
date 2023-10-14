#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO

#SBATCH --time=00-01:00:00
#SBATCH --job-name=TestAMGHunter
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=60G
#SBATCH --output=/mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/Job_Logs/%J.out
#SBATCH --error=/mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/Job_Logs/%J.err
#SBATCH --mail-type=All
#SBATCH --mail-user=fhernandes@icm.csic.es

##############################################
module load cesga/system miniconda3/22.11.1-1
conda activate /mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/Conda_Envs/basic

cd /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/Malaspina/Profiles/Prok_Viruses/
python3 /mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/repos/ICM_Code/AMG_Hunter.py --cds /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/Malaspina/Profiles/Prok_Viruses/Profiles_Malaspina_Viruses_CheckV_Trimmed_Genomes.faa --info_cds_output CDS_Info_Profiles_Malaspina_Viruses_CheckV_Trimmed.tsv --info_genome_output Scaffold_Annotation_Info_Profiles_Malaspina_Viruses_CheckV_Trimmed.tsv --annotate True --parse_only True --skip_uniref True

#Command dies due to memory error, even when using very reduced kegg and pfamoubput files. Possibly error during the pfam indexing, as skipin this stage leads to running succesfuly
