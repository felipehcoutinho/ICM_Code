#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=01-12:00:00
#SBATCH --job-name=CCfinder
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=15G
#SBATCH --output=/mnt/smart/scratch/vir/felipe/Job_Logs/jobLog_%J.out
#SBATCH --error=/mnt/smart/scratch/vir/felipe/Job_Logs/jobLog_%J.err
#SBATCH --mail-type=All
#SBATCH --mail-user=fhernandes@icm.csic.es

##############################################
module load miniconda
source /home/apps/miniconda/4.5.x/etc/profile.d/conda.sh

cd /mnt/smart/scratch/vir/felipe/Databases/Malaspina/DeepMAGs/

conda activate /mnt/smart/scratch/vir/felipe/envs/cctyper

for genome in *fasta ; do cctyper --no_plot --threads 23 $genome cctyper_output_${genome%.fasta} ; done ; 

# python3 /mnt/smart/scratch/vir/felipe/Repos/ICM_Code/Merge_Matrixes.py --matrixes /mnt/smart/scratch/vir/felipe/Databases/Malaspina/DeepMAGs/cctyper_output_mp-deep_mag-05*/*putative.tab  --output CCfinder_DeepMAGs_output.tsv --out_index Contig

# head -n 1 /mnt/smart/scratch/vir/felipe/Databases/Malaspina/DeepMAGs/cctyper_output_mp-deep_mag-0002/cas_operons_putative.tab > CCfinder_DeepMAGs_output.tsv

# cat /mnt/smart/scratch/vir/felipe/Databases/Malaspina/DeepMAGs/cctyper_out*/CRISPR_Cas.tab >> CCfinder_DeepMAGs_CRISPR_Cas_output.tsv

