#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=15-00:00:00
#SBATCH --job-name=OMDB_CCfinder
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=15G
#SBATCH --output=/mnt/smart/scratch/vir/felipe/Job_Logs/jobLog_%J.out
#SBATCH --error=/mnt/smart/scratch/vir/felipe/Job_Logs/jobLog_%J.err
#SBATCH --mail-type=All
#SBATCH --mail-user=fhernandes@icm.csic.es

##############################################
cd /mnt/smart/scratch/vir/felipe/Databases/OMDB/Individual_Genomes/

#split the sequences into individual files
module load python

#the equal sign and qutes are necessary so the arparser understands -scaffold_ is to be assigned to id_split_sep
python3 /mnt/smart/scratch/vir/felipe/Repos/ICM_Code/Get_Seqs.py --input_seq /mnt/smart/scratch/vir/felipe/Databases/OMDB/Renamed_OMDBv2.0_SC_G_R.fa --fetch_all True --id_split_pos 0 --group True --id_split_sep="-scaffold_"

# run cctyper, delete genome files as they are generated to save space
module load miniconda
source /home/apps/miniconda/4.5.x/etc/profile.d/conda.sh

conda activate /mnt/smart/scratch/vir/felipe/envs/cctyper

for genome in *fasta ; do cctyper --no_plot --threads 23 $genome cctyper_output_${genome%.fasta} ; rm -f $genome ;  done ; 

# python3 /mnt/smart/scratch/vir/felipe/Repos/ICM_Code/Merge_Matrixes.py --matrixes /mnt/smart/scratch/vir/felipe/Databases/Malaspina/DeepMAGs/cctyper_output_mp-deep_mag-05*/*putative.tab  --output CCfinder_DeepMAGs_output.tsv --out_index Contig

# head -n 1 /mnt/smart/scratch/vir/felipe/Databases/Malaspina/DeepMAGs/cctyper_output_mp-deep_mag-0002/cas_operons_putative.tab > CCfinder_DeepMAGs_output.tsv

# cat /mnt/smart/scratch/vir/felipe/Databases/Malaspina/DeepMAGs/cctyper_out*/CRISPR_Cas.tab >> CCfinder_DeepMAGs_CRISPR_Cas_output.tsv

