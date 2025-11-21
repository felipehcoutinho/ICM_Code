#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=01-01:00:00
#SBATCH --job-name=Bins_gvSAGs
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=50G
#SBATCH --output=/mnt/smart/scratch/vir/felipe/Job_Logs/jobLog_%J.out
#SBATCH --error=/mnt/smart/scratch/vir/felipe/Job_Logs/jobLog_%J.err
#SBATCH --mail-type=All
#SBATCH --mail-user=fhernandes@icm.csic.es

##############################################
conda activate /mnt/smart/scratch/vir/felipe/envs/vRhyme

# module load prodigal
# module load bowtie2
# module load samtools

cd /mnt/smart/scratch/vir/felipe/gvSAGs/vRhyme_Binning

module load miniconda
source /home/apps/miniconda/4.5.x/etc/profile.d/conda.sh

conda activate /mnt/smart/scratch/vir/felipe/envs/vRhyme

vRhyme -i /mnt/smart/scratch/emm2/xlopezal/01_SAGs/02_mapping_virfa_vs_SAGs/data/gvSAGs.fasta -g /mnt/smart/scratch/emm2/xlopezal/01_SAGs/02_mapping_virfa_vs_SAGs/data/gvSAGs.fna -p /mnt/smart/scratch/emm2/xlopezal/01_SAGs/02_mapping_virfa_vs_SAGs/data/gvSAGs.faa -b /mnt/smart/scratch/emm2/xlopezal/01_SAGs/02_mapping_virfa_vs_SAGs/*.bam -t 24 -o vRhyme_Output_R1/

vRhyme -i /mnt/smart/scratch/emm2/xlopezal/01_SAGs/02_mapping_virfa_vs_SAGs/data/gvSAGs.fasta -g /mnt/smart/scratch/emm2/xlopezal/01_SAGs/02_mapping_virfa_vs_SAGs/data/gvSAGs.fna -p /mnt/smart/scratch/emm2/xlopezal/01_SAGs/02_mapping_virfa_vs_SAGs/data/gvSAGs.faa -b /mnt/smart/scratch/emm2/xlopezal/01_SAGs/02_mapping_virfa_vs_SAGs/*.bam -t 24 --method longest -o vRhyme_Output_R2/

cd /mnt/smart/scratch/vir/felipe/gvSAGs/vRhyme_Binning/vRhyme_Output_R2

module load python 
python3 /mnt/smart/scratch/vir/felipe/Repos/ICM_Code/Concatenate_Sequences.py /mnt/smart/scratch/vir/felipe/gvSAGs/vRhyme_Binning/vRhyme_Output_R2/vRhyme_best_bins_fasta/*fasta

# Outputs:
# /mnt/smart/scratch/vir/felipe/gvSAGs/vRhyme_Binning/vRhyme_Output_R2/gvSAGs_Concatenated_Bins_vRhyme_R2.fasta
# /mnt/smart/scratch/vir/felipe/gvSAGs/vRhyme_Binning/vRhyme_Output_R2/Concatenated_Sequence_Info.tsv

###Run CheckV (CESGA)
module load cesga/2020 gcccore/system checkv/1.0.1
cd /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/gvSAGs/
checkv end_to_end gvSAGs_Concatenated_Bins_vRhyme_R2.fasta CheckV_Output -t 4

# Output moved to:
# /mnt/smart/scratch/vir/felipe/gvSAGs/vRhyme_Binning/vRhyme_Output_R2/CheckV_Output/