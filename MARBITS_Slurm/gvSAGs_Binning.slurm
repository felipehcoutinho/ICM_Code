#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=vir
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
# module load prodigal
# module load bowtie2
# module load samtools

cd /mnt/smart/scratch/vir/felipe/gvSAGs/vRhyme_Binning

module load miniconda
source /home/apps/miniconda/4.5.x/etc/profile.d/conda.sh

conda activate /mnt/smart/scratch/vir/felipe/envs/vRhyme

#with bams, default params
vRhyme -i /mnt/smart/scratch/emm2/xlopezal/01_SAGs/02_mapping_virfa_vs_SAGs/data/gvSAGs.fasta -g /mnt/smart/scratch/emm2/xlopezal/01_SAGs/02_mapping_virfa_vs_SAGs/data/gvSAGs.fna -p /mnt/smart/scratch/emm2/xlopezal/01_SAGs/02_mapping_virfa_vs_SAGs/data/gvSAGs.faa -b /mnt/smart/scratch/emm2/xlopezal/01_SAGs/02_mapping_virfa_vs_SAGs/*.bam -t 24 -o vRhyme_Output_R1/

#with bams, longest method
vRhyme -i /mnt/smart/scratch/emm2/xlopezal/01_SAGs/02_mapping_virfa_vs_SAGs/data/gvSAGs.fasta -g /mnt/smart/scratch/emm2/xlopezal/01_SAGs/02_mapping_virfa_vs_SAGs/data/gvSAGs.fna -p /mnt/smart/scratch/emm2/xlopezal/01_SAGs/02_mapping_virfa_vs_SAGs/data/gvSAGs.faa -b /mnt/smart/scratch/emm2/xlopezal/01_SAGs/02_mapping_virfa_vs_SAGs/*.bam -t 24 --method longest -o vRhyme_Output_R2/

#No bams, default params (does not work s vRhyme requires at least1 sample)
vRhyme -i /mnt/smart/scratch/emm2/xlopezal/01_SAGs/02_mapping_virfa_vs_SAGs/data/gvSAGs.fasta -g /mnt/smart/scratch/emm2/xlopezal/01_SAGs/02_mapping_virfa_vs_SAGs/data/gvSAGs.fna -p /mnt/smart/scratch/emm2/xlopezal/01_SAGs/02_mapping_virfa_vs_SAGs/data/gvSAGs.faa -t 24 -o vRhyme_Output_R3/

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


###metabat binning 
module load metabat
module load miniconda
source /home/apps/miniconda/4.5.x/etc/profile.d/conda.sh
conda activate metabat-2.15

cd /mnt/smart/scratch/vir/felipe/gvSAGs/MetaBat_Binning/

metabat --numThreads 8 --in /mnt/smart/scratch/emm2/xlopezal/01_SAGs/02_mapping_virfa_vs_SAGs/data/gvSAGs.fasta --minClsSize 100000 --minContig 10000 --minS 7 --maxEdges 75  --outFile gvSAGs_Aylward_MetaBins_R1 --saveCls

mv gvSAGs_Aylward_MetaBins_R1 Clusters_gvSAGs_Aylward_MetaBins_R1.tsv
mkdir R1_Bins
mv gvSAGs_Aylward_MetaBins_R1*fa R1_Bins/

#do the same usign the read mapping data
jgi_summarize_bam_contig_depths --outputDepth /mnt/smart/scratch/vir/felipe/gvSAGs/MetaBat_Binning/gvSAGS_jgi_depths_from_02_mapping_virfa_vs_SAGs.tsv /mnt/smart/scratch/emm2/xlopezal/01_SAGs/02_mapping_virfa_vs_SAGs/*.bam

metabat --abdFile /mnt/smart/scratch/vir/felipe/gvSAGs/MetaBat_Binning/gvSAGS_jgi_depths_from_02_mapping_virfa_vs_SAGs.tsv --numThreads 8 --in /mnt/smart/scratch/emm2/xlopezal/01_SAGs/02_mapping_virfa_vs_SAGs/data/gvSAGs.fasta --minClsSize 100000 --minContig 10000 --minS 7 --maxEdges 75  --outFile gvSAGs_Aylward_MetaBins_R2 --saveCls

mkdir R2_Bins
mv gvSAGs_Aylward_MetaBins_R2*fa R2_Bins/

#Louvain 67 (11 contigs):
grep -w "P4_62_K127_NODE_2_length_95293_cov_32.001923\|P3_51_K127_NODE_202_length_18609_cov_14.179634\|P4_35_K127_NODE_170_length_11938_cov_1.367454\|P4_68_K127_NODE_12_length_46384_cov_20.423806\|P3_40_K127_NODE_59_length_38315_cov_18.313056\|P3_66_K127_NODE_66_length_29713_cov_3.460353\|P4_31_K127_NODE_37_length_52353_cov_8.660476\|P4_69_K127_NODE_23_length_50714_cov_40.712021\|P4_49_K127_NODE_136_length_17678_cov_5.482878\|P4_39_K127_NODE_28_length_41147_cov_6.849196\|P4_02_K127_NODE_268_length_7529_cov_2.158741" R1_Bins/*fa

#returns R1_Bins/gvSAGs_Aylward_MetaBins_R1.24.fa (10 contigs)

#same search in R2 bins returns nothing. uneven contig coverages from MDA make it unsuitable for binning with metabat?

#vrheym appends string vRhyme_#__ to seq ids
grep "P4_62_K127_NODE_2_length_95293_cov_32.001923\|P3_51_K127_NODE_202_length_18609_cov_14.179634\|P4_35_K127_NODE_170_length_11938_cov_1.367454\|P4_68_K127_NODE_12_length_46384_cov_20.423806\|P3_40_K127_NODE_59_length_38315_cov_18.313056\|P3_66_K127_NODE_66_length_29713_cov_3.460353\|P4_31_K127_NODE_37_length_52353_cov_8.660476\|P4_69_K127_NODE_23_length_50714_cov_40.712021\|P4_49_K127_NODE_136_length_17678_cov_5.482878\|P4_39_K127_NODE_28_length_41147_cov_6.849196\|P4_02_K127_NODE_268_length_7529_cov_2.158741" /mnt/smart/scratch/vir/felipe/gvSAGs/vRhyme_Binning/vRhyme_Output_R1/vRhyme_best_bins_fasta/*.fasta

#only a single sequence binned: /mnt/smart/scratch/vir/felipe/gvSAGs/vRhyme_Binning/vRhyme_Output_R2/vRhyme_best_bins_fasta/vRhyme_bin_9.fasta:>vRhyme_9__P4_62_K127_NODE_2_length_95293_cov_32.001923

grep "P4_62_K127_NODE_2_length_95293_cov_32.001923\|P3_51_K127_NODE_202_length_18609_cov_14.179634\|P4_35_K127_NODE_170_length_11938_cov_1.367454\|P4_68_K127_NODE_12_length_46384_cov_20.423806\|P3_40_K127_NODE_59_length_38315_cov_18.313056\|P3_66_K127_NODE_66_length_29713_cov_3.460353\|P4_31_K127_NODE_37_length_52353_cov_8.660476\|P4_69_K127_NODE_23_length_50714_cov_40.712021\|P4_49_K127_NODE_136_length_17678_cov_5.482878\|P4_39_K127_NODE_28_length_41147_cov_6.849196\|P4_02_K127_NODE_268_length_7529_cov_2.158741" /mnt/smart/scratch/vir/felipe/gvSAGs/vRhyme_Binning/vRhyme_Output_R2/vRhyme_best_bins_fasta/*.fasta

#only a single sequence binned: /mnt/smart/scratch/vir/felipe/gvSAGs/vRhyme_Binning/vRhyme_Output_R1/vRhyme_best_bins_fasta/vRhyme_bin_2.fasta:>vRhyme_2__P4_62_K127_NODE_2_length_95293_cov_32.001923


