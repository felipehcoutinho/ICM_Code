#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=03-00:00:00
#SBATCH --job-name=Abundance_Profiles_Large_Fraction_MLProfiles_dRep_MAGs
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=150G
#SBATCH --output=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%A_%a.out
#SBATCH --error=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%A_%a.err
#SBATCH --array=1-6%3
#SBATCH --mail-type=ALL
#SBATCH --mail-user=fhernandes@icm.csic.es

##############################################
module load python/3.8.5
module load bowtie2/2.3.4.1
module load samtools

CHUNK_FILE=/mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/Abundance_MLProfiles_Large_Fraction/Chunk_${SLURM_ARRAY_TASK_ID}_Malaspina_Profiles_ALL_MG_Large_Fraction_3-20_Info.tsv

python3 /mnt/lustre/bio/users/fcoutinho/Repos/virathon/Virathon.py --bowtie_k 1 --genome /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/Abundance_dRep_MAGs/Merged_Sequences.fasta --bowtiedb /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/Abundance_dRep_MAGs/Merged_Sequences --abundance_table True --raw_read_table $CHUNK_FILE --threads 47

#python3 /mnt/lustre/bio/users/fcoutinho/Repos/virathon/Virathon.py --bowtie_k 1 --genome /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/Abundance_dRep_MAGs/Merged_Sequences.fasta --bowtiedb /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/Abundance_dRep_MAGs/Merged_Sequences --abundance_table True --raw_read_table /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Metadata/Kaust_Malaspina_Profiles_MG_Large_Fraction_3-20_Info.tsv --threads 47

#python3 /mnt/lustre/bio/users/fcoutinho/Repos/virathon/Virathon.py --genome /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/Abundance_dRep_MAGs/Merged_Sequences.fasta --bowtiedb /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/Abundance_dRep_MAGs/Merged_Sequences --abundance_table True --raw_read_table /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/Abundance_MLProfiles_Large_Fraction/Malaspina_Profiles+3_Unmatched_PA_Metagenome_Info.tsv --abundance_rpkm True --parse_only True

#python3 /mnt/lustre/bio/users/fcoutinho/Repos/virathon/Virathon.py --genome /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/Abundance_dRep_MAGs/Merged_Sequences.fasta --bowtiedb /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/Abundance_dRep_MAGs/Merged_Sequences --abundance_table True --raw_read_table /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Metadata/Kaust_Malaspina_Profiles_MG_Large_Fraction_3-20_Info.tsv --abundance_rpkm True  --threads 47

#python3 /mnt/lustre/bio/users/fcoutinho/Repos/virathon/Virathon.py --genome /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/Abundance_dRep_MAGs/Merged_Sequences.fasta --bowtiedb /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/Abundance_dRep_MAGs/Merged_Sequences --abundance_table True --raw_read_table /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/Abundance_MLProfiles_Large_Fraction/Malaspina_Profiles+3_Unmatched_PA_Metagenome_Info.tsv --abundance_rpkm True --parse_only True

#python3 /mnt/lustre/bio/users/fcoutinho/Repos/virathon/Virathon.py --genome /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/Abundance_dRep_MAGs/Merged_Sequences.fasta --bowtiedb /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/Abundance_dRep_MAGs/Merged_Sequences --abundance_table True --raw_read_table /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/Abundance_MLProfiles_Large_Fraction/MP0320_Info.tsv --threads 24

#python3 /mnt/lustre/bio/users/fcoutinho/Repos/virathon/Virathon.py --genome /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/Abundance_dRep_MAGs/Merged_Sequences.fasta --bowtiedb /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/Abundance_dRep_MAGs/Merged_Sequences --abundance_table True --raw_read_table /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/Abundance_MLProfiles_Large_Fraction/Malaspina_Profiles+3_Unmatched_PA_Metagenome_Info.tsv --threads 24

#python3 /mnt/lustre/bio/users/fcoutinho/Scripts/Virathon.py --genome /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/Abundance_dRep_MAGs/Merged_Sequences.fasta --bowtiedb /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/Abundance_dRep_MAGs/Merged_Sequences --abundance_table True --abundance_rpkm True --raw_read_table /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/Abundance_MLProfiles_Large_Fraction/Missing_Sample.tsv --threads 45

#python3 /mnt/lustre/bio/users/fcoutinho/Scripts/Virathon.py --genome /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/Abundance_dRep_MAGs/Merged_Sequences.fasta --bowtiedb /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/Abundance_dRep_MAGs/Merged_Sequences --abundance_table True --abundance_rpkm True --raw_read_table /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/Abundance_MLProfiles_Large_Fraction/Malaspina_Profiles+3_Unmatched_PA_Metagenome_Info.tsv --threads 45 --parse_only True
