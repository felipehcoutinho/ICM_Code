#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=vir
#SBATCH --time=15-00:00:00
#SBATCH --job-name=OctoMicroxENA
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=50G
#SBATCH --output=/mnt/smart/scratch/vir/felipe/Job_Logs/%J.out
#SBATCH --error=/mnt/smart/scratch/vir/felipe/Job_Logs/%J.err
#SBATCH --mail-type=All
#SBATCH --mail-user=fhernandes@icm.csic.es

##############################################
cd /mnt/smart/scratch/vir/felipe/OctoMicro/MAGs/MAGs_Abundance/ENA_MGs/ 

module load bowtie2/2.3.4.1
module load samtools

GENOMEFILE="/mnt/smart/scratch/vir/felipe/Something.fasta"
DBPREFIX="/mnt/smart/scratch/vir/felipe/OctoMicro/Viruses/Virus_Abundances/Test/All_Genomic"
MGTABLE="/mnt/smart/scratch/vir/felipe/OctoMicro/MAGs/MAGs_Abundance/ENA_MGs/Selected_Filtered_ENA_Illumina_PE_Query_results_read_run.tsv"
THREADS=4

bowtie2-build --threads $THREADS $GENOMEFILE $DBPREFIX

while read line; do sra_id=$(echo "$line" | cut -f 1 ); r1=$(echo "$line" | cut -f 2 ) ; r2=$(echo "$line" | cut -f 3) ; echo "Processing ${sra_id} R1 ${r1} R2 ${r2}"; bowtie2  -x $DBPREFIX -q -1 $r1 -2 $r2 -S ${sra_id}.sam --very-sensitive-local --no-unal --threads $THREADS ; samtools view -bS ${sra_id}.sam > ${sra_id}.bam ; samtools sort ${sra_id}.bam -o ${sra_id}.sorted.bam ; samtools index ${sra_id}.sorted.bam ; samtools idxstats ${sra_id}.sorted.bam > ${sra_id}.Counts.tsv ; rm -f ${sra_id}.sam ${sra_id}.bam ${sra_id}.sorted.bam ${sra_id}.sorted.bam.bai ; done < $MGTABLE

#Get the raw and rpwk read counts tables
module load python/3.8.5
python3 /mnt/smart/scratch/vir/felipe/Repos/ICM_Code/RPKM_from_Counts.py




