#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO

#SBATCH --time=00-04:00:00
#SBATCH --job-name=Antarctic_Lagoon_Clean_Reads
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=15G
#SBATCH --output=/mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/Job_Logs/%J.out
#SBATCH --error=/mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/Job_Logs/%J.err
#SBATCH --mail-type=All
#SBATCH --mail-user=fhernandes@icm.csic.es

##############################################
module load cutadapt fastqc

cd /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/Antarctic_Lagoon/Reads/Raw/

for r1 in *R1_001.fastq.gz; do r2=$(echo "$r1" | sed 's/_R1_/_R2_/') ; sid=$(echo "$r1" | sed 's/_S1_.*//g'); echo "Processing together Sample $sid R1:" $r1 "R2:" $r2 ;  grep -w -A 1 $sid  /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/Antarctic_Lagoon/Reads/indexes_merged.fasta > Sample_$sid.Sequences-to_Remove.fasta ; cat /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/Antarctic_Lagoon/Reads/primers_and_adapters.fasta >> Sample_$sid.Sequences-to_Remove.fasta ; cutadapt --cores 23 -g "file:Sample_$sid.Sequences-to_Remove.fasta" -G "file:Sample_$sid.Sequences-to_Remove.fasta" -a "file:Sample_$sid.Sequences-to_Remove.fasta" -A "file:Sample_$sid.Sequences-to_Remove.fasta" --times 3 --minimum-length 200  --output /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/Antarctic_Lagoon/Reads/Clean/Clean_$r1 --paired-output /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/Antarctic_Lagoon/Reads/Clean/Clean_$r2 $r1 $r2 > Cutadapt_Report_$r1.log ; done

cd /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/Antarctic_Lagoon/Reads/Clean/
mkdir Cutadapt_Reports
mv /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/Antarctic_Lagoon/Reads/Raw/*log Cutadapt_Reports/
mv /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/Antarctic_Lagoon/Reads/Raw/*Sequences-to_Remove.fasta Cutadapt_Reports/

#mkdir FastQC_Clean_Results
#mkdir fastqc_tmp
fastqc --quiet --threads 23 --format fastq --outdir FastQC_Clean_Results --dir fastqc_tmp Clean*fastq.gz

#Remove blank sample
rm -f Clean_BPCR2_*fastq.gz