#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=00-03:00:00
#SBATCH --job-name=Clean_WGS_Reads_Antarctic_Lagoons
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=15G
#SBATCH --output=/mnt/smart/users/fcoutinho/Job_Logs/jobLog_%J.out
#SBATCH --error=/mnt/smart/users/fcoutinho/Job_Logs/jobLog_%J.err
#SBATCH --mail-type=All
#SBATCH --mail-user=fhernandes@icm.csic.es

##############################################
#Run fastqc on raw reads
# cd /mnt/smart/users/fcoutinho/Antarctic_Lagoons/BMK_DATA_20241226153150_1/Data/rawdata/
# module load fastqc
# mkdir FastQC_Raw_Results
# mkdir fastqc_tmp

# # #Run FastQC to visualize the quality of the raw reads
# fastqc --quiet --threads 23 --format fastq --outdir FastQC_Raw_Results --dir fastqc_tmp *fq.gz

# # #Trim and Clean reads with cutadapt
# cd /mnt/smart/users/fcoutinho/Antarctic_Lagoons/BMK_DATA_20241226153150_1/Data/rawdata/
# module load cutadapt
# source activate
# # --discard-untrimmed  ommitted sicne these are not amplicons
# #for r1 in /mnt/smart/users/fcoutinho/OctoMicro/WGS_Raw_Reads/Unknown_BL327-0030001_1*gz
for r1 in /mnt/smart/users/fcoutinho/Antarctic_Lagoons/BMK_DATA_20241226153150_1/Data/rawdata/*_1.fq.gz; do r2=$(echo "$r1" | sed 's/_1.fq.gz/_2.fq.gz/') ; sid=$(echo "$r1" | sed 's/.*\///'| sed 's/_1.fq.gz//g'); echo "Processing together Sample: $sid R1:" $r1 "R2:" $r2 ;  cutadapt --cores 23 -b AGATCGGAAGAG -B AGATCGGAAGAG --quality-cutoff 30 --overlap 5 --cut 15 -U 15 --times 3 --minimum-length 100 --output /mnt/smart/users/fcoutinho/Antarctic_Lagoons/BMK_DATA_20241226153150_1/Data/Clean_Reads/Clean_${sid}_1.fq.gz --paired-output /mnt/smart/users/fcoutinho/Antarctic_Lagoons/BMK_DATA_20241226153150_1/Data/Clean_Reads/Clean_${sid}_2.fq.gz $r1 $r2 > Cutadapt_Report_$sid.log ; done

# deactivate

#Now re-run fastqc on the clean reads
cd /mnt/smart/users/fcoutinho/Antarctic_Lagoons/BMK_DATA_20241226153150_1/Data/Clean_Reads/
module load fastqc
mkdir FastQC_Clean_Results
mkdir fastqc_tmp

#fastqc --threads 23 --format fastq --outdir FastQC_Clean_Results --dir fastqc_tmp Clean_Unknown_BL327-0030001_2.fq.gz Clean_Unknown_BL327-0030001_1.fq.gz 
fastqc --quiet --threads 23 --format fastq --outdir FastQC_Clean_Results --dir fastqc_tmp *fq.gz


