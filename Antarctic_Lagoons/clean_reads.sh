#Generate QC report
module load biopython fastqc

cd /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/Antarctic_Lagoon/Reads/Raw/

mkdir FastQC_Raw_Results
mkdir fastqc_tmp

#Run fastQC for all files with the extension fastq.gz in the current directory
fastqc --quiet --threads 3 --format fastq --outdir FastQC_Raw_Results --dir fastqc_tmp *fastq.gz &

#Clean fastqc files with cuadapt
module load cutadapt
#--adapter
for r1 in *R1_001.fastq.gz; do r2=$(echo "$r1" | sed 's/_R1_/_R2_/') ; echo "Processing together R1:" $r1 "R2:" $r2 ; cutadapt --cores 3 -g GTGTCAGC -G CCGTCAAT --times 10 --minimum-length 200 --max-expected-errors 3 --output /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/Antarctic_Lagoon/Reads/Clean/Clean_$r1 --paired-output /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/Antarctic_Lagoon/Reads/Clean/Clean_$r2 $r1 $r2 ; done

#Run fastq again on the clean reads
cd /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/Antarctic_Lagoon/Reads/Clean/
mkdir FastQC_Clean_Results
mkdir fastqc_tmp
fastqc --quiet --threads 3 --format fastq --outdir FastQC_Clean_Results --dir fastqc_tmp Clean*fastq.gz &