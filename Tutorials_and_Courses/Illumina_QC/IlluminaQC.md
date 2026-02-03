# Illumina sequencing reads QC
## Summary
***

This tutorial will teach you how to perfomr quality control on Illumina paired-ends sequencing reads. 

***

## Introduction

[Illumina sequencing video](https://www.youtube.com/watch?v=fCd6B5HRaZ8)

- Illumina paired-end sequencing will generate two sqeuencing reads for every DNA molecule that was sucessfully sequenced. These reads will be split between two files (R1 and R2) that are analysed together for read-cleaning and most of the subsequent analysis steps.


## Software

- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [cutadapt](https://cutadapt.readthedocs.io/en/stable/)


## Commands
- Change directory to the location of your fastq files

### Run FastQC to visualize the quality of all the fastq files containing the raw reads
- In MARBITS, we will need to load the fastqc module to be able to run the software

`module load fastqc`

`mkdir FastQC_Raw_Results`

`fastqc --quiet --threads 4 --format fastq --outdir FastQC_Raw_Results --dir fastqc_tmp *fq.gz`

- The FastQC_Raw_Results file will contain an html report for each input fastq (i.e. each sample will generate two files)

- [Example FastQC Output from our lab](https://github.com/felipehcoutinho/ICM_Code/blob/main/Tutorials_and_Courses/Data/Unknown_BL327-0030008_1_fastqc.html)

- [FastQC Example of good illumina data](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/good_sequence_short_fastqc.html)

- [FastQC Example of bad illumina data](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/bad_sequence_fastqc.html)

### Trim and clean a sample paired reads with cutadapt
- In MARBITS, we will need to load the cutadapt module AND load its conda enviroment to be able to run the software

`module load cutadapt`

`source activate`

`for r1 *_1.fq.gz; do r2=$(echo "$r1" | sed 's/_1.fq.gz/_2.fq.gz/') ; sid=$(echo "$r1" | sed 's/.*\///'| sed 's/_1.fq.gz//g'); echo "Processing together Sample: $sid R1:" $r1 "R2:" $r2 ; cutadapt --cores 23 -b AGATCGGAAGAG -B AGATCGGAAGAG --quality-cutoff 30 --overlap 5 --cut 15 -U 15 --times 3 --minimum-length 100 --output Clean_${sid}_1.fq.gz --paired-output Clean_${sid}_2.fq.gz $r1 $r2 > Cutadapt_Report_$sid.log ; done `

`deactivate`

### Now re-run fastqc on the clean reads to determine if the cleaning was efficient 

`mkdir FastQC_Clean_Results`

`fastqc --quiet --threads 4 --format fastq --outdir FastQC_Clean_Results --dir fastqc_tmp Clean*fq.gz`


## References

- [FastQC]()
- [cudadapt](https://journal.embnet.org/index.php/embnetjournal/article/view/200)
