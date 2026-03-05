# Read mapping and abundance quantification

## Introduction

Read mapping is often used to calculate the relative abundances of sequences in a sample. In the context of metaomics, read mapping is necessary because the relative abundances of individual sequences will vary dramatially within and among samples. Also, from read mapping, the taxon/function abundance matrixes can be obtained and used for standard ecological analysis of diversity. In this case, only DNA sequences can be used as reference. Unlike alignment tools that work with protein references that usually can identify very low identity matches (< 30%), bowtie is used for high identity mapping at the nucleotide level (> 97%). Thus, this methodology is only suitable for quantifying reference sequences very similar to those in your short reads dataset. 

## Software

- [Bowtie2](https://github.com/BenLangmead/bowtie2)
- [SAMTools](https://github.com/samtools/samtools)
- [CoverM](https://github.com/wwood/CoverM)

## Commands

We first build a reference database that contains the long DNA sequences (e.g. genome fragments, individual genes, multi sequence genomes, etc). 

`module load bowtie2/2.3.4.1`

Below we build a Database named DB_Genomes from the Genomes.fasta file using 4 threads. This will generate multiple files with the .bt2 extension

`bowtie2-build --threads 4 Genomes.fasta DB_Genomes`

Now the short reads are mapped to the referene database:

`bowtie2 -x DB_Genomes -q -1 SampleA_R1.fastq -2 SampleA_R2.fastq -S SampleAxDB_Genomes.sam -k 1 -u 100000 --no-discordant --no-mixed --no-unal --threads 4`

- NOTE: bowtie will also accept gzip compressed input files (with the gz extension) passed to the -1 and -2 parameters

- Often reads can map to multiple sequences/regions in the same reference file. Depending on the intend use of the results, we can limit that to a maximum number of mapped reads by setting the value of the parameter -k
- We can limit the maximum number of short reads that bowtie will try to map by setting the value of the parameter -u. This is useuful for speeding up the analysis when you have a lot of input reads or to standardize the total reads mapped when using mutiple samples
- --no-unal: the output will not list the unaligned reads (useful for making the output smaller and easier to parse)
- --no-discordant: only consider as valid when both r1 and r2 reads match the sequence
- --no-mixed: do not allow r1 and r2 to match different sequences

After the alignment completes, bowtie will output a summary with basic stats:

    10000000 reads; of these:
    10000000 (100.00%) were paired; of these:
        9999950 (100.00%) aligned concordantly 0 times
        50 (0.00%) aligned concordantly exactly 1 time
        0 (0.00%) aligned concordantly >1 times        
    0.00% overall alignment rate

The generated .sam file (human-readable) needs to be converted to .bam (binary), sorted, and indexed. These steps are carried out with samtools. The contents of the sam file can be visualized. Explanation of each field can be found [here](https://samtools.github.io/hts-specs/SAMv1.pdf)

`module load samtools/1.8`

`samtools view -bS SampleAxDB_Genomes.sam > SampleAxDB_Genomes.bam`

`samtools sort SampleAxDB_Genomes.bam -o SampleAxDB_Genomes.sorted.bam`

`samtools index SampleAxDB_Genomes.sorted.bam`

Finally, we obtain a table with the number of INDIVIDUAL reads mapped to each sequence with idxstats

`samtools idxstats SampleAxDB_Genomes.sorted.bam > SampleAxDB_Genomes.Counts.tsv`

The first three columns of SampleAxDB_Genomes.Counts.tsv file are: Reference Sequence ID,  Reference Sequence length, and number of INDIVIDUAL reads mapped (in this case the sum of R1 and R2 reads)

if we have multiple counts files generated from mapping reads from multiple samples to the same reference database, we can put the number of PAIRED reads mapped together in a single table:

`module load python/3.8.5`

`python3 /mnt/smart/scratch/vir/felipe/Repos/ICM_Code/RPKM_from_Counts.py`

- The Raw_Abundance.tsv contains the number of PAIRED reads mapped to each sequence
- The RPKM.tsv contains the number of read  PAIRED reads mapped per Kilo base pair in the sequence, per each 1 million MAPPED pairs

## CoverM
`module load coverM/0.7.0`

use CoverM to count:

- Total numbr of mapped reads

`coverm contig --methods count --bam-files SampleAxDB_Genomes.bam --output-file CoverM_Counts_Report.tsv`

- Average number of aligned reads overlapping each position on the contig:

`coverm contig --methods mean --bam-files SampleAxDB_Genomes.bam --output-file CoverM_Mean_Report.tsv`

- Number of bases covered by 1 or more reads:

`coverm contig --methods covered_bases --bam-files SampleAxDB_Genomes.bam --output-file CoverM_Covered_Bases_Report.tsv`


## References

- [Bowtie2 Publication](http://www.nature.com/nmeth/journal/v9/n4/full/nmeth.1923.html)
- [Samtools Publication](https://doi.org/10.1093/gigascience/giab008)
- [CoverM Publication](https://doi.org/10.1093/bioinformatics/btaf147)


