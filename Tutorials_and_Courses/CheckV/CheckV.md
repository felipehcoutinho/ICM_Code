## Introduction

Tools to detect viral sequences often produce errors, either by incorrectely classifying cellular sequences as viral, or byincorrectely idenfiying the limits of integrated viral genomes. This lead to viral sequences that contain cellular contamination. Also, it is necessary to estimate to completeness of the viral genomes obtained from metagenomic datasets as a forma of quality control.

## Software

- [CheckV](https://bitbucket.org/berkeleylab/checkv/src/master/)

## Commands

Activate the environment

`conda activate /mnt/smart/scratch/vir/felipe/envs/checkv`

Run checkv on your set of putative viral genomes

`checkv end_to_end Viral_Genomes.fasta CheckV_Output -t 23`

CheckV will output viruses and integrated proviruses to separate files. We merge them for downstream analysis

`cat CheckV_Output/viruses.fna CheckV_Output/proviruses.fna > Trimmed_Viral_Genomes.fasta`

The CheckV_Output/quality_summary.tsv file is a table with the results of the CheckV QC. Be aware that CheckV will that trims regions of host contamination from the viral genomes, and the values in this table refer to the metrics in the original sequences, i.e. before this trimming step. If you want a final metric of contamination, completeness, length, etc for the trimmed sequences You should re-run CheckV a second time on Trimmed_Viral_Genomes.fasta

| contig_id | contig_length | provirus | proviral_length | gene_count | viral_genes | host_genes | checkv_quality | miuvig_quality | completeness | completeness_method | contamination | kmer_freq | warnings |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
Blanes_Viromes_ABR_Scaffold_1\|\|full | 157010 | No | NA | 168 | 39 | 12 | Low-quality | Genome-fragment | 43.37 | AAI-based (medium-confidence) | 0 | 1 | 
Blanes_Viromes_FEB_Scaffold_1\|\|full | 111929 | No | NA | 165 | 137 | 0 | High-quality | High-quality | 98.93 | AAI-based (high-confidence) | 0 | 1 | 
Blanes_Viromes_FEB_Scaffold_3\|\|full | 109941 | No | NA | 126 | 4 | 6 | Low-quality | Genome-fragment | 40.19 | HMM-based (lower-bound) | 0 | 1 | 
Blanes_Viromes_FEB_Scaffold_2\|\|full | 109892 | No | NA | 112 | 12 | 3 | High-quality | High-quality | 100 | AAI-based (high-confidence) | 0 | 1 | 
Blanes_Viromes_Sample_2409_S60_Scaffold_1\|\|full | 106895 | No | NA | 127 | 65 | 0 | Low-quality | Genome-fragment | 45.9 | AAI-based (high-confidence) | 0 | 1 | 
Blanes_Viromes_ABR_Scaffold_2\|\|full | 90296 | No | NA | 70 | 11 | 2 | High-quality | High-quality | 100 | AAI-based (medium-confidence) | 0 | 1 | 
Blanes_Viromes_FEB_Scaffold_4\|\|full | 87108 | Yes | 78302 | 119 | 71 | 6 | Low-quality | Genome-fragment | 41.45 | AAI-based (high-confidence) | 10.11 | 1 | 
Blanes_Viromes_Sample_2411_S62_Scaffold_1\|\|full | 81606 | No | NA | 98 | 24 | 2 | High-quality | High-quality | 100 | AAI-based (medium-confidence) | 0 | 1 | 
Blanes_Viromes_MAR_Scaffold_1\|\|full | 76474 | No | NA | 109 | 12 | 5 | Low-quality | Genome-fragment | 43.85 | AAI-based (high-confidence) | 0 | 1 | 
Blanes_Viromes_FEB_Scaffold_5\|\|full | 67511 | No | NA | 66 | 0 | 5 | Low-quality | Genome-fragment | 37.28 | AAI-based (high-confidence) | 0 | 1 | no viral genes detected
Blanes_Viromes_Sample_2411_S62_Scaffold_2\|\|full | 67164 | No | NA | 74 | 20 | 0 | High-quality | High-quality | 100 | AAI-based (medium-confidence) | 0 | 1 | 


## References

- [CheckV Publication](https://doi.org/10.1038/s41587-020-00774-7)