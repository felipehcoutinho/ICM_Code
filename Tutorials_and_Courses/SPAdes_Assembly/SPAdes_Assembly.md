### Assembling Illumina paired-end reads using SPAdes

## Introduction

[Basics of de Brujin Assembly explained by the amazing Rob Edwards](https://www.youtube.com/watch?v=OY9Q_rUCGDw)
[A more detailed explanation of de Brujin Assembly](https://www.youtube.com/watch?v=B3Oc_3zZvTg)
[Theoretical/Mathematical foundation of genome assembly using de Brujin graphs](https://www.nature.com/articles/nbt.2023)

## Software

[SPAdes](https://github.com/ablab/spades)

## Assembly Commands

`module load spades/4.2.0`

Assembling a single metagenome sample

`spades.py -1 /Clean_Reads/CJ961_1.fq.gz -2 /Clean_Reads/CJ961_1.fq.gz -o Assembly_CJ961 --threads 48 --memory 500 --meta`

NOTE: Omit the --meta flag when assembling isolate/pure culture genomes

Performing a co-assembly using multiple metagenome samples

`spades.py -1 /Clean_Reads/CJ961_1.fq.gz /Clean_Reads/CJ222_1.fq.gz /Clean_Reads/CJ333_1.fq.gz -2 /Clean_Reads/CJ961_2.fq.gz /Clean_Reads/CJ222_2.fq.gz /Clean_Reads/CJ333_2.fq.gz -o CoAssembly_CJ961_222_333 --threads 48 --memory 500 --meta`

## Assembly output files

## References

(SPAdes Citation)(https://doi.org/10.1002/cpbi.102)