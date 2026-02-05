### Assembling Illumina paired-end reads using SPAdes

## Introduction

- [Basics of de Brujin Assembly explained by the amazing Rob Edwards](https://www.youtube.com/watch?v=OY9Q_rUCGDw)
- [A more detailed explanation of de Brujin Assembly](https://www.youtube.com/watch?v=B3Oc_3zZvTg)
- [Theoretical/Mathematical foundation of genome assembly using de Brujin graphs](https://www.nature.com/articles/nbt.2023)

## Software

- [SPAdes](https://github.com/ablab/spades)

## Assembly Commands

`module load spades/4.2.0`

Assembling a single metagenome sample:

`spades.py -1 /Clean_Reads/CJ961_1.fq.gz -2 /Clean_Reads/CJ961_1.fq.gz -o Assembly_CJ961 --threads 48 --memory 500 --meta`

NOTE: Omit the --meta flag when assembling isolate/pure culture genomes

Performing a co-assembly using multiple metagenome samples:

`spades.py -1 /Clean_Reads/CJ961_1.fq.gz /Clean_Reads/CJ222_1.fq.gz /Clean_Reads/CJ333_1.fq.gz -2 /Clean_Reads/CJ961_2.fq.gz /Clean_Reads/CJ222_2.fq.gz /Clean_Reads/CJ333_2.fq.gz -o CoAssembly_CJ961_222_333 --threads 48 --memory 500 --meta`

## Assembly output files & directories

Each succesful assembly will generate the following in the specified output directory:

- assembly_graph_after_simplification.gfa
- assembly_graph.fastg
- assembly_graph_with_scaffolds.gfa
- before_rr.fasta
- contigs.fasta
- contigs.paths
- corrected
- dataset.info
- first_pe_contigs.fasta
- input_dataset.yaml
- K21
- K33
- K55
- misc
- params.txt
- pipeline_state
- run_spades.sh
- run_spades.yaml
- scaffolds.fasta
- scaffolds.paths
- spades.log
- strain_graph.gfa
- tmp

We are interested in the scaffolds.fasta files specifically.

First, we look at the basic stats of the assemblies

`module load perl/5.36`

`perl /mnt/smart/scratch/vir/felipe/Repos/ICM_Code/assembly_stats_table.pl  Assembly*/scaffolds.fasta  > Assembly_Stats.tsv`

SPAdes can often ouput short scaffolds hich need to be filtered out. Also, the naming scheme used for the scaffolds is nto specific for each sample. Next we filter all scaffolds < 10000 bp and rename them according to the samples from which they come from

`module load python/3.8.5`

`for asb in Assembly_*/scaffolds.fasta ; do nam=$(echo "$asb" | sed 's/Assembly_//' | sed 's/\/scaffolds.fasta//'); echo $nam ; python3 /mnt/smart/scratch/vir/felipe/Repos/ICM_Code/rename_sequences.py --input $asb --output Filtered_10Kbp_Renamed_${nam}.fasta --min_length 10000 --string_rename ${nam}_Scaffold_ ; done`

We check the status of the assemblies filtered and renamed files

`perl /mnt/smart/scratch/vir/felipe/Repos/ICM_Code/assembly_stats_table.pl  Filtered_10Kbp_Renamed*.fasta  > Filtered_Assembly_Stats.tsv`


## References

- [SPAdes Citation](https://doi.org/10.1002/cpbi.102)