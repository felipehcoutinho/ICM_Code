## Introduction

This tutorial will teach you how to identify viral sequences in metagenomic assemblies or integrated into host genomes 

## Software

- [VirSorter2](https://github.com/jiarong/VirSorter2)

## Commands

Load the module and activate the environment

`module load  virsorter/2.2.3`

`conda activate vs2`


If the dataset is small. It is feasible to merge all filtered assembly files of all samples into a single fasta and perform  single virsorter run. This can be done with the cat command, examples:

`cat Filtered_Assembly_Sample_A.fasta Filtered_Assembly_Sample_B.fasta Filtered_Assembly_Sample_C.fasta > ll_Filtered_Assebmlies.fasta`

OR

`cat Filtered_Assembly_Sample_*.fasta > All_Filtered_Assebmlies.fasta`


`virsorter run --db /mnt/smart/scratch/vir/felipe/Databases/VS2/ --seqfile All_Filtered_Assebmlies.fasta --exclude-lt2gene --jobs 8 --include-groups dsDNAphage,NCLDV,ssDNA,lavidaviridae --working-dir All_Filt_Asb_VS2_Output`

The sequences classified as viral will be in the All_Filt_Asb_VS2_Output/final-viral-combined.fa file

The output table All_Filt_Asb_VS2_Output/final-viral-score.tsv will be formatted as shown below:

seqname |  dsDNAphage |  NCLDV |  ssDNA |  lavidaviridae |  max_score |  max_score_group |  length |  hallmark |  viral |  cellular |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
Blanes_Viromes_ABR_Scaffold_1\|\|full |  1 |  0.387 |  0.12 |  0.247 |  1 |  dsDNAphage |  157010 |  12 |  71.4 |  0.6
Blanes_Viromes_FEB_Scaffold_1\|\|full |  1 |  0.207 |  0.147 |  0.127 |  1 |  dsDNAphage |  111929 |  14 |  90.9 |  0
Blanes_Viromes_FEB_Scaffold_3\|\|full |  0.993 |  0.247 |  0 |  0.173 |  0.993 |  dsDNAphage |  109941 |  1 |  31.7 |  1.6
Blanes_Viromes_FEB_Scaffold_2\|\|full |  0.987 |  0.233 |  0.053 |  0.16 |  0.987 |  dsDNAphage |  109892 |  4 |  64 |  0
Blanes_Viromes_FEB_Scaffold_80\|\|full |  0.84 |  0.3 |  0.313 |  0.873 |  0.873 |  lavidaviridae |  26191 |  0 |  78.1 |  0
Blanes_Viromes_ABR_Scaffold_2\|\|full |  0.987 |  0.32 |  0.273 |  0.4 |  0.987 |  dsDNAphage |  90296 |  7 |  76.3 |  1.3
Blanes_Viromes_FEB_Scaffold_4\|\|full |  1 |  0.38 |  0.193 |  0.653 |  1 |  dsDNAphage |  87108 |  9 |  95 |  0
Blanes_Viromes_Sample_2411_S62_Scaffold_1\|\|full |  1 |  0.233 |  0.233 |  0.62 |  1 |  dsDNAphage |  81606 |  11 |  67.7 |  1
Blanes_Viromes_MAR_Scaffold_1\|\|full |  1 |  0.913 |  0.327 |  0.587 |  1 |  dsDNAphage |  76474 |  2 |  83.5 |  0
Blanes_Viromes_Sample_2410_S61_Scaffold_1\|\|full |  0.973 |  0.187 |  0 |  0.027 |  0.973 |  dsDNAphage |  74784 |  6 |  28 |  1.2
Blanes_Viromes_ABR_Scaffold_3\|\|full |  0.993 |  0.14 |  0.007 |  0.013 |  0.993 |  dsDNAphage |  70907 |  2 |  33.7 |  2.3
Blanes_Viromes_Sample_2409_S60_Scaffold_2\|\|full |  0.993 |  0.467 |  0.187 |  0.607 |  0.993 |  dsDNAphage |  67998 |  11 |  64.8 |  0
Blanes_Viromes_FEB_Scaffold_6\|\|full |  1 |  0.913 |  0.293 |  0.62 |  1 |  dsDNAphage |  67539 |  6 |  88 |  1.3
Blanes_Viromes_FEB_Scaffold_5\|\|full |  1 |  0.227 |  0 |  0.053 |  1 |  dsDNAphage |  67511 |  2 |  57.6 |  4.5
Blanes_Viromes_Sample_2411_S62_Scaffold_2\|\|full |  1 |  0.427 |  0.287 |  0.327 |  1 |  dsDNAphage |  67164 |  8 |  52.1 |  0
Blanes_Viromes_FEB_Scaffold_7\|\|full |  0.947 |  1 |  0.073 |  0.307 |  1 |  NCLDV |  66952 |  8 |  77.3 |  5.2
Blanes_Viromes_Sample_2405_S56_Scaffold_1\|\|full |  1 |  0.26 |  0.5 |  0.873 |  1 |  dsDNAphage |  66178 |  2 |  47.6 |  0
Blanes_Viromes_Sample_2411_S62_Scaffold_3\|\|full |  0.993 |  0.293 |  0.087 |  0.187 |  0.993 |  dsDNAphage |  65688 |  7 |  43.2 |  0



## References

- [VirSorter2 Publication](https://link.springer.com/article/10.1186/s40168-020-00990-y)