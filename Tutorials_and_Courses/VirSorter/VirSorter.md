## Introduction

This tutorial will teach you how to identify viral sequences in metagenomic assemblies or integrated into host genomes. In metagenomic datasets there is usually a mix of organisms and their viruses. Even the viromes have some lelvel of contamination from celular organisms. Also, many viruses are capale of integrating their genomes into those of their host, even hen dealing with isolate genome data, these viral sequences might be present. Because of that, in metagenomic datasets we can never assume that all sequences in any datasets are all viral or all cellular. Often, these viral sequences will be fragmented, which makes their identification even harder. VirSorter2 is specially suitable for identifying viruses that infect Bacteria (i.e. bacteriophages) and Archaea, and diferrentiating between sequences that are entirely viral from those that are integrated into a prokaryote genome. Although it can also look for NCLDV whih ifnect eukaryotes, virophages (which infect the NCLDV), as well as ssDNA and RNA viruses that infect all organisms. 

## Software

- [VirSorter2](https://github.com/jiarong/VirSorter2)

## Commands

Load the module and activate the environment


`conda activate /mnt/smart/scratch/vir/felipe/envs/vs2`


The input files are assembled sequences (NOT sequencing reads). It is recommended to do a filtering to remove sequences < 5K bp becse those have little information to be confidentely classified as viral. Some protocols recommend removing sequences < 10 Kbp for even mroe accuracy. Keep in mind that dsDNA viruses of Bacteria have genomes in the range of (15 Kbp - 250 Kbp, with Jumbo phages reaching up to 750 Kbp), while the NCLDV can be > 1 Mbp

If the dataset is small. It is feasible to merge all filtered assembly files of all samples into a single fasta and perform  single virsorter run. This can be done with the cat command, examples:

`cat Filtered_Assembly_Sample_A.fasta Filtered_Assembly_Sample_B.fasta Filtered_Assembly_Sample_C.fasta > ll_Filtered_Assebmlies.fasta`

OR

`cat Filtered_Assembly_Sample_*.fasta > All_Filtered_Assebmlies.fasta`

We can speed up the run by setting virsorter to use more threads specifically during the hmmsearch step

`virsorter config --set HMMSEARCH_THREADS=8`

The command below will run VirSorter2, exclude sequences with less than 2 genes, as these cannot be reliabily classified as viruses. The program will look for viruses in the groups dsDNAphage, NCLDV, ssDNA, lavidaviridae. If we are interested in changing this behavior,we can run the command with instead with, e.g. `--include-groups dsDNAphage` or `--include-groups dsDNAphage,NCLDV`

`virsorter run --seqfile All_Filtered_Assebmlies.fasta --exclude-lt2gene --jobs 8 --include-groups dsDNAphage,NCLDV,ssDNA,lavidaviridae --db-dir /mnt/smart/scratch/vir/felipe/Databases/VS2/ --working-dir All_Filt_Asb_VS2_Output all`

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