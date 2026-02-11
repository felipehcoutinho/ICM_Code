## Introduction


This tutorial will teach you how to identify viral sequences in metagenomic assemblies or integrated into host genomes 

## Software

- [VirSorter2](https://github.com/jiarong/VirSorter2)

## Commands

Load the module and activate the environment

`module load virsorter/2.2.4`

module load  virsorter/2.2.3

conda activate vs2

`conda activate virsorter2-2.2.4`

`virsorter config --set LOCAL_SCRATCH=$TMPDIR`

`virsorter config --set HMMSEARCH_THREADS=23`

If the datatse is small. It is feasible to merge all filtered assembly files of all samples into a single fasta and perform  single virsorter run. This can be done with the cat command, examples:

`cat Filtered_Assembly_Sample_A.fasta Filtered_Assembly_Sample_B.fasta Filtered_Assembly_Sample_C.fasta > ll_Filtered_Assebmlies.fasta`

OR

`cat Filtered_Assembly_Sample_*.fasta > All_Filtered_Assebmlies.fasta`

`virsorter run --seqfile All_Filtered_Assebmlies.fasta --exclude-lt2gene --jobs 23 --include-groups dsDNAphage,NCLDV,ssDNA,lavidaviridae --working-dir All_Asb_VS2_Output`

virsorter setup --db-dir /mnt/smart/scratch/vir/felipe/Databases/VS2/ --jobs 12 --skip-deps-install

virsorter run --seqfile ../Assemblies/Assembly_QIK/scaffolds.fasta  --exclude-lt2gene --jobs 12 --include-groups dsDNAphage,NCLDV,ssDNA,lavidaviridae --working-dir QIK_VS2_Output