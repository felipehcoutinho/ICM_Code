# Gene calling and annotation of phage genomes


## Introduction

A thorough characterization of any genome requires identifying genes and assigning them functions. For each group of organisms, specific tools and databaes exist to achieve this.

Prodigal is the standar tool to identify open reading frames (ORFs) among prokaryote genomes and genomes of the viruses that infect them. Most Bacteria and Archaea use [genetic code 11](https://en.wikipedia.org/wiki/List_of_genetic_codes). But many of their viruses use alternative codes (most often 4 and 12), in which stop codons are repurposed to encode amino acids instead. Because this has a large impact on ORF identification, the use of tools that can accomodate this such as Prodigal-gv is warranted.

## Software

- [Prodigal-gv](https://github.com/apcamargo/prodigal-gv)
- [mmseqs2](https://github.com/soedinglab/MMseqs2)
- [hmmer3](http://hmmer.org/)

## Commands

A standard progigal-gv command in which the input file is a mix of viral genomes from a metagenomcis dataset:

`module load python/3.8.5`

`/mnt/smart/scratch/vir/felipe/Software/prodigal-gv -p meta -a Viruses_CDS_Prot_Seqs.faa -d Viruses_CDS_Nucl_Seqs.fna -f gff -i Viral_Genomic_Seqs_Trimmed.fasta -o Viruses_CDS_Coords.gff`

- The -p meta option is suitable for analysing genomic sequences derived from metagenomic datasets, in which the input sequences are devired from different genomes and many of them are incomplete genomes. This option can be omitted when the input file is derived from a single genome if the completeness is high
- The Viruses_CDS_Nucl_Seqs.fna file contains the nucleotide sequences of each individual coding dna sequence identified.
- The Viruses_CDS_Prot_Seqs.fna file contains the translated protein sequences of each individual coding dna sequence identified.
- The Viruses_CDS_Coords.gff file contains information about each of the ORFs regarding coordinates, presence of start and stop codons, frame, and other metrics related to the ORF identification process

When annotating genomes both searches with the nucleotide or amino acid sequences can be performed. Because multiple codons can encode for the same amino acid, two nucleotide sequences with little identity can still encode proteins with high identity at amino acid level. Thus, amino acid searches are more sensitive, which is of special relevance when the genomes being annotated have few close relatives with homologues in the reference database which is the case for viruses from metagenomic datasets

The standard tools for protein homology searches work by performing local alignments, i.e. finding similar regions between query (in our case, the sequences indentified in the viral genomes) and subject sequences (those in the reference database). For this purpose, BLAST, DIAMOND and mmseqs are the most often used tools

There are multiple reference databases available. The choice is important and dependes on the questions being asked. Also, the larger the database, the longer it will take to run the search, and the more likely it is we find false positives. For this exercise we will use a small reference database of curated protein sequences: SwissProt dataset from Uniprot

The fasta file containing the protein sequences is located in: /mnt/smart/scratch/vir/felipe/Databases/Uniprot_02-10-2024/uniprot_sprot.fasta

First, we set up th reference database from this files (this only needs to be done once than this databae cna be reused for additional searchs)

`module load mmseqs2`

`mmseqs createdb uniprot_sprot.fasta DB_SwissProt`

`mmseqs createindex DB_SwissProt tmp`

Now we query our protein sequences of interest against the reference database

`mmseqs easy-search Viruses_CDS_Prot_Seqs.faa DB_SwissProt Viruses_CDS_ProtxDB_SwissProt.m8 tmp --threads 4 --max-seqs 10 --min-seq-id 0.3 --min-aln-len 30 --format-mode 0 --format-output query,theader,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits`

The output file from this command will be in the standard BLASt m8 format, except that the second column will contain not only the identifier of the subjects equence but also its decription, which makes it more informative for this database because it gives us more than just an accession number:

| query | theader | pident | alnlen | mismatch | gapopen | qstart | qend | tstart | tend | evalue | bits |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
IMGVR_UViG_2504756050_000001\|\|2504756050\|\|2504761616_10 | sp\|\|O64240\|\|VG50_BPMD2 Putative adenosylcobalamin-dependent ribonucleoside-triphosphate reductase OS=Mycobacterium phage D29 OX=28369 GN=50 PE=3 SV=1 | 54.9 | 682 | 298 | 0 | 8 | 669 | 11 | 692 | 6.364E-230 | 724 
IMGVR_UViG_2524614843_000001\|\|2524614843\|\|2524708311_106 | sp\|\|O84834\|\|RIR1_CHLTR Ribonucleoside-diphosphate reductase subunit alpha OS=Chlamydia trachomatis serovar D (strain ATCC VR-885 / DSM 19411 / UW-3/Cx) OX=272561 GN=nrdA PE=3 SV=2 | 46.7 | 765 | 398 | 0 | 5 | 752 | 239 | 1003 | 3.145E-203 | 652
IMGVR_UViG_2507262019_000003\|\|2507262019\|\|2507271997_44 | sp\|\|P03708\|\|TERL_LAMBD Terminase, large subunit OS=Escherichia phage lambda OX=2681611 GN=A PE=1 SV=1 | 51.2 | 614 | 295 | 0 | 8 | 621 | 8 | 612 | 1.824E-199 | 636
IMGVR_UViG_2524614843_000001\|\|2524614843\|\|2524708311_128 | sp\|\|O68770\|\|DPO3A_YERPE DNA polymerase III subunit alpha OS=Yersinia pestis OX=632 GN=dnaE PE=3 SV=2 | 39.4 | 785 | 473 | 0 | 4 | 788 | 3 | 784 | 2.758E-164 | 539
IMGVR_UViG_2504756050_000001\|\|2504756050\|\|2504761616_14 | sp\|\|Q05254\|\|DPOL_BPML5 DNA polymerase OS=Mycobacterium phage L5 OX=31757 GN=44 PE=3 SV=1 | 44.3 | 590 | 328 | 0 | 1 | 590 | 1 | 590 | 2.3E-147 | 481
IMGVR_UViG_2519103125_000002\|\|2519103125\|\|2519118532_25 | sp\|\|P14892\|\|CWLA_BACSP N-acetylmuramoyl-L-alanine amidase CwlA OS=Bacillus sp. OX=1409 GN=cwlA PE=3 SV=1 | 91.2 | 251 | 22 | 0 | 1 | 251 | 1 | 251 | 1.707E-152 | 478
IMGVR_UViG_2515075000_000001\|\|2515075000\|\|2515075041_4 | sp\|\|P54308\|\|TERL_BPSPP Terminase, large subunit OS=Bacillus phage SPP1 OX=10724 GN=2 PE=1 SV=1 | 55.9 | 422 | 185 | 0 | 1 | 420 | 1 | 422 | 9.07E-148 | 474
IMGVR_UViG_2519103075_000001\|\|2519103075\|\|2519115393_2 | sp\|\|P45921\|\|YQBE_BACSU Putative prophage capsid protein YqbE OS=Bacillus subtilis (strain 168) OX=224308 GN=yqbE PE=3 SV=2 | 67.2 | 311 | 100 | 0 | 3 | 309 | 1 | 311 | 3.123E-131 | 420

Using hmmer for protein searches

`module load hmmer`

`hmmpress pfam.hmm`

For this exercise we will use a precomputed database:

`hmmsearch -o Viruses_CDS_ProtxPfam --tblout Viruses_CDS_ProtxPfam.tsv --noali --cpu 4 /mnt/smart/scratch/vir/felipe/Databases/PfamA/Pfam-A.hmm Viruses_CDS_Prot_Seqs.faa`

## References

- [Prodigal-gv Publication](https://www.nature.com/articles/s41587-023-01953-y)
- [MMseqs2 Publication](https://www.nature.com/articles/nbt.3988)
- [hmmer Publication](https://doi.org/10.1371/journal.pcbi.1002195)
- [SwissProt Database](https://www.expasy.org/resources/uniprotkb-swiss-prot)
