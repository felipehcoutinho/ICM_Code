# Gene calling and annotation of phage genomes


## Introduction

A thorough characterization of any genome requires identifying genes and assigning them functions. For each group of organisms, specific tools and databases exist to achieve this.

Prodigal is the standard tool to identifying open reading frames (ORFs) among prokaryote genomes and genomes of the viruses that infect them. Most Bacteria and Archaea use [genetic code 11](https://en.wikipedia.org/wiki/List_of_genetic_codes). But many of their viruses use alternative codes (most often 4 and 12), in which stop codons are repurposed to encode amino acids instead. Because this has a large impact on ORF identification, the use of tools that can identify genomes using alterntive codes automatically is warranted, for this exercise we will use Prodigal-gv for this.

## Software

- [Prodigal-gv](https://github.com/apcamargo/prodigal-gv)
- [mmseqs2](https://github.com/soedinglab/MMseqs2)
- [hmmer3](http://hmmer.org/)

## Commands

## Gene calling

A standard progigal-gv command in which the input file is a mix of viral genomes from a metagenomcis dataset:

`module load python/3.8.5`

`/mnt/smart/scratch/vir/felipe/Software/prodigal-gv -p meta -a Viruses_CDS_Prot_Seqs.faa -d Viruses_CDS_Nucl_Seqs.fna -f gff -i Viral_Genomic_Seqs_Trimmed.fasta -o Viruses_CDS_Coords.gff`

- The -p meta option is suitable for analysing genomic sequences derived from metagenomic datasets, in which the input sequences are devired from different genomes and many of them are incomplete genomes. This option can be omitted when the input file is derived from a single genome (if the completeness is high and the total sequences are of more than 200 Kbp)
- The Viruses_CDS_Nucl_Seqs.fna file contains the nucleotide sequences of each individual coding dna sequence identified.
- The Viruses_CDS_Prot_Seqs.fna file contains the translated amino acid sequences of each individual coding dna sequence identified.
- The Viruses_CDS_Coords.gff file contains information about each of the ORFs regarding coordinates, presence of start and stop codons, frame, and other metrics related to the ORF identification process

## Homology searches

When annotating genomes both searches at the nucleotide or amino acid sequence levels can be performed. Because multiple codons can encode for the same amino acid, two nucleotide sequences with little identity can still encode proteins with high identity at amino acid level. Thus, amino acid searches are more sensitive, which is of special relevance when the genomes being annottated have few close relatives with homologues in the reference database, which is the case for viruses from metagenomic datasets

The standard tools for protein homology searches work by performing local alignments, i.e. finding similar regions between query (in our case, the sequences indentified in the viral genomes) and subject sequences (those in the reference database). For this purpose, BLAST, DIAMOND and mmseqs are the most often used tools.

There are multiple reference databases available. The choice is important and dependes on the questions being asked. Also, the larger the database, the longer it will take to run the search, and the more likely it is we find false positives. For this exercise we will use a small reference database of curated protein sequences: SwissProt dataset from Uniprot

### mmseqs search

The fasta file containing the protein sequences is located in: /mnt/smart/scratch/vir/felipe/Databases/Uniprot_02-10-2024/uniprot_sprot.fasta

First, we set up a reference database by formatting the fasta file (this only needs to be done once than this database can be reused for additional searches with different queries)

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
IMGVR_UViG_2515075000_000001\|\|2515075000\|\|2515075041_4 | sp\|P45916\|YQAT_BACSU Uncharacterized protein YqaT OS=Bacillus subtilis (strain 168) OX=224308 GN=yqaT PE=4 SV=1 | 30.700 | 385 | 261 | 0 | 33 | 409 | 25 | 409 | 1.298E-33 | 140
IMGVR_UViG_2507262019_000003\|\|2507262019\|\|2507271997_44 | sp\|\|P03708\|\|TERL_LAMBD Terminase, large subunit OS=Escherichia phage lambda OX=2681611 GN=A PE=1 SV=1 | 51.2 | 614 | 295 | 0 | 8 | 621 | 8 | 612 | 1.824E-199 | 636
IMGVR_UViG_2524614843_000001\|\|2524614843\|\|2524708311_128 | sp\|\|O68770\|\|DPO3A_YERPE DNA polymerase III subunit alpha OS=Yersinia pestis OX=632 GN=dnaE PE=3 SV=2 | 39.4 | 785 | 473 | 0 | 4 | 788 | 3 | 784 | 2.758E-164 | 539
IMGVR_UViG_2504756050_000001\|\|2504756050\|\|2504761616_14 | sp\|\|Q05254\|\|DPOL_BPML5 DNA polymerase OS=Mycobacterium phage L5 OX=31757 GN=44 PE=3 SV=1 | 44.3 | 590 | 328 | 0 | 1 | 590 | 1 | 590 | 2.3E-147 | 481
IMGVR_UViG_2519103125_000002\|\|2519103125\|\|2519118532_25 | sp\|\|P14892\|\|CWLA_BACSP N-acetylmuramoyl-L-alanine amidase CwlA OS=Bacillus sp. OX=1409 GN=cwlA PE=3 SV=1 | 91.2 | 251 | 22 | 0 | 1 | 251 | 1 | 251 | 1.707E-152 | 478
IMGVR_UViG_2515075000_000001\|\|2515075000\|\|2515075041_4 | sp\|\|P54308\|\|TERL_BPSPP Terminase, large subunit OS=Bacillus phage SPP1 OX=10724 GN=2 PE=1 SV=1 | 55.9 | 422 | 185 | 0 | 1 | 420 | 1 | 422 | 9.07E-148 | 474
IMGVR_UViG_2519103075_000001\|\|2519103075\|\|2519115393_2 | sp\|\|P45921\|\|YQBE_BACSU Putative prophage capsid protein YqbE OS=Bacillus subtilis (strain 168) OX=224308 GN=yqbE PE=3 SV=2 | 67.2 | 311 | 100 | 0 | 3 | 309 | 1 | 311 | 3.123E-131 | 420 |


### hmmer search

Another strategy to annotate proteins relie son Hidden Markov Models. HMMs are generated from the concensus information obtaied from grouping together multiple homologous proteins from different organisms. Because of this, HMM searches are even more sentitive than protein searches performed with local alignment search tools such as BLAST and mmseqs. 

`module load hmmer`

First we set up a reference database from the hmm profiles

`hmmpress pfam.hmm`

For this exercise we will use a pre-computed database:

`hmmsearch -o Viruses_CDS_ProtxPfam --tblout Viruses_CDS_ProtxPfam.hmmer_tbl --noali --cpu 4 /mnt/smart/scratch/vir/felipe/Databases/PfamA/Pfam-A.hmm Viruses_CDS_Prot_Seqs.faa`

The generated Viruses_CDS_ProtxPfam.hmmer_tbl file will be formatted as shown below:

\#                                                               --- full sequence ---- --- best 1 domain ---- --- domain number estimation ----
\# target name        accession  query name           accession    E-value  score  bias   E-value  score  bias   exp reg clu  ov env dom rep inc description of target
\#------------------- ---------- -------------------- ---------- --------- ------ ----- --------- ------ -----   --- --- --- --- --- --- --- --- ---------------------

IMGVR_UViG_2519899781_000002\|2519899781\|2519959787_71 -          2-Hacid_dh_C         PF02826.24   9.4e-06   21.6   0.0   2.4e-05   20.3   0.0   1.6   2   0   0   2   2   2   1 \# 46352 \# 47455 \# -1 \# ID=17_71;partial=00;start_type=ATG;rbs_motif=AGGA/GGAG/GAGG;rbs_spacer=11-12bp;genetic_code=11;gc_cont=0.282
IMGVR_UViG_2524614843_000001\|2524614843\|2524708311_9 -          2Fe-2S_Ferredox      PF11591.13    0.0092   12.4   0.3     0.022   11.2   0.3   1.6   1   0   0   1   1   1   1 \# 7670 \# 8272 \# -1 \# ID=19_9;partial=00;start_type=ATG;rbs_motif=GGA/GAG/AGG;rbs_spacer=5-10bp;genetic_code=11;gc_cont=0.544
IMGVR_UViG_2519103147_000002\|2519103147\|2519120960_35 -          2TM                  PF13239.11     0.008   12.8   1.3    0.0086   12.7   0.8   1.4   1   1   0   1   1   1   1 \# 18332 \# 18556 \# -1 \# ID=13_35;partial=00;start_type=ATG;rbs_motif=AGGAGG;rbs_spacer=5-10bp;genetic_code=11;gc_cont=0.347
IMGVR_UViG_2519103127_000001\|2519103127\|2519118661_13 -          2TM                  PF13239.11    0.0081   12.8   1.3    0.0086   12.7   0.8   1.4   1   1   0   1   1   1   1 \# 9037 \# 9261 \# -1 \# ID=12_13;partial=00;start_type=ATG;rbs_motif=AGGAGG;rbs_spacer=5-10bp;genetic_code=11;gc_cont=0.382
IMGVR_UViG_2519103075_000001\|2519103075\|2519115393_20 -          2TM                  PF13239.11    0.0085   12.7   1.4    0.0093   12.6   0.9   1.4   1   1   0   1   1   1   1 \# 15064 \# 15288 \# 1 \# ID=9_20;partial=00;start_type=ATG;rbs_motif=AGGAGG;rbs_spacer=5-10bp;genetic_code=11;gc_cont=0.347
IMGVR_UViG_2519899585_000002\|2519899585\|2519929593_8 -          3keto-disac_hyd      PF06439.17    0.0011   15.5   0.0      0.21    8.1   0.0   2.6   2   0   0   2   2   2   2 \# 6280 \# 8655 \# -1 \# ID=14_8;partial=00;start_type=GTG;rbs_motif=AGGAG;rbs_spacer=5-10bp;genetic_code=11;gc_cont=0.641
IMGVR_UViG_2524614843_000001\|2524614843\|2524708311_123 -          5_3_exonuc           PF01367.25   2.3e-18   63.0   0.0   4.6e-18   62.0   0.0   1.5   1   0   0   1   1   1   1 \# 93439 \# 94353 \# -1 \# ID=19_123;partial=00;start_type=ATG;rbs_motif=GGAG/GAGG;rbs_spacer=5-10bp;genetic_code=11;gc_cont=0.589
IMGVR_UViG_2524614843_000001\|2524614843\|2524708311_123 -          5_3_exonuc_N         PF02739.21   1.8e-30  102.3   0.0   3.4e-30  101.4   0.0   1.4   1   0   0   1   1   1   1 \# 93439 \# 94353 \# -1 \# ID=19_123;partial=00;start_type=ATG;rbs_motif=GGAG/GAGG;rbs_spacer=5-10bp;genetic_code=11;gc_cont=0.589
IMGVR_UViG_2519103072_000001\|2519103072\|2519115199_9 -          60KD_IMP             PF02096.25     0.015   11.8   1.4      0.02   11.4   1.4   1.1   1   0   0   1   1   1   0 \# 6305 \# 6613 \# -1 \# ID=7_9;partial=00;start_type=ATG;rbs_motif=AGxAGG/AGGxGG;rbs_spacer=11-12bp;genetic_code=11;gc_cont=0.340
IMGVR_UViG_2519103075_000001\|2519103075\|2519115393_30 -          7TMR-DISM_7TM        PF07695.16      0.61    6.5  20.4     0.018   11.4   7.7   2.0   1   1   0   2   2   2   0 \# 22957 \# 23646 \# -1 \# ID=9_30;partial=00;start_type=ATG;rbs_motif=AGGAG;rbs_spacer=5-10bp;genetic_code=11;gc_cont=0.293
IMGVR_UViG_2515075000_000001\|2515075000\|2515075041_4 -          RPA_interact_N       PF14766.11    0.0096   12.2   0.1      0.02   11.1   0.1   1.6   1   0   0   1   1   1   1 \# 2629 \# 3906 \# -1 \# ID=5_4;partial=00;start_type=ATG;rbs_motif=None;rbs_spacer=None;genetic_code=11;gc_cont=0.357
IMGVR_UViG_2515075000_000001\|2515075000\|2515075041_4  -          Terminase_3          PF04466.18   5.4e-53  176.0   0.4   8.2e-53  175.4   0.4   1.2   1   0   0   1   1   1   1 \# 2629 \# 3906 \# -1 \# ID=5_4;partial=00;start_type=ATG;rbs_motif=None;rbs_spacer=None;genetic_code=11;gc_cont=0.357
IMGVR_UViG_2515075000_000001\|2515075000\|2515075041_4  -          Terminase_3C         PF17288.7    9.9e-17   58.2   0.7   1.7e-16   57.5   0.0   1.7   2   0   0   2   2   2   1 \# 2629 \# 3906 \# -1 \# ID=5_4;partial=00;start_type=ATG;rbs_motif=None;rbs_spacer=None;genetic_code=11;gc_cont=0.357
IMGVR_UViG_2515075000_000001\|2515075000\|2515075041_4  -          Terminase_6N         PF03237.20   1.8e-11   40.6   0.0   2.6e-11   40.1   0.0   1.2   1   0   0   1   1   1   1 \# 2629 \# 3906 \# -1 \# ID=5_4;partial=00;start_type=ATG;rbs_motif=None;rbs_spacer=None;genetic_code=11;gc_cont=0.357

### Conclusions

Lets take IMGVR_UViG_2515075000_000001|2515075000|2515075041_4 as an example. This protein produces hits in both searches. In both searches, it produces multiple hits. For this cases, we focus on the E-value (the lower the mroe reliable the hit) and Bitscore (the higher the more reliable the hit). For this case, our bets hit in SwissProt to [TERL_BPSPP](https://www.uniprot.org/uniprotkb/P54308/entry), a Terminase, large subunit associated with Bacillus phage SPP1. Likewise, the best match on Pfam is to [PF04466 Terminase_3](https://www.ebi.ac.uk/interpro/entry/pfam/PF04466/), which corroborates functional annotation obtained with the SwissProt search. But because we are dealing with an HMM, there is no taxonomic annotation for the hit. In this case, we can be confident of the functional annotation as a terminase large subuni, but because the the identity was of the SwissProt hit was low (only 55%), we cannot state this protein was derived from Bacillus phage SPP1.

## References

- [Prodigal-gv Publication](https://www.nature.com/articles/s41587-023-01953-y)
- [MMseqs2 Publication](https://www.nature.com/articles/nbt.3988)
- [hmmer Publication](https://doi.org/10.1371/journal.pcbi.1002195)
- [SwissProt Database](https://www.expasy.org/resources/uniprotkb-swiss-prot)
- [Pfam Database](http://pfam.xfam.org/)
