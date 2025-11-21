# Running KEGG Decoder

### References
- [Video tutorial](https://youtu.be/1v4UzjE7K2g?t=962)
- [KEGG](https://www.genome.jp/kegg/)
- [KEGG Decoder Pathway Definitions](https://github.com/bjtully/BioData/blob/master/KEGGDecoder/KOALA_definitions.txt)

### Downloading reference genomes from NCBI (or from any other FTP)
`wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/018/465/GCF_000018465.1_ASM1846v1/GCF_000018465.1_ASM1846v1_genomic.fna.gz`

`wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/015/665/GCF_000015665.1_ASM1566v1/GCF_000015665.1_ASM1566v1_genomic.fna.gz`

`wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/273/775/GCF_001273775.1_ASM127377v1/GCF_001273775.1_ASM127377v1_genomic.fna.gz`

`wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/028/281/495/GCA_028281495.1_ASM2828149v1/GCA_028281495.1_ASM2828149v1_genomic.fna.gz`

`wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/011/805/GCF_000011805.1_ASM1180v1/GCF_000011805.1_ASM1180v1_genomic.fna.gz`

### Decompress files if necessary
`gunzip *gz`

### Rename the sequence identifiers so they comply with Kegg-Decoder naming scheme (i.e. Genome ID is everything that comes before the first "_" character. The sequence IDs should NOT include characters "-", "_", or spaces)
- Note: you might need to change "fna" in the command below to match the extension of your own files (e.g. fasta, or fas)
 
`module lod python`
`for genome in *fna; do id=$(echo $genome | sed -e s'/_\+\|\W\+/./g' ); echo "Changing sequence IDs in" $genome; python3 rename_sequences.py --input $genome --string $id; done`

- Optionally you can rename the sequences in the genome files manually, ony by one, using some other criteria that you might prefer. e.g.:
`python3 rename_sequences.py --input GCF_000018465.1_ASM1846v1_genomic.fna --string Nitrosopumilus`
`python3 rename_sequences.py --input GCF_000015665.1_ASM1566v1_genomic.fna --string Prochlorococcus`
`python3 rename_sequences.py --input GCF_001273775.1_ASM127377v1_genomic.fna Nitrospira`

### Call protein encoding genes with prodigal
`module load prodigal`

- Note: Here also replace all the occurrences of "genome%.fna" by "genome%.fasta" if necessary to match the naming of your files

`for genome in Renamed*; do prodigal -q -p single -a "CDS_${genome%.fna}.faa" -d "Genes_${genome%.fna}.fna" -f gff -i $genome -o "${genome%.fna}.gff" ; done`

### Check if protein IDs follow the Kegg-Decoder notation (Everything that comes before the first _ is considered the genome identifier, evrything after it will be ignored. Proteins from the same genome MUST have the same genome identifier)
`for genome in CDS*faa; do grep ">" $genome | head ; done`

### Concatenate the protein sequence files into a single file

`cat CDS*faa > All_CDS.faa`

### Query proteins against the KOfam database using hmmsearch
`module load hmmer`

`hmmsearch --noali --tblout All_CDSxKOfam.tsv -o All_CDSxKOfam -E 0.001 --cpu 4 /mnt/smart/scratch/vir/felipe/Databases/KEGG/All_KOs.hmm All_CDS.faa`

### Parse the output of the search to contain matches within the desired e-value and bitscore cutoffs and to only print the best KO hit of each protein
`module load python`
`python3 /mnt/smart/scratch/vir/felipe/Repos/ICM_Code/AMG_Hunter.py --cds All_CDS.faa --parse_only True --kegg_annotate True --kegg_hits_file All_CDSxKOfam --kegg_min_score 50 --kegg_max_evalue 0.00001`

### The contents of CDS_Info.tsv will look s displayed below. This file contains only th best hit of each protein. For a full list of hits check KEGG_Info.tsv instead
Sequence |  Length |  KEGG_Best_Subject |  KEGG_Best_Subject_Score |  KEGG_Best_Subject_Function |  KEGG_Best_Subject_Pathways |  KEGG_Best_Subject_Modules |  Pfam_Subjects |
| --- | --- | --- | --- | --- | --- | --- | --- |
  GCA.028281495.1.ASM2828149v1.genomic.fna_Seq_1_1 |  321 |  NA |  0 |  NA |  set() |  set() |  set() |
  GCA.028281495.1.ASM2828149v1.genomic.fna_Seq_1_2 |  392 |  NA |  0 |  NA |  set() |  set() |  set() |
  GCA.028281495.1.ASM2828149v1.genomic.fna_Seq_1_3 |  284 |  K01769 |  335.2 |  E4.6.1.2; guanylate cyclase, other [EC:4.6.1.2] |  {'Purine metabolism'} |  {'Nucleotide metabolism'} |  set() |
  GCA.028281495.1.ASM2828149v1.genomic.fna_Seq_1_4 |  109 |  NA |  0 |  NA |  set() |  set() |  set() |
  GCA.028281495.1.ASM2828149v1.genomic.fna_Seq_1_5 |  119 |  NA |  0 |  NA |  set() |  set() |  set() |
  GCA.028281495.1.ASM2828149v1.genomic.fna_Seq_1_6 |  106 |  NA |  0 |  NA |  set() |  set() |  set() |
  GCA.028281495.1.ASM2828149v1.genomic.fna_Seq_1_7 |  440 |  K01768 |  125.4 |  E4.6.1.1; adenylate cyclase [EC:4.6.1.1] |  {'Meiosis - yeast', 'Biofilm formation - Pseudomonas aeruginosa', 'Purine metabolism'} |  {'Cellular community - prokaryotes', 'Cell growth and death', 'Nucleotide metabolism'} |  set() |
  GCA.028281495.1.ASM2828149v1.genomic.fna_Seq_1_8 |  514 |  NA |  0 |  NA |  set() |  set() |  set() |
  GCA.028281495.1.ASM2828149v1.genomic.fna_Seq_1_9 |  330 |  NA |  0 |  NA |  set() |  set() |  set() |
  GCA.028281495.1.ASM2828149v1.genomic.fna_Seq_1_10 |  276 |  NA |  0 |  NA |  set() |  set() |  set() |
  

### Edit the output of CDS_Info to have only the two columns used by KEGG-Decoder
`cut -f 1,3 CDS_Info.tsv | grep -E -w "K[0-9]+" > Parsed_CDS_Info.tsv`

### Run kegg decoder with the output generated by hmmsearch
`conda activate /mnt/smart/scratch/vir/felipe/envs/keggdecoder`

`KEGG-decoder --input Parsed_CDS_Info.tsv --output Parsed_CDS_Info_Decoder.tsv --vizoption static`

`conda deactivate`

### Visualize the results with your preferred tool
`module load R/4.2.2`

`Rscript /mnt/smart/scratch/vir/felipe/Repos/ICM_Code/Kegg_Decoder_Workshop/Make_Metabolic_Heatmap.R  Parsed_CDS_Info_Decoder.tsv`