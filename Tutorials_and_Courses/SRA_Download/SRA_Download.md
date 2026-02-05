module load python/3.8.5
module load bowtie2/2.3.4.1
module load samtools

#reads in: /mnt/biostore/data/fs/bio/malaspina/vertical-profiles-clean-reads/ cannot be accessed from compute nodes

cd /mnt/smart/scratch/vir/felipe/Databases/OceanDNA/

#build database

# bowtie2-build --threads 47 /mnt/smart/scratch/vir/felipe/Databases/OceanDNA/OceanDNA_All_Species_Rep_MAGs_Scaffolds.fasta OceanDNA_All_Species_Rep_MAGs_Scaffolds

###Map to malaspina deep

GENOMEFILE="/mnt/smart/scratch/vir/felipe/Databases/OceanDNA/OceanDNA_All_Species_Rep_MAGs_Scaffolds.fasta"
DBPREFIX="OceanDNA_All_Species_Rep_MAGs_Scaffolds"
MAXMATCHES=100
MAXREADS=10000000
THREADS=47
SRALISTFILE="/mnt/smart/scratch/vir/felipe/Metadata/Malaspina/Malaspina_Deep_SRA_IDs_bottom55.txt"
#"/mnt/smart/scratch/vir/felipe/Metadata/Malaspina/Malaspina_Deep_SRA_IDs.txt"
RAWREADTABLE="/mnt/smart/scratch/vir/felipe/Metadata/Malaspina/Malaspina_Deep_Metagenome_Info.tsv"
#"/mnt/smart/scratch/vir/felipe/Metadata/Malaspina/Malaspina_Deep_Metagenome_Info_Top58.tsv"
#

#Dont have the necessary glib modules to run version 3.3.0
# while read sra_id; do echo "Processing ${sra_id}"; /mnt/smart/scratch/vir/felipe/Software/sratoolkit.3.1.1-ubuntu64/bin/prefetch --output-directory . $sra_id; /mnt/smart/scratch/vir/felipe/Software/sratoolkit.3.1.1-ubuntu64/bin/fasterq-dump --threads $THREADS --split-files $sra_id ; bowtie2 -x $DBPREFIX -q -1 ${sra_id}_1.fastq -2 ${sra_id}_2.fastq -S ${sra_id}.sam -k $MAXMATCHES -u $MAXREADS --no-discordant --no-mixed --no-unal --threads $THREADS ; samtools view -bS ${sra_id}.sam > ${sra_id}.bam ; samtools sort ${sra_id}.bam -o ${sra_id}.sorted.bam ; samtools index ${sra_id}.sorted.bam ; samtools idxstats ${sra_id}.sorted.bam > ${sra_id}.Counts.tsv ; rm -fr $sra_id ; rm -f $sra_id*fastq ; rm -f ${sra_id}.sam ${sra_id}.bam ${sra_id}.sorted.bam ${sra_id}.sorted.bam.bai; done < $SRALISTFILE

# python3 /mnt/smart/scratch/vir/felipe/Repos/virathon/Virathon.py --genome $GENOMEFILE --abundance_table True --abundance_rpkm True --bowtie_k $MAXMATCHES --abundance_max_reads $MAXREADS --threads $THREADS --bowtiedb $DBPREFIX --raw_read_table $RAWREADTABLE --parse True

# rm -f /mnt/smart/scratch/vir/felipe/Databases/OceanDNA/All_Genomic.fasta

mv SRR* /Abundances/MalaspinaDeep/
mv *Abundance*tsv /Abundances/MalaspinaDeep/
###Map to malaspina profiles

#profile reads available locally:
# /mnt/smart/scratch/vir/felipe/Malaspina_Profiles/PostQC_Reads/

sed -i s"/\/mnt\/lustre\/scratch\/nlsas\/home\/csic\/eyg\/fhc\/Malaspina\/Profiles\/PostQC_Reads/\/mnt\/smart\/scratch\/vir\/felipe\/Malaspina_Profiles\/PostQC_Reads/g" /mnt/smart/scratch/vir/felipe/Metadata/Malaspina/Profiles_Malaspina_76_FL_Metagenome_Info.tsv

RAWREADTABLE="/mnt/smart/scratch/vir/felipe/Metadata/Malaspina/Profiles_Malaspina_76_FL_Metagenome_Info.tsv"

python3 /mnt/smart/scratch/vir/felipe/Repos/virathon/Virathon.py --genome $GENOMEFILE --abundance_table True --abundance_rpkm True --bowtie_k $MAXMATCHES --abundance_max_reads $MAXREADS --threads $THREADS --bowtiedb $DBPREFIX --raw_read_table $RAWREADTABLE

rm -f /mnt/smart/scratch/vir/felipe/Databases/OceanDNA/All_Genomic.fasta

mv MP* Abundances/MalaspinaProfiles/


###Download and keep mlp deep reads
cd /mnt/smart/scratch/vir/felipe/Databases/Malaspina/Reads/

THREADS=12
SRALISTFILE="/mnt/smart/scratch/vir/felipe/Metadata/Malaspina/Malaspina_Deep_SRA_IDs.txt"

while read sra_id; do echo "Processing ${sra_id}"; /mnt/smart/scratch/vir/felipe/Software/sratoolkit.3.1.1-ubuntu64/bin/prefetch --output-directory . $sra_id; /mnt/smart/scratch/vir/felipe/Software/sratoolkit.3.1.1-ubuntu64/bin/fasterq-dump --threads $THREADS --split-files $sra_id ; rm -fr $sra_id ; done < $SRALISTFILE

while read line; do mid=$(echo "$line" | cut -f 1 ); sid=$(echo "$line" | cut -f 2); echo "Renaming ${sid} to ${mid}"; mv ${sid}_1.fastq  ${mid}_1.fastq; mv ${sid}_2.fastq  ${mid}_2.fastq; done < /mnt/smart/scratch/vir/felipe/Databases/Malaspina/Reads/Malaspina_Deep_SRR_to_MP_IDs.txt


