#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=02-00:00:00
#SBATCH --job-name=Antartic_Lagoons
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=75G
#SBATCH --output=/mnt/smart/users/fcoutinho/Job_Logs/jobLog_%J.out
#SBATCH --error=/mnt/smart/users/fcoutinho/Job_Logs/jobLog_%J.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=fhernandes@icm.csic.es

##############################################
# module load python/3.8.5
# module load spades/4.0.0
# module load megahit

# cd /mnt/smart/users/fcoutinho/Antarctic_Lagoons/WGS_Assemblies/

# # python3 /mnt/smart/users/fcoutinho/Repos/virathon/Virathon.py --max_ram 730 --threads 48 --assemble True --raw_read_table /mnt/smart/users/fcoutinho/Antarctic_Lagoons/Metadata/WGS_Samples.tsv

# python3 /mnt/smart/users/fcoutinho/Repos/virathon/Virathon.py --max_ram 730 --threads 47 --assembler megahit --raw_read_table /mnt/smart/users/fcoutinho/Antarctic_Lagoons/Metadata/WGS_Samples.tsv

# perl /mnt/smart/users/fcoutinho/Repos/ICM_Code/assembly_stats_table.pl */*fa  > MG_Assembly_Stats.tsv

# #Filter and rename scaffolds
# for i in Assembly*/final.contigs.fa; do uid=$(echo "$i" | sed 's/Assembly_//g' | cut -d "/" -f 1); python3 /mnt/smart/users/fcoutinho/Repos/virathon/Virathon.py --genome_files $i --rename_seqs True --string_rename ${uid}_Seq_ --min_length 1000 ; mv All_Genomic.fasta ${uid}_Filtered_Scaffolds_1kbp.fasta; done

# perl /mnt/smart/users/fcoutinho/Repos/ICM_Code/assembly_stats_table.pl *Filtered_Scaffolds_1kbp.fasta  > Antarctic_Lagoons_Filtered_Scaffolds_1kbp_MG_Assembly_Stats.tsv

# cat *_Filtered_Scaffolds_1kbp.fasta > Antarctic_Lagoons_All_Filtered_Scaffolds_1kbp.fasta

# #Filter for viral identification

# python3 /mnt/smart/users/fcoutinho/Repos/virathon/Virathon.py --genome_files Antarctic_Lagoons_All_Filtered_Scaffolds_1kbp.fasta --min_length 5000

# mv All_Genomic.fasta Antarctic_Lagoons_All_Filtered_Scaffolds_5kbp.fasta

# #Filter for binning
# python3 /mnt/smart/users/fcoutinho/Repos/virathon/Virathon.py --genome_files Antarctic_Lagoons_All_Filtered_Scaffolds_1kbp.fasta --min_length 2500

# mv All_Genomic.fasta Antarctic_Lagoons_All_Filtered_Scaffolds_2.5kbp.fasta

###Contib abundances for binning
# module load python/3.8.5
# module load bowtie2/2.3.4.1
# module load samtools

# cd /mnt/smart/users/fcoutinho/Antarctic_Lagoons/Binning/

# #Notice command below max reads is set to 1M to speed up for a trial run, and bowtie_k set to 100 becaseu the contigs detabase includes all contigs from all 20 samples and redundancy in read mapping is expected
# # python3 /mnt/smart/users/fcoutinho/Repos/virathon/Virathon.py --genome_files /mnt/smart/users/fcoutinho/Antarctic_Lagoons/Binning/Antarctic_Lagoons_All_Filtered_Scaffolds_2.5kbp.fasta --abundance_table True --abundance_rpkm True --raw_read_table /mnt/smart/users/fcoutinho/Antarctic_Lagoons/Metadata/WGS_Samples.tsv --bowtie_k 100 --abundance_max_reads 1000000 --threads 47
# #using the command above, many scaffolds have abundance of 0 in their own sample. Bins are not formed using the abundance file derived from it either

# # python3 /mnt/smart/users/fcoutinho/Repos/virathon/Virathon.py --genome_files /mnt/smart/users/fcoutinho/Antarctic_Lagoons/Binning/Antarctic_Lagoons_All_Filtered_Scaffolds_2.5kbp.fasta --bowtiedb /mnt/smart/users/fcoutinho/Antarctic_Lagoons/Binning/All_Genomic --abundance_table True --abundance_rpkm True --raw_read_table /mnt/smart/users/fcoutinho/Antarctic_Lagoons/Metadata/WGS_Samples.tsv --bowtie_k 100 --abundance_max_reads 15000000 --threads 23

# python3 /mnt/smart/users/fcoutinho/Repos/virathon/Virathon.py --genome_files /mnt/smart/users/fcoutinho/Antarctic_Lagoons/Binning/Antarctic_Lagoons_All_Filtered_Scaffolds_2.5kbp.fasta --bowtiedb /mnt/smart/users/fcoutinho/Antarctic_Lagoons/Binning/All_Genomic --abundance_table True --abundance_rpkm True --raw_read_table /mnt/smart/users/fcoutinho/Antarctic_Lagoons/Metadata/Subset_WGS_Samples.tsv --bowtie_k 100 --abundance_max_reads 15000000 --threads 47

# #Split scaffolds per sample
# python3 /mnt/smart/users/fcoutinho/Repos/ICM_Code/Get_Seqs.py --input_sequences /mnt/smart/users/fcoutinho/Antarctic_Lagoons/Binning/Antarctic_Lagoons_All_Filtered_Scaffolds_2.5kbp.fasta --id_split_sep _ --id_split_pos 0 --fetch_all True --group True

# #Metabat binning
# module load metabat
# conda activate /mnt/smart/apps/condaEnvs/metabat-2.15

# #/mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/Software/metabat/bin/jgi_summarize_bam_contig_depths

# jgi_summarize_bam_contig_depths --outputDepth JGI_Depths_Antarctic_Lagoons_All_Filtered_Scaffolds_2.5kbp.tsv *bam

# #/mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/Software/metabat/bin/metabat2

# for i in C*.fasta; do smp=$(echo "$i" | sed 's/.fasta//g') ; echo $smp; metabat2 -i $i --abdFile  /mnt/smart/users/fcoutinho/Antarctic_Lagoons/Binning/JGI_Depths_Antarctic_Lagoons_All_Filtered_Scaffolds_2.5kbp.tsv -s 200000 -t 23 --outFile AL_R1_${smp}_Bin; done

# conda deactivate

# # rm -f C*fasta

# mkdir /mnt/smart/users/fcoutinho/Antarctic_Lagoons/Binning/AL_Bins_R1

# mv AL_R1_*Bin* /mnt/smart/users/fcoutinho/Antarctic_Lagoons/Binning/AL_Bins_R1/

# cd /mnt/smart/users/fcoutinho/Antarctic_Lagoons/Binning/AL_Bins_R1/

# module load perl

# #perl /mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/repos/ICM_Code/assembly_stats_table.pl fa > AL_R1_Metabat_Bins_Assembly_Stats.tsv

# perl  /mnt/smart/users/fcoutinho/Repos/ICM_Code/assembly_stats_table.pl fa > AL_R1_Metabat_Bins_Assembly_Stats.tsv
###CheckM2 Bin QC
# module load checkm2/1.0.1
# conda activate /mnt/smart/apps/condaEnvs/checkm2

###CHECKM2 with custom conda env
module load miniconda
source /home/apps/miniconda/4.5.x/etc/profile.d/conda.sh
conda activate /mnt/smart/users/fcoutinho/envs/checkm2
export CHECKM2DB="/mnt/smart/apps/bio/app-databases/checkm2/uniref100.KO.1.dmnd"

cd /mnt/smart/users/fcoutinho/Antarctic_Lagoons/Binning/

# export CHECKM2DB="/mnt/smart/apps/bio/app-databases/checkm2/uniref100.KO.1.dmnd"

# # checkm2 predict --input /mnt/smart/users/fcoutinho/Antarctic_Lagoons/Binning/AL_Bins_Subset/*fa --database_path /mnt/smart/apps/bio/app-databases/checkm2/uniref100.KO.1.dmnd --threads 23 --output-directory CheckM2_Bin_R1_Subset_QC

checkm2 predict --input /mnt/smart/users/fcoutinho/Antarctic_Lagoons/Binning/AL_Bins_R1/*fa --threads 47 --output-directory CheckM2_Bin_R1_QC --database_path /mnt/smart/apps/bio/app-databases/checkm2/uniref100.KO.1.dmnd

###Select MQ and HQ Bins
awk -F "\t" '{ if (($2 >= 50) && ($3 <= 5)) { print } }' /mnt/smart/users/fcoutinho/Antarctic_Lagoons/Binning/CheckM2_Bin_R1_QC/quality_report.tsv | cut -f 1 > AL_MQ_and_HQ_Bins_R1_List.txt

mkdir AL_MAGs_R1

while read i; do cp AL_Bins_R1/${i}.fa AL_MAGs_R1/${i}.fa ; done < AL_MQ_and_HQ_Bins_R1_List.txt

module load gtdbtk/2.4.0
source activate
export PYTHONNOUSERSITE=1

gtdbtk classify_wf --skip_ani_screen --extension fa --genome_dir AL_MAGs_R1/ --pplacer_cpus 12 --cpus 23 --prefix GTDBtk_Output_AL_MAGs_R1 --out_dir GTDBtk_Output_AL_MAGs_R1/ 

###dRep MAGs
module load drep
cd /mnt/smart/users/fcoutinho/Antarctic_Lagoons/Binning/

grep -w -f AL_MQ_and_HQ_Bins_R1_List.txt /mnt/smart/users/fcoutinho/Antarctic_Lagoons/Binning/CheckM2_Bin_R1_QC/quality_report.tsv > AL_MAGs_R1_CheckM_Info.tsv

#create a csv format in the specific format of drep, use regular expressions with sed to add extension to MAG names and convert to csv
echo "genome,completeness,contamination" > AL_MAGs_R1_Info_for_dRep.csv
cut -f 1,2,3 /mnt/smart/users/fcoutinho/Antarctic_Lagoons/Binning/AL_MAGs_R1_CheckM_Info.tsv| grep -v "Name" | sed 's/\t/,/g'  |  sed -E 's/(_Bin\.[0-9]+)/\1.fa/g' >> AL_MAGs_R1_Info_for_dRep.csv

dRep dereplicate /mnt/smart/users/fcoutinho/Antarctic_Lagoons/Binning/dRep_AL_MAGs_R1/ -p 24 -comp 50 -con 5 -l 200000 --genomeInfo AL_MAGs_R1_Info_for_dRep.csv -g /mnt/smart/users/fcoutinho/Antarctic_Lagoons/Binning/AL_MAGs_R1/*fa
