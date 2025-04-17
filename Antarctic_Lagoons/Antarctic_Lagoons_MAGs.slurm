#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=02-00:00:00
#SBATCH --job-name=Antartic_Lagoons
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=50G
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
module load python/3.8.5
module load bowtie2/2.3.4.1
module load samtools

cd /mnt/smart/users/fcoutinho/Antarctic_Lagoons/Binning/

#Notice command below max reads is set to 1M to speed up for a trial run, and bowtie_k set to 100 becaseu the contigs detabase includes all contigs from all 20 samples and redundancy in read mapping is expected
# python3 /mnt/smart/users/fcoutinho/Repos/virathon/Virathon.py --genome_files /mnt/smart/users/fcoutinho/Antarctic_Lagoons/Binning/Antarctic_Lagoons_All_Filtered_Scaffolds_2.5kbp.fasta --abundance_table True --abundance_rpkm True --raw_read_table /mnt/smart/users/fcoutinho/Antarctic_Lagoons/Metadata/WGS_Samples.tsv --bowtie_k 100 --abundance_max_reads 1000000 --threads 47
#using the command above, many scaffolds have abundance of 0 in their own sample. Bins are not formed using the abundance file derived from it either

python3 /mnt/smart/users/fcoutinho/Repos/virathon/Virathon.py --genome_files /mnt/smart/users/fcoutinho/Antarctic_Lagoons/Binning/Antarctic_Lagoons_All_Filtered_Scaffolds_2.5kbp.fasta --bowtiedb /mnt/smart/users/fcoutinho/Antarctic_Lagoons/Binning/All_Genomic --abundance_table True --abundance_rpkm True --raw_read_table /mnt/smart/users/fcoutinho/Antarctic_Lagoons/Metadata/WGS_Samples.tsv --bowtie_k 100 --abundance_max_reads 15000000 --threads 23

module load metabat
conda activate /mnt/smart/apps/condaEnvs/metabat-2.15

#/mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/Software/metabat/bin/jgi_summarize_bam_contig_depths

jgi_summarize_bam_contig_depths --outputDepth JGI_Depths_Antarctic_Lagoons_All_Filtered_Scaffolds_2.5kbp.tsv *bam

# module load python

#Split scaffolds per sample

python3 /mnt/smart/users/fcoutinho/Repos/ICM_Code/Get_Seqs.py --input_sequences /mnt/smart/users/fcoutinho/Antarctic_Lagoons/Binning/Antarctic_Lagoons_All_Filtered_Scaffolds_2.5kbp.fasta --id_split_sep _ --id_split_pos 0 --fetch_all True --group True

#/mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/Software/metabat/bin/metabat2

for i in C*.fasta; do smp=$(echo "$i" | sed 's/.fasta//g') ; echo $smp; metabat2 -i $i --abdFile  /mnt/smart/users/fcoutinho/Antarctic_Lagoons/Binning/JGI_Depths_Antarctic_Lagoons_All_Filtered_Scaffolds_2.5kbp.tsv -s 200000 -t 23 --outFile AL_R1_${smp}_Bin; done

conda deactivate

# rm -f C*fasta

mkdir /mnt/smart/users/fcoutinho/Antarctic_Lagoons/Binning/AL_Bins_R1

mv AL_R1_*Bin* /mnt/smart/users/fcoutinho/Antarctic_Lagoons/Binning/AL_Bins_R1/

cd /mnt/smart/users/fcoutinho/Antarctic_Lagoons/Binning/AL_Bins_R1/

module load perl

perl /mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/repos/ICM_Code/assembly_stats_table.pl fa > AL_R1_Metabat_Bins_Assembly_Stats.tsv

###CheckM2 Bin QC
# module load checkm2/1.0.1
# conda activate /mnt/smart/apps/condaEnvs/checkm2

# cd /mnt/smart/users/fcoutinho/Antarctic_Lagoons/Binning/AL_Bins_R1/

# checkm2 predict --input /mnt/smart/users/fcoutinho/OctoMicro/WGS_Assemblies/Metabat_Bins/*fa --threads 47 --output-directory CheckM2_Bin_QC

# cd  Metabat_Bins_R2/
# checkm2 predict --input *fa --threads 12 --output-directory CheckM2_Bin_QC

###CheckM1 R3 Bin QC
# module load cesga/2020 checkm/1.2.2-python-3.9.9

# cd /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/OCTOMicro/

# checkm lineage_wf -x fa -t 23 --pplacer_threads 12 --tab_table --file OctoMicro_Bins_R2_Checkm_Results.tsv Metabat_Bins_R3/ Checkm_Output_OctoMicro_Bins_R3/

# ###Select MQ and HQ Bins
# awk -F "\t" '{ if (($12 >= 50) && ($13 <= 5)) { print } }' OctoMicro_Bins_R2_Checkm_Results.tsv | cut -f 1 > OctoMicro_MQ_and_HQ_Bins_R3_List.txt

# mkdir OctoMicro_MAGs_R3

# while read i; do cp Metabat_Bins_R3/${i}.fa OctoMicro_MAGs_R3/${i}.fa ; done < OctoMicro_MQ_and_HQ_Bins_R3_List.txt

###GTDBtk
# module purge
# module load cesga/2020 gcccore/system gtdb-tk/2.3.2-python-3.9.9

# gtdbtk classify_wf --skip_ani_screen --extension fa --genome_dir OctoMicro_MAGs_R3/ --pplacer_cpus 12 --cpus 23 --prefix GTDBtk_Output_OctoMicro_MAGs_R3 --out_dir GTDBtk_Output_OctoMicro_MAGs_R3/ 


# ###Prokka mags R3
# module load cesga/2020  gcccore/system  prokka/1.14.5

# cd /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/OCTOMicro/OctoMicro_MAGs_R3/

# for file in *.fa; do prokka --kingdom Bacteria --cpus 4 --noanno --outdir "Prokka_${file%.fa}" --prefix "${file%.fa}" --locustag "${file%.fa}" "$file" ; done ;

# mkdir /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/OCTOMicro/Prokka_Output_OctoMicro_MAGs_R3

# mv Prokka_* /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/OCTOMicro/Prokka_Output_OctoMicro_MAGs_R3/

