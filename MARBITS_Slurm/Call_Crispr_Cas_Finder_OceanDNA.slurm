#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=07-00:00:00
#SBATCH --job-name=OceanDNA_CCfinder
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=15G
#SBATCH --output=/mnt/smart/scratch/vir/felipe/Job_Logs/jobLog_%J.out
#SBATCH --error=/mnt/smart/scratch/vir/felipe/Job_Logs/jobLog_%J.err
#SBATCH --mail-type=All
#SBATCH --mail-user=fhernandes@icm.csic.es

##############################################
cd /mnt/smart/scratch/vir/felipe/Databases/OceanDNA/Individual_Genomes/

#split the sequences into individual files
module load python

#the equal sign and qutes are necessary so the arparser understands -scaffold_ is to be assigned to id_split_sep
# python3 /mnt/smart/scratch/vir/felipe/Repos/ICM_Code/Get_Seqs.py --input_seq /mnt/smart/scratch/vir/felipe/Databases/OceanDNA/OceanDNA_All_Species_Rep_MAGs_Scaffolds.fasta --fetch_all True --id_split_pos 0 --group True --id_split_sep="_"

# run cctyper, delete genome files as they are generated to save space
module load miniconda
source /home/apps/miniconda/4.5.x/etc/profile.d/conda.sh

conda activate /mnt/smart/scratch/vir/felipe/envs/cctyper

for genome in *fasta ; do cctyper --no_plot --threads 47 $genome cctyper_output_${genome%.fasta} ; rm -f $genome ;  done ; 

# python3 /mnt/smart/scratch/vir/felipe/Repos/ICM_Code/Merge_Matrixes.py --matrixes /mnt/smart/scratch/vir/felipe/Databases/Malaspina/DeepMAGs/cctyper_output_mp-deep_mag-05*/*putative.tab  --output CCfinder_DeepMAGs_output.tsv --out_index Contig

# head -n 1 /mnt/smart/scratch/vir/felipe/Databases/Malaspina/DeepMAGs/cctyper_output_mp-deep_mag-0002/cas_operons_putative.tab > CCfinder_DeepMAGs_output.tsv

cat /mnt/smart/scratch/vir/felipe/Databases/OceanDNA/Individual_Genomes/cctyper_out*/CRISPR_Cas.tab | head -n 1 > /mnt/smart/scratch/vir/felipe/Databases/OceanDNA/CCfinder_OceanDNA_CRISPR_Cas_output.tsv

grep -v "Operon" cat /mnt/smart/scratch/vir/felipe/Databases/OceanDNA/Individual_Genomes/cctyper_out*/CRISPR_Cas.tab >> /mnt/smart/scratch/vir/felipe/Databases/OceanDNA/CCfinder_OceanDNA_CRISPR_Cas_output.tsv

cat /mnt/smart/scratch/vir/felipe/Databases/OceanDNA/Individual_Genomes/cctyper_out*/cas_operons.tab | head -n 1 > /mnt/smart/scratch/vir/felipe/Databases/OceanDNA/CCfinder_cas_operons.tsv

grep -v "Operon" /mnt/smart/scratch/vir/felipe/Databases/OceanDNA/Individual_Genomes/cctyper_out*/cas_operons.tab  >> /mnt/smart/scratch/vir/felipe/Databases/OceanDNA/CCfinder_cas_operons.tsv

#make a list of all proteins assigned to cas operons
cut -f 1,11 /mnt/smart/scratch/vir/felipe/Databases/OceanDNA/CCfinder_cas_operons.tsv | sed -E s'/(.*:)|( )|\[|]//g' | awk '{num=split($2,array,",");for(i=1;i<=num;i++){print $1,"_",array[i]}}' | sed -E s'/\s//g' | grep -v "Contig_Positions" > /mnt/smart/scratch/vir/felipe/Databases/OceanDNA/List_Ocean_DNA_cas_CDS.txt

echo "CDS,Class" > /mnt/smart/scratch/vir/felipe/Databases/OceanDNA/Ocean_DNA_cas_CDS_with_classification.csv

cut -f 1,10,11 /mnt/smart/scratch/vir/felipe/Databases/OceanDNA/CCfinder_cas_operons.tsv | sed -E s'/(.*:)|( )|\[|]//g' | awk '{num=split($2,class,",");numb=split($3,posits,",");for(i=1;i<=num;i++){print $1,"_",posits[i],",",class[i]}}' | sed -E s"/(\s)|(')//g"  | tail -n +1 > /mnt/smart/scratch/vir/felipe/Databases/OceanDNA/Ocean_DNA_cas_CDS_with_classification.csv

#Must call prodigal again because IDs in /mnt/smart/scratch/vir/felipe/Databases/OceanDNA/Sequences/OceanDNA_Species_Rep_MAGs_PEGs.faa do not match prodigal scheme and cctyper was not run with the option to keep temporary files

cd /mnt/smart/scratch/vir/felipe/Databases/OceanDNA/Prodigal/

module load python

#the equal sign and qutes are necessary so the arparser understands -scaffold_ is to be assigned to id_split_sep
python3 /mnt/smart/scratch/vir/felipe/Repos/ICM_Code/Get_Seqs.py --input_seq /mnt/smart/scratch/vir/felipe/Databases/OceanDNA/Sequences/OceanDNA_All_Species_Rep_MAGs_Scaffolds.fasta --fetch_all True --id_split_pos 0 --group True --id_split_sep="_"

# run cctyper, delete genome files as they are generated to save space
module load miniconda
source /home/apps/miniconda/4.5.x/etc/profile.d/conda.sh

conda activate /mnt/smart/scratch/vir/felipe/envs/cctyper

for genome in *fasta ; do prodigal -p single -i $genome -a ${genome%fasta}.faa -d ${genome%fasta}.ffn -f gff -i  $genome -o ${genome%fasta}.gff;  done ;

#Make a DB of all cas proteins from OceanDna genomes
cat /mnt/smart/scratch/vir/felipe/Databases/OceanDNA/Prodigal/*faa > /mnt/smart/scratch/vir/felipe/Databases/OceanDNA/OceanDNA_All_Species_Rep_MAGs_CDS.faa

#make subset of cas protein CDS
python3 /mnt/smart/scratch/vir/felipe/Repos/ICM_Code/Get_Seqs.py --input_seq /mnt/smart/scratch/vir/felipe/Databases/OceanDNA/OceanDNA_All_Species_Rep_MAGs_CDS.faa --list /mnt/smart/scratch/vir/felipe/Databases/OceanDNA/List_Ocean_DNA_cas_CDS.txt --matched_out /mnt/smart/scratch/vir/felipe/Databases/OceanDNA/OceanDNA_All_Species_Rep_MAGs_cas_only_CDS.faa

rm -f /mnt/smart/scratch/vir/felipe/Databases/OceanDNA/Unmatched_Sequences.fasta

###Build mmseqs DB
mkdir /mnt/smart/scratch/vir/felipe/Databases/OceanDNA/cas_Abundance
cd /mnt/smart/scratch/vir/felipe/Databases/OceanDNA/cas_Abundance/

module load mmseqs2

mmseqs createdb /mnt/smart/scratch/vir/felipe/Databases/OceanDNA/OceanDNA_All_Species_Rep_MAGs_cas_only_CDS.faa DB_OceanDNA_All_Species_Rep_MAGs_cas_only_CDS 

###Make list of R1 files to be mapped
#add profiles reads
grep "fastq" /mnt/smart/scratch/vir/felipe/Metadata/Malaspina/Profiles_Malaspina_76_FL_Metagenome_Info.tsv | cut -f 2 > /mnt/smart/scratch/vir/felipe/Databases/OceanDNA/cas_Abundance/Profiles_Malaspina_76_FL_Metagenome_R1_List.txt

#add deep reads
grep "fastq" /mnt/smart/scratch/vir/felipe/Databases/Malaspina/Malaspina_Deep_Read_Info.tsv | cut -f 2 >> /mnt/smart/scratch/vir/felipe/Databases/OceanDNA/cas_Abundance/Malaspina_Deep_R1_List.txt

#add kaust reads???

#iterate over all MLP profiles files and query them against the  cas1 DB using mmseqs2
sbatch /mnt/smart/scratch/vir/felipe/Repos/ICM_Code/MARBITS_Slurm/Calc_OceanDNA_cas_abd_MaLasPina_MGs.slurm


