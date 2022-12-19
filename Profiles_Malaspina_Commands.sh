#Generate metadata table
module load malaspina
echo -e "Sample\tStation\tOcean\tDepth\tFilter_Size\tZone\tBatch\tType\tGbp\tMillion_Reads\tComment\tDataset\tFraction" > Malaspina_Profiles_Sample_Info.tsv
mp.sh -e profiles -f FL -t metag -o table >> Malaspina_Profiles_Sample_Info.tsv

#Generate matching samples metadata
cut -f 1 Malaspina_Profiles_Sample_Info.tsv | grep -v "Sample" > List_Malaspina_Profiles_MG_IDs.txt
head -n 1 04_malaspina-160-metag-metadata.tsv > Selected_Malaspina_Profiles_Metadata.tsv
grep -f List_Malaspina_Profiles_MG_IDs.txt -w 04_malaspina-160-metag-metadata.tsv >> Selected_Malaspina_Profiles_Metadata.tsv

#DenfenseFinder
sbatch /mnt/lustre/bio/users/fcoutinho/Scripts/Call_DefenseFinder_Profiles_Malaspina_MAGs_Array.slurm

time python3 /mnt/lustre/bio/users/fcoutinho/Scripts/parse_defense_finder_output.py --output Profiles_Malaspina_MAGs_DenfenseFinder_Ordered_Results.tsv --input /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/Defense_Finder_Results/ &

#Generate assembly dictionary table fields
#R1
mp.sh -e profiles -f FL -t metag -o clean | grep "MP1527\|MP1682\|MP2235\|MP2817\|MP0323\|MP0888\|MP1421\|MP1517\|MP1857\|MP2243\|MP2819\|MP1154\|MP1411\|MP1521\|MP1409\|MP1523\|MP1849\|MP1851\|MP2231\|MP1417\|MP1672\|MP1684\|MP1855\|MP1529\|MP2241\|MP2815\|MP1674\|MP1676" | grep "pair1"

#R2
mp.sh -e profiles -f FL -t metag -o clean | grep "MP1527\|MP1682\|MP2235\|MP2817\|MP0323\|MP0888\|MP1421\|MP1517\|MP1857\|MP2243\|MP2819\|MP1154\|MP1411\|MP1521\|MP1409\|MP1523\|MP1849\|MP1851\|MP2231\|MP1417\|MP1672\|MP1684\|MP1855\|MP1529\|MP2241\|MP2815\|MP1674\|MP1676" | grep "pair2"

#Sample ID and Group
mp.sh -e profiles -f FL -t metag -o clean | grep "MP1527\|MP1682\|MP2235\|MP2817\|MP0323\|MP0888\|MP1421\|MP1517\|MP1857\|MP2243\|MP2819\|MP1154\|MP1411\|MP1521\|MP1409\|MP1523\|MP1849\|MP1851\|MP2231\|MP1417\|MP1672\|MP1684\|MP1855\|MP1529\|MP2241\|MP2815\|MP1674\|MP1676" | grep "pair1" | cut -d "/" -f 16 | cut -d "." -f 1


#Assemble with QueroBins
sbatch /mnt/lustre/bio/users/fcoutinho/Scripts/Call_QueroBins_Assembly_Malaspina_Profiles2_Part1.slurm
sbatch /mnt/lustre/bio/users/fcoutinho/Scripts/Call_QueroBins_Assembly_Malaspina_Profiles2_Part2.slurm


#Filter and Rename Assembled Scaffolds
python3 /mnt/lustre/bio/users/fcoutinho/Scripts/QueroBins2.py --min_scaffold_length 500 --assemblies_files ../Raw_Scaffolds/*fasta

#Annotate with CAT through QueroBins2
sbatch /mnt/lustre/bio/users/fcoutinho/Scripts/Call_QueroBins_CAT_Malaspina_Profiles.slurm

#Refine through QueroBins2
mv /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/CAT_Results/*classification.txt .
mv /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Metabat_Bins_Round_1/Bin_Filtered_Renamed_S*fa

python3 /mnt/lustre/bio/users/fcoutinho/Scripts/QueroBins2.py --refine_bins True --refinement_tax_level class --call_checkm True --parse_only True --assemblies /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Filtered_Renamed_S*fasta &

#VAMB Binning
module load metabat

jgi_summarize_bam_contig_depths --outputDepth Profiles_Malaspina_JGI_Depth.tsv *bam

#Corrected IDs in R generating Corrected_Profiles_Malaspina_JGI_Depth.tsv

module load vamb
source activate

concatenate.py --nozip -m 2500 Profiles_Malaspina_Scaffolds_Filtered_2500_bp.fasta /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Filtered_Renamed_Scafolds/*fasta

vamb --outdir /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/VAMB_Binning/VAMB_Bins_Round_1/ --fasta /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/VAMB_Binning/Profiles_Malaspina_Scaffolds_Filtered_2500_bp.fasta --jgi /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/VAMB_Binning/Corrected_Profiles_Malaspina_JGI_Depth.tsv -o C --minfasta 200000


#Metabat binning
sbatch /mnt/lustre/bio/users/fcoutinho/Scripts/Call_QueroBins_Metabat_Malaspina_Profiles.slurm

#MaxBin binning
module load perl
module load maxbin

i=2
while [ $i -ne 154 ];do  i=$(($i+2));  echo "$i"; sid=$(head -n 1 /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/VAMB_Binning/Corrected_Profiles_Malaspina_JGI_Depth.tsv | cut -f $i | cut -d "x" -f 1); echo "$sid"; cut -f 1,$i /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/VAMB_Binning/Corrected_Profiles_Malaspina_JGI_Depth.tsv | grep -v "contigName" > Abundance_$sid.tsv ; done;

ls | grep "Abundance_" > List_Abd_Files.txt

#Merge bin info
python /mnt/lustre/bio/users/fcoutinho/Scripts/Merge_Matrixes.py --matrixes /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Profiles_Malaspina_Metabat_MAGs_Round_1_GTDBtk_Output/Profiles_Malaspina_Metabat_MAGs_Round_1_GTDBtk.bac120.summary.tsv /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Profiles_Malaspina_Metabat_MAGs_Round_1_GTDBtk_Output/Profiles_Malaspina_Metabat_MAGs_Round_1_GTDBtk.ar122.summary.tsv --out Info_1.tsv --axis 0 --index user_genome

python /mnt/lustre/bio/users/fcoutinho/Scripts/Merge_Matrixes.py --matrixes /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Info_1.tsv /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/CheckM_Info_Metabat_Bins_Round_1.tsv /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Assembly_Stats_Metabat_MAGs_Round_1.tsv /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/MAG_to_Sample.tsv --axis 1 --index user_genome --out Info_2.tsv 

python /mnt/lustre/bio/users/fcoutinho/Scripts/Merge_Matrixes.py --matrixes /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Info_2.tsv /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Info_3.tsv --axis 1 --index user_genome --out Info_4.tsv

###GTDB Tree Decoration
#Bacteria Tree
perl /mnt/lustre/bio/users/fcoutinho/Scripts/Decorate_iTOL.pl --input_info /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/Phylogenomics/Info_GTDB_R202_Selected_Bac_Taxa.tsv --color_category Phylum --color_palette sasha --min_trait_count 50

perl /mnt/lustre/bio/users/fcoutinho/Scripts/Node_Colapse.pl --tree New_Labels_Profiles_Malaspina_MetaBat_Redo_MAGs_GTDBtk_Output.bac120.classify.newick --seq_info_file Info_GTDB_R202_Bac_Genomes.tsv --category Phylum --collapse_full ; mv Nodes_To_Colapse.txt Profiles_Malaspina_Bac_Nodes_To_Colapse.txt ; mv Node_New_IDS_Phylum.txt Profiles_Malaspina_Bac_Node_New_IDS_Phylum.txt ; grep "Node_\|GB_GCA_" Profiles_Malaspina_Bac_Nodes_To_Colapse.txt > List_Collapsed_Nodes_Profiles_Malaspina_Bac.txt ;

perl /mnt/lustre/bio/users/fcoutinho/Scripts/Decorate_iTOL.pl --color_palette zone --tree New_Labels_Profiles_Malaspina_MetaBat_Redo_MAGs_GTDBtk_Output.bac120.classify.newick --input_info ../Info_All_MAGs_Profiles_Malaspina_Redo.tsv --barplot_category Zone ;grep -v "Node_" Tree_Barplot_Zone.txt > Bact_Zone_Barplots_Selected.txt ; grep -f List_Collapsed_Nodes_Profiles_Malaspina_Bac.txt -w Tree_Barplot_Zone.txt >> Bact_Zone_Barplots_Selected.txt ;

#Archaea tree
perl /mnt/lustre/bio/users/fcoutinho/Scripts/Decorate_iTOL.pl --input_info /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/Phylogenomics/Info_GTDB_R202_Arc_Genomes.tsv --color_category Phylum --color_palette rainbow20 --min_trait_count 1 --output_info_file Decoration_Info_Profiles_Malaspina_Arch.tsv; mv Tree_Colors_By_Index_Number_Phylum-Palette_rainbow20.txt Profiles_malaspina_Arch_Tree_Colors_By_Index_Number_Phylum-Palette_rainbow20.txt ;

perl /mnt/lustre/bio/users/fcoutinho/Scripts/Give_Node_IDs.pl --input No_Labels_Profiles_Malaspina_MetaBat_Redo_MAGs_GTDBtk_Output.ar122.classify.newick --output New_Labels_Profiles_Malaspina_MetaBat_Redo_MAGs_GTDBtk_Output.ar122.classify.newick

perl /mnt/lustre/bio/users/fcoutinho/Scripts/Node_Colapse.pl --tree New_Labels_Profiles_Malaspina_MetaBat_Redo_MAGs_GTDBtk_Output.ar122.classify.newick --seq_info_file Info_GTDB_R202_Arc_Genomes.tsv --category Phylum --collapse_full ; mv Nodes_To_Colapse.txt Profiles_Malaspina_Arch_Nodes_To_Colapse.txt ; mv Node_New_IDS_Phylum.txt Profiles_Malaspina_Arch_Node_New_IDS_Phylum.txt ;

perl /mnt/lustre/bio/users/fcoutinho/Scripts/Decorate_iTOL.pl --color_palette zone --tree New_Labels_Profiles_Malaspina_MetaBat_Redo_MAGs_GTDBtk_Output.ar122.classify.newick --input_info ../Info_All_MAGs_Profiles_Malaspina_Redo.tsv --barplot_category Zone ; grep -v "Node_" Tree_Barplot_Zone.txt > Arch_Zone_Barplots_Selected.txt ; grep -f List_Collapsed_Arch_Nodes.txt -w Tree_Barplot_Zone.txt >> Arch_Zone_Barplots_Selected.txt ;

#Planctomycetota subtree
grep -w "Genome\|Planctomycetota" /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/Phylogenomics/Info_GTDB_R202_Bac_Genomes.tsv > Info_GTDB_R202_Planctomycetota_Genomes.tsv

grep -w "Genome\|Planctomycetes" /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/Phylogenomics/Info_GTDB_R202_Bac_Genomes.tsv > Info_GTDB_R202_Planctomycetes_Genomes.tsv

perl /mnt/lustre/bio/users/fcoutinho/Scripts/Decorate_iTOL.pl --input_info /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/Phylogenomics/Info_GTDB_R202_Planctomycetes_Genomes.tsv --color_category Order --color_palette sasha --min_trait_count 1

#Cyanobacteria subtree
grep -w "Genome\|Cyanobacteria" /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/Phylogenomics/Info_GTDB_R202_Bac_Genomes.tsv > Info_GTDB_R202_Cyanobacteria_Genomes.tsv

perl /mnt/lustre/bio/users/fcoutinho/Scripts/Decorate_iTOL.pl --input_info /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/Phylogenomics/Info_GTDB_R202_Cyanobacteria_Genomes.tsv --color_category Class --color_palette basic --min_trait_count 1

#Metabolic predictions manual curation
###Manual Curation

#IDS=$(grep -w $HOST /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Viruses/Profiles_Malaspina_Viruses_Full_Info.tsv | cut -f 1)
TAXON="Cyanobacteriia" #"SAR324"
IDS=$( grep -w $TAXON /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/Info_All_MAGs_Profiles_Malaspina_Redo.tsv | cut -f 1 )

touch List_MAGs_$TAXON.txt
for i in $IDS ; do echo $i >> List_MAGs_$TAXON.txt; done ;
head -n 1 Long_Format_No_Absent_Profiles_Malaspina_MAGs_METABOLIC_result_worksheet4.tsv > Long_Format_No_Absent_METABOLIC_result_worksheet4_$TAXON.tsv
grep -f List_MAGs_$TAXON.txt Long_Format_No_Absent_Profiles_Malaspina_MAGs_METABOLIC_result_worksheet4.tsv >> Long_Format_No_Absent_METABOLIC_result_worksheet4_$TAXON.tsv &


touch List_MAG_CDS_Prefix_$TAXON.txt
for i in $IDS ; do j="${i}_"; echo $j >> List_MAG_CDS_Prefix_$TAXON.txt; done ;
head -n 1 /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/Genes_and_CDS_All_MAGs/Split_All_CDS/Profiles_Malaspina_CDS_Info.tsv > CDS_Info_$TAXON.tsv
grep -f List_MAG_CDS_Prefix_$TAXON.txt /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/Genes_and_CDS_All_MAGs/Split_All_CDS/Profiles_Malaspina_CDS_Info.tsv >> CDS_Info_$TAXON.tsv &


TAXON="Cyanobacteriia" #"SAR324"
IDS=$( grep -w $TAXON /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/Info_All_MAGs_Profiles_Malaspina_Redo.tsv | cut -f 1 )
MOID="M00165"

grep "M00165+2" Long_Format_No_Absent_Profiles_Malaspina_MAGs_METABOLIC_result_worksheet4.tsv | grep -w "Cyanobacteriia" | cut -f 1 | sort | uniq > List_MAGs_Cyanobacteriia_Containing_Module_M00165+2.txt
touch List_Prefix_MAGs_Cyanobacteriia_Containing_Module_M00165+2.txt
while read i; do j="${i}_"; echo $j >> List_Prefix_MAGs_Cyanobacteriia_Containing_Module_M00165+2.txt; done  < List_MAGs_Cyanobacteriia_Containing_Module_M00165+2.txt
head -n 1 /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/Genes_and_CDS_All_MAGs/Split_All_CDS/Profiles_Malaspina_CDS_Info.tsv > CDS_Info_MAGs_Cyanobacteriia_Missing_Module_M00165+2.tsv
grep -f List_MAG_CDS_Prefix_$TAXON.txt /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/Genes_and_CDS_All_MAGs/Split_All_CDS/Profiles_Malaspina_CDS_Info.tsv | grep -v -f List_Prefix_MAGs_Cyanobacteriia_Containing_Module_M00165+2.txt >> CDS_Info_MAGs_Cyanobacteriia_Missing_Module_M00165+2.tsv

grep -w "UniRef_Best_Hit_ID\|K01601\|K01602" CDS_Info_MAGs_Cyanobacteriia_Missing_Module_M00165+2.tsv > short.tsv

MAGID=""
head -n 1 /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/Metabolic_Preds/Long_Format_No_Absent_Profiles_Malaspina_MAGs_METABOLIC_result_worksheet4.tsv > Long_Format_No_Absent_Profiles_Malaspina_METABOLIC_result_worksheet4_MAG_$MAGID.tsv
grep $MAGID /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/Metabolic_Preds/Long_Format_No_Absent_Profiles_Malaspina_MAGs_METABOLIC_result_worksheet4.tsv >> Long_Format_No_Absent_Profiles_Malaspina_METABOLIC_result_worksheet4_MAG_$MAGID.tsv  &


MAGID="Bin_S41.17"
MAGPREFIX="${MAGID}_"
head -n 1 /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/Genes_and_CDS_All_MAGs/Split_All_CDS/Profiles_Malaspina_CDS_Info.tsv > CDS_Info_MAG_$MAGID.tsv
grep $MAGPREFIX /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/Genes_and_CDS_All_MAGs/Split_All_CDS/Profiles_Malaspina_CDS_Info.tsv >> CDS_Info_MAG_$MAGID.tsv &


MOID="M00165"
TAXON="Cyanobacteriia"
head -n 1 /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/Metabolic_Preds/Long_Format_No_Absent_Profiles_Malaspina_MAGs_METABOLIC_result_worksheet4.tsv > Long_Format_No_Absent_Profiles_Malaspina_METABOLIC_result_worksheet4_MAG_$MAGID.$TAXON.tsv
grep $MOID /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/Metabolic_Preds/Long_Format_No_Absent_Profiles_Malaspina_MAGs_METABOLIC_result_worksheet4.tsv | grep $TAXON >> Long_Format_No_Absent_Profiles_Malaspina_METABOLIC_result_worksheet4_MAG_$MAGID.$TAXON.tsv 

head -n 1 /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/Metabolic_Preds/Profiles_Malaspina_MAGs_METABOLIC_Results/METABOLIC_result_each_spreadsheet/METABOLIC_result_worksheet4.tsv | cut -f 1-4 > Module_Info_$MOID.tsv
grep $MOID /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/Metabolic_Preds/Profiles_Malaspina_MAGs_METABOLIC_Results/METABOLIC_result_each_spreadsheet/METABOLIC_result_worksheet4.tsv | cut -f 1-4 >> Module_Info_$MOID.tsv &

#Module mean/max completeness Barplots
INFILE="/mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/Metabolic_Preds/Profiles_Malaspina_MAGs_METABOLIC_Module_Completeness_From_worksheet4.tsv"
TAXLEVEL="Phylum" ##  "Domain" #
TAXNAME="Marinisomatota" #"Bacteria" # "Actinobacteriota" #"Bacteroidota" #  "Archaea" #Chloroflexota"
GROUPLEVEL="Class" #"Phylum" "Class" "Order"
MINCOMP=0
BYZONE=FALSE

Rscript /mnt/lustre/bio/users/fcoutinho/Scripts/make_modulue_compl_barplot.R $INFILE $TAXLEVEL $TAXNAME $GROUPLEVEL $MINCOMP $BYZONE

#####Viruses
####Viruses genomes
cat /mnt/lustre/scratch/elopez/4_checkV_output/proviruses.fna /mnt/lustre/scratch/elopez/4_checkV_output/viruses.fna > Profiles_Malaspina_Viruses_CheckV_Trimmed_Genomes.fasta

python3 /mnt/lustre/bio/users/fcoutinho/Scripts/virathon_v0.1.py --filter True --min_length 1000 --genome_files Profiles_Malaspina_Viruses_CheckV_Trimmed_Genomes.fasta --cds Profiles_Malaspina_Viruses_CheckV_Trimmed_Genomes.faa --index True


###Viruses abundance
perl /mnt/lustre/bio/users/fcoutinho/Scripts/Index_Host_Taxonomy.pl --seq_info /mnt/lustre/scratch/elopez/RaFAH_malaspina__Seq_Info_Prediction.tsv --out Profiles_Malaspina_Virus_Host_Tax_Info.tsv

perl /mnt/lustre/bio/users/fcoutinho/Scripts/Summarize_Abundance_By_Traits.pl --info /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Viruses/Corrected_Abundance/Profiles_Malaspina_Virus_Host_Tax_Info.tsv --matrix_file /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Viruses/Corrected_Abundance/Profiles_Malaspina_Virus_RPKM.tsv --min_filter_val 0.14


###Match Viruses to MAGs for Prophage Counts
python3 /mnt/lustre/bio/users/fcoutinho/Scripts/Virathon.py --genome /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Prok_Viruses/Profiles_Malaspina_Viruses_CheckV_Trimmed_Genomes.fasta --index True --info_output Info_Profiles_Malaspina_Viruses_CheckV_Trimmed_Genomes.tsv

cut -f 6 /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/Profiles_Malaspina_All_Scaffolds_Info_Redo.tsv > ids
cut -f 2 /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/Profiles_Malaspina_All_Scaffolds_Info_Redo.tsv > samples
paste ids samples > /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/Prophage_Counts/Profiles_Malaspina_MAG_Scaffolds_Info_Redo_ID+Sample.tsv
grep -v -w "NA" Profiles_Malaspina_MAG_Scaffolds_Info_Redo_ID+Sample.tsv > No_NA_Profiles_Malaspina_MAG_Scaffolds_Info_Redo_ID+Sample.tsv

module load python blast

python3 /mnt/lustre/bio/users/fcoutinho/Scripts/Link_Viral_Seqs_to_MAGs.py --viral_genomes /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Prok_Viruses/Profiles_Malaspina_Viruses_CheckV_Trimmed_Genomes.fasta --mags_dir /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/MAGs_Original_and_Refined/ --mags_ext fa --threads 12 --scaffold_specific True --sample_specific True --host_scaffold_info /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/Prophage_Counts/No_NA_Profiles_Malaspina_MAG_Scaffolds_Info_Redo_ID+Sample.tsv --virus_scaffold_info /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/Prophage_Counts/Info_Profiles_Malaspina_Viruses_CheckV_Trimmed_Genomes+Sample.tsv --parse_only True &

###Match Viruses to MAGs for host pred
sbatch /mnt/lustre/bio/users/fcoutinho/Scripts/Match_Vir_MAGs_Profiles_Malaspina.slurm
 
python3 /mnt/lustre/bio/users/fcoutinho/Scripts/Link_Viral_Seqs_to_MAGs.py --viral_genomes Profiles_Malaspina_Viruses_CheckV_Trimmed_Genomes.fasta --mags_dir /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Metabat_MAGs_Round_1/ --mags_ext fasta --threads 48 --parse_only True --scaffold_specific True --remove_virus True &

head -n 1 PHIST_Predictions.csv > Positive_PHIST_Predictions.csv
grep "No_Virus_Bin_Filtered" PHIST_Predictions.csv >> Positive_PHIST_Predictions.csv


###Viruses PHIST Host Pred
mkdir Profiles_Malaspina_Viral_Genomes
cd Profiles_Malaspina_Viral_Genomes
python3 /mnt/lustre/bio/users/fcoutinho/Scripts/Explode_Fasta.py --assembly /mnt/lustre/scratch/elopez/4_checkV_output/all_pro_viruses_checkV.fna

python3 /mnt/lustre/bio/users/fcoutinho/PHIST/phist.py -t 24 /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Viruses/Profiles_Malaspina_Viral_Genomes/ /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Metabat_MAGs_Round_1/ PHIST_Output/PHIST_Kmers.csv /PHIST_Output/PHIST_Predictions.csv &


###Manual Curation
HOST="Marinisomatota"
IDS=$(grep -w $HOST /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Viruses/Profiles_Malaspina_Viruses_Full_Info.tsv | cut -f 1)
touch List_Host_$HOST.txt
for i in $IDS ; do j="${i}_"; echo $j >> List_Host_$HOST.txt; done ;
head -n 1 /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Viruses/CDS_Info.tsv > CDS_Info_$HOST.tsv
grep -f List_Host_$HOST.txt /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Viruses/CDS_Info.tsv >> CDS_Info_$HOST.tsv

grep -w "Virus_ID\|${HOST}" /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Viruses/Profiles_Malaspina_Viruses_Full_Info.tsv > Seq_Info_Host_$HOST.tsv

KO="K01698"
grep -w "Virus_ID\|${KO}" /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Viruses/Profiles_Malaspina_Viruses_Full_Info.tsv > Seq_Info_KO_$KO.tsv

SCAFFOLD="MP2241_2_"
head -n 1 /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Viruses/CDS_Info.tsv > CDS_Info_$SCAFFOLD.tsv
grep $SCAFFOLD /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Viruses/CDS_Info.tsv >> CDS_Info_$SCAFFOLD.tsv
