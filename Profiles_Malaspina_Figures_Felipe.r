###Libraries
library("vegan")
library(ggsignif)
require(maps)
library(data.table)
library(dplyr)
library(tidyr)

library(ggplot2)
library(RColorBrewer)
library(reshape2)
library(ggpubr, lib.loc="/mnt/lustre/bio/users/fcoutinho/Rlibs/")

###Aesthetics
watermass_cols<-brewer.pal(9,"PuBuGn")[c(3,6,8)]
names(watermass_cols)<-c("CDW","NADW","WSDW")

depth_col_pal<-brewer.pal(9,"Blues")[4:9]
depth_col_grad<-rev(colorRampPalette(depth_col_pal)(n=100))

o2_col_pal<-brewer.pal(11,"Spectral")
o2_col_grad<-rev(colorRampPalette(o2_col_pal)(n=299))


phylum_order<-c("Acidobacteriota","Actinobacteriota","Alphaproteobacteria","Gammaproteobacteria","Bacteroidota","Bdellovibrionota","Chlamydiota","Chloroflexota","Cyanobacteria","Desulfobacterota","Desulfobacterota_D","Eremiobacterota","Gemmatimonadota","Margulisbacteria","Marinisomatota","Myxococcota","Nitrospinota","Patescibacteria","Planctomycetota","Poribacteria","SAR324","UBP7","Verrucomicrobiota","Halobacteriota","Hydrothermarchaeota","Nanoarchaeota","Thermoplasmatota","Thermoproteota")

phylum_coloring<-c(brewer.pal(9,"Set1"),brewer.pal(4,"Pastel1"),brewer.pal(8,"Dark2")[8],brewer.pal(5,"Set2")[5],brewer.pal(8,"Accent"),brewer.pal(6,"Dark2"))
names(phylum_coloring)<-c("Actinobacteriota","Alphaproteobacteria","Cyanobacteria","Crenarchaeota","Bacteroidota","Nanoarchaeota","Verrucomicrobiota","Marinisomatota","Myxococcota","Patescibacteria","Gammaproteobacteria","Acidobacteriota","UBP7","Poribacteria","Eremiobacterota","Margulisbacteria","Gemmatimonadota","Desulfobacterota","Nitrospinota","Desulfobacterota_D","Bdellovibrionota","Hydrothermarchaeota","Chlamydiota","SAR324","Halobacteriota","Chloroflexota","Thermoplasmatota","Planctomycetota","Thermoproteota")


vir_fam_coloring<-brewer.pal(11,"Paired")
names(vir_fam_coloring)<-c("Globuloviridae","Myoviridae","Phycodnaviridae","Siphoviridae","Baculoviridae","Podoviridae","Inoviridae","Sphaerolipoviridae","Lavidaviridae","Mimiviridae","Microviridae")

custom_family_order<-c("Myoviridae","Siphoviridae","Podoviridae","Inoviridae","Microviridae","Phycodnaviridae","Mimiviridae")

metab_coloring<-c(brewer.pal(9,"Set1"),brewer.pal(12,"Set3"))
names(metab_coloring)<-c("Amino acid metabolism","Biosynthesis of other secondary metabolites","Carbohydrate metabolism","Nucleotide metabolism","Cellular community - eukaryotes","Energy metabolism","Glycan biosynthesis and metabolism","Folding, sorting and degradation","Cellular community - prokaryotes","Metabolism of cofactors and vitamins","Membrane transport","Lipid metabolism","Metabolism of other amino acids","Metabolism of terpenoids and polyketides","Cell growth and death","Replication and repair","Signal transduction","Transcription","Translation","Transport and catabolism","Xenobiotics biodegradation and metabolism")

custom_metab_order<-sort(names(metab_coloring))

fraction_coloring<-c("#B35806","#1F78B4")
names(fraction_coloring)<-c("Particle-Attached","Free-Living")

zone_coloring<-c(brewer.pal(9,"YlOrRd")[c(4)],brewer.pal(9,"YlGnBu")[c(5,9)])
names(zone_coloring)<-c("Epipelagic","Mesopelagic","Bathypelagic")


###Metadata
metadata<-read.table(file="/mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Metadata/Selected_Malaspina_Profiles_Metadata.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)
rownames(metadata)<-metadata$MPCode
metadata$Sample<-metadata$MPCode
metadata$Zone<-"Epipelagic"
metadata$Zone[which(metadata$Depth_m >= 200)] <- "Mesopelagic"
metadata$Zone[which(metadata$Depth_m >= 1000)] <- "Bathypelagic"

###Prokaryote Bin/MAG data
mag_data<-read.table(file="/mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/Info_All_MAGs_Profiles_Malaspina_Redo.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)
colnames(mag_data)[1]<-c("MAG")

#Most up-to-date viral scaffold data
full_vir_scaff_data<-read.table(file="/mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Viruses/Profiles_Malaspina_Viruses_Full_Info.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)
#summary(full_vir_scaff_data)

###Viral Scaffold Abundance Data
vir_scaff_abund<-read.table(file="/mnt/lustre/scratch/elopez/5_bowtie_results/After_checkV_output/RPKM.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)

###Summarize abundance by PHIST predicted host phylum (or class for the Proteobacteria) using the aggregate function
tax_level<-"PHIST_Host_2_Phylum"

full_vir_scaff_data$PHIST_Host_2_Phylum[which(full_vir_scaff_data$PHIST_Host_2_Phylum == "Proteobacteria")]<-full_vir_scaff_data$PHIST_Host_3_Class[which(full_vir_scaff_data$PHIST_Host_2_Phylum == "Proteobacteria")]

m_vir_scaff_abund<-melt(vir_scaff_abund,id="Sequence") # FPKM

colnames(m_vir_scaff_abund)<-c("Virus_ID","Sample","RPKM")

m_vir_scaff_abund<-merge(m_vir_scaff_abund,subset(full_vir_scaff_data,select=c("Virus_ID",tax_level)),by="Virus_ID",all.x=TRUE)

colnames(m_vir_scaff_abund)
table(m_vir_scaff_abund$Sample)
table(m_vir_scaff_abund[[tax_level]])

#backup_m_vir_scaff_abund<-m_vir_scaff_abund
#m_vir_scaff_abund<-backup_m_vir_scaff_abund

m_vir_scaff_abund<-m_vir_scaff_abund[!is.na(m_vir_scaff_abund[[tax_level]]),]

summary(m_vir_scaff_abund)

levels(m_vir_scaff_abund$Sample)[levels(m_vir_scaff_abund$Sample)=="MP2233"] <- "MP2233bis"

table(m_vir_scaff_abund$Sample)

m_vir_scaff_abund<-merge(m_vir_scaff_abund,metadata,by="Sample",all.x=TRUE)

vir_tax_counts<-as.data.frame(aggregate(m_vir_scaff_abund$RPKM, by=list(Taxon=m_vir_scaff_abund[[tax_level]],Sample=m_vir_scaff_abund$Sample), FUN=sum))

colnames(vir_tax_counts)[3]<-"Virus_Abundance"

m_vir_tax_counts<-merge(vir_tax_counts,subset(metadata,select=c("Sample","Station","Depth_m","Zone")),by="Sample",all.x=TRUE)

summary(m_vir_tax_counts)

write.table(m_vir_tax_counts,file=paste("Long_Format_Malaspina_Profiles_PHIST_",tax_level,"AbundancexSample+Metadata.tsv",sep=""),sep="\t",append=FALSE,row.names=FALSE,col.names=TRUE,quote=FALSE)

#m_vir_tax_counts<-read.table(file="/mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Viruses/Long_Format_Malaspina_Profiles_PHIST_Host_2_Phylum_RPKM_ABundancexSample+Metadata.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)

fm_vir_tax_counts<-m_vir_tax_counts[which(m_vir_tax_counts$Virus_Abundance >= 1000),]

summary(fm_vir_tax_counts)

for (taxon in unique(fm_vir_tax_counts$Taxon)) {
	if (taxon %in% names(phylum_coloring)) {
	} else {
		print(paste(taxon," missing from color palette!",sep=""))
	}
	
	if (taxon %in% phylum_order) {
	} else {
		print(paste(taxon," missing from custom taxon order!",sep=""))
	}
}

sub_phylum_coloring<-phylum_coloring[names(phylum_coloring) %in% as.vector(unique(fm_vir_tax_counts$Taxon))]

fm_vir_tax_counts$UID<-paste("St",fm_vir_tax_counts$Station,"|",fm_vir_tax_counts$Depth_m,"m",sep="")

fm_vir_tax_counts$UID<-factor(fm_vir_tax_counts$UID,levels=unique(fm_vir_tax_counts$UID[rev(order(fm_vir_tax_counts$Depth_m))]))

fm_vir_tax_counts$Zone<-factor(fm_vir_tax_counts$Zone,levels=c("Epipelagic","Mesopelagic","Bathypelagic"))

sub_phylum_order<-phylum_order[phylum_order %in% as.vector(unique(fm_vir_tax_counts$Taxon))]

fm_vir_tax_counts$Taxon<-factor(fm_vir_tax_counts$Taxon,levels=sub_phylum_order)

fig2B<-ggplot(fm_vir_tax_counts,aes(x=Virus_Abundance,y=UID,fill=Taxon))+geom_bar(position="stack",stat="identity",colour="black",alpha=0.9)+theme_bw()+guides(fill=guide_legend(ncol=1))+facet_grid(Zone ~ ., scales = "free", space = "free")+ylab("Sample")+xlab("RPKM")+scale_fill_manual(name="Taxon",values=sub_phylum_coloring,breaks=sub_phylum_order)
fig_name<-paste("Malaspina_Profiles_Virus_Abundance_by_",tax_level,"_Barplots.pdf")
ggsave(fig_name,plot=fig2B,device=cairo_pdf,width=10,height=12,pointsize=8)


###Host Taxon Abundance
host_scaff_data<-read.table(file="/mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/Profiles_Malaspina_MAG_Scaffolds_Info_Redo.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)

host_scaff_data<-host_scaff_data %>% separate(Bin, c("Bin_1","Bin_2"),sep=",") # creates 2 columns
host_scaff_data$Bin<-host_scaff_data$Bin_1
host_scaff_data$Bin[grepl("Refined",host_scaff_data$Bin_2)]<-host_scaff_data$Bin_2[grepl("Refined",host_scaff_data$Bin_2)]
host_scaff_data$Bin<-gsub("(\\[)|(\\])|(\\')","",host_scaff_data$Bin,perl=TRUE) # change format of column Bin
colnames(host_scaff_data)[6]<-c("Sequence")
#summary(host_scaff_data)


host_scaff_abund<-read.table(file="/mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/Abundance_dRep_MAGs/RPKM_Abundance_Merged_Sequences.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)

host_scaff_abund<-merge(host_scaff_abund,subset(host_scaff_data,select=c("Sequence","Bin")),by="Sequence",all.x=TRUE)
#summary(host_scaff_abund)
m_host_scaff_abund<-melt(host_scaff_abund,id=c("Sequence","Bin"))
m_host_scaff_abund<-m_host_scaff_abund[which(m_host_scaff_abund$value > 0),]
colnames(m_host_scaff_abund)[2]<-"MAG"
tax_level<-"Phylum"

r_mag_data<-mag_data
r_mag_data$Phylum[which(r_mag_data$Phylum == "Proteobacteria")]<-r_mag_data$Class[which(r_mag_data$Phylum == "Proteobacteria")]

m_host_scaff_abund<-merge(m_host_scaff_abund,subset(r_mag_data,select=c("MAG",tax_level)),by="MAG",all.x=TRUE)

summary(m_host_scaff_abund)

m_host_scaff_abund<-m_host_scaff_abund[!is.na(m_host_scaff_abund[[tax_level]]),]

host_tax_counts<-as.data.frame(aggregate(m_host_scaff_abund$value, by=list(Taxon=m_host_scaff_abund[[tax_level]],Sample=m_host_scaff_abund$variable), FUN=sum))

colnames(host_tax_counts)<-c("Taxon","Sample","Host_Abundance")

summary(host_tax_counts)


###PHIST Data
all_phist_preds<-read.table(file="/mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Viruses/New_PHIST/Positive_PHIST_Predictions.csv",sep=",",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)

all_phist_preds$MAG<-gsub("(No_Vir_)|(.fa)","",all_phist_preds$host,perl=TRUE)
all_phist_preds$Virus_ID<-gsub("\\.fasta","",all_phist_preds$phage,perl=TRUE)

all_phist_preds<-merge(all_phist_preds,mag_data,by=c("MAG"),all.x=TRUE)

filtered_phist_preds<-all_phist_preds[which(all_phist_preds$pvalue <= 2.384e-14),]

filtered_phist_preds<-as.data.table(filtered_phist_preds)

filtered_phist_preds<-filtered_phist_preds[filtered_phist_preds[, .I[which.min(pvalue)], by=Virus_ID]$V1]

filtered_phist_preds<-as.data.frame(filtered_phist_preds)

col_order <- c("Virus_ID","MAG","phage","host","X.common.kmers","pvalue","adj.pvalue","Domain","Phylum","Class","Order","Family","Genus","Species","Completeness","Contamination")

filtered_phist_preds<-filtered_phist_preds[,col_order]

colnames(filtered_phist_preds)<-c("Virus_ID","MAG","Virus_Genome_File","PHIST_Host_Genome_File","PHIST_Common_Kmers","PHIST_pvalue","PHIST_adj_pvalue","PHIST_Host_1_Domain","PHIST_Host_2_Phylum","PHIST_Host_3_Class","PHIST_Host_4_Order","PHIST_Host_5_Family","PHIST_Host_6_Genus","PHIST_Host_7_Species","PHIST_Host_Genome_Completeness","PHIST_Host_Genome_Contamination")

#summary(filtered_phist_preds)



###Viral Scaffold full data
vir_scaff_data<-read.table(file="/mnt/lustre/scratch/elopez/Merged_Matrixes.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)
summary(vir_scaff_data)
dim(vir_scaff_data)
colnames(vir_scaff_data)<-c("Virus_ID","RaFAH_Predicted_Host","RaFAH_Predicted_Host_1_Domain","RaFAH_Predicted_Host_2_Phylum","RaFAH_Predicted_Host_3_Class","RaFAH_Predicted_Host_4_Order","RaFAH_Predicted_Host_5_Family","RaFAH_Predicted_Host_6_Genus","RaFAH_Predicted_Host_7_Species","RaFAH_Predicted_Host_Score","VPFClass_Baltimore_Classification","VPFClass_Family_Classification","VPFClass_Genus_Classification","GC","Untrimmed_Length","Description","AMG_List","VIBRANT_Lifestyle_Prediction","VIBRANT_Quality","Original_file","Sample_Source_Zone")

sub_vir_scaff_data<-subset(vir_scaff_data,select=c("Virus_ID","RaFAH_Predicted_Host_1_Domain","RaFAH_Predicted_Host_2_Phylum","RaFAH_Predicted_Host_3_Class","RaFAH_Predicted_Host_4_Order","RaFAH_Predicted_Host_5_Family","RaFAH_Predicted_Host_6_Genus","RaFAH_Predicted_Host_Score","VPFClass_Baltimore_Classification","VPFClass_Family_Classification","VPFClass_Genus_Classification","GC","Untrimmed_Length","Description","AMG_List","VIBRANT_Lifestyle_Prediction","VIBRANT_Quality","Original_file"))

sub_vir_scaff_data$Original_Sample<-gsub("^(VIBRANT_filtered_renamed_)|(\\/VIBRANT_phages_filtered_renamed_MP(\\d)+\\/filtered_renamed_(.)+)$","",sub_vir_scaff_data$Original_file,perl=TRUE)
table(sub_vir_scaff_data$Original_Sample)

sub_vir_scaff_data<-merge(sub_vir_scaff_data,filtered_phist_preds,by="Virus_ID",all.x=TRUE)


checkv_scaff_data<-read.table(file="/mnt/lustre/scratch/elopez/3_checkV_output/quality_summary.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)
summary(checkv_scaff_data)
colnames(checkv_scaff_data)<-c("Virus_ID","Trimmed_Length","CheckV_Is_Provirus","CheckV_Proviral_Length","CheckV_Gene_Count","viral_genes","host_genes","CheckV_quality","MIUVIG_Quality","CheckV_Completeness","CheckV_Completeness_method","CheckV_Contamination","kmer_freq","CheckV_Warnings")
sub_vir_scaff_data<-merge(sub_vir_scaff_data,subset(checkv_scaff_data,select=c("Virus_ID","Trimmed_Length","CheckV_Is_Provirus","CheckV_Proviral_Length","CheckV_Gene_Count","CheckV_quality","MIUVIG_Quality","CheckV_Completeness","CheckV_Completeness_method","CheckV_Contamination","CheckV_Warnings")),by="Virus_ID",all.x=TRUE)
summary(sub_vir_scaff_data)

colnames(metadata)[1]<-"Original_Sample"
full_vir_scaff_data<-merge(sub_vir_scaff_data,subset(metadata,select=c("Original_Sample","Depth_m","Zone")),by="Original_Sample",all.x=TRUE)
col_order<-c("Virus_ID","RaFAH_Predicted_Host_1_Domain","RaFAH_Predicted_Host_2_Phylum","RaFAH_Predicted_Host_3_Class","RaFAH_Predicted_Host_4_Order","RaFAH_Predicted_Host_5_Family","RaFAH_Predicted_Host_6_Genus","RaFAH_Predicted_Host_Score","VPFClass_Baltimore_Classification","VPFClass_Family_Classification","VPFClass_Genus_Classification","GC","Untrimmed_Length","Description","AMG_List","VIBRANT_Lifestyle_Prediction","VIBRANT_Quality","Original_file","MAG","Virus_Genome_File","PHIST_Host_Genome_File","PHIST_Common_Kmers","PHIST_pvalue","PHIST_adj_pvalue","PHIST_Host_1_Domain","PHIST_Host_2_Phylum","PHIST_Host_3_Class","PHIST_Host_4_Order","PHIST_Host_5_Family","PHIST_Host_6_Genus","PHIST_Host_7_Species","PHIST_Host_Genome_Completeness","PHIST_Host_Genome_Contamination","Trimmed_Length","CheckV_Is_Provirus","CheckV_Proviral_Length","CheckV_Gene_Count","CheckV_quality","MIUVIG_Quality","CheckV_Completeness","CheckV_Completeness_method","CheckV_Contamination","CheckV_Warnings","Original_Sample","Depth_m","Zone")
full_vir_scaff_data<-full_vir_scaff_data[,col_order]
summary(full_vir_scaff_data)
dim(full_vir_scaff_data)

write.table(full_vir_scaff_data,file="/mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Viruses/Profiles_Malaspina_Viruses_Full_Info.tsv",sep="\t",append=FALSE,row.names=FALSE,col.names=TRUE,quote=FALSE)

# mine
#write.table(full_vir_scaff_data,file="/mnt/lustre/scratch/elopez/malaspina_virus_merged_2_PHISTVPFClass_completetaxonomy.tsv",sep="\t",append=FALSE,row.names=FALSE,col.names=TRUE,quote=FALSE)


###Viral Scaffold Sequence Data
vir_scaff_data<-read.table(file="/mnt/lustre/scratch/elopez/malaspina_virus_merged_2_RaFAHVPFClass_completetaxonomy.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)
colnames(vir_scaff_data)[1]<-"Virus_ID"
summary(vir_scaff_data)

###Viral x Host Phylum Level Abundance Scatterplots
all_tax_counts<-merge(host_tax_counts,vir_tax_counts,by=c("Taxon","Sample"),all.x=TRUE)
all_tax_counts[is.na(all_tax_counts)]<-0

all_tax_counts<-merge(all_tax_counts,metadata,by=c("Sample"),all.x=TRUE)

min_sample_prevalence<-10

taxon_prev_count<-table(all_tax_counts$Taxon[which(all_tax_counts$Virus_Abundance > 0 & all_tax_counts$Host_Abundance > 0)])

passed_taxa<-names(taxon_prev_count[which(taxon_prev_count >= min_sample_prevalence)])

f_tax_counts<-all_tax_counts[which(all_tax_counts$Taxon %in% passed_taxa),]

figS1<-ggplot(f_tax_counts,aes(y=Virus_Abundance,x=Host_Abundance))+geom_abline(alpha=0.8,colour="red")+geom_smooth(method = "loess",se=TRUE,  formula = y ~ x)+geom_point(aes(colour=Zone),alpha=0.7,shape=18,size=3)+theme_bw()+scale_colour_manual(name="Zone",values=zone_coloring)+facet_wrap(Taxon ~ .,scales="free")

ggsave("Malaspina_Profiles_HostxVirus_Abundance_by_Phylum_Scatterplots.pdf",plot=figS1,device=cairo_pdf,width=15,height=10,pointsize=8)

f_tax_counts$Zone<-factor(f_tax_counts$Zone,levels=c("Epipelagic","Mesopelagic","Bathypelagic"))

figS2<-ggplot(f_tax_counts,aes(y=Virus_Abundance,x=Host_Abundance))+geom_abline(alpha=0.8,colour="red")+geom_smooth(method = "loess",se=TRUE,  formula = y ~ x)+geom_point(aes(colour=Zone),alpha=0.7,shape=18,size=3)+theme_bw()+scale_colour_manual(name="Zone",values=zone_coloring)+facet_wrap(Taxon ~ Zone,scales="free",ncol=3)

ggsave("Malaspina_Profiles_HostxVirus_Abundance_by_PhylumxZone_Scatterplots.pdf",plot=figS2,device=cairo_pdf,width=8,height=45,pointsize=8)

