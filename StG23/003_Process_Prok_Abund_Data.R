args = commandArgs(trailingOnly=TRUE)
scaff_abund_file<-args[1]
mag_info_file<-args[2]
group_var<-args[3]
min_prev<-as.numeric(args[4])
output_file<-args[5]

library(dplyr)
library(tidyr)
library(tibble)
library(reshape2)

#scaff_abund_file<-"/mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/StG23/Prokaryotes/Abundance/RPKM_Abundance_OceanDNA_All_Species_Rep_MAGs_Scaffolds.tsv"
#mag_info_file<-"/mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/Databases/OceanDNA/Sub_OceanDNA_MAGs_Info.tsv"
#min_prev<-50 
#group_var<-"MAG_ID"
#output_file<-"/mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/StG23/Prokaryotes/Abundance/RPKM_Abundance_OceanDNA_All_Species_Rep_MAGs_by_Species_Cluster.tsv"

#scaff_abund_file<-"/mnt/lustre/scratch/fcoutinho/StG_23/Prokaryotes/Abundance/Short_RPKM.tsv"
#mag_info_file<-"/mnt/lustre/bio/users/fcoutinho/OceanDNA_MAGs/Sub_OceanDNA_MAGs_Info.tsv"
#min_prev<-50 
#output_file<-"/mnt/lustre/scratch/fcoutinho/StG_23/Prokaryotes/Abundance/RPKM_Abundance_OceanDNA_All_Species_Rep_MAGs_by_Species_Cluster.tsv"

#Read in scaffold abundance data frame
raw_response_df<-read.table(file=scaff_abund_file,sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE,row.names=1,check.names=FALSE)
print("Scaffold Abundance DF dimensions:")
print(dim(raw_response_df))

#For Ocean DNA Sscaffold ids follow the MAG_SCaffold_Number format. Thus we get MAG ID from scaffold IDs
raw_response_df$MAG_ID<-as.factor(gsub("_[0-9]+$","",rownames(raw_response_df),perl=TRUE))
print("Unique MAG IDs identified in Scaffold Abundance DF:")
print(length(unique(raw_response_df$MAG_ID)))

#Melt df and remove all scaffolds with 0 abundance to speed things up
m_raw_response_df<-melt(raw_response_df,id.vars=c("MAG_ID"),variable.name="Sample",value.name="Abundance")

f_m_raw_response_df<-m_raw_response_df[which(m_raw_response_df$Abundance > 0),]

#Calculate the sum of abundances for each MAG in each sample
m_mag_abd_sums_df<-aggregate(f_m_raw_response_df$Abundance,by=list(Sample=f_m_raw_response_df$Sample,MAG_ID=f_m_raw_response_df$MAG_ID),FUN=sum)

colnames(m_mag_abd_sums_df)[3]<-"Abundance"

# remove all MAGs with 0 abundance to speed things up
f_m_mag_abd_sums_df<-m_mag_abd_sums_df[which(m_mag_abd_sums_df$Abundance > 0),]

if (group_var != "MAG_ID") {
    mag_info_df<-read.table(file=mag_info_file,sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE,check.names=FALSE)

    f_m_mag_abd_sums_df<-merge(f_m_mag_abd_sums_df,mag_info_df[,c("genome","species_cluster")],by.x="MAG_ID",by.y="genome",all.x=TRUE)
}

#Calculate the sum of abundances for grouping var for  each sample
m_srep_abd_sums_df<-aggregate(f_m_mag_abd_sums_df$Abundance,by=list(Sample=f_m_mag_abd_sums_df$Sample,Taxon=f_m_mag_abd_sums_df[[group_var]]),FUN=sum)

colnames(m_srep_abd_sums_df)[3]<-"Abundance"

#Calculate the rpevalenc eof groupds across samples. Prev is absolute (i.e. number of samples with abundances above 0), not relative
prevalence<-table(m_srep_abd_sums_df$Taxon)

#remove the groups which do not meet the minimum prevalence cutoff
passed_prev_ids<-names(prevalence[which(prevalence >= min_prev)])

f_m_srep_abd_sums_df<-m_srep_abd_sums_df[which(m_srep_abd_sums_df$Taxon %in% passed_prev_ids),]

#Cast the data into a DF and print it to the output file
srep_adb_df<-dcast(f_m_srep_abd_sums_df,Taxon ~ Sample,value.var="Abundance",fun.aggregate=sum,fill=0)

#summary(srep_adb_df)

write.table(srep_adb_df,file=output_file,sep="\t",append=FALSE,row.names=FALSE,col.names=TRUE,quote=FALSE)