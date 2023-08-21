args = commandArgs(trailingOnly=TRUE)
scaff_abund_file<-args[1]
mag_info_file<-args[2]
min_prev<-as.numeric(args[3])
output_file<-args[4]

library(dplyr)
library(tidyr)
library(tibble)
library(reshape2)

#scaff_abund_file<-"/mnt/lustre/scratch/fcoutinho/StG_23/Prokaryotes/Abundance/Short_RPKM.tsv"
#mag_info_file<-"/mnt/lustre/bio/users/fcoutinho/OceanDNA_MAGs/Sub_OceanDNA_MAGs_Info.tsv"
#min_prev<-50 
#output_file<-"/mnt/lustre/scratch/fcoutinho/StG_23/Prokaryotes/Abundance/RPKM_Abundance_OceanDNA_All_Species_Rep_MAGs_by_Species_Cluster.tsv"

raw_response_df<-read.table(file=scaff_abund_file,sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE,row.names=1,check.names=FALSE)

#dim(raw_response_df)

raw_response_df$MAG_ID<-as.factor(gsub("_[0-9]+$","",rownames(raw_response_df),perl=TRUE))

m_raw_response_df<-melt(raw_response_df,id.vars=c("MAG_ID"),variable.name="Sample",value.name="Abundance")

f_m_raw_response_df<-m_raw_response_df[which(m_raw_response_df$Abundance > 0),]

m_mag_abd_sums_df<-aggregate(f_m_raw_response_df$Abundance,by=list(Sample=f_m_raw_response_df$Sample,MAG_ID=f_m_raw_response_df$MAG_ID),FUN=sum)

colnames(m_mag_abd_sums_df)[3]<-"Abundance"

f_m_mag_abd_sums_df<-m_mag_abd_sums_df[which(m_mag_abd_sums_df$Abundance > 0),]

mag_info_df<-read.table(file=mag_info_file,sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE,check.names=FALSE)

f_m_mag_abd_sums_df<-merge(f_m_mag_abd_sums_df,mag_info_df[,c("genome","species_cluster")],by.x="MAG_ID",by.y="genome",all.x=TRUE)

m_srep_abd_sums_df<-aggregate(f_m_mag_abd_sums_df$Abundance,by=list(Sample=f_m_mag_abd_sums_df$Sample,Species_Cluster_ID=f_m_mag_abd_sums_df$species_cluster),FUN=sum)

colnames(m_srep_abd_sums_df)[3]<-"Abundance"

prevalence<-table(m_srep_abd_sums_df$Species_Cluster_ID)

passed_prev_ids<-names(prevalence[which(prevalence >= min_prev)])

f_m_srep_abd_sums_df<-m_srep_abd_sums_df[which(m_srep_abd_sums_df$Species_Cluster_ID %in% passed_prev_ids),]

srep_adb_df<-dcast(f_m_srep_abd_sums_df,Species_Cluster_ID ~ Sample,value.var="Abundance",fun.aggregate=sum,fill=0)

#summary(srep_adb_df)

write.table(srep_adb_df,file=output_file,sep="\t",append=FALSE,row.names=FALSE,col.names=TRUE,quote=FALSE)