library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(ggpubr, lib.loc="/mnt/lustre/bio/users/fcoutinho/Rlibs/")
library(tidyr)

###Aesthetics
#depth_col_pal<-brewer.pal(9,"Blues")[2:9]
#depth_col_pal<-c(brewer.pal(9,"Blues")[2:9],brewer.pal(9,"Purples")[7:9])
depth_col_pal<-brewer.pal(9,"YlGnBu")[c(1,4,9)]
depth_col_grad<-rev(colorRampPalette(depth_col_pal)(n=299))

egs_pal<-brewer.pal(9,"RdBu")
EGS_Gradient<-rev(colorRampPalette(egs_pal)(n=299))
#zone_coloring<-c(brewer.pal(9,"YlOrRd")[c(4)],brewer.pal(9,"YlGnBu")[c(5,9)])
zone_coloring<-c(brewer.pal(9,"YlGnBu")[c(1,4,9)])
names(zone_coloring)<-c("Epipelagic","Mesopelagic","Bathypelagic")

phylum_coloring<-c(brewer.pal(9,"Set1"),brewer.pal(4,"Pastel1"),brewer.pal(8,"Dark2")[8],brewer.pal(5,"Set2")[5],brewer.pal(8,"Accent"),brewer.pal(5,"Dark2"))
names(phylum_coloring)<-c("Actinobacteriota","Alphaproteobacteria","Cyanobacteria","Crenarchaeota","Bacteroidota","Nanoarchaeota","Verrucomicrobiota","Marinisomatota","Myxococcota","Bdellovibrionota","Gammaproteobacteria","Acidobacteriota","UBP7","Poribacteria","Eremiobacterota","Margulisbacteria","Gemmatimonadota","Desulfobacterota","Nitrospinota","Desulfobacterota_D","Patescibacteria","Hydrothermarchaeota","Chlamydiota","SAR324","Halobacteriota","Chloroflexota","Thermoplasmatota","Planctomycetota")

custom_taxon_order<-c("Acidobacteriota","Actinobacteriota","Bacteroidota","Chlamydiota" ,"Chloroflexota","Cyanobacteria","Desulfobacterota_D","Gemmatimonadota","Marinisomatota","Nitrospinota","Patescibacteria"  ,   "Planctomycetota","Alphaproteobacteria","Gammaproteobacteria" ,    "SAR324","Thermoplasmatota","Thermoproteota","Verrucomicrobiota")

scaff_info_file<-"/mnt/lustre/scratch/fcoutinho/StG_23/Viruses/Sequence_Info/Full_Info_IMGVR_Scaffolds_Comp_75-100_Max_Conta_0_vOTU_Representatives.tsv"
scaff_info_df<-read.table(file=scaff_info_file,sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)
#head(scaff_info_df)
scaff_info_vars<-colnames(scaff_info_df)

#####Virus Abundance Analysis
group_var<-"Host_Phylum_PHIST"

raw_abd_file<-"/mnt/lustre/scratch/fcoutinho/StG_23/Viruses/Abundance/Raw_Abundance_IMGVR_Scaffolds_Comp_75-100_Max_Conta_0_vOTU_Representatives.tsv"
abd_df<-read.table(file=raw_abd_file,sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)
abd_df$X<-as.factor(gsub("\\|.*$","",abd_df$X,perl=TRUE))
head(abd_df$X)
colnames(abd_df)[1]<-"UVIG"
abd_df<-merge(abd_df,scaff_info_df[,c("UVIG",group_var)],by="UVIG",all.x=TRUE)
summary(abd_df)

m_abd_df<-melt(abd_df,id.vars=c("UVIG",group_var),variable.name="Sample",value.name="Abundance")
summary(m_abd_df)

f_m_abd_df<-m_abd_df[which(m_abd_df$Abundance > 0),]
summary(f_m_abd_df)

abd_sums_df<-aggregate(f_m_abd_df$Abundance,by=list(Sample=f_m_abd_df$Sample,Group=f_m_abd_df[[group_var]]),FUN=sum)

colnames(abd_sums_df)[3]<-"Abundance"

f_abd_sums_df<-abd_sums_df[which(abd_sums_df$Abundance > 1000),]

fig1<-ggplot(f_abd_sums_df,aes(x=Sample,y=Abundance,fill=Group))+geom_bar(position = "stack",stat="identity",colour="black",alpha=0.9)+theme_bw()+theme(axis.text.x = element_text(angle = 90,hjust = 1,size=9))

out_name<-paste("/mnt/lustre/scratch/fcoutinho/StG_23/Viruses/Abundance/StG_23_Viruses_",group_var,"_Abundance_Barplots.pdf",sep="")
ggsave(out_name,plot=fig1,device=cairo_pdf,width=20,height=7,pointsize=8)

#PHIST Host Pred Processing
library(data.table)
 
all_phist_preds<-read.table(file="/mnt/lustre/scratch/fcoutinho/StG_23/Viruses/PHIST_OceanDNA/PHIST_Output/Positive_PHIST_Predictions.csv",sep=",",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)

all_phist_preds$MAG<-gsub("\\.fa$","",all_phist_preds$host,perl=TRUE)
all_phist_preds$Virus_acc<-gsub(".fasta","",all_phist_preds$phage,perl=TRUE)

mags_tax<-read.table(file="/mnt/lustre/bio/users/fcoutinho/OceanDNA_MAGs/Sub_OceanDNA_MAGs_Info.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)

all_phist_preds<-merge(all_phist_preds,mags_tax,by.x="MAG",by.y="genome",all.x=TRUE)

summary(all_phist_preds)

filtered_phist_preds<-all_phist_preds[which(all_phist_preds$pvalue <= 2.384e-14),]

#Giant virus mags are often split into multiple contigs, so they should have gone in the same fata file. Since they are euk viruses, PHISt preds are not relevant, and these predictions are excluded for simplicity and avoiding having tor epate it all.
filtered_phist_preds<-filtered_phist_preds[!grepl("GVMAG",filtered_phist_preds$phage),]

summary(filtered_phist_preds)

filtered_phist_preds<-as.data.table(filtered_phist_preds)

filtered_phist_preds<-filtered_phist_preds[filtered_phist_preds[, .I[which.min(pvalue)], by=Virus_acc]$V1]

filtered_phist_preds<-as.data.frame(filtered_phist_preds)

col_order <- c("Virus_acc","MAG","phage","host","X.common.kmers","pvalue","adj.pvalue","Domain","Phylum","Class","Order","Family","Genus","Species")

filtered_phist_preds<-filtered_phist_preds[,col_order]

colnames(filtered_phist_preds)<-c("Virus_acc","MAG","phage_fasta_file","host_fasta_file","X.common.kmers","pvalue","adj.pvalue","Host_Domain","Host_Phylum","Host_Class","Host_Order","Host_Family","Host_Genus","Host_Species")

filtered_phist_preds$Virus_acc<-as.factor(filtered_phist_preds$Virus_acc)
filtered_phist_preds$MAG<-as.factor(filtered_phist_preds$MAG)

summary(filtered_phist_preds)

#There are differences in the IDS of virals equence sin the fasta files and the info tables that are accounted for in the line below.
filtered_phist_preds$Virus_acc<-as.factor(gsub("\\|.*$","",filtered_phist_preds$Virus_acc,perl=TRUE))

summary(filtered_phist_preds)

write.table(filtered_phist_preds,file="/mnt/lustre/bio/users/fcoutinho/Databases/IMG_VR_4.0/IMGVRxOcean_DNA_MAGs_Filtered_PHIST_Preds.tsv",sep="\t",append=FALSE,row.names=FALSE,col.names=TRUE,quote=FALSE)

scaff_info_df<-merge(scaff_info_df,filtered_phist_preds,by.x="UVIG",by.y="Virus_acc",all.x=TRUE,suffixes=c("_IMGVR","_PHIST"))

scaff_info_df$Host_Phylum_PHIST[which(scaff_info_df$Host_Phylum_PHIST == "Proteobacteria")]<-scaff_info_df$Host_Class_PHIST[which(scaff_info_df$Host_Phylum_PHIST == "Proteobacteria")]

summary(scaff_info_df)

write.table(scaff_info_df, file="/mnt/lustre/scratch/fcoutinho/StG_23/Viruses/Sequence_Info/Full_Info_IMGVR_Scaffolds_Comp_75-100_Max_Conta_0_vOTU_Representatives.tsv",sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)

