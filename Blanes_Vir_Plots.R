library(ggpubr)
library(reshape2)
library(ggplot2)
library(RColorBrewer)

debug<-TRUE
###Aesthetics
#depth_col_pal<-brewer.pal(9,"Blues")[2:9]
#depth_col_pal<-c(brewer.pal(9,"Blues")[2:9],brewer.pal(9,"Purples")[7:9])
depth_col_pal<-brewer.pal(9,"YlGnBu")
depth_col_grad<-rev(colorRampPalette(depth_col_pal)(n=299))

div_col_pal<-brewer.pal(9,"RdBu")[c(1:9)]
div_col_grad<-rev(colorRampPalette(div_col_pal)(n=99))

taxon_order<-c("Acidobacteriota","Actinobacteriota","Bacteroidota","Bdellovibrionota","Chlamydiota","Chloroflexota","Cyanobacteria","Desulfobacterota","DesulfobacterotaB","Eremiobacterota","Gemmatimonadota","Margulisbacteria","Marinisomatota","Myxococcota","Nitrospinota","Patescibacteria","Planctomycetota","Poribacteria","Alphaproteobacteria","Gammaproteobacteria","SAR324","Latescibacterota","Verrucomicrobiota","Halobacteriota","Hydrothermarchaeota","Nanoarchaeota","Thermoplasmatota","Thermoproteota")

phylum_coloring<-c(brewer.pal(9,"Set1"),brewer.pal(4,"Pastel1"),brewer.pal(8,"Dark2")[8],brewer.pal(5,"Set2")[5],brewer.pal(8,"Accent"),brewer.pal(6,"Dark2"))

names(phylum_coloring)<-c("Actinobacteriota","Alphaproteobacteria","Cyanobacteria","Crenarchaeota","Bacteroidota","Nanoarchaeota","Verrucomicrobiota","Marinisomatota","Myxococcota","Patescibacteria","Gammaproteobacteria","Acidobacteriota","Latescibacterota","Poribacteria","Eremiobacterota","Margulisbacteria","Gemmatimonadota","Desulfobacterota","Nitrospinota","DesulfobacterotaB","Bdellovibrionota","Hydrothermarchaeota","Chlamydiota","SAR324","Halobacteriota","Chloroflexota","Thermoplasmatota","Planctomycetota","Thermoproteota")

phylum_custom_selection<-c("Actinobacteriota","Alphaproteobacteria","Cyanobacteria","Bacteroidota","Verrucomicrobiota","Marinisomatota","Gammaproteobacteria","Nitrospinota","SAR324","Chloroflexota","Thermoplasmatota","Planctomycetota","Thermoproteota")

###Viral scaffold info
full_vir_scaff_data<-read.table(file="/mnt/smart/users/fcoutinho/Blanes_Vir_Nestor/Seq_Info/Info_VP_Rep_Propagated_Hosts_dsDNAphage_Blanes_virus.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)
summary(full_vir_scaff_data)



###Raw phist preds (FT3)
all_phist_preds<-read.table(file="/mnt/smart/users/fcoutinho/Blanes_Vir_Nestor/PHIST_Predictions/PHIST_Output_Clean_MAGs/Positive_PHIST_Predictions.csv",sep=",",header=TRUE,quote="",comment="",stringsAsFactors=TRUE,check.names=FALSE)
#all_phist_preds<-read.table(file="/mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/Blanes_Vir_Nestor/PHIST_Predictions/PHIST_Output/positive_predictions.csv",sep=",",header=TRUE,quote="",comment="",stringsAsFactors=TRUE,check.names=FALSE)

summary(all_phist_preds)

Host_ID_data<-read.table(file="/mnt/smart/users/fcoutinho/Blanes_Vir_Nestor/Seq_Info/Blanes_MAG_Info.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE,check.names=FALSE)
#Host_ID_data<-read.table(file="/mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/Blanes_Vir_Nestor/Blanes_MAG_Info.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE,check.names=FALSE)
colnames(Host_ID_data)[1]<-"Host_ID"

summary(Host_ID_data)

filter_phist_preds <- function(phist_df,host_df,max_p_value=2.384e-14) {
    library(data.table)

    colnames(all_phist_preds)[1]<-"Virus_ID"
    colnames(all_phist_preds)[2]<-"Host_ID"

    all_phist_preds$Host_ID<-as.factor(gsub("^(No_Vir_)|(\\.fasta)$|(\\.fa)$","",all_phist_preds$Host_ID,perl=TRUE))
    all_phist_preds$Virus_ID<-as.factor(gsub("(\\.fasta)$|(\\.fa)$","",all_phist_preds$Virus_ID,perl=TRUE))

    all_phist_preds<-merge(all_phist_preds,Host_ID_data,by=c("Host_ID"),all.x=TRUE)

    #summary(all_phist_preds)

    write.table(all_phist_preds,file="All_PHIST_Predictions_with_Host_Data.tsv",sep="\t",quote=FALSE,row.names=FALSE)

    filtered_phist_preds<-all_phist_preds[which(all_phist_preds$pvalue <= max_p_value),]

    filtered_phist_preds<-as.data.table(filtered_phist_preds)

    filtered_phist_preds<-filtered_phist_preds[filtered_phist_preds[, .I[which.min(pvalue)], by=Virus_ID]$V1]

    filtered_phist_preds<-as.data.frame(filtered_phist_preds)

    write.table(filtered_phist_preds,file="Best_Filtered_PHIST_Predictions_with_Host_Data.tsv",sep="\t",quote=FALSE,row.names=FALSE)

    #summary(filtered_phist_preds)

    return(filtered_phist_preds)
}

filtered_phist_preds<-filter_phist_preds(phist_df=all_phist_preds,host_df=Host_ID_data,max_p_value=2.384e-14)

#Remove viruses from the MAgs (redundant with viruses from the assemblies) 
filtered_phist_preds<-filtered_phist_preds[!grepl("bin",filtered_phist_preds$Virus_ID),]

summary(filtered_phist_preds)

write.table(filtered_phist_preds,file="/mnt/smart/users/fcoutinho/Blanes_Vir_Nestor/PHIST_Predictions/PHIST_Output_Clean_MAGs/Best_Filtered_PHIST_Predictions_with_Host_Data.tsv",sep="\t",quote=FALSE,row.names=FALSE)
###Parse Viral population and phist data (MARBITS)
vpop_data<-read.table(file="/mnt/smart/users/fcoutinho/Blanes_Vir_Nestor/Vpop/VPop_dsDNAphage_Blanes_virus.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)
summary(vpop_data)

phist_data<-read.table(file="/mnt/smart/users/fcoutinho/Blanes_Vir_Nestor/PHIST_Predictions/PHIST_Output_Clean_MAGs/Best_Filtered_PHIST_Predictions_with_Host_Data.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)
#phist_data<-read.table(file="/mnt/smart/users/fcoutinho/Blanes_Vir_Nestor/PHIST_Predictions/PHIST_Output/Best_Filtered_PHIST_Predictions_with_Host_Data.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)
summary(phist_data)
phist_data$GC<-NULL
phist_data$Length<-NULL

full_data<-merge(vpop_data,phist_data,by.x="Sequence", by.y="Virus_ID",all.x=TRUE)

write.table(full_data,file="/mnt/smart/users/fcoutinho/Blanes_Vir_Nestor/Info_Indiv_Hosts_dsDNAphage_Blanes_virus.tsv",sep="\t",append=FALSE,row.names=FALSE,col.names=TRUE,quote=FALSE)

vp_h_data<-full_data[which(full_data$Population_Representative == "True" & !is.na(full_data$Host_ID)),]

full_data2<-merge(vpop_data,vp_h_data[,c("Population","Host_ID","Domain","Phylum","Class","Order","Family", "Genus","Species","Complete.GTDB.taxonomy")],by="Population",all.x=TRUE)

# full_data2<-full_vir_scaff_data

full_data2$Phylum<-as.character(full_data2$Phylum)
full_data2$Class<-as.character(full_data2$Class)

table(full_data2$Phylum)

full_data2$Phylum[which(full_data2$Phylum == "Proteobacteria")]<-full_data2$Class[which(full_data2$Phylum == "Proteobacteria")]

full_data2$Phylum<-as.factor(full_data2$Phylum)

summary(full_data2$Phylum)

write.table(full_data2,file="/mnt/smart/users/fcoutinho/Blanes_Vir_Nestor/Info_VP_Rep_Propagated_Hosts_dsDNAphage_Blanes_virus.tsv",sep="\t",append=FALSE,row.names=FALSE,col.names=TRUE,quote=FALSE)



full_vir_scaff_data<-merge(vir_scaff_data,filtered_phist_preds,by="Virus_ID",all.x=TRUE)



###Parse metabolic output
for (i in 10:10) {
    metabolic_data<-read.table(file=paste0("/mnt/smart/users/fcoutinho/Blanes_Vir_Nestor/Metabolic_Output/Metabolic_Outputs_Batch_",i,"/METABOLIC_result_each_spreadsheet/METABOLIC_result_worksheet4.tsv"),sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)
    head(metabolic_data)
}

func_data<-read.table(file="/mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/Metabolic_Preds/Profiles_Malaspina_MAGs_METABOLIC_Results/METABOLIC_result_each_spreadsheet/METABOLIC_result_worksheet4.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)
dim(func_data)

long_func_data<-melt(func_data,id=c("Module.step","Module","KO.id","Module.Category"))

colnames(long_func_data)<-c("Module_Step","Module_Name","KO_IDs","Module_Category","MAG","Module_Step_Presence")

long_func_data$MAG<-gsub(".Module.step.presence$","",long_func_data$MAG,perl=TRUE)


###Abundance patterns
#Read in the viral abundance data
vir_scaff_abund<-read.table(file="/mnt/smart/users/fcoutinho/Blanes_Vir_Nestor/Abundance/RPKM_Abundance_All_Genomic.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)
dim(vir_scaff_abund)

rownames<-vir_scaff_abund$Sequence

vir_scaff_abund<-merge(vir_scaff_abund,subset(full_vir_scaff_data,select=c(Sequence,Population)),by="Sequence",all.x=TRUE)

m_vir_scaff_abund<-melt(vir_scaff_abund,id=c("Sequence","Population"),value.name="Abundance",variable.name="Sample") # FPKM

head(m_vir_scaff_abund)

fm_vir_scaff_abund<-m_vir_scaff_abund[which(m_vir_scaff_abund$Abundance > 0),]

vir_tax_counts<-as.data.frame(aggregate(fm_vir_scaff_abund$Abundance, by=list(Taxon=fm_vir_scaff_abund$Population,Sample=fm_vir_scaff_abund$Sample), FUN=sum))


colnames(vir_tax_counts)[3]<-"Abundance"

summary(vir_tax_counts)

m_vir_tax_counts<-vir_tax_counts[which(vir_tax_counts$Abundance >= 0),]

summary(m_vir_tax_counts)

compl_col_grad<-rev(colorRampPalette(brewer.pal(11,"Spectral"))(n=100))

fig2B<-ggplot(m_vir_tax_counts,aes(fill=log10(Abundance),x=Sample,y=Taxon))+
geom_tile()+
theme_bw()+
scale_fill_gradientn(colours = compl_col_grad)+
scale_y_discrete(labels = NULL, breaks = NULL)+
theme(axis.text.x = element_text(angle = 90,hjust = 1,size=8))+
guides(fill=guide_legend(ncol=1))+xlab("Sample")+ylab("VP")

ggsave("Virus_Abundance_by_VP_Tileplots.pdf",plot=fig2B,device=cairo_pdf,width=25,height=10,pointsize=8)

#melt df, remove abundances of zero and unassinged sequences
m_vir_scaff_abund<-melt(vir_scaff_abund,id=c("Sequence","Phylum"),value.name="Abundance",variable.name="Sample") # FPKM

head(m_vir_scaff_abund)

fm_vir_scaff_abund<-m_vir_scaff_abund[which((m_vir_scaff_abund$Abundance > 0) & !is.na(m_vir_scaff_abund$Phylum)),]

head(fm_vir_scaff_abund)

tax_level<-"Phylum"

vir_tax_counts<-as.data.frame(aggregate(fm_vir_scaff_abund$Abundance, by=list(Taxon=fm_vir_scaff_abund[[tax_level]],Sample=fm_vir_scaff_abund$Sample), FUN=sum))

colnames(vir_tax_counts)[3]<-"Abundance"


summary(vir_tax_counts)

m_vir_tax_counts<-vir_tax_counts[which(vir_tax_counts$Abundance >= 100),]

for (taxon in unique(m_vir_tax_counts$Taxon)) {
	if (taxon %in% names(phylum_coloring)) {
	} else {
		print(paste(taxon," missing from color palette!",sep=""))
	}
	
	if (taxon %in% taxon_order) {
	} else {
		print(paste(taxon," missing from custom taxon order!",sep=""))
	}
}


 [1] Actinobacteriota    Alphaproteobacteria Bacteroidota
 [4] Chloroflexota       Cyanobacteria       Gammaproteobacteria
 [7] Gemmatimonadota     Margulisbacteria    Marinisomatota
[10] Myxococcota         Planctomycetota     Thermoplasmatota
[13] Thermoproteota      Verrucomicrobiota   SAR324
[16] DesulfobacterotaB   Bdellovibrionota    Latescibacterota


#sub_phylum_coloring<-phylum_coloring[unique(m_vir_tax_counts$Taxon)]

#sub_custom_taxon_order<-custom_taxon_order[custom_taxon_order %in% unique(m_vir_tax_counts$Taxon)]

m_vir_tax_counts$Date<-as.Date(gsub("BL","",m_vir_tax_counts$Sample),format="%d%m%y")

m_vir_tax_counts$Sample<-factor(m_vir_tax_counts$Sample,levels=sort(unique(m_vir_tax_counts$Sample)))

summary(m_vir_tax_counts)

fig2B<-ggplot(m_vir_tax_counts,aes(y=Abundance,x=Sample,fill=Taxon))+
geom_bar(position="stack",stat="identity",colour="black",alpha=0.9)+
theme_bw()+
#scale_fill_manual(name="Taxon",values=phylum_coloring)+ #,breaks=sub_custom_taxon_order
theme(axis.text.x = element_text(angle = 90,hjust = 1,size=8))+
guides(fill=guide_legend(ncol=1))+xlab("Sample")+ylab("RPKM")
#scale_fill_manual(name="Taxon",values=sub_phylum_coloring,breaks=sub_custom_taxon_order)+

ggsave("Virus_Abundance_by_PHIST_Phylum_Barplots.pdf",plot=fig2B,device=cairo_pdf,width=25,height=10,pointsize=8)
