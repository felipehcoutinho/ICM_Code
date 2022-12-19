library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(ggpubr, lib.loc="/mnt/lustre/bio/users/fcoutinho/Rlibs/")

###Aesthetics
#depth_col_pal<-brewer.pal(9,"Blues")[2:9]
#depth_col_pal<-c(brewer.pal(9,"Blues")[2:9],brewer.pal(9,"Purples")[7:9])
depth_col_pal<-brewer.pal(9,"YlGnBu")
depth_col_grad<-rev(colorRampPalette(depth_col_pal)(n=299))

#zone_coloring<-c(brewer.pal(9,"YlOrRd")[c(4)],brewer.pal(9,"YlGnBu")[c(5,9)])
zone_coloring<-c(brewer.pal(9,"YlGnBu")[c(1,4,9)])
names(zone_coloring)<-c("Epipelagic","Mesopelagic","Bathypelagic")

taxon_order<-c("Acidobacteriota","Actinobacteriota","Bacteroidota","Bdellovibrionota","Chlamydiota","Chloroflexota","Cyanobacteria","Desulfobacterota","Desulfobacterota_D","Eremiobacterota","Gemmatimonadota","Margulisbacteria","Marinisomatota","Myxococcota","Nitrospinota","Patescibacteria","Planctomycetota","Poribacteria","Alphaproteobacteria","Gammaproteobacteria","SAR324","UBP7","Verrucomicrobiota","Halobacteriota","Hydrothermarchaeota","Nanoarchaeota","Thermoplasmatota","Thermoproteota")

phylum_coloring<-c(brewer.pal(9,"Set1"),brewer.pal(4,"Pastel1"),brewer.pal(8,"Dark2")[8],brewer.pal(5,"Set2")[5],brewer.pal(8,"Accent"),brewer.pal(6,"Dark2"))
names(phylum_coloring)<-c("Actinobacteriota","Alphaproteobacteria","Cyanobacteria","Crenarchaeota","Bacteroidota","Nanoarchaeota","Verrucomicrobiota","Marinisomatota","Myxococcota","Patescibacteria","Gammaproteobacteria","Acidobacteriota","UBP7","Poribacteria","Eremiobacterota","Margulisbacteria","Gemmatimonadota","Desulfobacterota","Nitrospinota","Desulfobacterota_D","Bdellovibrionota","Hydrothermarchaeota","Chlamydiota","SAR324","Halobacteriota","Chloroflexota","Thermoplasmatota","Planctomycetota","Thermoproteota")


###Metadata
metadata<-read.table(file="/mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Metadata/Selected_Malaspina_Profiles_Metadata.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)
rownames(metadata)<-metadata$MPCode
metadata$Sample<-metadata$MPCode
metadata$Zone<-"Epipelagic"
metadata$Zone[which(metadata$Depth_m >= 200)] <- "Mesopelagic"
metadata$Zone[which(metadata$Depth_m >= 1000)] <- "Bathypelagic"
metadata$UID<-paste("St",metadata$Station,"|",metadata$Depth_m,"m",sep="")

###Updated MAG data
mag_data<-read.table(file="/mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/Info_All_MAGs_Profiles_Malaspina_Redo.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)

mag_data$Zone<-factor(mag_data$Zone,levels=c("Epipelagic","Mesopelagic","Bathypelagic"))

#MAG_Abundance
scaff_data<-read.table(file="/mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/Profiles_Malaspina_MAG_Scaffolds_Info_Redo.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)

scaff_data<-scaff_data %>% separate(Bin, c("Bin_1","Bin_2"),sep=",")
scaff_data$Bin<-scaff_data$Bin_1
scaff_data$Bin[grepl("Refined",scaff_data$Bin_2)]<-scaff_data$Bin_2[grepl("Refined",scaff_data$Bin_2)]

scaff_data$Bin<-gsub("(\\[)|(\\])|(\\')","",scaff_data$Bin,perl=TRUE)

colnames(scaff_data)[6]<-c("Sequence")

summary(scaff_data)

scaff_abund<-read.table(file="/mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/Abundance_dRep_MAGs/RPKM_Abundance_Merged_Sequences.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)


scaff_abund<-merge(scaff_abund,subset(scaff_data,select=c("Sequence","Bin")),by="Sequence",all.x=TRUE)
#scaff_abund<-merge(scaff_abund,subset(mag_data,select=c("Bin",tax_level)),by="Bin",all.x=TRUE)

summary(scaff_abund)

m_scaff_abund<-melt(scaff_abund,id=c("Sequence","Bin"))

all_levels_tax_counts<-as.data.frame(aggregate(m_scaff_abund$value, by=list(Taxon=m_scaff_abund$Bin,Sample=m_scaff_abund$variable), FUN=sum))

colnames(all_levels_tax_counts)<-c("MAG","Sample","Abundance")

all_levels_tax_counts$MAG<-gsub(".*Original_and_Refined_Bins.","",all_levels_tax_counts$MAG,perl=TRUE)

head(all_levels_tax_counts)

write.table(all_levels_tax_counts,file="/mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/Abundance_dRep_MAGs/Profiles_Malaspina_MAG_RPKM_Abundance.tsv",sep="\t",append=FALSE,row.names=FALSE,col.names=TRUE,quote=FALSE)


###Key KO
keyko_data<-read.table(file="/mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/Metabolic_Preds/Key_KOs_Long_Format_No_Absent_Profiles_Malaspina_MAGs_METABOLIC_result_worksheet4.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)

keyko_data$Full_Taxonomy<-paste(keyko_data$Domain,keyko_data$Phylum,keyko_data$Class,keyko_data$Order,keyko_data$Family,keyko_data$Genus,keyko_data$MAG,sep="|")

#keyko_data$Zone<-factor(keyko_data$Zone,levels=c("Epipelagic","Mesopelagic","Bathypelagic"))

keyko_data$Full_Taxonomy<-factor(keyko_data$Full_Taxonomy,levels=sort(unique(keyko_data$Full_Taxonomy),decreasing=TRUE))

keyko_data$Module_Category<-str_wrap(keyko_data$Module_Category,width=20)
keyko_data$Module_Name<-str_wrap(keyko_data$Module_Name,width=20)

raster_fig<-ggplot(keyko_data,aes(x=Module_Step,y=Full_Taxonomy,fill=Zone))+geom_tile()+scale_fill_manual(name="Zone",values=zone_coloring)+theme_bw()+theme(axis.text.x = element_text(angle = 90,hjust = 1,size=6),axis.text.y = element_text(size=1),strip.text.x.top = element_text(angle = 90),strip.text.y.right = element_text(angle = 0,size=1),legend.position="top")+facet_grid(Phylum ~ Module_Category + Module_Name, space = "free", scales = "free")#+scale_fill_gradientn(colours =depth_col_grad)

fig_name<-paste("MLP_All_Genera_Key_KOs_Metab_Raster_Plot.pdf",sep="_")
ggsave(fig_name,plot=raster_fig,device=cairo_pdf,width=15,height=45,pointsize=8)
	
selected_genera<-c("Sulfitobacter","MED-G11","Alcanivorax","Marinisoma","UBA1096","Thalassarchaeum","Nitrosopelagicus","Nitrosopumilus")

sub_keyko_data<-keyko_data[which(keyko_data$Genus %in% selected_genera),]

sub_keyko_data$Zone<-factor(sub_keyko_data$Zone,levels=c("Epipelagic","Mesopelagic","Bathypelagic"))

raster_fig<-ggplot(sub_keyko_data,aes(x=Module_Step,y=Full_Taxonomy,fill=Zone))+geom_tile()+scale_fill_manual(name="Zone",values=zone_coloring)+theme_bw()+theme(axis.text.x = element_text(angle = 90,hjust = 1,size=6),axis.text.y = element_text(size=3),strip.text.x.top = element_text(angle = 90),strip.text.y.right = element_text(angle = 0,size=3),legend.position="top")+facet_grid(Phylum ~ Module_Category + Module_Name, space = "free", scales = "free")#+scale_fill_gradientn(colours =depth_col_grad)

fig_name<-paste("MLP_Selected_Genera_Key_KOs_Metab_Raster_Plot.pdf",sep="_")
ggsave(fig_name,plot=raster_fig,device=cairo_pdf,width=10,height=20,pointsize=8)