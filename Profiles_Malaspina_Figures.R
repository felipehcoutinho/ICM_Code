###Libraries
library("GGally")
options(bitmapType='cairo')
library("fpc")
require(maps)
library("vegan")
library(ggsignif)
library(dplyr)
library(tidyr)
library(tibble)

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


###Prophage prevalence plots
pro_data<-read.table(file="/mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/Prophage_Counts/MAG_Info.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)

pro_data$MAG<-gsub("((.)+\\/)|((.)fa$)","",pro_data$MAG,perl=TRUE)

head(pro_data)

pro_mag_data<-merge(mag_data,subset(pro_data,select=c("MAG","Virus_List","Virus_Regions","Prophage_Count"
,"Prophage_Base_Pairs")),by="MAG",all.x=TRUE)

pro_mag_data$Zone<-factor(pro_mag_data$Zone,levels=c("Epipelagic","Mesopelagic","Bathypelagic"))
pro_mag_data$Prophage_Perc<-(pro_mag_data$Prophage_Base_Pairs/pro_mag_data$Bases)*100

summary(pro_mag_data)

fig_S4A<-ggplot(pro_mag_data,aes(y=log2(Prophage_Count+1),x=Phylum,fill=Zone))+geom_boxplot(position = position_dodge(preserve = "single"),lwd=0.3,outlier.size = 0.3)+xlim(taxon_order)+theme_bw()+scale_fill_manual(name="MAG Source Zone",values=zone_coloring)+theme(legend.position='top',text = element_text(size = 9),axis.text.x = element_text(angle = 45,hjust = 1))+labs(x = NULL,y="log2(Prophage Count +1)")
ggsave("Malaspina_Profiles_MAGs_PhylumxZone_Prophage_Count_Boxplot.pdf",plot=fig_S4A,device=cairo_pdf,width=7,height=4,pointsize=8)

fig_S4B<-ggplot(pro_mag_data,aes(y=log10(Prophage_Perc+1),x=Phylum,fill=Zone))+geom_boxplot(position = position_dodge(preserve = "single"),lwd=0.3,outlier.size = 0.3)+xlim(taxon_order)+theme_bw()+scale_fill_manual(name="MAG Source Zone",values=zone_coloring)+theme(legend.position='top',text = element_text(size = 9),axis.text.x = element_text(angle = 45,hjust = 1))+labs(x = NULL,y="Log10(% of Genome identified as Prophage(s))")
ggsave("Malaspina_Profiles_MAGs_PhylumxZone_Prophage_Perc_Boxplot.pdf",plot=fig_S4B,device=cairo_pdf,width=7,height=4,pointsize=8)

fig_S4B<-ggplot(pro_mag_data,aes(y=log10(Prophage_Perc+1),x=Phylum,fill=Zone))+geom_boxplot(position = position_dodge(preserve = "single"),lwd=0.3,outlier.size = 0.3)+xlim(taxon_order)+theme_bw()+scale_fill_manual(name="MAG Source Zone",values=zone_coloring)+theme(legend.position='top',text = element_text(size = 9),axis.text.x = element_text(angle = 45,hjust = 1))+labs(x = NULL,y="Log10(% of Genome identified as Prophage(s))")
ggsave("Malaspina_Profiles_MAGs_PhylumxZone_Prophage_Perc_Boxplot.pdf",plot=fig_S4B,device=cairo_pdf,width=7,height=4,pointsize=8)

fig_S4C<-ggplot(pro_mag_data,aes(y=log10(Prophage_Perc+1),x=Bases))+geom_point(aes(colour=Phylum),alpha=0.9,shape=18,size=3)+geom_smooth(method = "loess",colour="black",se=TRUE,  formula = y ~ x)+theme_bw()+scale_colour_manual(name="Taxon",values=phylum_coloring)+theme(legend.position='top')#+facet_wrap(Phylum ~ .,scales="free",ncol=7)
ggsave("Malaspina_Profiles_MAGs_PhylumxZone_Genome_LengthxProphage_Count_Scatterplot.pdf",plot=fig_S4C,device=cairo_pdf,width=9,height=7,pointsize=8)

min_mag_count<-10
phylum_mag_count<-table(mag_data$Phylum)
passed_phyla<-names(phylum_mag_count[which(phylum_mag_count >= min_mag_count)])
sub_pro_mag_data<-pro_mag_data[which(pro_mag_data$Phylum %in% passed_phyla),]
summary(sub_pro_mag_data)

fig_S4D<-ggplot(sub_pro_mag_data,aes(y=log10(Prophage_Count+1),x=Bases,group=Zone,fill=Zone),colour="black")+geom_point(alpha=0.9,shape=23,size=5)+geom_smooth(method = "lm",se=TRUE,  formula = y ~ x,aes(colour=Zone))+theme_bw()+scale_fill_manual(name="Zone",values=zone_coloring)+scale_colour_manual(name="Zone",values=zone_coloring)+theme(legend.position='top',axis.text.x = element_text(angle = 45,hjust = 1))+scale_x_continuous(labels = function(x) format(x, scientific = TRUE))+labs(x="Assembly length (bp)")+facet_wrap(Phylum ~ .,scales="free",ncol=6)
ggsave("Malaspina_Profiles_MAGs_PhylumxZone_Genome_LengthxProphage_Count_lm_by_PhylumScatterplot.pdf",plot=fig_S4D,device=cairo_pdf,width=21,height=14,pointsize=8)

fig_S4E<-ggplot(sub_pro_mag_data,aes(y=log10(Prophage_Perc+1),x=Bases,group=Zone,fill=Zone),colour="black")+geom_point(alpha=0.9,shape=23,size=5)+geom_smooth(method = "lm",se=TRUE,  formula = y ~ x,aes(colour=Zone))+theme_bw()+scale_fill_manual(name="Zone",values=zone_coloring)+scale_colour_manual(name="Zone",values=zone_coloring)+theme(legend.position='top',axis.text.x = element_text(angle = 45,hjust = 1))+scale_x_continuous(labels = function(x) format(x, scientific = TRUE))+labs(x="Assembly length (bp)")+facet_wrap(Phylum ~ .,scales="free",ncol=6)
ggsave("Malaspina_Profiles_MAGs_PhylumxZone_Genome_LengthxProphage_Perc_lm_by_PhylumScatterplot.pdf",plot=fig_S4E,device=cairo_pdf,width=21,height=14,pointsize=8)


###DefenseFinder Figures
dfinder_data<-read.table(file="/mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/Profiles_Malaspina_MAGs_DenfenseFinder_Ordered_Results.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)

mag_dfinder_data<-merge(subset(mag_data,select=c(MAG,Depth_m,Zone,Estimated_Genome_Size,Cluster_Representative,Total_CDS,Total_KOs,Unique_KOs,TKO_Density,UKO_Density,CDS_Density,Quality,Completeness,Contamination,Strain.heterogeneity,Domain,Phylum,Class,Order,Family,Genus,Species,Contigs,Bases,Max,N50,N90,Sample)),subset(dfinder_data,select=c(MAG,Count_of_Unique_Defense_Systems,Count_of_Unique_Genes_in_Systems,List_of_Defense_SubSystems_in_Genome,List_of_Unique_Genes_in_Systems)),by="MAG",all.x=TRUE)

mag_dfinder_data[is.na(mag_dfinder_data)]<-0

mag_dfinder_data$DS_Prevalence_per_Mbp<-(mag_dfinder_data$Count_of_Unique_Defense_Systems/ mag_dfinder_data$Bases) * 1000000
mag_dfinder_data$DS_Gene_Prevalence_per_PEG_Perc<-(mag_dfinder_data$Count_of_Unique_Genes_in_Systems / mag_dfinder_data$Total_CDS) * 100

summary(mag_dfinder_data)

write.table(mag_dfinder_data,file="/mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/Info_All_MAGs_Profiles_Malaspina_Redo_DF_Summary+MAG_Metadata.tsv",sep="\t",append=FALSE,row.names=FALSE,col.names=TRUE,quote=FALSE)

###
mag_dfinder_data<-read.table(file="/mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/Info_All_MAGs_Profiles_Malaspina_Redo_DF_Summary+MAG_Metadata.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)
summary(mag_dfinder_data)

pro_mag_df_data<-merge(pro_mag_data,subset(mag_dfinder_data,select=c("MAG","Count_of_Unique_Defense_Systems","Count_of_Unique_Genes_in_Systems","List_of_Defense_SubSystems_in_Genome","List_of_Unique_Genes_in_Systems","DS_Prevalence_per_Mbp","DS_Gene_Prevalence_per_PEG_Perc")),by="MAG",all.x=TRUE)

pro_mag_df_data$Zone<-factor(pro_mag_df_data$Zone,levels=c("Epipelagic","Mesopelagic","Bathypelagic"))

summary(pro_mag_df_data)


#pro_mag_df_data$Depth_m<-as.factor(pro_mag_df_data$Depth_m)
pdf("Malaspina_Profiles_Metadata_MAG_DS_Prophages_Correlogram_All_Zones.pdf",width=18,height=18,pointsize=8)
plot<-ggpairs(pro_mag_df_data,columns=c("Bases","Estimated_Genome_Size","Prophage_Count","Prophage_Base_Pairs","Prophage_Perc","Count_of_Unique_Defense_Systems","Count_of_Unique_Genes_in_Systems","DS_Prevalence_per_Mbp","DS_Gene_Prevalence_per_PEG_Perc","Depth_m"),lower = list(continuous = scatter_fn),upper = list(continuous = cor_fn),aes(alpha = 0.5))
print(plot)
dev.off()

pdf("Malaspina_Profiles_Metadata_MAG_DS_Prophages_Correlogram_by_Zone.pdf",width=18,height=18,pointsize=8)
plot<-ggpairs(pro_mag_df_data,columns=c("Bases","Estimated_Genome_Size","Prophage_Count","Prophage_Base_Pairs","Prophage_Perc","Count_of_Unique_Defense_Systems","Count_of_Unique_Genes_in_Systems","DS_Prevalence_per_Mbp","DS_Gene_Prevalence_per_PEG_Perc"),aes(color=Zone,alpha = 0.5),lower = list(continuous = scatter_fn_no_se),upper = list(continuous = cor_fn))+scale_fill_manual(name="Zone",values=zone_coloring)+scale_colour_manual(name="Zone",values=zone_coloring)
print(plot)
dev.off()
###
pro_mag_df_data$Zone<-factor(pro_mag_df_data$Zone,levels=rev(c("Epipelagic","Mesopelagic","Bathypelagic")))
fig5D<-ggplot(pro_mag_df_data,aes(y=Prophage_Count,x=Phylum,fill=Zone))+geom_boxplot(position = position_dodge(preserve = "single"),lwd=0.3,outlier.size = 0.3)+theme_bw()+scale_fill_manual(name="MAG Source Zone",values=zone_coloring)+theme(legend.position='right',text = element_text(size = 9),axis.text.x = element_text(angle = 45,hjust = 1))+labs(x = NULL,y="Prophages in MAGs")+xlim(taxon_order)
fig_name<-"Profiles_Malaspina_MAGs_Defense_Systems_Count_Bar_Plots_by_Zone.pdf"
ggsave(fig_name,plot=fig5D,device=cairo_pdf,width=7,height=5,pointsize=8)

fig5E<-ggplot(pro_mag_df_data,aes(Prophage_Count,fill=Zone))+geom_density(alpha=0.9)+theme_bw()+scale_fill_manual(name="MAG Source Zone",values=zone_coloring)+theme(legend.position='right',text = element_text(size = 9))+labs(y="Prophages in MAGs")
fig_name<-"Profiles_Malaspina_Prophage_Count_Density_Plots_by_Zone.pdf"
ggsave(fig_name,plot=fig5E,device=cairo_pdf,width=7,height=5,pointsize=8)

pro_mag_df_data$Zone<-factor(pro_mag_df_data$Zone,levels=c("Epipelagic","Mesopelagic","Bathypelagic"))
fig5F<-ggplot(pro_mag_df_data,aes(x=Zone,y=Count_of_Unique_Defense_Systems,fill=Zone))+geom_violin(alpha=0.9)+geom_signif(comparisons = list(c("Epipelagic","Mesopelagic"),c("Mesopelagic","Bathypelagic"),c("Epipelagic","Bathypelagic")),map_signif_level = TRUE, textsize=2.5,test = "wilcox.test",step_increase=0.05,tip_length = 0.01,size=0.1,lwd=0.1)+theme_bw()+scale_fill_manual(name="MAG Source Zone",values=zone_coloring)+theme(legend.position='right',text = element_text(size = 9))+labs(y="Defense Systems in MAGs")
fig_name<-"Profiles_Malaspina_MAGs_Defense_Systems_Count_Violin_Plots_by_Zone.pdf"
ggsave(fig_name,plot=fig5F,device=cairo_pdf,width=7,height=5,pointsize=8)

fig5G<-ggplot(pro_mag_df_data,aes(x=Zone,y=Prophage_Count,fill=Zone))+geom_violin(alpha=0.9)+geom_signif(comparisons = list(c("Epipelagic","Mesopelagic"),c("Mesopelagic","Bathypelagic"),c("Epipelagic","Bathypelagic")),map_signif_level = TRUE, textsize=2.5,test = "wilcox.test",step_increase=0.05,tip_length = 0.01,size=0.1,lwd=0.1)+theme_bw()+scale_fill_manual(name="MAG Source Zone",values=zone_coloring)+theme(legend.position='right',text = element_text(size = 9))+labs(y="Prophages in MAGs")
fig_name<-"Profiles_Malaspina_MAGs_Prophage_Count_Violin_Plots_by_Zone.pdf"
ggsave(fig_name,plot=fig5G,device=cairo_pdf,width=7,height=5,pointsize=8)

#f_pro_mag_df_data<-pro_mag_df_data[which(pro_mag_df_data$Prophage_Count > 0 | pro_mag_df_data$Count_of_Unique_Defense_Systems > 0 ),]
#summary(f_pro_mag_df_data)

fig5H<-ggplot(pro_mag_df_data,aes(colour=Zone,x=Prophage_Count,y=Count_of_Unique_Defense_Systems))+geom_jitter(alpha=0.9,shape=18,width=0.25, height=0.25)+geom_smooth(method = "lm",se=TRUE,  formula = y ~ x, color = "black", size=0.5,alpha = 0.25)+theme_bw()+scale_colour_manual(name="MAG Source Zone",values=zone_coloring)+theme(legend.position='right',text = element_text(size = 9))+labs(y="Count of Defense Systems in MAG",x="Prophages in MAG")
fig_name<-"Profiles_Malaspina_MAGs_Prophage_CountxDS_Count_Scatter_Plot_by_Zone.pdf"
ggsave(fig_name,plot=fig5H,device=cairo_pdf,width=7,height=5,pointsize=8)

###
fig5A<-ggplot(mag_dfinder_data,aes(y=Count_of_Unique_Defense_Systems,x=Phylum,fill=Zone))+geom_boxplot(position = position_dodge(preserve = "single"),lwd=0.3,outlier.size = 0.3)+theme_bw()+scale_fill_manual(name="MAG Source Zone",values=zone_coloring)+theme(legend.position='right',text = element_text(size = 9),axis.text.x = element_text(angle = 45,hjust = 1))+labs(x = NULL,y="Defense Systems in MAGs")+xlim(taxon_order)

fig_name<-"Profiles_Malaspina_MAGs_Defense_Systems_Count_Boxplots_by_Phylum_by_Zone.pdf"
ggsave(fig_name,plot=fig5A,device=cairo_pdf,width=7,height=5,pointsize=8)


fig5B<-ggplot(mag_dfinder_data,aes(y=DS_Prevalence_per_Mbp,x=Phylum,fill=Zone))+geom_boxplot(position = position_dodge(preserve = "single"),lwd=0.3,outlier.size = 0.3)+theme_bw()+scale_fill_manual(name="MAG Source Zone",values=zone_coloring)+theme(legend.position='right',text = element_text(size = 9),axis.text.x = element_text(angle = 45,hjust = 1))+labs(x = NULL,y="# DS per Mbp")+xlim(taxon_order)

fig_name<-"Profiles_Malaspina_MAGs_Defense_Systems_Prevalence_per_Mbp_Boxplots_by_Phylum_by_Zone.pdf"
ggsave(fig_name,plot=fig5B,device=cairo_pdf,width=7,height=5,pointsize=8)

fig5C<-ggplot(mag_dfinder_data,aes(y=DS_Gene_Prevalence_per_PEG_Perc,x=Phylum,fill=Zone))+geom_boxplot(position = position_dodge(preserve = "single"),lwd=0.3,outlier.size = 0.3)+theme_bw()+scale_fill_manual(name="MAG Source Zone",values=zone_coloring)+theme(legend.position='right',text = element_text(size = 9),axis.text.x = element_text(angle = 45,hjust = 1))+labs(x = NULL,y="% of CDS assigned to Defense Systems")+xlim(taxon_order)

fig_name<-"Profiles_Malaspina_MAGs_Defense_System_Gene_Percentage_Boxplots_by_Phylum_by_Zone.pdf"
ggsave(fig_name,plot=fig5C,device=cairo_pdf,width=7,height=5,pointsize=8)


#Arrange and print Figure 5
fig_5pt1<-ggarrange(fig5D, fig5A, labels=c("A","B"), ncol=1, nrow=2, common.legend=TRUE)
fig_5pt2<-ggarrange(fig5G, fig5F, fig5H, labels=c("C","D","E"), ncol=3, nrow=1, common.legend=TRUE)

fig_5<-ggarrange(fig_5pt1,fig_5pt2, ncol=1, nrow=2, heights=c(2,1))
ggsave(filename="Malaspina_Profiles_Figure_5.pdf",plot=fig_5,device=cairo_pdf,width=13,height=11,pointsize=8)


###Genome metrics Correlogram
cor_fn <- function(data, mapping, method="p", use="pairwise", ...){

              # grab data
              x <- eval_data_col(data, mapping$x)
              y <- eval_data_col(data, mapping$y)

              # calculate correlation
			  corr_val<-0
			  if ((length(x) > 5) & (length(y) > 5)) {
				corr_val <- cor(x, y, method=method, use="pairwise.complete.obs")
				}
              # calculate colour based on correlation value
              # Here I have set a correlation of minus one to blue, 
              # zero to white, and one to red 
              # Change this to suit: possibly extend to add as an argument of `my_fn`
              colFn <- colorRampPalette(rev(brewer.pal(9,"RdBu")), interpolate ='spline')
              fill <- colFn(100)[findInterval(corr_val, seq(-1, 1, length=100))]

              ggally_cor(data = data, size = 4, colour="black", mapping = mapping, ...) + 
                theme_void() +
                theme(panel.background = element_rect(fill=fill))
            }
			
scatter_fn <- function(data, mapping, method="p", use="pairwise", ...){
	
              # grab data
              x <- eval_data_col(data, mapping$x)
              y <- eval_data_col(data, mapping$y)

	ggally_smooth_loess(data = data, mapping=mapping, ...)+theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

scatter_fn_no_se <- function(data, mapping, method="p", use="pairwise", ...){
	
              # grab data
              x <- eval_data_col(data, mapping$x)
              y <- eval_data_col(data, mapping$y)

	ggally_smooth_loess(data = data, mapping=mapping, se=FALSE, ...)+theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

#text = element_text(size = 25, color="black"),

pdf("Malaspina_Profiles_Metadata_MAG_Genome_MetricsxZone_Correlogram.pdf",width=15,height=15,pointsize=8)
plot<-ggpairs(mag_data,columns=c("Bases","Estimated_Genome_Size","Total_CDS","CDS_Density","Total_KOs","Unique_KOs","TKO_Density","UKO_Density"),aes(color=Zone,alpha = 0.5),lower = list(continuous = scatter_fn),upper = list(continuous = cor_fn))+scale_fill_manual(name="Zone",values=zone_coloring)+scale_colour_manual(name="Zone",values=zone_coloring)
print(plot)
dev.off()

#MAG Subset Correlogram by Phylum
mag_meta<-merge(mag_data,subset(metadata,select=c("Sample","Temp","Conductivity","Fluo","PAR","SPAR","Turb_FTU","Sal_PSU","Salinity_WOA13","NO3_WOA13","PO4_WOA13","SiO4_WOA13","percentPAR","MLD","Oxygen","sigma","O2_umol_kg","O2_corr_umol_kg","O2_sat","AOU_corr_umol_kg","Chla_ugl","Fmax1_resp_prok","Fmax2_resp_euk","Fmax3_tirosina","Fmax4_triptofano","TEP","POC_uM","Turb","pmol_leu","SE","LNA","HNA","All_BT","percentHNA","cell_size","Bacterial_cell_C","Biomass","ugC_l_d","d_1","HNF","low_virus","medium_virus" ,"high_virus","all_virus","VBR")),by="Sample",all.x=TRUE)

sub_vars<-c("Estimated_Genome_Size","CDS_Density","Total_KOs","TKO_Density","Depth_m","Temp","Sal_PSU","NO3_WOA13","PO4_WOA13","SiO4_WOA13","O2_umol_kg")

min_mag_count<-10
phylum_mag_count<-table(mag_data$Phylum)
passed_phyla<-names(phylum_mag_count[which(phylum_mag_count >= min_mag_count)])

sub_mag_meta<-mag_meta[which(mag_meta$Phylum %in% passed_phyla),]

sub_mag_meta$PHIST_Host_2_Phylum[which(sub_mag_meta$PHIST_Host_2_Phylum == "Proteobacteria")]<-sub_mag_meta$PHIST_Host_3_Class[which(sub_mag_meta$PHIST_Host_2_Phylum == "Proteobacteria")]

cor_fn_small <- function(data, mapping, method="p", use="pairwise", ...){

              # grab data
              x <- eval_data_col(data, mapping$x)
              y <- eval_data_col(data, mapping$y)

              # calculate correlation
			  corr_val<-0
			  if ((length(x) > 5) & (length(y) > 5)) {
				corr_val <- cor(x, y, method=method, use="pairwise.complete.obs")
				}
              # calculate colour based on correlation value
              # Here I have set a correlation of minus one to blue, 
              # zero to white, and one to red 
              # Change this to suit: possibly extend to add as an argument of `my_fn`
              colFn <- colorRampPalette(rev(brewer.pal(9,"RdBu")), interpolate ='spline')
              fill <- colFn(100)[findInterval(corr_val, seq(-1, 1, length=100))]

              ggally_cor(data = data, size = 2, colour="black", mapping = mapping, ...) + 
                theme_void() +
                theme(panel.background = element_rect(fill=fill))
            }
			

scatter_fn_no_se <- function(data, mapping, method="p", use="pairwise", ...){
	
              # grab data
              x <- eval_data_col(data, mapping$x)
              y <- eval_data_col(data, mapping$y)

	ggally_smooth_loess(data = data, mapping=mapping, se=FALSE, ...)+theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

pdf("Malaspina_Profiles_Metadata_MAG_Genome_MetricsxMetadataxPhylum_Correlogram.pdf",width=18,height=18,pointsize=8)
plot<-ggpairs(sub_mag_meta,columns=sub_vars,aes(color=Phylum,alpha = 0.5),lower = list(continuous = scatter_fn_no_se),upper = list(continuous = cor_fn_small))+scale_fill_manual(name="Taxon",values=phylum_coloring)+scale_colour_manual(name="Taxon",values=phylum_coloring)
print(plot)
dev.off()


###MAG Info Assembly Size/Estimated Genome Size x Depth Plots
subset_taxon_level<-"Phylum"
subset_taxon_name<-"Proteobacteria"
group_taxon_level<-"Order"

levels(mag_data$Phylum) <- c(levels(mag_data$Phylum),"Proteobacteria")
mag_data$Phylum[which(grepl("proteobacteria",mag_data$Phylum))] <- "Proteobacteria"

sub_mag_data<-mag_data[which(mag_data[[subset_taxon_level]] == subset_taxon_name),]

min_mag_count<-5

taxa_mag_count<-table(sub_mag_data[[group_taxon_level]])

print(taxa_mag_count)

passed_taxa<-names(taxa_mag_count[which(taxa_mag_count >= min_mag_count)])

print(passed_taxa)

f_mag_data<-sub_mag_data[which(sub_mag_data[[group_taxon_level]] %in% passed_taxa),]

summary(f_mag_data)

fig_2A_subpanel<-ggplot(f_mag_data,aes(y=Estimated_Genome_Size,x=Depth_m))+geom_point(aes(fill=Zone),shape=21)+geom_smooth(method = "loess",se=TRUE,  formula = y ~ x,, color = "black", size=0.5,alpha = 0.25)+theme_bw()+scale_fill_manual(name="Zone",values=zone_coloring)+theme(legend.position='none',text = element_text(size = 7))+labs(x = "Depth (m)",y="Estimated Genome Size (bp)")+scale_y_continuous(labels = function(x) format(x, scientific = TRUE))+facet_wrap(Class ~ Order, scales="free_y", nrow=5)

ggsave("Malaspina_Profiles_Figure_2A_subpanel.pdf",plot=fig_2A_subpanel,device=cairo_pdf,width=10,height=10,pointsize=8)


###CDS, UKO, TKO counts and Density plots
fig_2B_subpanel<-ggplot(f_mag_data,aes(y=CDS_Density,x=Order,fill=Zone))+geom_boxplot(position = position_dodge(preserve = "single"),lwd=0.3,outlier.size = 0.3)+theme_bw()+scale_fill_manual(name="MAG Source Zone",values=zone_coloring)+theme(legend.position='right',text = element_text(size = 9),axis.text.x = element_text(angle = 45,hjust = 1))+labs(x = NULL,y="Density (PEGs/Mbp)")+facet_wrap(Class ~ ., scales="free_x", nrow=5)

ggsave("Malaspina_Profiles_Figure_2B_subpanel.pdf",plot=fig_2B_subpanel,device=cairo_pdf,width=10,height=10,pointsize=8)

fig_2C_subpanel<-ggplot(f_mag_data,aes(y=Total_KOs,x=Order,fill=Zone))+geom_boxplot(position = position_dodge(preserve = "single"),lwd=0.3,outlier.size = 0.3)+theme_bw()+scale_fill_manual(name="MAG Source Zone",values=zone_coloring)+theme(legend.position='right',text = element_text(size = 9),axis.text.x = element_text(angle = 45,hjust = 1))+labs(x = NULL,y="Total KOs")+facet_wrap(Class ~ ., scales="free_x", nrow=5)

ggsave("Malaspina_Profiles_Figure_2C_subpanel.pdf",plot=fig_2C_subpanel,device=cairo_pdf,width=10,height=10,pointsize=8)


###MAG Info Assembly Size/Estimated Genome Size x Depth Plots 
mag_data$Phylum[which(mag_data$Phylum == "Proteobacteria")]<-mag_data$Class[which(mag_data$Phylum == "Proteobacteria")]

min_mag_count<-10

phylum_mag_count<-table(mag_data$Phylum)

passed_phyla<-names(phylum_mag_count[which(phylum_mag_count >= min_mag_count)])

f_mag_data<-mag_data[which(mag_data$Phylum %in% passed_phyla),]

sub_phylum_coloring<-phylum_coloring[which(names(phylum_coloring) %in% passed_phyla)]

sub_taxon_order<-taxon_order[which(taxon_order %in% passed_phyla)]

f_mag_data$Phylum<-factor(f_mag_data$Phylum,levels=sub_taxon_order)


fig_2A<-ggplot(f_mag_data,aes(y=Estimated_Genome_Size,x=Depth_m))+geom_point(aes(fill=Zone),shape=21)+geom_smooth(method = "loess",se=TRUE,  formula = y ~ x,, color = "black", size=0.5,alpha = 0.25)+theme_bw()+scale_fill_manual(name="Zone",values=zone_coloring)+theme(legend.position='none',text = element_text(size = 7))+labs(x = "Depth (m)",y="Estimated Genome Size (bp)")+scale_y_continuous(labels = function(x) format(x, scientific = TRUE))+facet_wrap(Phylum ~ .,scales="free_y", nrow=3)

ggsave("Malaspina_Profiles_Figure_2A_alt_2.pdf",plot=fig_2A,device=cairo_pdf,width=10,height=5,pointsize=8)


fig_2A_alt_2<-ggplot(f_mag_data,aes(y=Estimated_Genome_Size,x=Depth_m))+geom_point(aes(fill=Phylum),shape=21)+geom_smooth(method = "loess",se=TRUE,  formula = y ~ x,, color = "black", size=0.5,alpha = 0.25)+theme_bw()+scale_fill_manual(name="Taxon",values=sub_phylum_coloring,breaks=sub_taxon_order)+theme(legend.position='none',text = element_text(size = 7))+labs(x = "Depth (m)",y="Estimated Genome Size (bp)")+scale_y_continuous(labels = function(x) format(x, scientific = TRUE))+facet_wrap(Phylum ~ .,scales="free_y", nrow=3)

ggsave("Malaspina_Profiles_Figure_2A.pdf",plot=fig_2A_alt_2,device=cairo_pdf,width=10,height=5,pointsize=8)

fig_2A_alt_1<-ggplot(f_mag_data,aes(y=Estimated_Genome_Size,x=Depth_m,group=Phylum,colour=Phylum))+geom_point(alpha=0.7)+geom_smooth(method = "loess",se=FALSE,  formula = y ~ x)+theme_bw()+labs(x = "Depth (m)",y="Estimated Genome Size (bp)")+scale_x_reverse(limits = c(4000, 0))+scale_y_continuous(trans='log2',labels = function(x) format(x, scientific = TRUE))+theme(legend.position='right')+guides(col = guide_legend(ncol = 1))+coord_flip()+scale_colour_manual(name="Taxon",values=sub_phylum_coloring,breaks=sub_taxon_order)#+theme(text = element_text(size = 15))+scale_y_continuous(labels = function(x) format(x, scientific = TRUE))

ggsave("Malaspina_Profiles_Figure_2A_alt_1.pdf",plot=fig_2A_alt_1,device=cairo_pdf,width=7,height=6,pointsize=8)

fig_S<-ggplot(mag_data,aes(y=Estimated_Genome_Size,x=Phylum,fill=Zone))+geom_boxplot(position = position_dodge(preserve = "single"),lwd=0.3,outlier.size = 0.3)+theme_bw()+scale_fill_manual(name="MAG Source Zone",values=zone_coloring)+theme(legend.position='top',text = element_text(size = 9),axis.text.x = element_text(angle = 45,hjust = 1))+labs(x = NULL,y="Estimated Genome Size")

ggsave("Malaspina_Profiles_EGSxPhylumxZone_Boxplots.pdf",plot=fig_S,device=cairo_pdf,width=7,height=5,pointsize=8)


###CDS, UKO, TKO counts and Density plots
#mag_data<-read.table(file="/mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/Info_All_MAGs_Profiles_Malaspina_Redo.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)

mag_data$Phylum[which(mag_data$Phylum == "Proteobacteria")]<-mag_data$Class[which(mag_data$Phylum == "Proteobacteria")]

min_mag_count<-10

phylum_mag_count<-table(mag_data$Phylum)

passed_phyla<-names(phylum_mag_count[which(phylum_mag_count >= min_mag_count)])

mag_data$Phylum<-factor(mag_data$Phylum,levels=sub_taxon_order)

mag_data$Zone<-factor(mag_data$Zone,levels=c("Epipelagic","Mesopelagic","Bathypelagic"))

f_mag_data<-mag_data[which(mag_data$Phylum %in% passed_phyla),]

sub_phylum_coloring<-phylum_coloring[which(names(phylum_coloring) %in% passed_phyla)]

sub_taxon_order<-taxon_order[which(taxon_order %in% passed_phyla)]

fig_2B<-ggplot(f_mag_data,aes(y=CDS_Density,x=Phylum,fill=Zone))+geom_boxplot(position = position_dodge(preserve = "single"),lwd=0.3,outlier.size = 0.3)+theme_bw()+scale_fill_manual(name="MAG Source Zone",values=zone_coloring)+theme(legend.position='right',text = element_text(size = 9),axis.text.x = element_text(angle = 45,hjust = 1))+labs(x = NULL,y="Density (PEGs/Mbp)")

ggsave("Malaspina_Profiles_Figure_2B.pdf",plot=fig_2B,device=cairo_pdf,width=7,height=4,pointsize=8)

fig_2C<-ggplot(f_mag_data,aes(y=Total_KOs,x=Phylum,fill=Zone))+geom_boxplot(position = position_dodge(preserve = "single"),lwd=0.3,outlier.size = 0.3)+theme_bw()+scale_fill_manual(name="MAG Source Zone",values=zone_coloring)+theme(legend.position='right',text = element_text(size = 9),axis.text.x = element_text(angle = 45,hjust = 1))+labs(x = NULL,y="Total KOs")

ggsave("Malaspina_Profiles_Figure_2C.pdf",plot=fig_2C,device=cairo_pdf,width=7,height=4,pointsize=8)

#Fig2
fig_2<-ggarrange(fig_2A, ggarrange(fig_2B, fig_2C, labels=c("B","C"), ncol=2, nrow=1, common.legend=TRUE), labels = c("A"), ncol=1, nrow=2, common.legend=FALSE, heights=c(1.2,1))

ggsave(filename="Malaspina_Profiles_Figure_2.pdf",plot=fig_2,device=cairo_pdf,width=10,height=8,pointsize=8)


###Generate Bin/MAG data
bin_data<-read.table(file="/mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/Info_All_Bins_Profiles_Malaspina_Redo.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)

rownames(bin_data)<-bin_data[,1]

bin_data$Quality<-(bin_data$Completeness - (3*(bin_data$Contamination)))

mag_data<-bin_data[which(bin_data$Completeness >= 50 & bin_data$Contamination <= 5),]

mag_data$Estimated_Genome_Size<-(mag_data$Bases/mag_data$Completeness)*100

mag_data$classification<-gsub("((d__)|(p__)|(c__)|(o__)|(f__)|(g__)|(s__))","",mag_data$classification,perl=TRUE)

mag_data<-mag_data %>% separate(classification, c("Domain","Phylum","Class","Order","Family","Genus","Species"),sep=";")

mag_data$MAG<-mag_data$user_genome
mag_data$Bin<-mag_data$user_genome

nr_mags<-read.table(file="/mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/List_dREP_genomes.txt",sep="\t",header=FALSE,quote="",comment="",stringsAsFactors=TRUE)
colnames(nr_mags)[1]<-"MAG"
nr_mags$MAG<-gsub(".fa$","",nr_mags$MAG,perl=TRUE)
nr_mags$Cluster_Representative<-TRUE

mag_data<-merge(mag_data,nr_mags,by="MAG",all.x=TRUE)
mag_data$Cluster_Representative[is.na(mag_data$Cluster_Representative)]<-FALSE

summary(mag_data)
#write.table(mag_data,file="/mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/Info_All_MAGs_Profiles_Malaspina_Redo.tsv",sep="\t",append=FALSE,row.names=FALSE,col.names=TRUE,quote=FALSE)

###Taxon x Zone barplot with SD
scaff_data<-read.table(file="/mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/Profiles_Malaspina_MAG_Scaffolds_Info_Redo.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)

#scaff_data<-read.table(file="/mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/short.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)

dim(scaff_data)
#summary(scaff_data)
summary(is.na(scaff_data$Bin))
unique(scaff_data$Bin)[1:10]

scaff_data<-scaff_data %>% separate(Bin, c("Bin_1","Bin_2"),sep=",", extra = "drop", fill = "right")
#scaff_data<-scaff_data %>% separate(Bin, c("Bin_1","Bin_2"),sep=",")
scaff_data$Bin<-NULL
scaff_data$Bin<-scaff_data$Bin_1
scaff_data$Bin[grepl("Refined",scaff_data$Bin_2)]<-scaff_data$Bin_2[grepl("Refined",scaff_data$Bin_2)]

scaff_data$Bin<-gsub("(\\[)|(\\])|(\\')|(\\s)+","",scaff_data$Bin,perl=TRUE)

colnames(scaff_data)[6]<-c("Sequence")

#summary(scaff_data)
summary(is.na(scaff_data$Bin))
unique(scaff_data$Bin)[1:10]

scaff_abund<-read.table(file="/mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/Abundance_dRep_MAGs/RPKM_Abundance_Merged_Sequences.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)

#summary(scaff_abund)

scaff_bin_abund<-merge(scaff_abund,subset(scaff_data,select=c("Sequence","Bin")),by="Sequence",all.x=TRUE)

#summary(scaff_bin_abund)
summary(is.na(scaff_bin_abund$Bin))
unique(scaff_bin_abund$Bin)[1:10]

tax_level<-"Phylum"
e_mag_data<-mag_data
if (tax_level == "Phylum") {
	e_mag_data<-mag_data
	e_mag_data$Phylum[which(e_mag_data$Phylum == "Proteobacteria")]<-e_mag_data$Class[which(e_mag_data$Phylum == "Proteobacteria")]
}

scaff_bin_tax_abund<-merge(scaff_bin_abund,subset(e_mag_data,select=c("Bin",tax_level)),by="Bin",all.x=TRUE)

summary(scaff_bin_tax_abund)

m_scaff_bin_tax_abund<-melt(scaff_bin_tax_abund,id=c("Sequence","Bin",tax_level))
colnames(m_scaff_bin_tax_abund)[4]<-"Sample"
colnames(m_scaff_bin_tax_abund)[5]<-"Abundance"

summary(m_scaff_bin_tax_abund)

fm_scaff_bin_tax_abund<-m_scaff_bin_tax_abund[which(m_scaff_bin_tax_abund$Abundance > 0),]

dim(fm_scaff_bin_tax_abund)
summary(fm_scaff_bin_tax_abund)

miss_data<-fm_scaff_bin_tax_abund[is.na(fm_scaff_bin_tax_abund[[tax_level]]),]
dim(miss_data)

fm_scaff_bin_tax_meta_abund<-fm_scaff_bin_tax_abund#merge(fm_scaff_bin_tax_abund,subset(metadata,select=c("Sample")),by="Sample",all.x=TRUE)

dim(fm_scaff_bin_tax_meta_abund)
head(fm_scaff_bin_tax_meta_abund)

sum_fm_scaff_bin_tax_meta_abund<-as.data.frame(aggregate(fm_scaff_bin_tax_meta_abund$Abundance, by=list(Taxon=fm_scaff_bin_tax_meta_abund[[tax_level]],Sample=fm_scaff_bin_tax_meta_abund$Sample), FUN=sum))

colnames(sum_fm_scaff_bin_tax_meta_abund)[3]<-"Abundance"
dim(sum_fm_scaff_bin_tax_meta_abund)
summary(sum_fm_scaff_bin_tax_meta_abund)
head(sum_fm_scaff_bin_tax_meta_abund)

sum_fm_scaff_bin_tax_abund_meta<-merge(sum_fm_scaff_bin_tax_meta_abund,subset(metadata,select=c("Sample","Zone")),by="Sample",all.x=TRUE)
dim(sum_fm_scaff_bin_tax_abund_meta)
summary(sum_fm_scaff_bin_tax_abund_meta)
head(sum_fm_scaff_bin_tax_abund_meta)

mn_tax_ab<-as.data.frame(aggregate(sum_fm_scaff_bin_tax_abund_meta$Abundance, by=list(Taxon=sum_fm_scaff_bin_tax_abund_meta$Taxon,Zone=sum_fm_scaff_bin_tax_abund_meta$Zone), FUN=mean))
colnames(mn_tax_ab)[3]<-"Mean_Abundance"
summary(mn_tax_ab)

sd_tax_ab<-as.data.frame(aggregate(sum_fm_scaff_bin_tax_abund_meta$Abundance, by=list(Taxon=sum_fm_scaff_bin_tax_abund_meta$Taxon,Zone=sum_fm_scaff_bin_tax_abund_meta$Zone), FUN=sd))
colnames(sd_tax_ab)[3]<-"SD_Abundance"
summary(sd_tax_ab)

se_tax_ab<-as.data.frame(aggregate(sum_fm_scaff_bin_tax_abund_meta$Abundance, by=list(Taxon=sum_fm_scaff_bin_tax_abund_meta$Taxon,Zone=sum_fm_scaff_bin_tax_abund_meta$Zone), FUN=std.error))
colnames(se_tax_ab)[3]<-"SE_Abundance"
summary(se_tax_ab)

full_tax_ab<-merge(mn_tax_ab,sd_tax_ab,by=c("Taxon","Zone"),all.x=TRUE)
#full_tax_ab<-merge(mn_tax_ab,merge(se_tax_ab,sd_tax_ab,by=c("Taxon","Zone"),all.x=TRUE),by=c("Taxon","Zone"),all.x=TRUE)

dim(full_tax_ab)
summary(full_tax_ab)
write.table(full_tax_ab,file="/mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/Abundance_dRep_MAGs/Abundance_Stats_by_Phylum_from_RPKM_Abundance_Merged_Sequences.tsv",sep="\t",append=FALSE,row.names=FALSE,col.names=TRUE,quote=FALSE)


full_tax_ab<-read.table(file="/mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/Abundance_dRep_MAGs/Abundance_Stats_by_Phylum_from_RPKM_Abundance_Merged_Sequences.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)

full_tax_ab$Zone<-factor(full_tax_ab$Zone,levels=c("Epipelagic","Mesopelagic","Bathypelagic"))

sub_taxon_order<-taxon_order[which(taxon_order %in% unique(full_tax_ab$Taxon))]

#sub_tax_ab<-full_tax_ab[which(full_tax_ab$Mean_Abundance > 0),]

#+scale_y_log10(limits=c(0.01,100000),breaks=c(0.01,0.1,1,10,100,1000,10000,100000),minor_breaks=c(1 %o% 10^(-2:4)),oob=scales::squish((full_tax_ab$Mean_Abundance, range = c(0.01, 100000), only.finite = TRUE))
#+scale_y_log10(limits=c(0.01,100000),breaks=c(0.01,0.1,1,10,100,1000,10000,100000))+scale_y_continuous(trans = pseudo_log_trans(base = 10))

fig_S3<-ggplot(full_tax_ab,aes(x=Taxon,y=Mean_Abundance,fill=Zone))+theme_bw()+geom_bar(position="dodge",stat="identity",colour="black",alpha=0.9,size=0.1)+scale_fill_manual(name="Zone",values=zone_coloring)+geom_errorbar(position="dodge",stat = "identity",aes(ymin=Mean_Abundance-SD_Abundance,ymax=Mean_Abundance+SD_Abundance),na.rm=TRUE,size=0.4,colour="#999999",alpha=0.5)+xlim(sub_taxon_order)+theme(axis.text.x = element_text(angle = 45,hjust = 1,size=10),legend.position="top")+labs(y="Mean abundance ± s.d. (RPKM)",x=NULL)

ggsave("/mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/Abundance_dRep_MAGs/Malaspina_Profiles_Figure_S3.pdf",plot=fig_S3,device=cairo_pdf,width=10,height=5,pointsize=8)

full_tax_ab$Taxon<-factor(full_tax_ab$Taxon,levels=sub_taxon_order)

fig_S3_alt1<-ggplot(full_tax_ab,aes(x=Taxon,y=Mean_Abundance,fill=Zone))+theme_bw()+geom_bar(position="dodge",stat="identity",colour="black",alpha=0.9)+scale_fill_manual(name="Zone",values=zone_coloring)+geom_errorbar(position="dodge",stat = "identity",aes(ymin=Mean_Abundance-SD_Abundance,ymax=Mean_Abundance+SD_Abundance),na.rm=TRUE,size=0.4,colour="#999999")+theme(axis.text.x=element_blank(),legend.position="top")+labs(y="Mean abundance ± s.d. (RPKM)",x=NULL)+facet_wrap(Taxon ~ .,scales="free",ncol=7)

ggsave("/mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/Abundance_dRep_MAGs/Malaspina_Profiles_Figure_S3_alt1.pdf",plot=fig_S3_alt1,device=cairo_pdf,width=14,height=9,pointsize=8)


###Metabolic module completeness x MAG calculation
func_data<-read.table(file="/mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/Metabolic_Preds/Profiles_Malaspina_MAGs_METABOLIC_Results/METABOLIC_result_each_spreadsheet/METABOLIC_result_worksheet4.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)
dim(func_data)

long_func_data<-melt(func_data,id=c("Module.step","Module","KO.id","Module.Category"))

colnames(long_func_data)<-c("Module_Step","Module_Name","KO_IDs","Module_Category","MAG","Module_Step_Presence")

long_func_data$MAG<-gsub(".Module.step.presence$","",long_func_data$MAG,perl=TRUE)

long_func_data<-merge(long_func_data,subset(mag_data,select=c("MAG","Domain","Phylum","Class","Order","Family","Genus","Species","Sample","Depth_m","Zone")),id=c("MAG"),all.x=TRUE)

summary(long_func_data)

write.table(long_func_data,file="/mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/Metabolic_Preds/Long_Format_Profiles_Malaspina_MAGs_METABOLIC_result_worksheet4.tsv",sep="\t",append=FALSE,row.names=FALSE,col.names=TRUE,quote=FALSE)

f_long_func_data<-long_func_data[which(long_func_data$Module_Step_Presence == "Present"),]

summary(f_long_func_data)

write.table(f_long_func_data,file="/mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/Metabolic_Preds/Long_Format_No_Absent_Profiles_Malaspina_MAGs_METABOLIC_result_worksheet4.tsv",sep="\t",append=FALSE,row.names=FALSE,col.names=TRUE,quote=FALSE)

colnames(func_data)<-gsub(".Module.step.presence","",colnames(func_data),perl=TRUE)

func_data[,c(5:ncol(func_data))]<-data.frame(apply(func_data[,c(5:ncol(func_data))],2,function(x) as.numeric(sub("Absent",0,sub("Present",1,x)))))

func_data$Module_UID<-gsub("\\+(\\d)+$","",func_data$Module.step,perl=TRUE)

m_func_data<-melt(func_data,id=c("Module.step","Module","KO.id","Module.Category","Module_UID"))

dim(m_func_data)
summary(m_func_data)

mag_step_counts<-as.data.frame(aggregate(m_func_data$value, by=list(MAG=m_func_data$variable,Module_UID=m_func_data$Module_UID,Module=m_func_data$Module,Module_Category=m_func_data$Module.Category), FUN=sum))
colnames(mag_step_counts)[5]<-"Step_Count"#c("Module_UID","Total_Steps")

mag_step_counts$Module_Category[which(m_func_data$Module_Category == "")]<-NA

dim(mag_step_counts)
summary(mag_step_counts)

module_total_steps<-as.data.frame(table(func_data$Module_UID))
colnames(module_total_steps)<-c("Module_UID","Total_Steps")

mag_step_counts<-merge(mag_step_counts,module_total_steps,by=c("Module_UID"),all.x=TRUE)

mag_step_counts$Module_Completeness<-(mag_step_counts$Step_Count/mag_step_counts$Total_Steps)*100

mag_step_counts<-merge(mag_step_counts,subset(mag_data,select=c("MAG","Domain","Phylum","Class","Order","Family","Genus","Species","Sample","Depth_m","Zone")),by=c("MAG"),all.x=TRUE)

summary(mag_step_counts)

write.table(mag_step_counts,file="/mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/Metabolic_Preds/Profiles_Malaspina_MAGs_METABOLIC_Module_Completeness_From_worksheet4.tsv",sep="\t",append=FALSE,row.names=FALSE,col.names=TRUE,quote=FALSE)

f_mag_step_counts<-mag_step_counts[which(mag_step_counts$Module_Completeness > 0),]

write.table(f_mag_step_counts,file="/mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/Metabolic_Preds/Profiles_Malaspina_MAGs_METABOLIC_Module_Non_Zero_Completeness_From_worksheet4.tsv",sep="\t",append=FALSE,row.names=FALSE,col.names=TRUE,quote=FALSE)

#
tax_compl_metrics<-as.data.frame(aggregate(mag_step_counts$Module_Completeness, by=list(Taxon=mag_step_counts$Phylum,Module_UID=mag_step_counts$Module_UID,Module=mag_step_counts$Module,Module_Category=mag_step_counts$Module_Category, Zone=mag_step_counts$Zone), FUN=mean))

colnames(tax_compl_metrics)[6]<-"Mean_Module_Completeness"

summary(tax_compl_metrics)

write.table(tax_compl_metrics,file="/mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/Metabolic_Preds/Profiles_Malaspina_PhylumxZone_METABOLIC_Module_Completeness_Mean_From_worksheet4.tsv",sep="\t",append=FALSE,row.names=FALSE,col.names=TRUE,quote=FALSE)

f_tax_compl_metrics<-tax_compl_metrics[which(tax_compl_metrics$Mean_Module_Completeness > 0),]

#f_tax_compl_metrics<-tax_compl_metrics[which(tax_compl_metrics$Module_Category == "Cofactor and vitamin metabolism" & tax_compl_metrics$Mean_Module_Completeness > 0),]

zone_coloring<-c(brewer.pal(9,"YlOrRd")[c(4)],brewer.pal(9,"YlGnBu")[c(5,9)])
names(zone_coloring)<-c("Epipelagic","Mesopelagic","Bathypelagic")

f_tax_compl_metrics$Zone<-factor(f_tax_compl_metrics$Zone,levels=rev(c("Epipelagic","Mesopelagic","Bathypelagic")))

col_pal<-brewer.pal(11,"Spectral")
depth_col_grad<-rev(colorRampPalette(col_pal)(n=299))

fig2<-ggplot(f_tax_compl_metrics,aes(x=Module,y=Zone,fill=Mean_Module_Completeness))+geom_tile()+theme_bw()+theme(axis.text.x = element_text(angle = 90,hjust = 1,size=9),strip.text.y.right = element_text(angle = 0),strip.text.x.top = element_text(angle = 90))+scale_fill_gradientn(colours =depth_col_grad)+facet_grid(Taxon ~ Module_Category, space = "free", scales = "free")

ggsave("Malaspina_Profiles_PhylumxModule_CompletenessxDepth_Raster.pdf",plot=fig2,device=cairo_pdf,width=40,height=25,pointsize=8)


fig1<-ggplot(f_tax_compl_metrics,aes(x=Module,y=Mean_Module_Completeness,fill=Zone))+geom_bar(position = position_dodge(preserve = "single"),stat="identity",colour="black",alpha=0.9)+theme_bw()+theme(axis.text.x = element_text(angle = 90,hjust = 1,size=9),strip.text.y.right = element_text(angle = 0),strip.text.x.top = element_text(angle = 90))+scale_fill_manual(name="Zone",values=zone_coloring)+facet_grid(Taxon ~ Module_Category, space = "free", scales = "free_x")

ggsave("Malaspina_Profiles_PhylumxModule_CompletenessxDepth_Barplots.pdf",plot=fig1,device=cairo_pdf,width=30,height=30,pointsize=8)

fig3<-ggplot(f_tax_compl_metrics,aes(x=Module,y=Mean_Module_Completeness,fill=Zone))+geom_bar(position = position_dodge(preserve = "single"),stat="identity",colour="black",alpha=0.9)+theme_bw()+theme(axis.text.x = element_text(angle = 90,hjust = 1,size=9),strip.text.y.right = element_text(angle = 0),strip.text.x.top = element_text(angle = 90))+scale_fill_manual(name="Zone",values=zone_coloring)+facet_wrap(Module_Category ~ Taxon, scales = "free_x")

ggsave("Malaspina_Profiles_PhylumxModule_CompletenessxDepth_Barplots_Wrap.pdf",plot=fig3,device=cairo_pdf,width=30,height=30,pointsize=8)


mag_step_counts$Zone<-factor(mag_step_counts$Zone,levels=c("Epipelagic","Mesopelagic","Bathypelagic"))

#f_mag_step_counts<-mag_step_counts[which(mag_step_counts$Module_Completeness > 0),]

pdf("Profiles_Malaspina_MAGs_Module_Completeness_Boxplot.pdf",width=60,height=60,pointsize=8)
my_plot<-ggplot(f_mag_step_counts)+geom_boxplot(aes(y=Module_Completeness,x=Module,fill=Zone),position = position_dodge(preserve = "single"))+theme_bw()+theme(axis.text.x = element_text(angle = 90,hjust = 1,size=9),axis.text.y = element_text(size=3),strip.text.y.right = element_text(angle = 0),strip.text.x.top = element_text(angle = 90))+scale_fill_manual(name="MAG_Source_Zone",values=zone_coloring)+facet_grid(Phylum ~ Module_Category, space = "free", scales = "free_x")
#+geom_boxplot(position = position_dodge(preserve = "single")) #+geom_signif(comparisons = list(c("Epipelagic","Mesopelagic"),c("Mesopelagic","Bathypelagic"),c("Epipelagic","Bathypelagic")),map_signif_level = FALSE, textsize=3,test = "wilcox.test")
print(my_plot)
dev.off()


#Module Completeness x MAG Tile plot 
ideal_height<-max(15,min(49,length(unique(f_func_data$MAG))/10))

compl_col_grad<-rev(colorRampPalette(brewer.pal(11,"Spectral"))(n=100))

fig2<-ggplot(f_func_data)+geom_tile(aes(x=Module,y=MAG,fill=Module_Completeness,group=Module),stat = "identity",position = "identity", color = "black")+theme_bw()+theme(axis.text.x = element_text(angle = 90,hjust = 1,size=5),axis.text.y = element_text(size=5),strip.text.y.right = element_text(angle = 0),strip.text.x.top = element_text(angle = 90))+scale_fill_gradientn(colours = compl_col_grad)+facet_grid(Family ~ Module_Category, space = "free", scales = "free")

fig_name<-paste("Malaspina_Profiles",tax_level,tax_name,"xModule_CompletenessxDepth_Tile.pdf",sep="_")

ggsave(fig_name,plot=fig2,device=cairo_pdf,width=30,height=ideal_height,pointsize=8)


###Metabolic Heatmap
func_data<-read.table(file="/mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/Metabolic_Preds/Profiles_Malaspina_MAGs_METABOLIC_Results/METABOLIC_result_each_spreadsheet/METABOLIC_result_worksheet3.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)

colnames(func_data)<-gsub(".Module.presence","",colnames(func_data),perl=TRUE)

m_func_data<-melt(func_data,id=c("Module.ID","Module","Module.Category"))

dim(m_func_data)

m_func_data$Module.Category[which(m_func_data$Module.Category == "")]<-NA

dim(m_func_data)

colnames(m_func_data)<-c("Module_ID","Module","Module_Category","MAG","Module_Presence")

m_func_data<-m_func_data[which(m_func_data$Module_Presence == "Present"),]

dim(m_func_data)

bin_data<-read.table(file="/mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/Info_All_Bins_Profiles_Malaspina_Redo.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)
colnames(bin_data)[1]<-"MAG"

m_func_data<-merge(m_func_data,subset(bin_data,select=c("MAG","classification","Sample","Depth_m","Zone")),by=c("MAG"),all.x=TRUE)

dim(m_func_data)

m_func_data<-m_func_data %>% separate(classification, c("Domain","Phylum","Class","Order","Family","Genus","Species"),sep=";")

dim(m_func_data)

m_func_data$MAG_Classification<-paste(m_func_data$MAG,m_func_data$Domain,m_func_data$Phylum,m_func_data$Class,m_func_data$Order,m_func_data$Family,m_func_data$Genus,m_func_data$Species,sep="|")

dim(m_func_data)

m_func_data$MAG_Taxonomy<-paste(m_func_data$Domain,m_func_data$Phylum,m_func_data$Class,m_func_data$Order,m_func_data$Family,m_func_data$Genus,m_func_data$Species,sep="|")

dim(m_func_data)

m_func_data$MAG_Classification<-factor(m_func_data$MAG_Classification,levels=unique(m_func_data$MAG_Classification[order(m_func_data$Zone)]))

dim(m_func_data)
m_func_data$Zone<-factor(m_func_data$Zone,levels=c("Epipelagic","Mesopelagic","Bathypelagic"))

dim(m_func_data)

zone_coloring<-c(brewer.pal(9,"YlOrRd")[c(4)],brewer.pal(9,"YlGnBu")[c(5,9)])
names(zone_coloring)<-c("Epipelagic","Mesopelagic","Bathypelagic")

mag_order<-rev(unique(as.vector(m_func_data$MAG_Classification[order(m_func_data$MAG_Taxonomy,m_func_data$Zone)])))

pdf("Profiles_Malaspina_MAGs_Pathway_Presence_Heatmap_Zone_Ordered.pdf",width=60,height=170,pointsize=8)
my_plot<-ggplot(m_func_data,aes(fill=Zone,y=MAG_Classification,x=Module))+geom_tile()+theme_bw()+theme(axis.text.x = element_text(angle = 90,hjust = 1,size=9),axis.text.y = element_text(size=3),strip.text.y.right = element_text(angle = 0),strip.text.x.top = element_text(angle = 90))+scale_fill_manual(name="MAG_Source_Zone",values=zone_coloring)+facet_grid(Phylum ~ Module_Category, scales = "free", space = "free")#+scale_y_discrete(limits=mag_order)
print(my_plot)
dev.off()

m_func_data$MAG_Classification<-factor(m_func_data$MAG_Classification,levels=unique(m_func_data$MAG_Classification[order(m_func_data$MAG_Taxonomy)]))

m_func_data$Zone<-factor(m_func_data$Zone,levels=c("Epipelagic","Mesopelagic","Bathypelagic"))

#dim(m_func_data)
#m_func_data[c(1:10000),]

pdf("Profiles_Malaspina_MAGs_Pathway_Presence_Heatmap_Taxonomy_Ordered.pdf",width=60,height=170,pointsize=8)
my_plot<-ggplot(m_func_data,aes(fill=Zone,y=MAG_Classification,x=Module))+geom_tile()+theme_bw()+theme(axis.text.x = element_text(angle = 90,hjust = 1,size=9),axis.text.y = element_text(size=3),strip.text.y.right = element_text(angle = 0),strip.text.x.top = element_text(angle = 90))+scale_fill_manual(name="MAG_Source_Zone",values=zone_coloring)+facet_grid(Domain + Phylum ~ Module_Category, scales = "free", space = "free")
print(my_plot)
dev.off()


fdata<-m_func_data[which(m_func_data$Module_Presence == "Present"),]

write.table(fdata,file="/mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/Metabolic_Preds/Profiles_Malaspina_MAGs_Filtered_and_Annotated_METABOLIC_worksheet3.tsv",sep="\t",append=FALSE,row.names=FALSE,col.names=TRUE,quote=FALSE)

fdata<-m_func_data[which(m_func_data$Module_Category == "Aromatics degradation" & m_func_data$Phylum == "p__Bacteroidota" & m_func_data$Module_Presence == "Present"),]

fdata<-m_func_data[which(m_func_data$Module_Category == "Carbon fixation" & m_func_data$Phylum == "p__Patescibacteria" & m_func_data$Module_Presence == "Present"),]


###Vir X Host Abundance Plots
#Viral Scaffold Sequence Data
vir_scaff_data<-read.table(file="/mnt/lustre/scratch/elopez/malaspina_virus_merged_2_RaFAHVPFClass_completetaxonomy.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)
colnames(vir_scaff_data)[1]<-"Virus_ID"
summary(vir_scaff_data)

#Viral Scaffold Raw Abundance Data
raw_vir_scaff_abund<-read.table(file="/mnt/lustre/scratch/elopez/5_bowtie_results/After_checkV_output/counts.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)

rownames(raw_vir_scaff_abund)<-raw_vir_scaff_abund[,1]

perc_vir_scaff_abund<-as.data.frame(t((t(raw_vir_scaff_abund[,-1])/colSums(raw_vir_scaff_abund[,-1]))*100))

summary(perc_vir_scaff_abund)

perc_vir_scaff_abund$Virus_ID<-rownames(perc_vir_scaff_abund)

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

colnames(filtered_phist_preds)<-c("Virus_ID","MAG","Virus_Genome_File","Host_Genome_File","Common_Kmers","PHIST_pvalue","PHIST_adj_pvalue","PHIST_Host_1_Domain","PHIST_Host_2_Phylum","PHIST_Host_3_Class","PHIST_Host_4_Order","PHIST_Host_5_Family","PHIST_Host_6_Genus","PHIST_Host_7_Species","PHIST_Host_Genome_Completeness","PHIST_Host_Genome_Contamination")

#summary(filtered_phist_preds)

full_vir_scaff_data<-merge(vir_scaff_data,filtered_phist_preds,by="Virus_ID",all.x=TRUE)

#summary(full_vir_scaff_data)
#write.table(full_vir_scaff_data,file="/mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Viruses/Profiles_Malaspina_Viruses_Info+PHIST.tsv",sep="\t",append=FALSE,row.names=FALSE,col.names=TRUE,quote=FALSE)

#Updated full vir scaff data
full_vir_scaff_data<-read.table(file="/mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Viruses/Profiles_Malaspina_Viruses_Info+PHIST.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)

###Summarize abundance by PHIST predicted host phylum (or class for the Proteobacteria) using the aggregate function
full_vir_scaff_data$PHIST_Host_2_Phylum[which(full_vir_scaff_data$PHIST_Host_2_Phylum == "Proteobacteria")]<-full_vir_scaff_data$PHIST_Host_3_Class[which(full_vir_scaff_data$PHIST_Host_2_Phylum == "Proteobacteria")]

m_vir_scaff_abund<-melt(vir_scaff_abund,id="Sequence") # FPKM

colnames(m_vir_scaff_abund)<-c("Virus_ID","Sample","RPKM")

m_vir_scaff_abund<-merge(m_vir_scaff_abund,full_vir_scaff_data,by="Virus_ID",all.x=TRUE)

colnames(m_vir_scaff_abund)
summary(m_vir_scaff_abund)
table(m_vir_scaff_abund$Sample)

#backup_m_vir_scaff_abund<-m_vir_scaff_abund
#m_vir_scaff_abund<-backup_m_vir_scaff_abund

m_vir_scaff_abund<-m_vir_scaff_abund[!is.na(m_vir_scaff_abund$PHIST_Host_2_Phylum),]

#m_vir_scaff_abund$Sample[which(m_vir_scaff_abund$Sample == "MP2233")]<-"MP2233bis"
#m_vir_scaff_abund$Sample<-replace(m_vir_scaff_abund$Sample, m_vir_scaff_abund$Sample == "MP2233", "MP2233bis")

levels(m_vir_scaff_abund$Sample)[levels(m_vir_scaff_abund$Sample)=="MP2233"] <- "MP2233bis"

table(m_vir_scaff_abund$Sample)

m_vir_scaff_abund<-merge(m_vir_scaff_abund,metadata,by="Sample",all.x=TRUE)

tax_level<-"PHIST_Host_2_Phylum"


vir_tax_counts<-as.data.frame(aggregate(m_vir_scaff_abund$RPKM, by=list(Taxon=m_vir_scaff_abund[[tax_level]],Sample=m_vir_scaff_abund$Sample), FUN=sum))

colnames(vir_tax_counts)[3]<-"Virus_Abundance"

m_vir_tax_counts<-merge(vir_tax_counts,metadata,by="Sample",all.x=TRUE)

summary(m_vir_tax_counts)

m_vir_tax_counts$Zone<-factor(m_vir_tax_counts$Zone,levels=c("Epipelagic","Mesopelagic","Bathypelagic"))

for (taxon in unique(m_vir_tax_counts$Taxon)) {
	if (taxon %in% names(phylum_coloring)) {
	} else {
		print(paste(taxon," missing from color palette!",sep=""))
	}
	
	if (taxon %in% custom_taxon_order) {
	} else {
		print(paste(taxon," missing from custom taxon order!",sep=""))
	}
}

m_vir_tax_counts<-m_vir_tax_counts[which(vir_tax_counts$Virus_Abundance >= 500),]

sub_phylum_coloring<-phylum_coloring[unique(m_vir_tax_counts$Taxon)]

sub_custom_taxon_order<-custom_taxon_order[custom_taxon_order %in% unique(m_vir_tax_counts$Taxon)]

m_vir_tax_counts$UID<-paste("St",m_vir_tax_counts$Station,"|",m_vir_tax_counts$Depth_m,"m",sep="")

m_vir_tax_counts$UID<-factor(m_vir_tax_counts$UID,levels=unique(m_vir_tax_counts$UID[rev(order(m_vir_tax_counts$Depth_m))]))


fig2B<-ggplot(m_vir_tax_counts,aes(x=Virus_Abundance,y=UID,fill=Taxon))+geom_bar(position="stack",stat="identity",colour="black",alpha=0.9)+theme_bw()+scale_fill_manual(name="Taxon",values=sub_phylum_coloring,breaks=sub_custom_taxon_order)+guides(fill=guide_legend(ncol=1))+facet_grid(Zone ~ ., scales = "free", space = "free")+ylab("Sample")+xlab("RPKM")

ggsave("Malaspina_Profiles_Virus_Abundance_by_PHIST_Phylum_Barplots.pdf",plot=fig2B,device=cairo_pdf,width=10,height=12,pointsize=8)

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



###CDS, UKO, TKO counts and Density Data prep
cds_data<-read.table(file="/mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/Func_Richness/Reduced_Profiles_Malaspina_MAGs_CDS_Info.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)

cds_data$MAG<-gsub("_(\\d)+$","",cds_data$Sequence,perl=TRUE)


count_unique_ko<-function (x) {
	ukc<-length(na.omit(unique(x)))
	return(ukc)
}

mag_ukcs<-as.data.frame(aggregate(cds_data$KEGG_Best_Subject, by=list(MAG=cds_data$MAG), FUN=count_unique_ko))
colnames(mag_ukcs)[2]<-"Unique_KOs"

count_total_ko<-function (x) {
	tkc<-length(na.omit(x))
	return(tkc)
}


mag_tkcs<-as.data.frame(aggregate(cds_data$KEGG_Best_Subject, by=list(MAG=cds_data$MAG), FUN=count_total_ko))
colnames(mag_tkcs)[2]<-"Total_KOs"


count_cds<-function (x) {
	cds_count<-length(x)
	return(cds_count)
}

mag_cds_count<-as.data.frame(aggregate(cds_data$Sequence, by=list(MAG=cds_data$MAG), FUN=count_cds))
colnames(mag_cds_count)[2]<-"Total_CDS"


mag_data<-merge(merge(mag_tkcs,mag_ukcs,by="MAG",all.x=TRUE),merge(mag_data,mag_cds_count,by="MAG",all.x=TRUE),by="MAG",all.x=TRUE)

mag_data$TKO_Density<-(mag_data$Total_KOs/(mag_data$Bases/100000))

mag_data$UKO_Density<-(mag_data$Unique_KOs/(mag_data$Bases/100000))

mag_data$CDS_Density<-(mag_data$Total_CDS/(mag_data$Bases/100000))

write.table(mag_data,file="/mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/Info_All_MAGs_Profiles_Malaspina_Redo.tsv",sep="\t",append=FALSE,row.names=FALSE,col.names=TRUE,quote=FALSE)

#MAG_Abundance
scaff_data<-read.table(file="/mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/Profiles_Malaspina_MAG_Scaffolds_Info_Redo.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)

scaff_data<-scaff_data %>% separate(Bin, c("Bin_1","Bin_2"),sep=",")
scaff_data$Bin<-scaff_data$Bin_1
scaff_data$Bin[grepl("Refined",scaff_data$Bin_2)]<-scaff_data$Bin_2[grepl("Refined",scaff_data$Bin_2)]

scaff_data$Bin<-gsub("(\\[)|(\\])|(\\')","",scaff_data$Bin,perl=TRUE)

colnames(scaff_data)[6]<-c("Sequence")

summary(scaff_data)

scaff_abund<-read.table(file="/mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/Abundance_dRep_MAGs/RPKM_Abundance_Merged_Sequences.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)

tax_level<-"Phylum"

scaff_abund<-merge(scaff_abund,subset(scaff_data,select=c("Sequence","Bin")),by="Sequence",all.x=TRUE)
scaff_abund<-merge(scaff_abund,subset(mag_data,select=c("Bin",tax_level)),by="Bin",all.x=TRUE)

summary(scaff_abund)

m_scaff_abund<-melt(scaff_abund,id=c("Sequence","Bin",tax_level))

m_scaff_abund<-m_scaff_abund[which(m_scaff_abund$value > 0),]

all_levels_tax_counts<-as.data.frame(aggregate(m_scaff_abund$value, by=list(Taxon=m_scaff_abund[[tax_level]],Sample=m_scaff_abund$variable), FUN=sum))

colnames(all_levels_tax_counts)<-c("Taxon","Sample","Abundance")

sub_tax_counts<-all_levels_tax_counts[which(all_levels_tax_counts$Abundance >= 1000),]

sub_tax_counts<-merge(sub_tax_counts,metadata,by="Sample",all.x=TRUE)

sub_tax_counts$Zone<-factor(sub_tax_counts$Zone,levels=c("Epipelagic","Mesopelagic","Bathypelagic"))

sub_tax_counts$Station<-factor(sub_tax_counts$Station,levels=sort(unique(sub_tax_counts$Station),decreasing=FALSE))

sub_tax_counts$Depth_m<-factor(sub_tax_counts$Depth_m,levels=sort(unique(sub_tax_counts$Depth_m),decreasing=FALSE))

sub_tax_counts$UID<-paste("St",sub_tax_counts$Station,"|",sub_tax_counts$Depth_m,"m",sep="")

sub_tax_counts$UID<-factor(sub_tax_counts$UID,levels=unique(sub_tax_counts$UID[order(sub_tax_counts$Depth_m,sub_tax_counts$Station,decreasing=TRUE)]))

passed_phyla<-unique(sub_tax_counts$Taxon)
sub_phylum_coloring<-phylum_coloring[which(names(phylum_coloring) %in% passed_phyla)]
sub_taxon_order<-taxon_order[which(taxon_order %in% passed_phyla)]

for (taxon in unique(sub_tax_counts$Taxon)) {
	if (taxon %in% names(phylum_coloring)) {
	} else {
		print(paste(taxon," missing from color palette!",sep=""))
	}
}

fig_S2<-ggplot(sub_tax_counts,aes(x=Abundance,y=UID,fill=Taxon))+geom_bar(position="stack",stat="identity",colour="black",alpha=0.9)+theme_bw()+scale_fill_manual(name="Taxon",values=sub_phylum_coloring,breaks=sub_taxon_order)+guides(fill=guide_legend(ncol=1))+facet_grid(Zone ~ ., scales = "free", space="free")

ggsave("Malaspina_Profiles_Figure_S2.pdf",plot=fig_S2,device=cairo_pdf,width=9,height=9,pointsize=8)

###CDS length boxplot
cds_data<-read.table(file="/mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/Func_Richness/Reduced_Profiles_Malaspina_MAGs_CDS_Info.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)

cds_data$MAG<-gsub("_(\\d)+$","",cds_data$Sequence,perl=TRUE)

cds_data<-merge(cds_data,subset(mag_data,select=c("MAG","Domain","Phylum","Class","Order","Family","Genus","Species","Sample","Depth_m","Zone")),by="MAG",all.x=TRUE)

cds_data$Zone<-factor(cds_data$Zone,levels=c("Epipelagic","Mesopelagic","Bathypelagic"))

my_plot_A<-ggplot(cds_data,aes(y=log10(Length),x=Zone,fill=Zone))+geom_boxplot(position = position_dodge(preserve = "single"))+geom_signif(comparisons = list(c("Epipelagic","Mesopelagic"),c("Mesopelagic","Bathypelagic"),c("Epipelagic","Bathypelagic")),map_signif_level = FALSE, textsize=3,test = "wilcox.test", tip_length = 0.01, step_increase=0.1,size=0.7)+theme_bw()+theme(legend.position='none')+ylab("Log10(CDS Length)")+scale_fill_manual(name="MAG Source Zone",values=zone_coloring)

ggsave("Profiles_Malaspina_MAG_CDS_LengthxZone_Boxplot.pdf",plot=my_plot,device=cairo_pdf,width=5,height=5,pointsize=8)

my_plot_B<-ggplot(cds_data,aes(y=log10(Length),x=Phylum,fill=Zone))+geom_boxplot(position = position_dodge(preserve = "single"))+theme_bw()+theme(legend.position='right',axis.text.x = element_text(angle = 45,hjust = 1))+ylab("Log10(CDS Length)")+scale_fill_manual(name="MAG Source Zone",values=zone_coloring)

ggsave("Profiles_Malaspina_MAG_CDS_LengthxPhylumxZone_Boxplot.pdf",plot=my_plot_B,device=cairo_pdf,width=8,height=5,pointsize=8)


###Sampling site map
#Load world map data
world_map <- map_data("world")

nr_st_metadata<-subset(metadata[!duplicated(metadata[,c('Station')]),],select=c("Station","Lat","Long"))

fig_1A<-ggplot(world_map, aes(x = long, y = lat, group = group))+geom_polygon(fill="lightgray", colour = "lightgray")+geom_point(metadata,mapping=aes(y=Lat,x=Long,group=MPCode,colour=Depth_m,size=Depth_m),shape=18)+coord_cartesian(xlim=c(-180,180),ylim=c(-45,45))+geom_text(data=nr_st_metadata,size=4,aes(y=Lat+4,x=Long,label=Station,group=Station))+theme_bw()+xlab("Longitude")+ylab("Latitude")+guides(size = "none")+scale_colour_gradientn(colours = depth_col_grad,limits = c(4000, 0),name="Depth",trans = 'reverse')

ggsave("Malaspina_Profiles_Figure_1A.pdf",plot=fig_1A,device=cairo_pdf,width=8,height=3,pointsize=8)


###PCA
pca_metadata<-subset(metadata,select=c(MPCode,Temp,Sal_PSU,NO3_WOA13,PO4_WOA13,SiO4_WOA13,Oxygen))
rownames(pca_metadata)<-pca_metadata$MPCode

pca_metadata<-na.omit(subset(metadata,select=c(Temp,Sal_PSU,NO3_WOA13,PO4_WOA13,SiO4_WOA13,Oxygen)))

sdata<-as.data.frame(scale(pca_metadata,center=TRUE,scale=TRUE))

set.seed(666)
rdadata<-rda(sdata)

xlabel<-round((summary(rdadata)$cont$importance[2,1]*100),digits=1)
xlabel<-paste("PC1 (",xlabel,"% explained)",sep="")

ylabel<-round((summary(rdadata)$cont$importance[2,2]*100),digits=1)
ylabel<-paste("PC2 (",ylabel,"% explained)",sep="")

uscores <- data.frame(rdadata$CA$u)
#uscores1 <- inner_join(rownames_to_column(pca_metadata), rownames_to_column(data.frame(uscores)), type = "right", by = "rowname")
uscores1<-merge(uscores,metadata,by="row.names",all.x=TRUE)

vscores <- data.frame(rdadata$CA$v)

uscores1$Zone<-factor(uscores1$Zone,levels=c("Epipelagic","Mesopelagic","Bathypelagic"))

rownames(vscores)[rownames(vscores) == 'Sal_PSU'] <- 'Sal'
rownames(vscores)[rownames(vscores) == 'NO3_WOA13'] <- 'NO3'
rownames(vscores)[rownames(vscores) == 'PO4_WOA13'] <- 'PO4'
rownames(vscores)[rownames(vscores) == 'SiO4_WOA13'] <- 'SiO4'
rownames(vscores)[rownames(vscores) == 'Oxygen'] <- 'Oxy'

fig_1B<-ggplot(uscores1)+geom_segment(data = vscores, aes(x = 0, y = 0, xend = PC1, yend = PC2), arrow=arrow(length=unit(0.2,"cm")),alpha = 0.5, color = 'darkgreen')+geom_label(data = vscores, aes(x = PC1, y = PC2, label = rownames(vscores)), fill="grey",col = 'black',size=3,label.padding = unit(0.05, "lines"))+geom_point(data=uscores1,aes(x = PC1, y = PC2, fill=Depth_m),shape=23,alpha=1,size=4)+scale_fill_gradientn(colours =depth_col_grad,limits = c(4000, 0),name="Depth",trans = 'reverse')+labs(y=ylabel,x=xlabel)+theme_bw()+theme(legend.position='none')#+geom_text(data=uscores1,aes(x = PC1, y = PC2, label=Station),size=1.2)

#position=position_jitter(width=0.005,height=0.005)

ggsave("Malaspina_Profiles_Figure_1B.pdf",plot=fig_1B,device=cairo_pdf,width=7,height=5,pointsize=8)



###mTAGs OTU abundance NMDS
all_otu_abund<-read.table(file="/mnt/lustre/bio/shared/malaspina/cnag-metagenomes/malaspina.mTags/output/2019-09-30.1_malaspina-160-metagenomes-cnag-SILVA132-0.99-merged-tables/Merged_table_OTU.txt",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)

colnames(all_otu_abund)<-gsub("\\.(.)+$","",colnames(all_otu_abund),perl=TRUE)
rownames(all_otu_abund)<-all_otu_abund$X
all_otu_abund$X<-NULL

prok_otu_abund<-all_otu_abund[grepl("Bacteria|Archaea", rownames(all_otu_abund)),]

prok_otu_abund<-as.data.frame(t((t(prok_otu_abund)/colSums(prok_otu_abund))*100))

summary(prok_otu_abund)

prok_otu_abund<-as.data.frame(t(prok_otu_abund))

mdata<-merge(prok_otu_abund,subset(metadata,select=c("layer")),by="row.names",all.x=TRUE)

mdata<-mdata[which(mdata$layer != "NA"),]
mdata$layer<-NULL

rownames(mdata)<-mdata$Row.names
mdata$Row.names<-NULL

mdata<-as.data.frame(t(mdata))

summary(mdata)

dist_metric<-"bray"

dists<-vegdist(t(mdata), method = dist_metric)
set.seed(666)
mdsresult<-metaMDS(dists,distance = dist_metric,k = 2,maxit = 999)
data.scores<-as.data.frame(scores(mdsresult))
data.scores$Sample<-rownames(data.scores)

mdata<-merge(data.scores,metadata,by="Sample",all.x=TRUE)

mdata$Zone<-factor(mdata$Zone,levels=c("Epipelagic","Mesopelagic","Bathypelagic"))

summary(mdata)

fig_1C<-ggplot(mdata, aes(x=NMDS1,y=NMDS2))+geom_point(shape=23,size=3.5, stroke = 0.2,aes(fill=Depth_m))+scale_fill_gradientn(colours = depth_col_grad,limits = c(4000, 0),name="Depth",trans = 'reverse')+theme_bw()+theme(text=element_text(size=16))+theme(legend.position='none')

ggsave("Malaspina_Profiles_Figure_1C.pdf",plot=fig_1C,device=cairo_pdf,width=6,height=5,pointsize=8)


fmdata<-subset(mdata[which(mdata$Zone == "Epipelagic"),],select=c(NMDS1,NMDS2,Depth_m,Temp,Sal_PSU,NO3_WOA13,PO4_WOA13,SiO4_WOA13,Oxygen))

pdf("Malaspina_Profiles_nMDS_Epipelagic_Correlogram.pdf",width=15,height=15,pointsize=8)
plot<-ggpairs(fmdata)
print(plot)
dev.off()


fmdata<-subset(mdata,select=c(NMDS1,NMDS2,Depth_m,Temp,Sal_PSU,NO3_WOA13,PO4_WOA13,SiO4_WOA13,Oxygen))

pdf("Malaspina_Profiles_nMDS_All_Samples_Correlogram.pdf",width=15,height=15,pointsize=8)
plot<-ggpairs(fmdata)
print(plot)
dev.off()

fmdata<-subset(mdata[which(mdata$Zone != "Epipelagic"),],select=c(NMDS1,NMDS2,Depth_m,Temp,Sal_PSU,NO3_WOA13,PO4_WOA13,SiO4_WOA13,Oxygen))

pdf("Malaspina_Profiles_nMDS_Meso_and_Bathy_Correlogram.pdf",width=15,height=15,pointsize=8)
plot<-ggpairs(fmdata)
print(plot)
dev.off()


###mTAGs diversity index
all_otu_abund<-read.table(file="/mnt/lustre/bio/shared/malaspina/cnag-metagenomes/malaspina.mTags/output/2019-09-30.1_malaspina-160-metagenomes-cnag-SILVA132-0.99-merged-tables/Merged_table_OTU.txt",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)

colnames(all_otu_abund)<-gsub("\\.(.)+$","",colnames(all_otu_abund),perl=TRUE)
rownames(all_otu_abund)<-all_otu_abund$X
all_otu_abund$X<-NULL

prok_otu_abund<-all_otu_abund[grepl("Bacteria|Archaea", rownames(all_otu_abund)),]

summary(prok_otu_abund)

shanon_div<-as.data.frame(diversity(t(prok_otu_abund), index = "shannon"))
shanon_div$Sample<-rownames(shanon_div)
colnames(shanon_div)[1]<-"Shanon_Diversity_Index"

metadata<-merge(metadata,shanon_div,by="Sample",all.x=TRUE)


fig_1D<-ggplot(metadata,aes(y=Shanon_Diversity_Index,x=Depth_m))+geom_smooth(method = "loess",se=TRUE,  formula = y ~ x)+scale_x_reverse(limits = c(4000, 0))+geom_point(aes(fill=Depth_m),shape=23,size=3.5, stroke = 0.2)+coord_flip()+theme_bw()+scale_fill_gradientn(colours = depth_col_grad,limits = c(4000, 0),name="Depth",trans = 'reverse')+labs(x="Depth",y="Shanon Diversity Index")+theme(legend.position='none')

ggsave("Malaspina_Profiles_Figure_1D.pdf",plot=fig_1D,device=cairo_pdf,width=6,height=5,pointsize=8)

pdf("Profiles_Malaspina_mTAGs_OTUs_Shanon_DiversityxDepth_Boxplot.pdf",width=7,height=7,pointsize=8)
my_plot<-ggplot(metadata)+geom_boxplot(aes(y=Shanon_Diversity_Index,x=Zone,fill=Zone),position = position_dodge(preserve = "single"))+theme_bw()+theme(axis.text.x = element_text(angle = 45,hjust = 1))+scale_fill_manual(name="Zone",values=zone_coloring)
print(my_plot)
dev.off()

###Figure 1 Merged
fig_1<-ggarrange(fig_1A, ggarrange(fig_1B, fig_1C, fig_1D, labels = c("B","C","D"), nrow=1,ncol=3, widths=c(1.2,1,1), common.legend = TRUE,  legend="none", align="h"), nrow=2,ncol=1,labels=c("A"),heights=c(1,1),common.legend = TRUE,  legend="right")

ggsave("Malaspina_Profiles_Figure_1.pdf",plot=fig_1,device=cairo_pdf,width=10,height=7,pointsize=8)

png(filename="Malaspina_Profiles_Figure_1.png",width=34,height=23.8,units="cm",pointsize=8,res=300)
plot(fig_1)
dev.off()

###MAG CDS Density
cds_data<-read.table(file="/mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/Func_Richness/Reduced_Profiles_Malaspina_MAGs_CDS_Info.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)

cds_data$MAG<-gsub("_(\\d)+$","",cds_data$Sequence,perl=TRUE)

mag_data<-read.table(file="/mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/Info_All_Bins_Profiles_Malaspina_Redo.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)

colnames(mag_data)[1]<-c("MAG")

mag_data$classification<-gsub("((d__)|(p__)|(c__)|(o__)|(f__)|(g__)|(s__))","",mag_data$classification,perl=TRUE)

mag_data<-mag_data %>% separate(classification, c("Domain","Phylum","Class","Order","Family","Genus","Species"),sep=";")

count_cds<-function (x) {
	cds_count<-length(x)
	return(cds_count)
}

mag_cds_count<-as.data.frame(aggregate(cds_data$Sequence, by=list(MAG=cds_data$MAG), FUN=count_cds))
colnames(mag_cds_count)[2]<-"Total_CDS"

mag_data<-merge(mag_cds_count,mag_data,by="MAG",all.x=TRUE)

mag_data$CDS_Density<-(mag_data$Total_CDS/(mag_data$Bases/100000))

mag_data$Phylum[which(mag_data$Phylum == "Proteobacteria")]<-mag_data$Class[which(mag_data$Phylum == "Proteobacteria")]

mag_data$Zone<-factor(mag_data$Zone,levels=c("Epipelagic","Mesopelagic","Bathypelagic"))

fig3A<-ggplot(mag_data,aes(y=CDS_Density,x=Zone))+geom_boxplot(aes(fill=Zone))+geom_signif(comparisons = list(c("Epipelagic","Mesopelagic"),c("Mesopelagic","Bathypelagic")),map_signif_level = FALSE, textsize=3,test = "wilcox.test", tip_length = 0.05)+geom_signif(comparisons = list(c("Epipelagic","Bathypelagic")),y_position=max(mag_data$CDS_Density)*1.05,vjust=0.2,map_signif_level = FALSE, textsize=3,test = "wilcox.test", tip_length = 0.05)+theme_bw()+scale_fill_manual(name="MAG_Source_Zone",values=zone_coloring)+theme(legend.position='none',text = element_text(size = 15))+labs(x = "Zone",y="Density (PEGs/Mbp)")

ggsave("Malaspina_Profiles_MAG_Density_by_Zone_Boxplot.pdf",plot=fig3A,device=cairo_pdf,width=5,height=5,pointsize=8)


fig3B<-ggplot(f_mag_data,aes(y=CDS_Density,x=Zone))+geom_boxplot(aes(fill=Zone))+geom_signif(comparisons = list(c("Epipelagic","Mesopelagic"),c("Mesopelagic","Bathypelagic"),c("Epipelagic","Bathypelagic")),map_signif_level = FALSE, textsize=3,test = "wilcox.test", tip_length = 0.01, step_increase=0.05,size=0.7)+theme_bw()+scale_fill_manual(name="MAG_Source_Zone",values=zone_coloring)+theme(legend.position='none',text = element_text(size = 15),axis.text.x = element_text(angle = 45,hjust = 1))+labs(x = "Zone",y="Density (PEGs/Mbp)")+facet_wrap(Phylum ~ ., ncol=6)

ggsave("Malaspina_Profiles_MAG_Density_by_ZonexPhylum_Boxplot.pdf",plot=fig3B,device=cairo_pdf,width=15,height=10,pointsize=8)


fig3C<-ggplot(f_mag_data,aes(y=CDS_Density,x=Depth_m))+geom_point(aes(colour=Phylum))+geom_smooth(method = "loess",se=TRUE,  formula = y ~ x)+theme_bw()+scale_colour_manual(name="Taxon",values=phylum_coloring,breaks=taxon_order)+theme(legend.position='none',text = element_text(size = 15))+labs(x = "Depth (m)",y="Density (PEGs/Mbp)")+facet_wrap(Phylum ~ .,scales="free_y", nrow=3)

ggsave("Malaspina_Profiles_MAG_DensityxDepthxPhylum_Scatterplot.pdf",plot=fig3C,device=cairo_pdf,width=18,height=10,pointsize=8)

pdf("Profiles_Malaspina_MAG_CDS_Density_TaxonxZone_Boxplot.pdf",width=15,height=8,pointsize=8)
my_plot<-ggplot(mag_data)+geom_boxplot(aes(y=CDS_Density,x=Phylum,fill=Zone),position = position_dodge(preserve = "single"))+theme_bw()+theme(axis.text.x = element_text(angle = 45,hjust = 1))+ylab("# CDS/Mbp")+scale_fill_manual(name="MAG_Source_Zone",values=zone_coloring)
print(my_plot)
dev.off()




#####MAG novelty
bin_data<-read.table(file="/mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/Info_All_Bins_Profiles_Malaspina_Redo.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)

mag_data<-bin_data[which(bin_data$Completeness >= 50 & bin_data$Contamination <= 5),]
#colnames(mag_data)[1]<-c("Bin")
mag_data$classification<-gsub("((d__)|(p__)|(c__)|(o__)|(f__)|(g__)|(s__))","",mag_data$classification,perl=TRUE)

mag_data<-mag_data %>% separate(classification, c("Domain","Phylum","Class","Order","Family","Genus","Species"),sep=";")

#summary(mag_data)

mag_data$Source<-"Original"
mag_data$Source[grepl("Refined",mag_data$user_genome)]<-"Refined"

split_level<-"Domain"
level_taxa<-unique(mag_data[[split_level]])

for (taxon in level_taxa) {
	m_mag_data<-melt(subset(mag_data[which(mag_data[[split_level]] == taxon),],select=c("user_genome","Sample","Source","Domain","Phylum","Class","Order","Family","Genus","Species")),id=c("user_genome","Sample","Source"))

	colnames(m_mag_data)<-c("Genome","Sample","Source","Level","Taxon")

	m_mag_data$Novelty<-"Known"
	m_mag_data$Novelty[which(m_mag_data$Taxon == "")]<-"Unknown"
	m_mag_data$Novelty<-as.factor(m_mag_data$Novelty)

	#summary(m_mag_data)

	sub_m_mag_data<-m_mag_data
	set_counts<-as.data.frame(aggregate(sub_m_mag_data$Novelty, by=list(Source=sub_m_mag_data$Source,Level=sub_m_mag_data$Level), FUN=table))
	set_counts<-cbind(subset(set_counts,select=c("Source","Level")),set_counts$x)
	set_counts$Total_MAGs<-rowSums(set_counts[,c(3:4)])

	m_set_counts<-melt(set_counts,id=c("Source","Total_MAGs","Level"))

	m_set_counts$Percentage<-(m_set_counts$value/m_set_counts$Total_MAGs)*100

	colnames(m_set_counts)<-c("Source","Total_MAGs","Level","Novelty","Count","Percentage")

	#summary(m_set_counts)

	m_set_counts$Level<-factor(m_set_counts$Level,levels=rev(c("Domain","Phylum","Class","Order","Family","Genus","Species")))

	fig1A<-ggplot(m_set_counts,aes(y=Count,x=Level,fill=Novelty))+geom_bar(position="stack",stat="identity",colour="black",alpha=0.9)+coord_flip()+theme_bw()+facet_wrap(Source ~ .,scales="free")

	out_name<-paste("Malaspina_Profiles_",taxon,"_MAGs_NoveltyxDataset.pdf",sep="")
	ggsave(out_name,plot=fig1A,device=cairo_pdf,width=12,height=7,pointsize=8)

}

###VIRUS PHIST Prediction Info and Filtering
all_phist_preds<-read.table(file="/mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Viruses/PHIST_Output/Positive_PHIST_Predictions.csv",sep=",",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)

colnames(all_phist_preds)<-c("Virus","MAG","Common_Kmers","pvalue","adj_pvalue")

all_phist_preds$MAG<-gsub("(No_Virus_)|(.fa)","",all_phist_preds$MAG,perl=TRUE)
all_phist_preds$Virus<-gsub(".fasta","",all_phist_preds$Virus,perl=TRUE)

sorted_data<-all_phist_preds[order(all_phist_preds$Virus,all_phist_preds$pvalue,decreasing=FALSE),]

nr_sorted_data<-sorted_data[ !duplicated(sorted_data$Virus), ]

f_nr_sorted_data<-nr_sorted_data[which(nr_sorted_data$pvalue <= 2.384e-14),]

bin_data<-read.table(file="/mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/Profiles_Malaspina_Redo_Bins_Original_and_Refined_CheckM_GTDBtk_Assembly_Stats_Info.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)

mag_data<-bin_data[which(bin_data$Completeness >= 50 & bin_data$Contamination <= 5),]
colnames(mag_data)[1]<-c("MAG")

mag_data<-mag_data %>% separate(classification, c("Domain","Phylum","Class","Order","Family","Genus","Species"),sep=";")

f_nr_sorted_data_info<-merge(f_nr_sorted_data,mag_data,by=c("MAG"),all.x=TRUE)

summary(f_nr_sorted_data_info)


cols_to_keep<-c("MAG","Virus","Common_Kmers","pvalue","adj_pvalue","Domain","Phylum","Class","Order","Family","Genus","Species")

f_nr_sorted_data_info<-f_nr_sorted_data_info[,cols_to_keep]

#File below needs to be checked as it might be pre checkV data
vir_data<-read.table(file="/mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Viruses/Corrected_Abundance/Profiles_Malaspina_Virus_Host_Tax_Info.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)
colnames(vir_data)[1]<-"Virus"

vir_data<-merge(vir_data,f_nr_sorted_data_info,by=c("Virus"),all.x=TRUE)

write.table(vir_data,file="/mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Viruses/PHIST_Output/Profiles_Malaspina_Viruses_Info+Filtered_PHIST_Predictions.tsv",sep="\t",append=FALSE,row.names=FALSE,col.names=TRUE,quote=FALSE)


fdata<-vir_data[which(!is.na(vir_data$Domain) & vir_data$CDS_Count >= 30),]


###mTAGs Class to Phylum abundance
all_class_abund<-read.table(file="/mnt/lustre/bio/shared/malaspina/cnag-metagenomes/malaspina.mTags/output/2019-09-30.1_malaspina-160-metagenomes-cnag-SILVA132-0.99-merged-tables/Merged_table_class.txt",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)

colnames(all_class_abund)<-gsub("\\.(.)+$","",colnames(all_class_abund),perl=TRUE)
rownames(all_class_abund)<-all_class_abund$X
all_class_abund$X<-NULL

prok_class_abund<-all_class_abund[grepl("Bacteria|Archaea", rownames(all_class_abund)),]

prok_class_abund<-as.data.frame(t((t(prok_class_abund)/colSums(prok_class_abund))*100))

summary(prok_class_abund)

taxon<-data.frame(do.call('rbind', strsplit(as.character(rownames(prok_class_abund)),';',fixed=TRUE)))

taxon$X2[which(taxon$X2 == "Proteobacteria")]<-taxon$X3[which(taxon$X2 == "Proteobacteria")]

prok_class_abund$Phylum<-taxon$X2

m_prok_class_abund<-melt(prok_class_abund,id=c("Phylum"))

colnames(m_prok_class_abund)<-c("Phylum","Sample","Abundance")

tax_sample_counts<-aggregate(m_prok_class_abund$Abundance, by=list(Phylum = m_prok_class_abund$Phylum, Sample = m_prok_class_abund$Sample), FUN=sum)

colnames(tax_sample_counts)<-c("Phylum","Sample","Abundance")

metadata<-read.table(file="/mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Metadata/Selected_Malaspina_Profiles_Metadata.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)
metadata$Sample<-metadata$MPCode

mdata<-merge(tax_sample_counts,subset(metadata,select=c("Sample","Station","Depth_m","layer")),by="Sample",all.x=TRUE)

colnames(mdata)<-c("Sample","Taxon","Abundance","Station","Depth_m","Zone")

fmdata<-mdata[which(mdata$Abundance >= 1 & mdata$Zone != "NA"),]

summary(fmdata)

host_coloring<-c(brewer.pal(9,"Set1"),brewer.pal(4,"Pastel1"),brewer.pal(8,"Dark2")[8],brewer.pal(9,"Set1")[2],brewer.pal(8,"Accent"))
names(host_coloring)<-c("Actinobacteria","Alphaproteobacteria","Cyanobacteria","Euryarchaeota","Bacteroidetes","Deltaproteobacteria","Verrucomicrobia","Firmicutes","Fusobacteria","Planctomycetes","Gammaproteobacteria","Betaproteobacteria","Crenarchaeota","Unknown","Proteobacteria","Thaumarchaeota","Chloroflexi","Gemmatimonadetes","Nitrospinae","WPS-2","Patescibacteria","Acidobacteria","Chlamydiae")

for (taxon in unique(fmdata$Taxon)) {
	if (taxon %in% names(host_coloring)) {
	} else {
		print(paste(taxon," missing from color palette!",sep=""))
	}
}

host_coloring<-host_coloring[c(as.vector(unique(fmdata$Taxon)))]

summary(fmdata)

fmdata$Zone<-gsub("Epi","Epipelagic",fmdata$Zone)
fmdata$Zone<-gsub("Meso","Mesopelagic",fmdata$Zone)
fmdata$Zone<-gsub("Bathy","Bathypelagic",fmdata$Zone)

fmdata$Zone<-factor(fmdata$Zone,levels=c("Epipelagic","Mesopelagic","Bathypelagic"))

fmdata$UID<-paste("St_",fmdata$Station,"|",fmdata$Depth_m,"m","|", fmdata$Sample,sep="")

fmdata$UID<-factor(fmdata$UID,levels=unique(fmdata$UID[order(fmdata$Depth_m,fmdata$Station,decreasing=TRUE)]))


fig1A<-ggplot(fmdata,aes(x=Abundance,y=UID,fill=Taxon))+geom_bar(position="stack",stat="identity",colour="black",alpha=0.9)+theme_bw()+scale_fill_manual(name="Taxon",values=host_coloring)+facet_grid(Zone ~ ., scales="free_y")

ggsave("Malaspina_Profiles_mTAGs_phylum_proteobacteria_class_AbundancexSample.pdf",plot=fig1A,device=cairo_pdf,width=12,height=9,pointsize=8)


###mTAGs Phylum abundance
all_phylum_abund<-read.table(file="/mnt/lustre/bio/shared/malaspina/cnag-metagenomes/malaspina.mTags/output/2019-09-30.1_malaspina-160-metagenomes-cnag-SILVA132-0.99-merged-tables/Merged_table_phylum.txt",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)

colnames(all_phylum_abund)<-gsub("\\.(.)+$","",colnames(all_phylum_abund),perl=TRUE)
rownames(all_phylum_abund)<-all_phylum_abund$X
all_phylum_abund$X<-NULL

prok_phylum_abund<-all_phylum_abund[grepl("Bacteria|Archaea", rownames(all_phylum_abund)),]

prok_phylum_abund<-as.data.frame(t((t(prok_phylum_abund)/colSums(prok_phylum_abund))*100))

summary(prok_phylum_abund)

taxon<-data.frame(do.call('rbind', strsplit(as.character(rownames(prok_phylum_abund)),';',fixed=TRUE)))

rownames(prok_phylum_abund)<-taxon$X2

prok_phylum_abund<-t(prok_phylum_abund)

summary(prok_phylum_abund)

metadata<-read.table(file="/mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Metadata/Selected_Malaspina_Profiles_Metadata.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)
rownames(metadata)<-metadata$MPCode


mdata<-merge(prok_phylum_abund,subset(metadata,select=c("Station","Depth_m","layer")),by="row.names",all.x=TRUE)


mdata<-melt(mdata,id=c("Row.names","Station","Depth_m","layer"))

colnames(mdata)<-c("Sample","Station","Depth_m","Zone","Taxon","Abundance")

mdata$Sample<-as.character(mdata$Sample)

fmdata<-mdata[which(mdata$Abundance >= 1 & mdata$Zone != "NA"),]

summary(fmdata)

host_coloring<-c(brewer.pal(9,"Set1"),brewer.pal(4,"Pastel1"),brewer.pal(8,"Dark2")[8],brewer.pal(9,"Set1")[2],brewer.pal(8,"Accent"))
names(host_coloring)<-c("Actinobacteria","Alphaproteobacteria","Cyanobacteria","Euryarchaeota","Bacteroidetes","Deinococcus−Thermus","Verrucomicrobia","Firmicutes","Fusobacteria","Planctomycetes","Gammaproteobacteria","Betaproteobacteria","Crenarchaeota","Unknown","Proteobacteria","Thaumarchaeota","Chloroflexi","Gemmatimonadetes","Nitrospinae","WPS-2","Patescibacteria","Acidobacteria","Chlamydiae")

for (taxon in unique(fmdata$Taxon)) {
	if (taxon %in% names(host_coloring)) {
	} else {
		print(paste(taxon," missing from color palette!",sep=""))
	}
}

host_coloring<-host_coloring[c(as.vector(unique(fmdata$Taxon)))]

summary(fmdata)

fmdata$Zone<-gsub("Epi","Epipelagic",fmdata$Zone)
fmdata$Zone<-gsub("Meso","Mesopelagic",fmdata$Zone)
fmdata$Zone<-gsub("Bathy","Bathypelagic",fmdata$Zone)

fmdata$Zone<-factor(fmdata$Zone,levels=c("Epipelagic","Mesopelagic","Bathypelagic"))

fmdata$UID<-paste("St_",fmdata$Station,"|",fmdata$Depth_m,"m","|", fmdata$Sample,sep="")

fmdata$UID<-factor(fmdata$UID,levels=unique(fmdata$UID[order(fmdata$Depth_m,fmdata$Station,decreasing=TRUE)]))


fig1A<-ggplot(fmdata,aes(x=Abundance,y=UID,fill=Taxon))+geom_bar(position="stack",stat="identity",colour="black",alpha=0.9)+theme_bw()+scale_fill_manual(name="Taxon",values=host_coloring)+facet_grid(Zone ~ ., scales="free_y")

ggsave("Malaspina_Profiles_mTAGs_Phylum_AbundancexSample.pdf",plot=fig1A,device=cairo_pdf,width=10,height=12,pointsize=8)


###Virus NMDS
vir_scaff_abund<-read.table(file="/mnt/lustre/scratch/elopez/5_bowtie_results/After_checkV_output/RPKM.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)

fdata<-vir_scaff_abund
rownames(fdata)<-fdata[,1]
colnames(fdata)[which(colnames(fdata) == "MP2233")]<-"MP2233bis"
fdata<-fdata[,-1]

dist_metric<-"bray"

dists<-vegdist(t(fdata), method = dist_metric)
set.seed(666)
mdsresult<-metaMDS(dists,distance = dist_metric,k = 2,maxit = 999)
data.scores<-as.data.frame(scores(mdsresult))
data.scores$Sample<-rownames(data.scores)

mdata<-merge(data.scores,metadata,by="Sample",all.x=TRUE)

mdata$Zone<-"Epipelagic"
mdata$Zone[which(mdata$Depth_m > 200)] <- "Mesopelagic"
mdata$Zone[which(mdata$Depth_m > 1000)] <- "Bathypelagic"



summary(mdata)

zone_coloring<-c(brewer.pal(9,"YlOrRd")[c(4)],brewer.pal(9,"YlGnBu")[c(5,9)])
names(zone_coloring)<-c("Epipelagic","Mesopelagic","Bathypelagic")

fig_nmds<-ggplot(mdata, aes(x=NMDS1,y=NMDS2))+geom_point(shape=18,size=2,alpha=0.75,aes(color=Zone))+theme_bw()+theme(text=element_text(size=16))+scale_colour_manual(name="Zone",values=zone_coloring)

ggsave("Malaspina_Profiles_Viruses_NMDS.pdf",plot=fig_nmds,device=cairo_pdf,width=9,height=7,pointsize=8)

###Virus abundance
#Host Phylum
vir_host_phylum_abund<-read.table(file="/mnt/lustre/scratch/fcoutinho/elopez/Filtered_by_Predicted_Host_Score-Predicted_Host_2_Phylum-Original_IDs_RPKM.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)

mdata<-melt(vir_host_phylum_abund)
colnames(mdata)<-c("Taxon","Sample","Abundance")
summary(mdata)

mdata<-merge(mdata,metadata,by="Sample",all.x=TRUE)

mdata$Zone<-"Epipelagic"
mdata$Zone[which(mdata$Depth_m > 200)] <- "Mesopelagic"
mdata$Zone[which(mdata$Depth_m > 1000)] <- "Bathypelagic"

mdata$Zone<-factor(mdata$Zone,levels=c("Epipelagic","Mesopelagic","Bathypelagic"))

fmdata<-mdata[which(mdata$Abundance >= 1000),]

host_coloring<-c(brewer.pal(9,"Set1"),brewer.pal(4,"Pastel1"),brewer.pal(8,"Dark2")[8])
names(host_coloring)<-c("Actinobacteria","Alphaproteobacteria","Cyanobacteria","Euryarchaeota","Bacteroidetes","Deinococcus−Thermus","Verrucomicrobia","Firmicutes","Fusobacteria","Planctomycetes","Gammaproteobacteria","Betaproteobacteria","Crenarchaeota","Unknown")

host_coloring<-host_coloring[c(as.vector(unique(fmdata$Taxon)))]

fig1A<-ggplot(fmdata,aes(y=Abundance,x=Sample,fill=Taxon))+geom_bar(position="stack",stat="identity",colour="black",alpha=0.9)+coord_flip()+theme_bw()+theme(axis.text.x = element_text(angle = 45,hjust = 1,size=9))+scale_fill_manual(name="Taxon",values=host_coloring)+facet_grid(Zone ~ ., scales="free_y")

ggsave("Malaspina_Profiles_Virus_RaFAH_Filtered_0.14_Host_Phylum_AbundancexSample.pdf",plot=fig1A,device=cairo_pdf,width=10,height=12,pointsize=8)

#Family
vir_fam_abund<-read.table(file="/mnt/lustre/scratch/fcoutinho/elopez/family-Original_IDs_RPKM.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)

mdata<-melt(vir_fam_abund)
colnames(mdata)<-c("Taxon","Sample","Abundance")
summary(mdata)

mdata<-merge(mdata,metadata,by="Sample",all.x=TRUE)

mdata$Zone<-"Epipelagic"
mdata$Zone[which(mdata$Depth_m > 200)] <- "Mesopelagic"
mdata$Zone[which(mdata$Depth_m > 1000)] <- "Bathypelagic"
mdata$Zone<-factor(mdata$Zone,levels=c("Epipelagic","Mesopelagic","Bathypelagic"))

fmdata<-mdata[which(mdata$Abundance >= 1000),]


vir_fam_coloring<-c(brewer.pal(11,"Spectral"),brewer.pal(8,"Dark2")[8])
names(vir_fam_coloring)<-c("Lavidaviridae","Myoviridae","Microviridae","Inoviridae","Baculoviridae","Globuloviridae","Phycodnaviridae","Sphaerolipoviridae","Podoviridae","Mimiviridae","Siphoviridae","Unknown")

vir_fam_coloring<-vir_fam_coloring[sort(as.vector(unique(fmdata$Taxon)))]


fig1B<-ggplot(fmdata,aes(y=Abundance,x=Sample,fill=Taxon))+geom_bar(position="stack",stat="identity",colour="black",alpha=0.9)+coord_flip()+theme_bw()+theme(axis.text.x = element_text(angle = 45,hjust = 1,size=9))+facet_grid(Zone ~ ., scales="free_y")+scale_fill_manual(name="Taxon",values=vir_fam_coloring)

ggsave("Malaspina_Profiles_Virus_VPFClass_Family_AbundancexSample.pdf",plot=fig1B,device=cairo_pdf,width=10,height=12,pointsize=8)


#Genus
vir_genus_abund<-read.table(file="/mnt/lustre/scratch/fcoutinho/elopez/genus-Original_IDs_RPKM.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)

mdata<-melt(vir_genus_abund)
colnames(mdata)<-c("Taxon","Sample","Abundance")
summary(mdata)

mdata<-merge(mdata,metadata,by="Sample",all.x=TRUE)

mdata$Zone<-"Epipelagic"
mdata$Zone[which(mdata$Depth_m > 200)] <- "Mesopelagic"
mdata$Zone[which(mdata$Depth_m > 1000)] <- "Bathypelagic"
mdata$Zone<-factor(mdata$Zone,levels=c("Epipelagic","Mesopelagic","Bathypelagic"))

fmdata<-mdata[which(mdata$Abundance >= 100),]

#fig1C<-ggplot(fmdata,aes(y=Abundance,x=Sample,fill=Taxon))+geom_bar(position="stack",stat="identity",colour="black",alpha=0.9)+coord_flip()+theme_bw()+theme(axis.text.x = element_text(angle = 45,hjust = 1,size=9))+facet_grid(Zone ~ ., scales="free_y")

col_pal<-brewer.pal(11,"Spectral")
depth_col_grad<-rev(colorRampPalette(col_pal)(n=299))

fig1C<-ggplot(fmdata,aes(fill=log10(Abundance),y=Sample,x=Taxon))+geom_raster()+theme_bw()+theme(axis.text.x = element_text(angle = 45,hjust = 1,size=9))+facet_grid(Zone ~ ., scales="free")+scale_fill_gradientn(colours =depth_col_grad)

ggsave("Malaspina_Profiles_Virus_VPFClass_Genus_AbundancexSample.pdf",plot=fig1C,device=cairo_pdf,width=25,height=15,pointsize=8)

###NMDS
#Scaffold


###Virus genomes 
vir_data<-read.table(file="/mnt/lustre/scratch/elopez/malaspina_virus_merged_3.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)

fvir_data<-vir_data[which(vir_data$Completeness >= 50 & vir_data$Contamination <= 5),]

vir_data<-read.table(file="/mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Viruses/PHIST_Output/Positive_PHIST_Predictions+Host_Taxonomy.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)

vir_data$Sample<-gsub("_[0-9]+.fasta","",vir_data$phage,perl=TRUE)

for (i in 1:nrow(vir_data)) {
	sample_id<-vir_data$Sample[i]
	vir_data$Sample_Depth[i]<-metadata[sample_id[1],"Depth_m"]
}

summary(vir_data)


vir_data$Count<-1
tax_sample_counts<-aggregate(vir_data$Count, by=list(Phylum = vir_data$Host_Phylum, Sample =vir_data$Sample, Depth=vir_data$Sample_Depth), FUN=sum)
colnames(tax_sample_counts)[4]<-"Count"
summary(tax_sample_counts)

tax_sample_counts$Zone<-"Epipelagic"
tax_sample_counts$Zone[which(tax_sample_counts$Depth > 200)] <- "Mesopelagic"
tax_sample_counts$Zone[which(tax_sample_counts$Depth > 1000)] <- "Bathypelagic"

tax_sample_counts$Zone<-factor(tax_sample_counts$Zone,levels=c("Epipelagic","Mesopelagic","Bathypelagic"))

f_tax_sample_counts<-tax_sample_counts[which(tax_sample_counts$Count >= 10),]

fig1A<-ggplot(f_tax_sample_counts,aes(y=Count,x=Sample,fill=Phylum))+geom_bar(position="stack",stat="identity",colour="black",alpha=0.9)+coord_flip()+theme_bw()+theme(axis.text.x = element_text(angle = 45,hjust = 1,size=9))+facet_grid(Zone ~ ., scales="free_y")

ggsave("Malaspina_Profiles_Virus_PHIST_Host_Phylum_CountsxSample.pdf",plot=fig1A,device=cairo_pdf,width=13,height=9,pointsize=8)




###MAG Info Virus Prevalence  x Depth Plots 
bin_data<-read.table(file="/mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Viruses/Profiles_Malaspina_Metabat_Bins_Round_1_Info+Virus_Info.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)

rownames(bin_data)<-bin_data[,1]

mag_data<-bin_data[which(bin_data$Completeness >= 50 & bin_data$Contamination <= 5),]

mag_data$Sample_Depth<-NA

for (i in 1:nrow(mag_data)) {
	sample_id<-mag_data$Sample[i]
	mag_data$Sample_Depth[i]<-metadata[sample_id[1],"Depth_m"]
}

fig1A<-ggplot(mag_data,aes(y=Prophage_Count,x=Sample_Depth))+geom_point(aes(colour=Phylum))+geom_smooth(method = "loess",se=TRUE,  formula = y ~ x)+theme_bw()

ggsave("Malaspina_Profiles_MAG_Virus_CountxDepth.pdf",plot=fig1A,device=cairo_pdf,width=12,height=7,pointsize=8)

fig2A<-ggplot(mag_data,aes(y=Prophage_Count,x=Sample_Depth))+geom_point(aes(colour=Phylum))+geom_smooth(method = "loess",se=TRUE,  formula = y ~ x)+theme_bw()+facet_wrap(Phylum ~ .,scales="free")

ggsave("Malaspina_Profiles_MAG_Virus_CountxDepth_by_Phylum.pdf",plot=fig2A,device=cairo_pdf,width=20,height=15,pointsize=8)

fig3A<-ggplot(mag_data,aes(y=Prophage_Count,x=Phylum))+geom_boxplot(aes(fill=Phylum))+theme(axis.text.x = element_text(angle = 45,hjust = 1,size=9))

ggsave("Malaspina_Profiles_MAG_Virus_CountxPhylum_Boxplot.pdf",plot=fig3A,device=cairo_pdf,width=15,height=10,pointsize=8)

####Tax sample counts
bin_data<-read.table(file="/mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Profiles_Malaspina_Metabat_Bins_Round_1_Info.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)

rownames(bin_data)<-bin_data[,1]

bin_data$Quality<-(bin_data$Completeness - (3*(bin_data$Contamination)))

mag_data<-bin_data[which(bin_data$Completeness >= 50 & bin_data$Contamination <= 5),]

mag_data$Estimated_Genome_Size<-(mag_data$Bases/mag_data$Completeness)*100

mag_data$Sample_Depth<-NA

for (i in 1:nrow(mag_data)) {
	sample_id<-mag_data$Sample[i]
	mag_data$Sample_Depth[i]<-metadata[sample_id[1],"Depth_m"]
}

summary(mag_data)


mag_data$Count<-1
tax_sample_counts<-aggregate(mag_data$Count, by=list(Phylum = mag_data$Phylum, Sample =mag_data$Sample, Depth=mag_data$Sample_Depth), FUN=sum)
colnames(tax_sample_counts)[4]<-"Count"
summary(tax_sample_counts)


tax_sample_counts$Zone<-"Epipelagic"
tax_sample_counts$Zone[which(tax_sample_counts$Depth > 200)]<-"Mesopelagic"
tax_sample_counts$Zone[which(tax_sample_counts$Depth > 1000)]<-"Bathypelagic"

tax_sample_counts$Zone<-factor(tax_sample_counts$Zone,levels=c("Epipelagic","Mesopelagic","Bathypelagic"))

f_tax_sample_counts<-tax_sample_counts[which(tax_sample_counts$Count >= 3),]
fig1A<-ggplot(f_tax_sample_counts,aes(y=Count,x=Sample,fill=Phylum))+geom_bar(position="stack",stat="identity",colour="black",alpha=0.9)+coord_flip()+theme_bw()+theme(axis.text.x = element_text(angle = 45,hjust = 1,size=9))+facet_grid(Zone ~ ., scales="free_y")

ggsave("Malaspina_Profiles_MAG_Phylum_CountsxSample.pdf",plot=fig1A,device=cairo_pdf,width=13,height=7,pointsize=8)



###Correlogram
library("GGally")

cgram_metadata<-subset(metadata,select=c(Depth_m,Temp,Sal_PSU,NO3_WOA13,PO4_WOA13,SiO4_WOA13,Oxygen,All_BT,all_virus))

pdf("Malaspina_Profiles_Metadata_Correlogram.pdf",width=25,height=25,pointsize=8)
plot<-ggpairs(cgram_metadata)
print(plot)
dev.off()


###Virus Host Filtered (0.14) abund barplot

predhost_abund<-read.table(file="/mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Viruses/Corrected_Abundance/Predicted_Host_2_Phylum-Profiles_Malaspina_Virus_RPKM.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)

fdata<-predhost_abund
colnames(fdata)[which(colnames(fdata) == "MP2233")]<-"MP2233bis"
rownames(fdata)<-fdata[,1]

fdata<-merge(t(fdata[,-1]),subset(metadata,select=c("Depth_m","layer","Station")),by="row.names",all.x=TRUE)

mdata<-melt(fdata,id=c("Row.names","Depth_m","layer","Station"))
colnames(mdata)<-c("Sample","Depth","Layer","Station","Taxon","Abundance")
summary(mdata)

#mdata<-na.omit(mdata)
#mdata$Abundance<-as.numeric(mdata$Abundance)

min_abund<-1000
mdata<-mdata[which(mdata$Abundance >= min_abund),]
mdata<-mdata[which(mdata$Taxon != "Unknown"),]

summary(mdata)

mdata$Sample_ID<-paste("St",mdata$Station,"_",mdata$Depth,"m",sep="")

summary(mdata)

mdata$Layer<-factor(mdata$Layer,levels=c("Epi","Meso","Bathy"))

mdata$Sample_ID<-factor(mdata$Sample_ID,levels=unique(mdata$Sample_ID[rev(order(mdata$Depth))]))

host_coloring<-c(brewer.pal(9,"Set1"),brewer.pal(4,"Pastel1"))
names(host_coloring)<-c("Actinobacteria","Alphaproteobacteria","Cyanobacteria","Euryarchaeota","Bacteroidetes","Deinococcus−Thermus","Verrucomicrobia","Firmicutes","Fusobacteria","Planctomycetes","Gammaproteobacteria","Betaproteobacteria","Crenarchaeota")

host_coloring<-host_coloring[c(as.vector(unique(mdata$Taxon)))]

fig_2B<-ggplot(mdata, aes(x=Sample_ID, y=Abundance, fill=Taxon))+geom_bar(position="stack",stat="identity",colour="black",alpha=0.9)+coord_flip()+theme_bw()+theme(legend.position="right",text = element_text(size = 18),axis.text.x = element_text(angle = 45,hjust = 1,size=9),legend.key.size = unit(0.5, "cm"))+scale_fill_manual(name="Taxon",values=host_coloring)+facet_grid(Layer ~ . , scales="free_y",space="free")

ggsave("Malaspina_Profiles_Viruses_Host_Abund.pdf",plot=fig_2B,device=cairo_pdf,width=10,height=15,pointsize=8)


###DBSCAN
kmer_data<-read.table(file="/mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Kmer_Counts/Kmer_Counts.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)
rownames(kmer_data)<-kmer_data[,1]
kmer_data<-kmer_data[,-1]

perc_kmer_data<-t(t(kmer_data) / colSums(kmer_data))

dist_metric<-"bray"#"euclidean"#"manhattan"#

dists<-vegdist(t(perc_kmer_data), method = dist_metric)

set.seed(999)
eps_val<-0.02
db <- fpc::dbscan(dists, eps = eps_val, MinPts = 2, method="dist")

all_samples_dbclusters<-db$cluster
names(all_samples_dbclusters)<-colnames(perc_kmer_data)

metadata<-read.table(file="/mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Metadata/Malaspina_Profiles_Sample_Info.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)

rownames(metadata)<-metadata[,1]
metadata<-metadata[,-1]

metadata$Station<-as.character(metadata$Station)

metadata$ZoneB<-"Epipelagic"
metadata$ZoneB[which(metadata$Depth > 200)]<-"Mesopelagic"
metadata$ZoneB[which(metadata$Depth > 1000)]<-"Bathypelagic"

set.seed(999)
mdsresult<-metaMDS(dists,distance = dist_metric,k = 2,maxit = 999)
data.scores<-as.data.frame(scores(mdsresult))
data.scores$Sample<-rownames(data.scores)
mdata<-data.scores
stress<-round(mdsresult$stress,digits=4)

for (i in 1:ncol(metadata)) {
	var_name<-colnames(metadata)[i]
	mdata[[var_name]]<-NA
	for (j in 1:nrow(mdata)) {
		sample_name<-rownames(mdata)[j]
		mdata[sample_name,var_name]<-as.vector(metadata[sample_name,var_name])
	}
}

mdata$DBSCAN_Cluster<-as.character(db$cluster)

title_string<-paste("EPS:",eps_val,"Stress:",stress,dist_metric,sep=" ")
pdf("Malaspina_Profiles_Kmer_NMDS_by_ZoneB_DBSCAN_CLusters.pdf",width=30,height=28,pointsize=8)
plot<-ggplot(mdata, aes(x=NMDS1,y=NMDS2, label=Sample))+geom_point(size=6,alpha=0.75,aes(color=DBSCAN_Cluster,shape=ZoneB))+geom_text()+theme_bw()+theme(text=element_text(size=16))
print(plot)
dev.off()

out_table_file<-paste("Malaspina_Profiles_DBSCAN",eps_val,"Kmer_Results.tsv",sep="_")
write.table(mdata,file=out_table_file,sep="\t",append=FALSE,row.names=FALSE,col.names=TRUE,quote=FALSE)

###Kmer Counts NMDS
library(ggplot2)
library(vegan)
library(RColorBrewer)

kmer_data<-read.table(file="/mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Kmer_Counts/Kmer_Counts.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)

rownames(kmer_data)<-kmer_data[,1]
kmer_data<-kmer_data[,-1]

metadata<-read.table(file="/mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Metadata/Malaspina_Profiles_Sample_Info.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)

rownames(metadata)<-metadata[,1]
metadata<-metadata[,-1]

metadata$Station<-as.character(metadata$Station)

metadata$ZoneB<-"Epipelagic"
metadata$ZoneB[which(metadata$Depth > 200)]<-"Mesopelagic"
metadata$ZoneB[which(metadata$Depth > 1000)]<-"Bathypelagic"


metadata2<-read.table(file="/mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Metadata/04_malaspina-160-metag-metadata.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)

rownames(metadata2)<-metadata2[,1]
metadata2<-metadata2[,-1]

metadata2<-subset(metadata2,select=c(Temp,Conductivity,Fluo,PAR,SPAR,Turb_FTU,Sal_PSU,Salinity_WOA13,NO3_WOA13,PO4_WOA13,SiO4_WOA13,MLD,Oxygen,sigma,O2_umol_kg,O2_corr_umol_kg,O2_sat,AOU_corr_umol_kg,Chla_ugl,Fmax1_resp_prok,Fmax2_resp_euk,Fmax3_tirosina,Fmax4_triptofano,TEP,POC_uM,Turb ,pmol_leu,SE,LNA,HNA,All_BT,percentHNA,cell_size,Bacterial_cell_C,Biomass,ugC_l_d,d_1,turnover_days,HNF,low_virus,medium_virus,high_virus,all_virus,VBR))

perc_kmer_data<-t(t(kmer_data) / colSums(kmer_data))

dist_metric<-"bray"#"euclidean"#"manhattan"#

dists<-vegdist(t(perc_kmer_data), method = dist_metric)

mdsresult<-metaMDS(dists,distance = dist_metric,k = 2,maxit = 999)
data.scores<-as.data.frame(scores(mdsresult))
data.scores$Sample<-rownames(data.scores)
mdata<-data.scores
stress<-round(mdsresult$stress,digits=4)

for (i in 1:ncol(metadata)) {
	var_name<-colnames(metadata)[i]
	mdata[[var_name]]<-NA
	for (j in 1:nrow(mdata)) {
		sample_name<-rownames(mdata)[j]
		mdata[sample_name,var_name]<-as.vector(metadata[sample_name,var_name])
	}
}

for (i in 1:ncol(metadata2)) {
	var_name<-colnames(metadata2)[i]
	mdata[[var_name]]<-NA
	for (j in 1:nrow(mdata)) {
		sample_name<-rownames(mdata)[j]
		mdata[sample_name,var_name]<-as.vector(metadata2[sample_name,var_name])
	}
}

sdata<-subset(mdata,select=c(NMDS1,NMDS2,Temp,Conductivity,Fluo,PAR,SPAR,Turb_FTU,Sal_PSU,Salinity_WOA13,NO3_WOA13,PO4_WOA13,SiO4_WOA13,Oxygen,sigma,O2_umol_kg,O2_corr_umol_kg,O2_sat,AOU_corr_umol_kg,Fmax1_resp_prok,Fmax2_resp_euk,Fmax3_tirosina,Fmax4_triptofano,TEP,Turb ,pmol_leu,SE,LNA,HNA,All_BT,percentHNA,cell_size,Bacterial_cell_C,Biomass,ugC_l_d,d_1,turnover_days,HNF,low_virus,medium_virus,high_virus,all_virus,VBR))

pdf("Malaspina_Profiles_NMDS_Metadata_Correlogram_All_Samples.pdf",width=25,height=25,pointsize=8)
plot<-ggpairs(sdata)
print(plot)
dev.off()


title_string<-paste("Stress:",stress,dist_metric,sep=" ")
pdf("Malaspina_Profiles_Kmer_NMDS_by_Zone.pdf",width=7,height=5,pointsize=8)
plot<-ggplot(mdata, aes(x=NMDS1,y=NMDS2))+geom_point(size=2,alpha=0.75,aes(color=Zone))+ggtitle(title_string)+theme_bw()+theme(text=element_text(size=16))
print(plot)
dev.off()


title_string<-paste("Stress:",stress,dist_metric,sep=" ")
pdf("Malaspina_Profiles_Kmer_NMDS_by_ZoneB.pdf",width=7,height=5,pointsize=8)
plot<-ggplot(mdata, aes(x=NMDS1,y=NMDS2))+geom_point(size=2,alpha=0.75,aes(color=ZoneB))+ggtitle(title_string)+theme_bw()+theme(text=element_text(size=16))
print(plot)
dev.off()

col_pal<-brewer.pal(9,"BuPu")
depth_col_grad<-colorRampPalette(col_pal[4:9])(n=299)

depth_correl<-cor.test(mdata$NMDS1,mdata$Depth)
title_string<-paste("Stress:",stress,"PCC NMDS1 x Depth:",round(depth_correl$estimate,digits=2),dist_metric,sep=" ")
pdf("Malaspina_Profiles_Kmer_NMDS_by_Depth.pdf",width=7,height=5,pointsize=8)
plot<-ggplot(mdata, aes(x=NMDS1,y=NMDS2))+geom_point(size=2,alpha=0.75,aes(color=Depth))+scale_colour_gradientn(colours =depth_col_grad)+ggtitle(title_string)+theme_bw()+theme(text=element_text(size=16))
print(plot)
dev.off()

##NMDS splitting by Zone
library("GGally")

kmer_data<-read.table(file="/mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Kmer_Counts/Kmer_Counts.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)

rownames(kmer_data)<-kmer_data[,1]
kmer_data<-kmer_data[,-1]

metadata<-read.table(file="/mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Metadata/Malaspina_Profiles_Sample_Info.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)

rownames(metadata)<-metadata[,1]
metadata<-metadata[,-1]

metadata$Station<-as.character(metadata$Station)

metadata$ZoneB<-"Epipelagic"
metadata$ZoneB[which(metadata$Depth > 200)]<-"Mesopelagic"
metadata$ZoneB[which(metadata$Depth > 1000)]<-"Bathypelagic"


metadata2<-read.table(file="/mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Metadata/04_malaspina-160-metag-metadata.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)

rownames(metadata2)<-metadata2[,1]
metadata2<-metadata2[,-1]

metadata2<-subset(metadata2,select=c(Temp,Conductivity,Fluo,PAR,SPAR,Turb_FTU,Sal_PSU,Salinity_WOA13,NO3_WOA13,PO4_WOA13,SiO4_WOA13,MLD,Oxygen,sigma,O2_umol_kg,O2_corr_umol_kg,O2_sat,AOU_corr_umol_kg,Chla_ugl,Fmax1_resp_prok,Fmax2_resp_euk,Fmax3_tirosina,Fmax4_triptofano,TEP,POC_uM,Turb ,pmol_leu,SE,LNA,HNA,All_BT,percentHNA,cell_size,Bacterial_cell_C,Biomass,ugC_l_d,d_1,turnover_days,HNF,low_virus,medium_virus,high_virus,all_virus,VBR))

perc_kmer_data<-t(t(kmer_data) / colSums(kmer_data))

dist_metric<-"bray"#"euclidean"#"manhattan"#

dists<-vegdist(t(perc_kmer_data), method = dist_metric)


all_zones<-unique(metadata$ZoneB)
col_pal<-brewer.pal(11,"Spectral")
depth_col_grad<-rev(colorRampPalette(col_pal)(n=299))

make_correlograms<-FALSE

for (zone in all_zones) {
	sample_list<-rownames(metadata)[which(metadata$ZoneB == zone)]
	subdata<-subset(perc_kmer_data,select=sample_list)
	
	dists<-vegdist(t(subdata), method = dist_metric)
	mdsresult<-metaMDS(dists,distance = dist_metric,k = 2,maxit = 999)
	data.scores<-as.data.frame(scores(mdsresult))
	data.scores$Sample<-rownames(data.scores)
	mdata<-data.scores
	stress<-round(mdsresult$stress,digits=4)
		
	for (i in 1:ncol(metadata)) {
		var_name<-colnames(metadata)[i]
		mdata[[var_name]]<-NA
		for (j in 1:nrow(mdata)) {
			sample_name<-rownames(mdata)[j]
			mdata[sample_name,var_name]<-as.vector(metadata[sample_name,var_name])
		}
	}

	for (i in 1:ncol(metadata2)) {
		var_name<-colnames(metadata2)[i]
		mdata[[var_name]]<-NA
		for (j in 1:nrow(mdata)) {
			sample_name<-rownames(mdata)[j]
			mdata[sample_name,var_name]<-as.vector(metadata2[sample_name,var_name])
		}
	}


	mdata$DBSCAN_Cluster<-NA
	for (j in 1:nrow(mdata)) {
		sample_name<-rownames(mdata)[j]
		mdata[sample_name,"DBSCAN_Cluster"]<-as.character(all_samples_dbclusters[sample_name])
	}
	
	fig_name<-paste("Malaspina_Profiles_Kmer_NMDS_",zone,".pdf",sep="")
	title_string<-paste("Stress:",stress,dist_metric,sep=" ")
	pdf(fig_name,width=7,height=5,pointsize=8)
	plot<-ggplot(mdata, aes(x=NMDS1,y=NMDS2))+geom_point(size=2,alpha=0.75,aes(color=pmol_leu))+scale_colour_gradientn(colours =depth_col_grad)+ggtitle(title_string)+theme_bw()+theme(text=element_text(size=16))
	print(plot)
	dev.off()
	
	fig_name<-paste("Malaspina_Profiles_Kmer_NMDS_",zone,"_DBSCAN_Clusters.pdf",sep="")
	title_string<-paste("Stress:",stress,dist_metric,sep=" ")
	pdf(fig_name,width=9,height=7,pointsize=8)
	plot<-ggplot(mdata, aes(x=NMDS1,y=NMDS2, label=Station))+geom_text()+geom_point(size=2,alpha=0.75,aes(color=DBSCAN_Cluster))+ggtitle(title_string)+theme_bw()+theme(text=element_text(size=16))
	print(plot)
	dev.off()
	
	if (make_correlograms == TRUE) {
		sdata<-subset(mdata,select=c(NMDS1,NMDS2,Depth,Temp,Conductivity,Fluo,PAR,SPAR,Turb_FTU,Sal_PSU,Salinity_WOA13,NO3_WOA13,PO4_WOA13,SiO4_WOA13,Oxygen,sigma,O2_umol_kg,O2_corr_umol_kg,O2_sat,AOU_corr_umol_kg,Fmax1_resp_prok,Fmax2_resp_euk,Fmax3_tirosina,Fmax4_triptofano,TEP,Turb,pmol_leu,SE,LNA,HNA,All_BT,percentHNA,cell_size,Bacterial_cell_C,Biomass,ugC_l_d,d_1,turnover_days,low_virus,medium_virus,high_virus,all_virus,VBR))
		
		fig_name<-paste("Malaspina_Profiles_NMDS_Metadata_Correlogram_",zone,".pdf",sep="")
		pdf(fig_name,width=30,height=30,pointsize=8)
		plot<-ggpairs(sdata)
		print(plot)
		dev.off()
	}

}

###Assembly Info
library(ggpubr, lib.loc="/mnt/lustre/bio/users/fcoutinho/Rlibs/")
library(ggplot2)
library(RColorBrewer)
library(gridExtra)

spades_data<-read.table(file="/mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Profiles_Malaspina_Filtered_Assembly_Stats.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)
spades_data$Assembler<-"SPAdes"

megahit_data<-read.table(file="/mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Profiles_Malaspina_MegaHIT_Full_Assembly_Stats.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)
megahit_data$Assembler<-"MEGAHIT"

merged_data<-rbind(spades_data,megahit_data)

merged_data$Average_Contig_Length<-merged_data$Bases/merged_data$Contigs

asb_coloring<-c(brewer.pal(3,"Paired")[2],brewer.pal(3,"Dark2")[2])
names(asb_coloring)<-unique(merged_data$Assembler)

contigs_plot<-ggplot(merged_data,aes(y=Contigs,x=Assembler))+geom_boxplot(aes(fill=Assembler))+geom_signif(comparisons = list(c("SPAdes","MEGAHIT")),map_signif_level = FALSE, textsize=3,test = "wilcox.test")+scale_fill_manual(name="Assembler",values=asb_coloring)

bases_plot<-ggplot(merged_data,aes(y=Bases,x=Assembler))+geom_boxplot(aes(fill=Assembler))+geom_signif(comparisons = list(c("SPAdes","MEGAHIT")),map_signif_level = FALSE, textsize=3,test = "wilcox.test")+scale_fill_manual(name="Assembler",values=asb_coloring)

n50_plot<-ggplot(merged_data,aes(y=N50,x=Assembler))+geom_boxplot(aes(fill=Assembler))+geom_signif(comparisons = list(c("SPAdes","MEGAHIT")),map_signif_level = FALSE, textsize=3,test = "wilcox.test")+scale_fill_manual(name="Assembler",values=asb_coloring)

avg_len_plot<-ggplot(merged_data,aes(y=Average_Contig_Length,x=Assembler))+geom_boxplot(aes(fill=Assembler))+geom_signif(comparisons = list(c("SPAdes","MEGAHIT")),map_signif_level = FALSE, textsize=3,test = "wilcox.test")+scale_fill_manual(name="Assembler",values=asb_coloring)

max_len_plot<-ggplot(merged_data,aes(y=Max,x=Assembler))+geom_boxplot(aes(fill=Assembler))+geom_signif(comparisons = list(c("SPAdes","MEGAHIT")),map_signif_level = FALSE, textsize=3,test = "wilcox.test")+scale_fill_manual(name="Assembler",values=asb_coloring)

ggarrange(contigs_plot,bases_plot,n50_plot,avg_len_plot,max_len_plot,common.legend = TRUE, legend="right",labels="AUTO")
ggsave("Malaspina_Profiles_Assembly_Comparison.pdf",width=14,height=12,pointsize=8)

###Bin Info
library(ggpubr, lib.loc="/mnt/lustre/bio/users/fcoutinho/Rlibs/")
library(ggplot2)
library(gridExtra)
library(reshape2)
library(ggsignif)
library(RColorBrewer)



mtb_data<-read.table(file="/mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/CheckM_Info_Metabat_Bins_Round_1.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)
mtb_data$Method<-"MetaBat2"
colnames(mtb_data)[1]<-"Bin.Id"
cat_data<-read.table(file="/mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/CheckM_Info_Metabat_CAT_Bins_Round_1.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)
cat_data$Method<-"MetaBat2+CAT"

mxb_data<-read.table(file="/mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/MaxBin_Binning/Profiles_Malaspina_CheckM_MaxBin_Bins_Round_1.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)
mxb_data$Method<-"MaxBin"

bin_data<-rbind(mtb_data,cat_data,mxb_data)

bin_data$Quality<-(bin_data$Completeness - (3*(bin_data$Contamination)))


fig_comp_bin<-ggplot(bin_data,aes(y=Completeness,x=Method))+geom_boxplot(aes(fill=Method))#+geom_signif(comparisons = list(c("MetaBat2","MetaBat2+CAT")),map_signif_level = FALSE, textsize=3,test = "wilcox.test")

fig_conta_bin<-ggplot(bin_data,aes(y=Contamination,x=Method))+geom_boxplot(aes(fill=Method))#+geom_signif(comparisons = list(c("MetaBat2","MetaBat2+CAT")),map_signif_level = FALSE, textsize=3,test = "wilcox.test")

fig_qual_bin<-ggplot(bin_data,aes(y=Quality,x=Method))+geom_boxplot(aes(fill=Method))#+geom_signif(comparisons = list(c("MetaBat2","MetaBat2+CAT")),map_signif_level = FALSE, textsize=3,test = "wilcox.test")


fig_bin<-ggarrange(fig_comp_bin, fig_conta_bin, fig_qual_bin, nrow=1, labels = "AUTO") 

ggsave("Profiles_Malaspina_CheckM_Bins.pdf",plot=fig_bin,device=cairo_pdf,width=15,height=7,pointsize=8)


mag_data<-bin_data[which(bin_data$Completeness >= 50 & bin_data$Contamination <= 5),]
summary(mag_data)

fig_comp_mag<-ggplot(mag_data,aes(y=Completeness,x=Method))+geom_boxplot(aes(fill=Method))#+geom_signif(comparisons = list(c("MetaBat2","MetaBat2+CAT")),map_signif_level = FALSE, textsize=3,test = "wilcox.test")

fig_conta_mag<-ggplot(mag_data,aes(y=Contamination,x=Method))+geom_boxplot(aes(fill=Method))#+geom_signif(comparisons = list(c("MetaBat2","MetaBat2+CAT")),map_signif_level = FALSE, textsize=3,test = "wilcox.test")

fig_qual_mag<-ggplot(mag_data,aes(y=Quality,x=Method))+geom_boxplot(aes(fill=Method))#+geom_signif(comparisons = list(c("MetaBat2","MetaBat2+CAT")),map_signif_level = FALSE, textsize=3,test = "wilcox.test")


fig_bin<-ggarrange(fig_comp_mag, fig_conta_mag, fig_qual_mag, nrow=1, labels = "AUTO") 

ggsave("Profiles_Malaspina_CheckM_MAGs.pdf",plot=fig_bin,device=cairo_pdf,width=15,height=7,pointsize=8)

#MAG novelty plots

bin_data<-read.table(file="/mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/All_MAGs/All_MAGs_Full_Info.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)

mag_data<-bin_data[which(bin_data$Completeness >= 50 & bin_data$Contamination <= 5),]

mag_data$Source<-"Original"
mag_data$Source[grepl("Refined",mag_data$user_genome)]<-"Refined"
mag_data$Source<-as.factor(mag_data$Source)


#MAG ANI Comps TOPC Malaspina Deep
 ani_data<-read.table(file="/mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/All_MAGs/ANI_Comps/fastANI_Result_Profiles_Malaspina_dRep_MAGsxTOPC_Deep_Malaspina_MAGs.tsv",sep="\t",header=FALSE,quote="",comment="",stringsAsFactors=TRUE)
 colnames(ani_data)<-c("Query_Genome","Reference_Genome","ANI","BFV","TQF")

ani_data$Reference_Source<-NA

ani_data$Reference_Source[grepl("TARA_Polar",ani_data$Reference_Genome)]<-"TOPC"
ani_data$Reference_Source[grepl("bathypelagic",ani_data$Reference_Genome)]<-"Malaspina_Deep"

ani_data$Query_Genome<-gsub(".fasta$","",ani_data$Query_Genome,perl=TRUE)
ani_data$Query_Genome<-gsub("(.)*\\/","",ani_data$Query_Genome,perl=TRUE)

bin_data<-read.table(file="/mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/All_MAGs/All_MAGs_Full_Info.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)

mag_data<-bin_data[which(bin_data$Completeness >= 50 & bin_data$Contamination <= 5),]
colnames(mag_data)[1]<-c("Query_Genome")

ani_data<-merge(ani_data,mag_data,by="Query_Genome",all.x=TRUE)

fig3A<-ggplot(ani_data,aes(y=ANI,x=Reference_Source))+geom_boxplot(aes(fill=Phylum))+theme_bw()

ggsave("Malaspina_Profiles_MAG_Virus_CountxPhylum_Boxplot.pdf",plot=fig3A,device=cairo_pdf,width=15,height=10,pointsize=8)


 