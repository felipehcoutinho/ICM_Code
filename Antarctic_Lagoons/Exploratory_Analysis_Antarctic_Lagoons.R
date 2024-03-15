# install.packages("indispecies",lib="/mnt/lustre/bio/users/fcoutinho/Rlibs/",dep=TRUE) # library not available for current R verison in marbits
library(vegan)
library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)

source("/mnt/smart/users/fcoutinho/ICM_Code/Microbiome_Analysis.R")

#"/mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/Antarctic_Lagoon/Metadata/Sample_Metadata.tsv"
#"/mnt/smart/users/fcoutinho/Antarctic_Lagoons/Metadata/Sample_Metadata.tsv"
metadata_file<-"/mnt/smart/users/fcoutinho/Antarctic_Lagoons/Metadata/Sample_Metadata.tsv"
metadata_df<-read.table(file = metadata_file, sep = "\t", header=TRUE)
rownames(metadata_df)<-metadata_df$Sample_UID

#metadata_df$Depth <- metadata_df$Depth - 0.5

# metadata_df$MG_ID<-rownames(metadata_df)
# metadata_df$MG_ID<-gsub("Core(\\s)+","",metadata_df$MG_ID,perl=TRUE)
# metadata_df$MG_ID<-gsub("Cm(\\s)","",metadata_df$MG_ID,perl=TRUE)
# metadata_df$MG_ID<-as.factor(paste(metadata_df$MG_ID,"_S1",sep=""))

# core_id<-data.frame(do.call('rbind', strsplit(as.character(metadata_df$MG_ID),'_',fixed=TRUE)))
# colnames(core_id)<-c("Core","Depth_cm","Replicate")
# core_id$Depth_cm<-as.numeric(as.character(core_id$Depth_cm))
# core_id$MG_ID<-metadata_df$MG_ID
# summary(core_id)
# metadata_df<-merge(metadata_df,core_id,by="MG_ID",all.x=TRUE)
# metadata_df$Depth<-NULL
# metadata_df$Core_ID<-NULL
# summary(metadata_df)
#write.table(metadata_df,file="/mnt/smart/users/fcoutinho/Antarctic_Lagoons/Metadata/Sample_Metadata.tsv",sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)

depth_col_pal<-brewer.pal(9,"BrBG")#[c(1,4,9)]
depth_col_grad<-colorRampPalette(depth_col_pal)(n=99)

year_col_pal<-brewer.pal(11,"Spectral")
year_col_grad<-rev(colorRampPalette(year_col_pal)(n=299))

phylum_coloring<-c(brewer.pal(9,"Set1"),brewer.pal(4,"Pastel1"),brewer.pal(8,"Dark2")[8],brewer.pal(5,"Set2")[5],brewer.pal(8,"Accent"),brewer.pal(5,"Dark2"))
names(phylum_coloring)<-c("Actinobacteriota","Alphaproteobacteria","Synergistota","Crenarchaeota","Bacteroidota","Nanoarchaeota","Verrucomicrobiota","Marinisomatota","Others","Bdellovibrionota","Gammaproteobacteria","Acidobacteriota","UBP7","Cloacimonadota","Firmicutes","Margulisbacteria","Gemmatimonadota","Desulfobacterota","Nitrospinota","Desulfobacterota_D","Patescibacteria","Hydrothermarchaeota","Chlamydiota","SAR324","Halobacterota","Chloroflexi","Thermoplasmatota","Myxococcota")

#ASV data
#/mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/Antarctic_Lagoon/DADA/DADA2_ASVs_Percentage_Abundances_Antarctic_Lagoons.tsv
#/mnt/smart/users/fcoutinho/Antarctic_Lagoons/DADA_Output/DADA2_ASVs_Percentage_Abundances_Antarctic_Lagoons.tsv
asv_abd_file<-"/mnt/smart/users/fcoutinho/Antarctic_Lagoons/DADA_Output/DADA2_ASVs_Percentage_Abundances_Antarctic_Lagoons.tsv"
asv_perc_abd_df<-read.table(file = asv_abd_file, sep = "\t", header=TRUE)
rownames(asv_perc_abd_df)<-asv_perc_abd_df[,1]
colnames(asv_perc_abd_df)[1]<-"MG_ID"
asv_raw_abd_file<-"/mnt/smart/users/fcoutinho/Antarctic_Lagoons/DADA_Output/DADA2_ASVs_Raw_Abundances_Antarctic_Lagoons.tsv"
asv_raw_abd_df<-read.table(file = asv_raw_abd_file, sep = "\t", header=TRUE)
rownames(asv_raw_abd_df)<-asv_raw_abd_df[,1]
colnames(asv_raw_abd_df)[1]<-"MG_ID"

#/mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/Antarctic_Lagoon/DADA/DADA2_ASVs_Info_Antarctic_Lagoons.tsv
#/mnt/smart/users/fcoutinho/Antarctic_Lagoons/DADA_Output/DADA2_ASVs_Info_Antarctic_Lagoons.tsv
asv_info_file<-"/mnt/smart/users/fcoutinho/Antarctic_Lagoons/DADA_Output/DADA2_ASVs_Info_Antarctic_Lagoons.tsv"
asv_info_df<-read.table(file = asv_info_file, sep = "\t", header=TRUE)


###Sum abundances of matching vs non matching ASVs per sample
group_abd_df<-calc_group_sums(abd_df=asv_perc_abd_df,info_df=asv_info_df,first_group_var="Matches_Isolate")
group_abd_df$MG_ID<-as.factor(group_abd_df$Sample_UID)
group_abd_df$Sample_UID<-NULL
group_abd_df<-merge(group_abd_df,metadata_df,by="MG_ID",all.x=TRUE)
#colnames(group_abd_df)[1]<-"MG_ID"
colnames(group_abd_df)[2]<-"NonMatching_Abundance_Sum"
colnames(group_abd_df)[3]<-"Matching_Abundance_Sum"
summary(group_abd_df)


figX<-ggplot(group_abd_df,aes(fill=Depth,y=Matching_Abundance_Sum,x=Year))+
scale_fill_gradientn(colours = depth_col_grad)+ 
geom_bar(position= "dodge2",stat="identity",alpha=0.9,colour="black")+
coord_flip()+
theme_bw()+
xlab("Sample year")+
ylab("Sum of % Abundance of ASVs matching Isolates")+
facet_wrap(Core_ID ~.)

ggsave("/mnt/smart/users/fcoutinho/Antarctic_Lagoons/DADA_Output/Antarctic_Lagoons_Matched_ASV_Abundance_Sums_by_Year.pdf",plot=figX,device=cairo_pdf,width=12,height=8,pointsize=8)

figX<-ggplot(group_abd_df,aes(x=Depth,y=Matching_Abundance_Sum,fill=Year))+
scale_x_reverse()+
scale_fill_gradientn(colours = year_col_grad)+ 
geom_bar(position="stack",stat="identity",alpha=0.9,colour="black")+
coord_flip()+
theme_bw()+
xlab("Sample Depth (cm)")+
ylab("Sum of % Abundance of ASVs matching Isolates")+
facet_wrap(Core_ID ~.)

ggsave("/mnt/smart/users/fcoutinho/Antarctic_Lagoons/DADA_Output/Antarctic_Lagoons_Matched_ASV_Abundance_Sums_by_Depth.pdf",plot=figX,device=cairo_pdf,width=12,height=8,pointsize=8)

###% of total sample ASVs with non zero abundances matching Isolates 


m_asv_perc_abd_df<-melt(asv_perc_abd_df,id="MG_ID",value.name="Abundance",variable.name="ASV_UID")
f_m_asv_perc_abd_df<-m_asv_perc_abd_df[which(m_asv_perc_abd_df$Abundance > 0),]
f_m_asv_perc_abd_df<-merge(f_m_asv_perc_abd_df,asv_info_df,by="ASV_UID",all.x=TRUE)

f_m_asv_perc_abd_df$Dummy<-1

smp_asv_tot_count<-aggregate(f_m_asv_perc_abd_df$Dummy, by=list(Sample=f_m_asv_perc_abd_df$MG_ID), FUN=sum)
colnames(smp_asv_tot_count)<-c("MG_ID","ASV_Richness")

smp_asv_match_count<-aggregate(f_m_asv_perc_abd_df$Matches_Isolate, by=list(Sample=f_m_asv_perc_abd_df$MG_ID), FUN=sum)

smp_asv_tot_count$Matching_ASVs<-smp_asv_match_count$x
smp_asv_tot_count$Percentage_Matching_ASVs<-smp_asv_tot_count$Matching_ASVs/smp_asv_tot_count$ASV_Richness*100

group_abd_df<-merge(smp_asv_tot_count,metadata_df,by="MG_ID",all.x=TRUE)



figX<-ggplot(group_abd_df,aes(x=Depth,y=Percentage_Matching_ASVs,fill=Year))+
scale_x_reverse()+
scale_fill_gradientn(colours = year_col_grad)+ 
geom_bar(position="stack",stat="identity",alpha=0.9,colour="black")+
coord_flip()+
theme_bw()+
xlab("Sample Depth (cm)")+
ylab("% of Sample ASVs that match Isolates")+
facet_wrap(Core_ID ~.)

ggsave("/mnt/smart/users/fcoutinho/Antarctic_Lagoons/DADA_Output/Antarctic_Lagoons_Percentage_of_ASVs_Matching_Isolates.pdf",plot=figX,device=cairo_pdf,width=12,height=8,pointsize=8)

###ASV x isolates blast
blast_df<-read.table(file = "/mnt/smart/users/fcoutinho/Antarctic_Lagoons/All_Isolates-x-All_ASVs.blastn", sep = "\t", header=FALSE, comment="",stringsAsFactors=TRUE)
colnames(blast_df)<-c("ASV_UID","Isolate_ID","Identity","Alignment_Length","Mismatches","Gap_Openings","Q_Start","Q_End","S_Start","S_End","E_Value","Bit_Score")
summary(blast_df)

f_blast_df<-blast_df[which(blast_df$Identity == 100 & blast_df$Alignment_Length >= 300),]

f_blast_df<-merge(f_blast_df,asv_info_df,by="ASV_UID",all.x=TRUE)

library(data.table)

iso_f_blast_df<-f_blast_df[f_blast_df[, .I[which.max(Bit_Score)], by=Isolate_ID]$V1]

summary(iso_f_blast_df)

write.table(iso_f_blast_df, file = "/mnt/smart/users/fcoutinho/Antarctic_Lagoons/Antarctic_Lagoons_ASVsxIsolates_ID_100_Ali_300bp_Best_Isolate_Matches.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

asv_f_blast_df<-f_blast_df[f_blast_df[, .I[which.max(Bit_Score)], by=ASV_UID]$V1]

summary(asv_f_blast_df)

write.table(asv_f_blast_df, file = "/mnt/smart/users/fcoutinho/Antarctic_Lagoons/Antarctic_Lagoons_ASVsxIsolates_ID_100_Ali_300bp_Best_ASV_Matches.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

write.table(asv_info_df, file = "/mnt/smart/users/fcoutinho/Antarctic_Lagoons/DADA_Output/DADA2_ASVs_Info_Antarctic_Lagoons.tsv", sep = "\t", quote = FALSE, row.names = FALSE)


###Abundance stats
stats_df<-calc_abund_stats(abd_df=asv_perc_abd_df,exclude_cols=c("MG_ID"))
stats_df<-merge(stats_df,asv_info_df,by.x="UID",by.y="ASV_UID",all.x=TRUE)
summary(stats_df)
#Make scatterplot of mean abundance (y-axis) vs prevalence (x-axis) colouring dots by phylum and highlighting which ones match the isolates
figX<-ggplot(stats_df,aes(y=Mean_Abundance,x=Prevalence,fill=Phylum))+
geom_point(aes(size=Matches_Isolate),shape=23, stroke = 0.2,alpha=0.8)+ # 
#scale_y_continuous(trans = log_trans())+
scale_y_log10()+
theme_bw()+
theme(legend.position='top')+
#scale_fill_manual(name="Phylum",values=phylum_coloring)+#, limits = c(35, 0), name="Depth", trans="reverse")+
labs(x="ASV Prevalence",y="Mean Abundance")

ggsave("/mnt/smart/users/fcoutinho/Antarctic_Lagoons/DADA_Output/Antarctic_Lagoons_ASV_PrevalencexAbundance_Scatterplot.pdf",plot=figX,device=cairo_pdf,width=15,height=10,pointsize=8)

figX<-ggplot(stats_df[which(stats_df$Matches_Isolate==TRUE),],aes(y=Mean_Abundance,x=Prevalence))+
#geom_point(shape=23, stroke = 0.2,alpha=0.9,size=3)+ # 
#scale_y_continuous(trans = log_trans())+
geom_label(aes(label=UID,y=Mean_Abundance,x=Prevalence,fill=Genus),size=2,inherit.aes=FALSE,position=position_jitter(width=12,height=0))+
scale_y_log10()+
theme_bw()+
theme(legend.position='top')+
#scale_fill_manual(name="Phylum",values=phylum_coloring)+#, limits = c(35, 0), name="Depth", trans="reverse")+
labs(x="ASV Prevalence",y="Mean Abundance")

ggsave("/mnt/smart/users/fcoutinho/Antarctic_Lagoons/DADA_Output/Antarctic_Lagoons_ASV_PrevalencexAbundance_Scatterplot_Matching_ASVs_Only.pdf",plot=figX,device=cairo_pdf,width=7,height=7,pointsize=8)


###Species accumulation curves by sample
#Create a vector or sample sizes
sample_vec<-seq(from=100, to=150000, by=500)
#Caculate rarefactioin (exclude first column as it is only the sample ID)
raref_result<-as.data.frame(rarefy(asv_raw_abd_df[,-1],sample=sample_vec))
#Conver the output to a dataframe, then melt
raref_result$MG_ID<-as.factor(asv_raw_abd_df$MG_ID)
m_raref_result<-melt(raref_result,id="MG_ID",value.name="Richness",variable.name="Sample_Size")
#Convert the Sample_Sixe column to numeric
m_raref_result$Sample_Size<-as.numeric(gsub("N","",m_raref_result$Sample_Size))
#Merge with metadata and make scatterplot (use geom_point)
m_raref_result<-merge(m_raref_result,metadata_df,by="MG_ID",all.x=TRUE)

figX<-ggplot(m_raref_result,aes(y=Richness,x=Sample_Size,fill=Phylum))+
geom_point(aes(colour=Year,shape=Core_ID),size=1.5, stroke = 0.2,alpha=0.8)+
theme_bw()+
scale_fill_manual(name="Phylum",values=phylum_coloring)+#, limits = c(35, 0), name="Depth", trans="reverse")+
labs(x="Sample Size",y="ASV Richness")

ggsave("/mnt/smart/users/fcoutinho/Antarctic_Lagoons/DADA_Output/Antarctic_Lagoons_Rarefaction_Curves.pdf",plot=figX,device=cairo_pdf,width=11,height=7,pointsize=8)

# spec_acc_result<-specaccum(asv_raw_abd_df[,-1], method = "rarefaction", permutations = 10)
# #raref_result<-rarecurve(asv_raw_abd_df[,c(2:201)],step=100,sample=1000,tidy=TRUE)
# #raref_result<-rarecurve(asv_raw_abd_df[,-1],step=100,sample=10000,tidy=TRUE)
# spec_acc_result<-specaccum(asv_raw_abd_df[,2:2000], method = "rarefaction", permutations = 3)

# rich_df<-as.data.frame(cbind(spec_acc_result$sites,spec_acc_result$richness,spec_acc_result$individuals))
# colnames(rich_df)<-c("Sites","Richness","Individuals")

###ASV x Matched x Abundance
m_asv_perc_abd_df<-melt(asv_perc_abd_df,id="MG_ID",value.name="Abundance",variable.name="ASV_UID")

colnames(m_asv_perc_abd_df)[1]<-"MG_ID"

m_asv_perc_abd_df<-merge(m_asv_perc_abd_df,asv_info_df[,c("ASV_UID","Matches_Isolate")],by="ASV_UID",all.x=TRUE)

m_asv_perc_abd_df<-merge(m_asv_perc_abd_df,subset(metadata_df,select=c(MG_ID,Ti,V,Cr,Mn,Fe,Co,Ni,Cu,Zn,Mo,Cd,Pb,Year,Depth_cm),by="MG_ID",all.x=TRUE))

summary(m_asv_perc_abd_df[which(m_asv_perc_abd_df$Matches_Isolate == TRUE),"Abundance"])
summary(m_asv_perc_abd_df[which(m_asv_perc_abd_df$Matches_Isolate == FALSE),"Abundance"])

figX<-ggplot(m_asv_perc_abd_df,aes(x=Matches_Isolate,y=log10(Abundance)))+
geom_boxplot(alpha=0.9)+
geom_signif(comparisons = list(c("TRUE","FALSE")),map_signif_level = TRUE, textsize=2.5,test = "wilcox.test",step_increase=0.05,tip_length = 0.01,size=0.1,lwd=0.1)+ #,c("Epipelagic","Bathypelagic")
theme_bw()+
#scale_fill_manual(name="Zone",values=zone_coloring)+
theme(legend.position='top',text = element_text(size = 9))
#+labs(y="Defense Systems in MAGs")

ggsave("test.pdf",plot=figX,device=cairo_pdf,width=7,height=5,pointsize=8)


#ASV x Isolate blast
f_blast_df<-read.table(file ="/mnt/smart/users/fcoutinho/Antarctic_Lagoons/Antarctic_Lagoons_ASVsxIsolates_ID_100_Ali_300bp.tsv", sep = "\t", header=TRUE, comment="",stringsAsFactors=TRUE)

unq_f_blast_df<-f_blast_df[!duplicated(f_blast_df$ASV_UID),]

unq_f_blast_df$Matches_Isolate<-TRUE

asv_info_df<-merge(asv_info_df,unq_f_blast_df[,c("ASV_UID","Matches_Isolate")],by="ASV_UID",all.x=TRUE)

###ASV x Metadata correls
abd_meta_df<-merge(asv_perc_abd_df,subset(metadata_df,select=c(MG_ID,Ti,V,Cr,Mn,Fe,Co,Ni,Cu,Zn,Mo,Cd,Pb,Year,Depth_cm)),by="MG_ID",all.x=TRUE)

cor_df<-cor(abd_meta_df[complete.cases(abd_meta_df),c(2:8638)],abd_meta_df[complete.cases(abd_meta_df),c(8639:8652)], method="spearman",use="everything")
m_cor_df<-melt(cor_df,value.name="Spearman_Rho")

colnames(m_cor_df)<-c("ASV_UID","Variable","Spearman_Rho")
summary(m_cor_df)

#Correl heatmap
fm_cor_df<-m_cor_df[!is.na(m_cor_df$Spearman_Rho),]
summary(fm_cor_df)

fm_cor_df<-fm_cor_df[which((fm_cor_df$Spearman_Rho >= 0.5) | (fm_cor_df$Spearman_Rho <= -0.5)),] # 
summary(fm_cor_df)

fm_cor_df<-merge(fm_cor_df,asv_info_df,by="ASV_UID",all.x=TRUE)

fm_cor_df<-merge(fm_cor_df,unq_f_blast_df[,c("ASV_UID","Matches_Isolate")],by="ASV_UID",all.x=TRUE)
fm_cor_df$Alias<-""
fm_cor_df$Alias[which(!is.na(fm_cor_df$Matches_Isolate))]<-"***" # Add asterisk to ASVs with isolate matches
fm_cor_df$Alias<-as.factor(paste(fm_cor_df$Alias,fm_cor_df$ASV_UID,fm_cor_df$Genus,sep="|"))
fm_cor_df$Alias<-factor(fm_cor_df$Alias,levels=unique(fm_cor_df$Alias[order(fm_cor_df$Genus,decreasing=TRUE)]))
summary(fm_cor_df)

compl_col_grad<-rev(colorRampPalette(brewer.pal(11,"Spectral"))(n=100))
figX<-ggplot(fm_cor_df,aes(fill=Spearman_Rho,y=Alias,x=Variable))+
geom_tile(colour="black")+
scale_fill_gradientn(colours = compl_col_grad, limits=c(-0.75,0.75))+
theme_bw()+
theme(axis.text.x = element_text(size=10, angle = 45, hjust=1), axis.text.y = element_text(size=5), strip.text.y = element_text(angle = 0), legend.position="top")+
facet_grid(Phylum ~ . , scales="free_y", space="free")
ggsave("ASVxEnv_Var_Heatmap_Antarctic_Lagoons.pdf",plot=figX,device=cairo_pdf,width=7,height=17,pointsize=8)

###ASV x PhytoREf blast
phyto_df<-read.table(file = "/mnt/smart/users/fcoutinho/Antarctic_Lagoons/PhytoRef/Antarctic_Lagoons_ASVsxPhytoRef", sep = "\t", header=FALSE, comment="",stringsAsFactors=TRUE)
colnames(phyto_df)<-c("ASV_UID","PhytoRef","Identity","Alignment_Length","Mismatches","Gap_Openings","Q_Start","Q_End","S_Start","S_End","E_Value","Bit_Score")
f_phyto_df<-phyto_df[which(phyto_df$Identity >= 97 & phyto_df$Alignment_Length >= 300),]
summary(f_phyto_df)

f_phyto_df$PhytoRef<-as.factor(gsub("#","|",f_phyto_df$PhytoRef,perl=TRUE))

fu_phyto_df<-f_phyto_df#[!duplicated(f_phyto_df$ASV_UID),]

taxon<-data.frame(do.call('rbind', strsplit(as.character(fu_phyto_df$PhytoRef),'|',fixed=TRUE)))

fu_phyto_df<-as.data.frame(cbind(fu_phyto_df,taxon))

write.table(fu_phyto_df, file = "/mnt/smart/users/fcoutinho/Antarctic_Lagoons/PhytoRef/Antarctic_Lagoons_ASVsxPhytoRef_97_Ali_300bp.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

###NMDS
dist_metric<-"bray"
dists<-vegdist(asv_perc_abd_df[,-1], method = dist_metric)

set.seed(666)
mdsresult<-metaMDS(dists,distance = dist_metric,k = 2,maxit = 999)
data.scores<-as.data.frame(scores(mdsresult))
data.scores$Sample<-rownames(data.scores)

mdata<-data.scores

core_id<-data.frame(do.call('rbind', strsplit(as.character(mdata$Sample),'_',fixed=TRUE)))
colnames(core_id)<-c("Core","Depth_cm","Replicate")
mdata<-as.data.frame(cbind(mdata,core_id))
mdata$Depth_cm<-as.numeric(as.character(mdata$Depth_cm))
mdata<-merge(mdata,metadata_df,by.x="Sample",by.y="MG_ID",all.x=TRUE)
summary(mdata)

###NMDS x Env Var Correlations
clean_data<-subset(mdata,select=c(NMDS1,NMDS2,Depth,Ti,V,Cr,Mn,Fe,Co,Ni,Cu,Zn,Mo,Cd,Pb,Year))
clean_data<-clean_data[complete.cases(clean_data),]
cor_df<-cor(clean_data)

write.table(cor_df,file="NMDS_Env_Var_Correlations_Antarctic_Lagoons.tsv",sep="\t",quote=FALSE,row.names=TRUE,col.names=NA)

####NMDS Year
figX<-ggplot(mdata, aes(x=NMDS1,y=NMDS2))+
geom_point(size=3.5, aes(colour=Year, shape=Core))+
theme_bw()+theme(text=element_text(size=16))+
scale_colour_gradientn(colours = year_col_grad, name="Year",trans = 'reverse')
#+scale_fill_gradientn(colours = depth_col_grad, name="Depth",trans = 'reverse') ,limits = c(4000, 0),
ggsave("NMDS_Antarctic_Lagoons_Year.pdf",plot=figX,device=cairo_pdf,width=6,height=5,pointsize=8)

###NMDS Depth
depth_col_pal<-brewer.pal(9,"RdBu")
depth_col_grad<-rev(colorRampPalette(depth_col_pal)(n=299))
figX<-ggplot(mdata, aes(x=NMDS1,y=NMDS2))+geom_point(size=3.5, aes(colour=Depth_cm, shape=Core))+theme_bw()+theme(text=element_text(size=16))+scale_colour_gradientn(colours = depth_col_grad, name="Depth",trans = 'reverse')
#+scale_fill_gradientn(colours = depth_col_grad, name="Depth",trans = 'reverse') ,limits = c(4000, 0),
ggsave("NMDS_Antarctic_Lagoons_Depth.pdf",plot=figX,device=cairo_pdf,width=6,height=5,pointsize=8)


###Rank Abundance
mean_abd_df<-as.data.frame(cbind(colnames(asv_perc_abd_df[,-1]),colMeans(asv_perc_abd_df[,-1])))

colnames(mean_abd_df)<-c("ASV","Mean_Abundance")

mean_abd_df$Mean_Abundance<-as.numeric(as.character(mean_abd_df$Mean_Abundance))

mean_abd_df$Rank<-order(mean_abd_df$Mean_Abundance,decreasing=TRUE)

summary(mean_abd_df)
head(mean_abd_df)

write.table(mean_abd_df,file="/mnt/smart/users/fcoutinho/Antarctic_Lagoons/DADA_Output/DADA2_ASVs_Rank_Abundance_Antarctic_Lagoons.tsv",sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)

###ANOSIM
full_data<-merge(asv_perc_abd_df,metadata_df,by="MG_ID")

anosim_result<-anosim(full_data[,c(2:ncol(asv_perc_abd_df))], full_data$Core_ID, distance="bray",permutations=999)

summary(anosim_result)

###Vriation partitioning (requires at least 2 tables)
vpart_result<-varpart(full_data[,c(2:ncol(asv_perc_abd_df))], subset(full_data,select=c(Depth,Ti,V,Cr,Mn,Fe,Co,Ni,Cu,Zn,Mo,Cd,Pb)))

###taxon abundance plots
asv_info_df$Phylum[which(asv_info_df$Phylum == "Proteobacteria")]<-asv_info_df$Class[which(asv_info_df$Phylum == "Proteobacteria")]

summary(asv_info_df)

colnames(asv_perc_abd_df)[1]<-"Sample"
m_asv_perc_abd_df<-melt(asv_perc_abd_df,id="Sample",value.name="Abundance",variable.name="ASV_UID")
m_asv_perc_abd_df<-merge(m_asv_perc_abd_df,asv_info_df,by="ASV_UID")

summary(m_asv_perc_abd_df)

tax_level<-"Genus"
tax_abd_df<-aggregate(m_asv_perc_abd_df$Abundance, by=list(Sample=m_asv_perc_abd_df$Sample, Taxon=m_asv_perc_abd_df[[tax_level]]), FUN=sum)
colnames(tax_abd_df)<-c("Sample","Taxon","Abundance")

summary(tax_abd_df)

#Calculate mean abundances per taxon
tax_mean_abds<-aggregate(tax_abd_df$Abundance, by=list(Taxon=tax_abd_df$Taxon), FUN=mean)
summary(tax_mean_abds)

valid_tax<-tax_mean_abds$Taxon[which(tax_mean_abds$x >= 0.01)]

m_asv_perc_abd_df[[tax_level]]<-as.character(m_asv_perc_abd_df[[tax_level]])
m_asv_perc_abd_df[[tax_level]][!(m_asv_perc_abd_df[[tax_level]] %in% valid_tax)]<-"Others"
m_asv_perc_abd_df[[tax_level]]<-as.factor(m_asv_perc_abd_df[[tax_level]])
summary(m_asv_perc_abd_df)

tax_abd_df<-aggregate(m_asv_perc_abd_df$Abundance, by=list(Sample=m_asv_perc_abd_df$Sample, Taxon=m_asv_perc_abd_df[[tax_level]]), FUN=sum)
colnames(tax_abd_df)<-c("Sample","Taxon","Abundance")
summary(tax_abd_df)

#Now Add core info to df
core_id<-data.frame(do.call('rbind', strsplit(as.character(tax_abd_df$Sample),'_',fixed=TRUE)))
colnames(core_id)<-c("Core","Depth_cm","Replicate")
tax_abd_df<-as.data.frame(cbind(tax_abd_df,core_id))
tax_abd_df$Depth_cm<-as.numeric(tax_abd_df$Depth_cm)*-1
summary(tax_abd_df)

###Heatmap
compl_col_grad<-rev(colorRampPalette(brewer.pal(11,"Spectral"))(n=100))

figX<-ggplot(tax_abd_df,aes(fill=log10(Abundance),y=Depth_cm,x=Taxon))+
geom_tile()+
scale_fill_gradientn(colours = compl_col_grad)+
theme_bw()+
theme(axis.text.x = element_text(angle = 45,hjust = 1,size=4))+
facet_wrap(Core ~.)

ggsave("Genus_Heatmap_Antarctic_Lagoons.pdf",plot=figX,device=cairo_pdf,width=35,height=8,pointsize=8)

#Calculate medin abundances per taxon
tax_median_abds<-aggregate(tax_abd_df$Abundance, by=list(Taxon=tax_abd_df$Taxon), FUN=median)
summary(tax_median_abds)
valid_tax<-tax_median_abds$Taxon[which(tax_median_abds$x >= 0.5)]

#Recalculate taxon abundances sums grouping small abundance phyla into "others"
tax_level<-"Phylum"
m_asv_perc_abd_df[[tax_level]]<-as.character(m_asv_perc_abd_df[[tax_level]])
m_asv_perc_abd_df[[tax_level]][!(m_asv_perc_abd_df[[tax_level]] %in% valid_tax)]<-"Others"
m_asv_perc_abd_df[[tax_level]]<-as.factor(m_asv_perc_abd_df[[tax_level]])
summary(m_asv_perc_abd_df)

tax_abd_df<-aggregate(m_asv_perc_abd_df$Abundance, by=list(Sample=m_asv_perc_abd_df$Sample, Taxon=m_asv_perc_abd_df$Phylum), FUN=sum)
colnames(tax_abd_df)<-c("Sample","Taxon","Abundance")

summary(tax_abd_df)

#Now Add core info to df
core_id<-data.frame(do.call('rbind', strsplit(as.character(tax_abd_df$Sample),'_',fixed=TRUE)))
colnames(core_id)<-c("Core","Depth_cm","Replicate")

tax_abd_df<-as.data.frame(cbind(tax_abd_df,core_id))

tax_abd_df$Depth_cm<-as.numeric(tax_abd_df$Depth_cm)*-1

summary(tax_abd_df)

#Make Phylum barplot

fig2B<-ggplot(tax_abd_df,aes(y=Abundance,x=Depth_cm,fill=Taxon))+
geom_bar(position="stack",stat="identity",alpha=0.9,colour="black")+
scale_fill_manual(name="Phylum",values=phylum_coloring)+#,breaks=custom_metab_order
theme(legend.position="top")+
theme_bw()+
coord_flip()+facet_wrap(Core ~.)

ggsave("Phylum_Barplots_Antarctic_Lagoons.pdf",plot=fig2B,device=cairo_pdf,width=15,height=10,pointsize=8)




#Perform PCA ordination
pca_metadata<-subset(metadata_df,select=c(Sample_UID,Ti,V,Cr,Mn,Fe,Co,Ni,Cu,Zn,Mo,Cd,Pb))
rownames(pca_metadata)<-pca_metadata$Sample_UID
pca_metadata$Sample_UID<-NULL
sdata<-as.data.frame(scale(pca_metadata,center=TRUE,scale=TRUE))


set.seed(666)
rdadata<-rda(sdata)

#Get variance explained values
xlabel<-round((summary(rdadata)$cont$importance[2,1]*100),digits=1)
xlabel<-paste("PC1 (",xlabel,"% explained)",sep="")

ylabel<-round((summary(rdadata)$cont$importance[2,2]*100),digits=1)
ylabel<-paste("PC2 (",ylabel,"% explained)",sep="")

#Convert output to DF compatible with ggplot2
uscores <- data.frame(rdadata$CA$u)

uscores1<-merge(uscores,metadata_df,by="row.names",all.x=TRUE)

rownames(uscores1)[1]<-"Sample_UID"

vscores <- data.frame(rdadata$CA$v)


#Make plots
figX<-ggplot(uscores1)+
geom_segment(data = vscores, aes(x = 0, y = 0, xend = PC1, yend = PC2), arrow=arrow(length=unit(0.2,"cm")),alpha = 0.5, color = 'darkgreen')+
geom_label(data = vscores, aes(x = PC1, y = PC2, label = rownames(vscores)), fill="grey",col = 'black',size=3,label.padding = unit(0.05, "lines"))+
geom_point(data=uscores1,aes(x = PC1, y = PC2, colour=Year ,shape=Core_ID),alpha=1,size=4)+
scale_colour_gradientn(colours = year_col_grad)+ #,limits = c(35, 0),name="Depth",
labs(y=ylabel,x=xlabel)+
theme_bw()+
theme(legend.position='right')#+geom_text(data=uscores1,aes(x = PC1, y = PC2, label=Station),size=1.2)
ggsave("Antarctic_Lagoons_RDA_Year.pdf",plot=figX,device=cairo_pdf,width=7,height=5,pointsize=8)


#Make plots
figX<-ggplot(uscores1)+
geom_segment(data = vscores, aes(x = 0, y = 0, xend = PC1, yend = PC2), arrow=arrow(length=unit(0.2,"cm")),alpha = 0.5, color = 'darkgreen')+
geom_label(data = vscores, aes(x = PC1, y = PC2, label = rownames(vscores)), fill="grey",col = 'black',size=3,label.padding = unit(0.05, "lines"))+
geom_point(data=uscores1,aes(x = PC1, y = PC2, colour=Depth ,shape=Core_ID),alpha=1,size=4)+
scale_colour_gradientn(colours = depth_col_grad)+ #,limits = c(35, 0),name="Depth",
labs(y=ylabel,x=xlabel)+
theme_bw()+
theme(legend.position='right')#+geom_text(data=uscores1,aes(x = PC1, y = PC2, label=Station),size=1.2)
ggsave("Antarctic_Lagoons_RDA_Depth.pdf",plot=figX,device=cairo_pdf,width=7,height=5,pointsize=8)



####Alpha diversity
asv_perc_abd_df<-read.table(file = "/mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/Antarctic_Lagoon/DADA/DADA2_ASVs_Percentage_Abundances_Antarctic_Lagoons.tsv", row.names=1, sep = "\t", header=TRUE)

shanon_div<-as.data.frame(diversity(asv_perc_abd_df, index = "shannon"))
shanon_div$Sample_UID<-rownames(shanon_div)
colnames(shanon_div)[1]<-"Shanon_Diversity_Index"

shanon_div<-merge(shanon_div,metadata_df,by.x="Sample_UID",by.y="MG_ID",all.x=TRUE)

fig_1D<-ggplot(shanon_div,aes(y=Shanon_Diversity_Index,x=Depth))+
geom_smooth(method = "loess",se=TRUE,  formula = y ~ x)+
scale_x_reverse(limits = c(35, 0))+
geom_point(aes(colour=Depth,shape=Core_ID),size=3.5, stroke = 0.2)+
coord_flip()+
theme_bw()+
scale_colour_gradientn(colours = depth_col_grad)+#, limits = c(35, 0), name="Depth", trans="reverse")+
labs(x="Depth",y="Shanon Diversity Index")+
theme(legend.position='none')

ggsave("Antarctic_Lagoons_DiversityxDepth_Scatterplots.pdf",plot=fig_1D,device=cairo_pdf,width=5,height=7,pointsize=8)
