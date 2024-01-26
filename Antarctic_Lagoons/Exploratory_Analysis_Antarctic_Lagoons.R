library(vegan)
library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)

metadata_df<-read.table(file = "/mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/Antarctic_Lagoon/Metadata/Sample_Metadata.tsv", sep = "\t", header=TRUE)
rownames(metadata_df)<-metadata_df$Sample_UID

metadata_df$Depth <- metadata_df$Depth - 0.5

metadata_df$MG_ID<-rownames(metadata_df)
metadata_df$MG_ID<-gsub("Core(\\s)+","",metadata_df$MG_ID,perl=TRUE)
metadata_df$MG_ID<-gsub("Cm(\\s)","",metadata_df$MG_ID,perl=TRUE)
metadata_df$MG_ID<-as.factor(paste(metadata_df$MG_ID,"_S1",sep=""))

pca_metadata<-subset(metadata_df,select=c(Sample_UID,Ti,V,Cr,Mn,Fe,Co,Ni,Cu,Zn,Mo,Cd,Pb))
rownames(pca_metadata)<-pca_metadata$Sample_UID
pca_metadata$Sample_UID<-NULL
sdata<-as.data.frame(scale(pca_metadata,center=TRUE,scale=TRUE))

###Phylum barplot
asv_perc_abd_df<-read.table(file = "/mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/Antarctic_Lagoon/DADA/DADA2_ASVs_Percentage_Abundances_Antarctic_Lagoons.tsv", sep = "\t", header=TRUE)

asv_info_df<-read.table(file = "/mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/Antarctic_Lagoon/DADA/DADA2_ASVs_Info_Antarctic_Lagoons.tsv", sep = "\t", header=TRUE)

asv_info_df$Phylum[which(asv_info_df$Phylum == "Proteobacteria")]<-asv_info_df$Class[which(asv_info_df$Phylum == "Proteobacteria")]

summary(asv_info_df)

colnames(asv_perc_abd_df)[1]<-"Sample"
m_asv_perc_abd_df<-melt(asv_perc_abd_df,id="Sample",value.name="Abundance",variable.name="ASV_UID")
m_asv_perc_abd_df<-merge(m_asv_perc_abd_df,asv_info_df,by="ASV_UID")

summary(m_asv_perc_abd_df)

tax_abd_df<-aggregate(m_asv_perc_abd_df$Abundance, by=list(Sample=m_asv_perc_abd_df$Sample, Taxon=m_asv_perc_abd_df$Phylum), FUN=sum)
colnames(tax_abd_df)<-c("Sample","Taxon","Abundance")

summary(tax_abd_df)

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
phylum_coloring<-c(brewer.pal(9,"Set1"),brewer.pal(4,"Pastel1"),brewer.pal(8,"Dark2")[8],brewer.pal(5,"Set2")[5],brewer.pal(8,"Accent"),brewer.pal(5,"Dark2"))
names(phylum_coloring)<-c("Actinobacteriota","Alphaproteobacteria","Synergistota","Crenarchaeota","Bacteroidota","Nanoarchaeota","Verrucomicrobiota","Marinisomatota","Others","Bdellovibrionota","Gammaproteobacteria","Acidobacteriota","UBP7","Cloacimonadota","Firmicutes","Margulisbacteria","Gemmatimonadota","Desulfobacterota","Nitrospinota","Desulfobacterota_D","Patescibacteria","Hydrothermarchaeota","Chlamydiota","SAR324","Halobacterota","Chloroflexi","Thermoplasmatota","Myxococcota")

fig2B<-ggplot(tax_abd_df,aes(y=Abundance,x=Depth_cm,fill=Taxon))+
geom_bar(position="stack",stat="identity",alpha=0.9,colour="black")+
scale_fill_manual(name="Phylum",values=phylum_coloring)+#,breaks=custom_metab_order
theme_bw()+
coord_flip()+facet_wrap(Core ~.)

ggsave("Phylum_Barplots_Antarctic_Lagoons.pdf",plot=fig2B,device=cairo_pdf,width=15,height=10,pointsize=8)

###NMDS
dist_metric<-"bray"
dists<-vegdist(asv_perc_abd_df, method = dist_metric)

set.seed(666)
mdsresult<-metaMDS(dists,distance = dist_metric,k = 2,maxit = 999)
data.scores<-as.data.frame(scores(mdsresult))
data.scores$Sample<-rownames(data.scores)

mdata<-data.scores

core_id<-data.frame(do.call('rbind', strsplit(as.character(mdata$Sample),'_',fixed=TRUE)))
colnames(core_id)<-c("Core","Depth_cm","Replicate")
mdata<-as.data.frame(cbind(mdata,core_id))
summary(mdata)

library(RColorBrewer)
depth_col_pal<-brewer.pal(9,"RdBu")
depth_col_grad<-rev(colorRampPalette(depth_col_pal)(n=299))

mdata$Depth_cm<-as.numeric(as.character(mdata$Depth_cm))

figX<-ggplot(mdata, aes(x=NMDS1,y=NMDS2))+geom_point(size=3.5, aes(colour=Depth_cm, shape=Core))+theme_bw()+theme(text=element_text(size=16))+scale_colour_gradientn(colours = depth_col_grad, name="Depth",trans = 'reverse')
#+scale_fill_gradientn(colours = depth_col_grad, name="Depth",trans = 'reverse') ,limits = c(4000, 0),

ggsave("NMDS_Antarctic_Lagoons.pdf",plot=figX,device=cairo_pdf,width=6,height=5,pointsize=8)



#Perform PCA ordination
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

depth_col_pal<-brewer.pal(9,"BrBG")#[c(1,4,9)]
depth_col_grad<-colorRampPalette(depth_col_pal)(n=99)

#Make plots
figX<-ggplot(uscores1)+
geom_segment(data = vscores, aes(x = 0, y = 0, xend = PC1, yend = PC2), arrow=arrow(length=unit(0.2,"cm")),alpha = 0.5, color = 'darkgreen')+
geom_label(data = vscores, aes(x = PC1, y = PC2, label = rownames(vscores)), fill="grey",col = 'black',size=3,label.padding = unit(0.05, "lines"))+
geom_point(data=uscores1,aes(x = PC1, y = PC2, colour=Depth ,shape=Core_ID),alpha=1,size=4)+
scale_colour_gradientn(colours = depth_col_grad)+ #,limits = c(35, 0),name="Depth",
labs(y=ylabel,x=xlabel)+
theme_bw()+
theme(legend.position='right')#+geom_text(data=uscores1,aes(x = PC1, y = PC2, label=Station),size=1.2)
ggsave("Antarctic_Lagoons_RDA.pdf",plot=figX,device=cairo_pdf,width=7,height=5,pointsize=8)


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

###NMDS
dist_metric<-"bray"
dists<-vegdist(asv_perc_abd_df, method = dist_metric)

set.seed(666)
mdsresult<-metaMDS(dists,distance = dist_metric,k = 2,maxit = 999)
data.scores<-as.data.frame(scores(mdsresult))
data.scores$Sample<-rownames(data.scores)

mdata<-data.scores

core_id<-data.frame(do.call('rbind', strsplit(as.character(mdata$Sample),'_',fixed=TRUE)))
colnames(core_id)<-c("Core","Depth_cm","Replicate")
mdata<-as.data.frame(cbind(mdata,core_id))
summary(mdata)

library(RColorBrewer)
depth_col_pal<-brewer.pal(9,"RdBu")
depth_col_grad<-rev(colorRampPalette(depth_col_pal)(n=299))

mdata$Depth_cm<-as.numeric(as.character(mdata$Depth_cm))

figX<-ggplot(mdata, aes(x=NMDS1,y=NMDS2))+geom_point(size=3.5, aes(colour=Depth_cm, shape=Core))+theme_bw()+theme(text=element_text(size=16))+scale_colour_gradientn(colours = depth_col_grad, name="Depth",trans = 'reverse')
#+scale_fill_gradientn(colours = depth_col_grad, name="Depth",trans = 'reverse') ,limits = c(4000, 0),

ggsave("NMDS_Antarctic_Lagoons.pdf",plot=figX,device=cairo_pdf,width=6,height=5,pointsize=8)

