library(vegan)
library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)

metadata_df<-read.table(file = "/mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/Antarctic_Lagoon/Metadata/Sample_Metadata.tsv", sep = "\t", header=TRUE)
rownames(metadata_df)<-metadata_df$Sample_UID

metadata_df$MG_ID<-rownames(metadata_df)
metadata_df$MG_ID<-gsub("Core(\\s)+","",metadata_df$MG_ID,perl=TRUE)
metadata_df$MG_ID<-gsub("Cm(\\s)","",metadata_df$MG_ID,perl=TRUE)
metadata_df$MG_ID<-as.factor(paste(metadata_df$MG_ID,"_S1",sep=""))

pca_metadata<-subset(metadata_df,select=c(Sample_UID,Ti,V,Cr,Mn,Fe	,Co,Ni,Cu,Zn,Mo,Cd,Pb))
rownames(pca_metadata)<-pca_metadata$Sample_UID
pca_metadata$Sample_UID<-NULL

sdata<-as.data.frame(scale(pca_metadata,center=TRUE,scale=TRUE))

set.seed(666)
rdadata<-rda(sdata)

xlabel<-round((summary(rdadata)$cont$importance[2,1]*100),digits=1)
xlabel<-paste("PC1 (",xlabel,"% explained)",sep="")

ylabel<-round((summary(rdadata)$cont$importance[2,2]*100),digits=1)
ylabel<-paste("PC2 (",ylabel,"% explained)",sep="")

uscores <- data.frame(rdadata$CA$u)

uscores1<-merge(uscores,metadata_df,by="row.names",all.x=TRUE)

rownames(uscores1)[1]<-"Sample_UID"

vscores <- data.frame(rdadata$CA$v)

depth_col_pal<-brewer.pal(9,"BrBG")#[c(1,4,9)]
depth_col_grad<-colorRampPalette(depth_col_pal)(n=99)

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
scale_colour_gradientn(colours = depth_col_grad,limits = c(35, 0),name="Depth",trans = 'reverse')+
labs(x="Depth",y="Shanon Diversity Index")+
theme(legend.position='none')

ggsave("Antarctic_Lagoons_DiversityxDepth_Scatterplots.pdf",plot=fig_1D,device=cairo_pdf,width=6,height=5,pointsize=8)

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

