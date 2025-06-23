library(reshape2)
library(ggplot2)
library(RColorBrewer)

metab_df<-read.table(file = "/mnt/smart/scratch/vir/felipe/Kegg_Decoder_Workshop/All_CDSxKOfam_Decoder.tsv", sep = "\t", header=TRUE, comment="",stringsAsFactors=TRUE,check.names=FALSE)

colnames(metab_df)[1]<-"Genome"

mdata<-melt(metab_df, id.vars="Genome", variable.name="Pathway", value.name="Completeness")

col_pal<-brewer.pal(9,"GnBu")
col_grad<-colorRampPalette(col_pal)(n=100)

figa<-ggplot(mdata,aes(fill=Completeness,y=Genome,x=Pathway))+geom_tile(colour = "grey50")+theme_bw()+theme(axis.text.x = element_text(angle = 90,hjust = 1))+scale_fill_gradientn(name="% Completeness",colours =col_grad,limits = c(0,1))+theme(legend.position="top")
ggsave("Pathway_Completenes_Heatmap.pdf",plot=figa,device=cairo_pdf,width=35,height=10)


figb<-ggplot(mdata[mdata$Completeness>=0.5,],aes(fill=Completeness,y=Genome,x=Pathway))+geom_tile(colour = "grey50")+theme_bw()+theme(axis.text.x = element_text(angle = 90,hjust = 1))+geom_text(size=2,aes(label=round(Completeness,digits=1)))+scale_fill_gradientn(name="% Completeness",colours =col_grad,limits = c(0,1))+theme(legend.position="top")
ggsave("Pathway_Completenes_Heatmap_Min_50.pdf",plot=figb,device=cairo_pdf,width=35,height=8)

