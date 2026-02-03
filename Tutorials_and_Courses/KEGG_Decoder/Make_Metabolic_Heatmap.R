library(reshape2)
library(ggplot2)
library(RColorBrewer)

args <- commandArgs(trailingOnly = TRUE)

metab_df<-read.table(file = args[1], sep = "\t", header=TRUE, comment="",stringsAsFactors=TRUE,check.names=FALSE)

colnames(metab_df)[1]<-"Genome"

mdata<-melt(metab_df, id.vars="Genome", variable.name="Pathway", value.name="Completeness")

# col_pal<-brewer.pal(9,"GnBu")
# col_pal<-brewer.pal(11,"BrBG")[-6]
col_pal<-brewer.pal(9,"Greens")
col_grad<-colorRampPalette(col_pal)(n=100)

#Cap the values of the completeness to 1
mdata$Completeness[mdata$Completeness>1]<-1

mdata$Group<-"All"

if (args[2] == "edit") {
    mdata<-mdata[grepl("TRUE",mdata$Genome),]
    library(tidyr)
    genomebkp<-mdata$Genome
    mdata<-mdata %>% separate(Genome, c("Cluster","Taxonomy","is_Cluster_Rep","OID"),sep="\\|")
    mdata$Genome<-genomebkp
    mdata$Group<-mdata$Cluster

}


figa<-ggplot(mdata,aes(fill=Completeness,y=Genome,x=Pathway))+geom_tile(colour = "grey50")+theme_bw()+theme(axis.text.x = element_text(angle = 90,hjust = 1))+scale_fill_gradientn(name="% Completeness",colours =col_grad,limits = c(0,1))+theme(legend.position="top")
ggsave("Pathway_Completenes_Heatmap.pdf",plot=figa,device=cairo_pdf,width=35,height=10)+facet_grid(Group ~ ., scales="free", space="free")


figb<-ggplot(mdata[mdata$Completeness>=0.5,],aes(fill=Completeness,y=Genome,x=Pathway))+geom_tile(colour = "grey50")+theme_bw()+theme(axis.text.x = element_text(angle = 90,hjust = 1))+geom_text(size=2,aes(label=round(Completeness,digits=1)))+scale_fill_gradientn(name="% Completeness",colours =col_grad,limits = c(0,1))+theme(legend.position="top")+facet_grid(Group ~ ., scales="free", space="free")
ggsave("Pathway_Completenes_Heatmap_Min_50.pdf",plot=figb,device=cairo_pdf,width=35,height=12)

