library(readxl)
# library(tidyverse)
# library(tibble)


bta_df<-as.data.frame(read_xlsx("/mnt/smart/users/fcoutinho/EAIMM/ARB_Quantification/EAIMM_AntibioticExp.xlsx",  
sheet = "BARCELONETA"))

summary(bta_df)

bla_df<-as.data.frame(read_xlsx("/mnt/smart/users/fcoutinho/EAIMM/ARB_Quantification/EAIMM_AntibioticExp.xlsx",  
sheet = "BLANES"))

summary(bla_df)

full_df<-as.data.frame(rbind(bta_df[,c(1:10)], bla_df[c(129:nrow(bla_df)),c(1:10)]))

#convert colums 1, 3, 4 , 6 ,7 and 8 of full_df to factors
for (i in c(1, 3, 4 , 6 ,7, 8)) {
  full_df[,i] <- as.factor(full_df[,i])
}

full_df$Temperature<-as.numeric(gsub("[a-z]|\\W","",full_df$Temperature,perl=TRUE,ignore=TRUE))


full_df<-full_df[!is.na(full_df$Treatment),]

colnames(full_df)[10]<-"CFU_per_ml"

full_df$Site<-as.factor(gsub("_.*$","",full_df$ID,perl=TRUE))

full_df$DAY<-as.character(full_df$DAY)

summary(full_df)

library(ggplot2)
library(reshape2)

sub_df<-full_df[which(full_df$Type=="Sea water" & full_df$Temperature==21),]


figA<-ggplot(sub_df,aes(x=DAY,y=CFU_per_ml,fill=Replicate))+
#geom_point(alpha=0.9,aes(fill=Treatment,shape=Replicate))+
geom_bar(position="dodge",stat="identity",alpha=0.9)+
#geom_boxplot(alpha=0.9,aes(fill=Treatment))+
#geom_bar(alpha=0.9,shape=23,aes(fill=Treatment))+
theme_bw()+
#scale_y_continuous(limits = c(0, 0.025))+
#labs(x = "Abundance of Archaea (RPKM)",y="Cas1 Abundance (% Mapped Reads)")+
#scale_fill_manual(name="Zone",values=zone_coloring)+
theme(legend.position='top')+
facet_grid( Treatment ~ Site, scales="free_y")

ggsave("/mnt/smart/users/fcoutinho/EAIMM/ARB_Quantification/ARB_Quantification_BlanesxBarceloneta_Comp.pdf",plot=figA,width=5,height=9,pointsize=8)

figA<-ggplot(full_df,aes(x=DAY,y=CFU_per_ml,fill=Treatment))+
#geom_point(alpha=0.9,aes(fill=Treatment,shape=Replicate))+
#geom_bar(position="dodge",stat="identity",alpha=0.9,aes(fill=Treatment))+
geom_boxplot(alpha=0.9,aes(fill=Treatment))+
#geom_bar(alpha=0.9,shape=23,aes(fill=Treatment))+
theme_bw()+
#scale_y_continuous(limits = c(0, 0.025))+
#labs(x = "Abundance of Archaea (RPKM)",y="Cas1 Abundance (% Mapped Reads)")+
#scale_fill_manual(name="Zone",values=zone_coloring)+
theme(legend.position='top')+
facet_grid(Type + Temperature ~ Site, scales="free")

ggsave("/mnt/smart/users/fcoutinho/EAIMM/ARB_Quantification/ARB_Quantification_DotPlots.pdf",plot=figA,width=8,height=5,pointsize=8)

bl_ctrl_df<-full_df[which(full_df$Treatment=="Control" & full_df$Site=="BL2"),c("CFU_per_ml","Temperature","DAY")]
colnames(bl_ctrl_df)[1]<-"Ctrl_CFU_per_ml"
summary(bl_ctrl_df)


bla_tmp_df<-full_df[which(full_df$Site=="BL2"),]

bla_tmp_df<-merge(bla_tmp_df,bl_ctrl_df,by=c("DAY","Temperature"),all.x=TRUE)

bla_tmp_df$ARB_Ratio<-bla_tmp_df$CFU_per_ml/bla_tmp_df$Ctrl_CFU_per_ml

summary(bla_tmp_df)

figA<-ggplot(bla_tmp_df,aes(x=DAY,y=CFU_per_ml,fill=Treatment))+
#geom_point(alpha=0.9,aes(fill=Treatment,shape=Replicate))+
geom_bar(position="dodge",stat="identity",alpha=0.9)+
#geom_boxplot(alpha=0.9,aes(fill=Treatment))+
#geom_bar(alpha=0.9,shape=23,aes(fill=Treatment))+
theme_bw()+
#scale_y_continuous(limits = c(0, 0.025))+
#labs(x = "Abundance of Archaea (RPKM)",y="Cas1 Abundance (% Mapped Reads)")+
#scale_fill_manual(name="Zone",values=zone_coloring)+
theme(legend.position='top')+
facet_grid( Treatment ~ Temperature, scales="free_y")

ggsave("/mnt/smart/users/fcoutinho/EAIMM/ARB_Quantification/ARB_Quantification_Blanes_Temp_Absolute.pdf",plot=figA,width=5,height=9,pointsize=8)

figA<-ggplot(bla_tmp_df,aes(x=DAY,y=ARB_Ratio,fill=Treatment))+
#geom_point(alpha=0.9,aes(fill=Treatment,shape=Replicate))+
geom_bar(position="dodge",stat="identity",alpha=0.9)+
#geom_boxplot(alpha=0.9,aes(fill=Treatment))+
#geom_bar(alpha=0.9,shape=23,aes(fill=Treatment))+
theme_bw()+
#scale_y_continuous(limits = c(0, 0.025))+
#labs(x = "Abundance of Archaea (RPKM)",y="Cas1 Abundance (% Mapped Reads)")+
#scale_fill_manual(name="Zone",values=zone_coloring)+
theme(legend.position='top')+
facet_grid( Treatment ~ Temperature, scales="free_y")

ggsave("/mnt/smart/users/fcoutinho/EAIMM/ARB_Quantification/ARB_Quantification_Blanes_Temp_Ratio.pdf",plot=figA,width=5,height=9,pointsize=8)

