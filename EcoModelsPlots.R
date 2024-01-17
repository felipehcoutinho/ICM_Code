library(ggpubr)
library(reshape2)
library(ggplot2)
library(RColorBrewer)

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

phylum_custom_selection<-c("Actinobacteriota","Alphaproteobacteria","Cyanobacteria","Bacteroidota","Verrucomicrobiota","Marinisomatota","Gammaproteobacteria","Nitrospinota","SAR324","Chloroflexota","Thermoplasmatota","Planctomycetota","Thermoproteota")

env_var_coloring<-c(brewer.pal(11,"RdYlBu")[c(2,10)],brewer.pal(9,"BuPu")[7],brewer.pal(11,"RdYlGn")[c(8,10)],brewer.pal(8,"Dark2")[c(6,7,8)])
#names(env_var_coloring)<-c("Temperature","Salinity","Oxygen","Chlorophyll_A","NPP.8d.VGPM..mgC.m2.day.","Ammonium_5m","Mean.Flux.at.150m","Iron_5m")
names(env_var_coloring)<-c("Temperature","Salinity","Oxygen","ChlorophyllA","NPP.8d.VGPM..mgC.m2.day.","Ammonium.5m","Mean.Flux.at.150m","Iron.5m")

###Predictor Taxon x Response Env Var (from RF) Importance plots
#Euk files
metadata_file<-"/mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/StG23/Metadata/TARA_MGs_Paired_With_Delmont_Data_Info+Metadata.tsv"
#Virus files
mes_vir_abd_file<-"/mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/StG23/Viruses/Abundance/RPKM_Abundance_NA.tsv"
importance_file<-"/mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/ICCOMM/VirRF_Predictor_Importance.tsv"
pred_info_file<-"/mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/StG23/Viruses/Info_Sequences/Full_Info_IMGVR_Scaffolds_Comp_75-100_Max_Conta_0_vOTU_Representatives.tsv"
group_var<-"Host_Phylum_PHIST"#"Host_Phylum_IMGVR"#
prefix<-"RF_Vir"

#Prok files
importance_file<-"/mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/ICCOMM/Prokaryotes/Niche_Models/ProkZANN_Predictor_SHAP_Importance.tsv"#/mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/ICCOMM/TestProkZANN_Predictor_SHAP_Relative_Importance.tsv"#"/mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/ICCOMM/ProkZANN_Predictor_SHAP_Relative_Importance.tsv"#"/mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/ICCOMM/ProkNoHptRF_Predictor_Importance.tsv"
mes_prok_abd_file<-"/mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/StG23/Prokaryotes/Abundance/Transposed_RPKM_Abundance_OceanDNA_All_Species_Rep_MAGs_by_MAG_ID.tsv"
pred_prok_abd_file<-"/mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/ICCOMM/Prokaryotes/Niche_Models/ProkZANNTraining_Models_Predicted_Values.tsv"
pred_info_file<-"/mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/Databases/OceanDNA/Sub_OceanDNA_MAGs_Info_Species_Rep_only.tsv"
metadata_file<-"/mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/StG23/Metadata/TARA_MGs_Paired_With_Salazar_Data_Info+Metadata_Updated_with_Guidi.tsv"
group_var<-"Phylum"
prefix<-"ANN_Prok"

###Calc abundance stats summary
abd_df<-read.table(file=pred_prok_abd_file,sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE,check.names=FALSE,row.names=1)

outname<-"/mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/ICCOMM/Prokaryotes/Niche_Models/Stats_ProkZANNTraining_Models_Predicted_Values.tsv"

calc_abund_stats<-function(abd_df=NA,outfile=NA) {
    abd_df<-abd_df[,which((colnames(abd_df) != "Dataset") & (colnames(abd_df) != "Data_Type"))]
    #calculate how many non zero values occur in each column of abd_df
    col_nz_count<-apply(abd_df,2,function(x) length(which(x > 0)))
    #calculate median values in each column of abd_df
    col_medians<-apply(abd_df,2,function(x) median(x))
    col_sd<-apply(abd_df,2,function(x) sd(x))
    col_means<-colMeans(abd_df,na.rm=TRUE)
    stats_df<-cbind(colnames(abd_df),col_nz_count,col_means,col_medians,col_sd)
    colnames(stats_df)<-c("UID","Prevalence","Mean_Abundance","Median_Abundance","SD_Abundance")
    if (!is.na(outfile)) {
        write.table(stats_df,file=outfile,sep="\t",quote=FALSE,row.names=FALSE)
    }
    return(as.data.frame(stats_df))
}

stats_df<-calc_abund_stats(abd_df=abd_df,outfile=outname)

###Parsing importance data and making SHAP importance violin plots
imp_df<-read.table(file=importance_file,sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE,check.names=FALSE)
colnames(imp_df)[1]<-"UID"
summary(imp_df)

#Divide the values in each row of imp_df by the max value of each row in abs_imp_df
abs_imp_df<-imp_df
abs_imp_df[,2:7]<-abs(abs_imp_df[,2:7])
for (i in 1:nrow(imp_df)) {
    imp_df[i,2:7]<-imp_df[i,2:7]/max(abs_imp_df[i,2:7])
}
summary(imp_df)


pred_info_df<-read.table(file=pred_info_file,sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE,check.names=FALSE)
colnames(pred_info_df)[1]<-"UID"

metadata_df<-read.table(file=metadata_file,sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE,check.names=FALSE)
summary(metadata_df)

full_df<-melt(imp_df,id=c("UID","Response_Variable"),variable.name="Predictor",value.name="Importance")
full_df<-merge(full_df,subset(metadata_df,select=c(Group, Layer_from_Table_W4, Ocean.region)),by.x="UID",by.y="Group",all.x=TRUE,suffixes=c("_Importance","_Raw_Value"))
full_df<-merge(full_df,subset(pred_info_df,select=c("UID","Domain","Phylum","Class","Order","Family","Genus","Species")),by.x="Response_Variable",by.y="UID",all.x=TRUE)

#head(full_df)
#m_full_df<-melt(full_df,id=c("Response_Variable","Domain","UID","Phylum","Class","Order","Family","Genus","Species","Layer_from_Table_W4","Ocean.region"),variable.name="Predictor",value.name="Importance")

###Shap Importance violin plots
make_shap_imp_violinplots<-function(data_df=NA,valid_groups=c(),first_group_var=NA,second_group_var=NA) {
    
    if (length(valid_groups) > 0) {
        data_df<-data_df[which(data_df[[first_group_var]] %in% valid_groups),]
    }

    if (is.na(first_group_var)) {
        data_df$Group<-"All"
    } else {
        #rename first_group_var column to Group
        colnames(data_df)[colnames(data_df) == first_group_var]<-"Group_1"
        colnames(data_df)[colnames(data_df) == second_group_var]<-"Group_2"
    }

    #summary(data_df)

    figX<-ggplot(data=data_df)+
    geom_violin(aes(x=Predictor, y=Importance, fill=Predictor),alpha=0.9)+
    #scale_y_continuous(limits=c(-100,100))+
    theme_bw()+
    theme(legend.position="top", axis.text.x = element_text(angle = 90,hjust = 1,size=8))+ #, strip.text.y.right = element_text(angle = 0, size =9),legend.title=element_text(size=9)
    scale_fill_manual(name="Predictor",values=env_var_coloring)+
    facet_wrap(. ~ Group_1, ncol=4)

    outname<-paste(prefix,"_Importance_by_",first_group_var,"_All_Samples_Violinlots.pdf",sep="")
    ggsave(outname,plot=figX,device=cairo_pdf,width=12,height=12,pointsize=8)

    if (!is.na(second_group_var)) {
        figY<-ggplot(data=data_df)+
        geom_violin(aes(x=Predictor, y=Importance, fill=Predictor),alpha=0.9)+
        #scale_y_continuous(limits=c(-100,100))+
        theme_bw()+
        theme(legend.position="top", axis.text.x = element_text(angle = 90,hjust = 1,size=8))+ #, strip.text.y.right = element_text(angle = 0, size =9),legend.title=element_text(size=9)
        scale_fill_manual(name="Predictor",values=env_var_coloring)+
        facet_grid(Group_2 ~ Group_1)

        outname<-paste(prefix,"_Importance_by_",first_group_var,"_and_",second_group_var,"_Violinlots.pdf",sep="")
        ggsave(outname,plot=figY,device=cairo_pdf,width=22,height=22,pointsize=8)
    }
}

make_shap_imp_violinplots(data_df=full_df,valid_groups=phylum_custom_selection,first_group_var="Phylum",second_group_var="Ocean.region")


###Current vs Projected abundance plots for single OTU
# resp_mes_file<-"/mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/StG23/Prokaryotes/Abundance/Transposed_RPKM_Abundance_OceanDNA_All_Species_Rep_MAGs_by_MAG_ID.tsv"

# resp_mes_df<-read.table(file=resp_mes_file,sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)

# m_resp_mes<-melt(resp_mes_df,id=c("Sample_UID"),variable.name="TaxonID",value.name="Abundance")

# m_resp_mes_df$Dataset<-"Training"
# m_resp_mes_df$Dat_Type<-"Measured"

resp_train_pred_file<-"/mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/ICCOMM/ProkZANNTraining_Models_Predicted_Values.tsv"

resp_train_pred_df<-read.table(file=resp_train_pred_file,sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)

colnames(resp_train_pred_df)[1]<-"Sample_UID"

m_resp_train_pred<-melt(resp_train_pred_df,id=c("Sample_UID","Data_Type","Dataset"),variable.name="TaxonID",value.name="Abundance")

resp_test_pred_file<-"/mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/ICCOMM/ProkZANNTesting_Models_Predicted_Values.tsv"

resp_test_pred_df<-read.table(file=resp_test_pred_file,sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)

colnames(resp_test_pred_df)[1]<-"Sample_UID"

m_resp_test_pred<-melt(resp_test_pred_df,id=c("Sample_UID","Data_Type","Dataset"),variable.name="TaxonID",value.name="Abundance")

#Best Cyano predictor of NPP: "OceanDNA.b16232"
#Best Cyano predictor of Chrlorphyll: "OceanDNA.b16700"
tid<-"OceanDNA.b16700"

m_resp_pred<-merge(m_resp_train_pred[,c("Sample_UID","TaxonID","Abundance")],m_resp_test_pred[,c("Sample_UID","TaxonID","Abundance")],suffixes=c("_Train","_Test"),by=c("Sample_UID","TaxonID"),all=TRUE)

sm_resp_pred<-m_resp_pred[which(m_resp_pred$TaxonID==tid),]

summary(sm_resp_pred)

#replace negative values in Abundance_Train and Abundance_Test with 0
sm_resp_pred$Abundance_Train[which(sm_resp_pred$Abundance_Train < 0)]<-0
sm_resp_pred$Abundance_Test[which(sm_resp_pred$Abundance_Test < 0)]<-0
summary(sm_resp_pred)

metadata<-read.table(file="/mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/StG23/Metadata/TARA_MGs_Paired_With_Salazar_Data_Info+Metadata_Updated_with_Guidi.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE,check.names=FALSE)

sm_resp_pred<-merge(sm_resp_pred,metadata[,c("Sample_UID","Latitude_from_Table_W4","Longitude_from_Table_W4","Layer_from_Table_W1","Depth.nominal","Ocean.region")],by.x="Sample_UID",all.x=TRUE)

colnames(sm_resp_pred)[colnames(sm_resp_pred) == "Latitude_from_Table_W4"]<-"Latitude"
colnames(sm_resp_pred)[colnames(sm_resp_pred) == "Longitude_from_Table_W4"]<-"Longitude"
colnames(sm_resp_pred)[colnames(sm_resp_pred) == "Layer_from_Table_W1"]<-"Layer"
colnames(sm_resp_pred)[colnames(sm_resp_pred) == "Depth.nominal"]<-"Depth"
colnames(sm_resp_pred)[colnames(sm_resp_pred) == "Ocean.region"]<-"Ocean_Region"

sm_resp_pred$Layer<-factor(sm_resp_pred$Layer,levels=c("SRF","DCM","MIX","MES"))

figX<-ggplot(sm_resp_pred)+
geom_point(aes(x=Abundance_Train, y=Abundance_Test, fill=Ocean_Region), alpha=0.9, shape=23)+
geom_abline(intercept=0,slope=1,colour="black",linetype="dashed")+
theme_bw()+
xlab("Current Predicted Abundance")+
ylab("Projected Predicted Abundance")+
theme(legend.position="top",legend.text=element_text(size=4))+
facet_wrap(. ~ Layer, ncol=2, nrow=2)

ggsave("Single_OTU_Test_TrainedxTest_Scatter_Plot.pdf",plot=figX,device=cairo_pdf,width=10,height=12,pointsize=8)


figX<-ggplot(sm_resp_pred)+
geom_density(aes(Abundance_Train),alpha=0.9,fill="green")+
geom_density(aes(Abundance_Test),alpha=0.9,"blue")+
theme_bw()+
theme(legend.position="top")+
facet_wrap(. ~ Layer, ncol=2, nrow=2)

ggsave("Single_OTU_Test_TrainedxTest_Density_Plot.pdf",plot=figX,device=cairo_pdf,width=10,height=10,pointsize=8)



###Trained x Test Predicted Values for Community diversity
#Training set
library(vegan, lib.loc="/mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/Rlibs/")

train_set_preds_df<-read.table(file="ANNProkZANNTraining_Models_Predicted_Values.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE,check.names=FALSE)

colnames(train_set_preds_df)[1]<-"Sample_UID"

m_tr_set<-melt(train_set_preds_df,id=c("Sample_UID", "Data_Type","Dataset"),variable.name="TaxonID",value.name="Abundance")

#Replace all negative values in training set with 0
m_tr_set$Abundance[which(m_tr_set$Abundance < 0)]<-0

summary(m_tr_set)

train_set_preds_df<-dcast(data=m_tr_set, formula = Sample_UID ~ TaxonID, fun.aggregate = NULL, fill=0, value.var = "Abundance") 
rownames(train_set_preds_df)<-train_set_preds_df$Sample_UID
train_set_preds_df$Sample_UID<-NULL

tr_div<-as.data.frame(diversity(train_set_preds_df, index = "shannon"))
tr_div$Sample_UID<-as.factor(rownames(tr_div))
colnames(tr_div)[1]<-"Shanon_Diversity_Index"

#Testing set
test_set_preds_df<-read.table(file="ANNProkZANNTesting_Models_Predicted_Values.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE,check.names=FALSE)

colnames(test_set_preds_df)[1]<-"Sample_UID"

m_ts_set<-melt(test_set_preds_df,id=c("Sample_UID", "Data_Type","Dataset"),variable.name="TaxonID",value.name="Abundance")

#Replace all negative values in testing set with 0
m_ts_set$Abundance[which(m_ts_set$Abundance < 0)]<-0

summary(m_ts_set)

test_set_preds_df<-dcast(data=m_ts_set, formula = Sample_UID ~ TaxonID, fun.aggregate = NULL, fill=0, value.var = "Abundance") 
rownames(test_set_preds_df)<-test_set_preds_df$Sample_UID
test_set_preds_df$Sample_UID<-NULL

ts_div<-as.data.frame(diversity(test_set_preds_df, index = "shannon"))
ts_div$Sample_UID<-as.factor(rownames(ts_div))
colnames(ts_div)[1]<-"Shanon_Diversity_Index"
summary(ts_div)

all_div<-merge(tr_div,ts_div,by="Sample_UID",all.xy=TRUE,suffixes=c("_Train","_Test"))

all_div$Shannon_L2_Fold<-log2(all_div$Shanon_Diversity_Index_Test/all_div$Shanon_Diversity_Index_Train)

figX<-ggplot(data=all_div,aes(y=Shannon_L2_Fold, x=Sample_UID))+
geom_bar(position="dodge",stat="identity",colour="black",alpha=0.9)+ #,group=f_sum_imp_df[[group_var]]
#scale_fill_manual(name="Taxon",values=phylum_coloring)+
theme_bw()+
theme(legend.position="top", axis.text.x = element_text(angle = 45,hjust = 1,size=5))+ #, strip.text.y.right = element_text(angle = 0, size =9),legend.title=element_text(size=9)
#facet_grid(. ~ Response,scales="free_x")+
#guides(col = guide_legend(nrow = 1))+
xlab("Sample")+
ylab("Log2(Fold Change) of Shannon Diversity Index (Test/Train)")

brp_outname<-paste("ANN_Prok_Warming+3_Scenario_Sample_Shannon_Diversity_Barplots.pdf",sep="")
ggsave(brp_outname,plot=figX,device=cairo_pdf,width=15,height=7,pointsize=8)

#Map of log2 shannon fold
metadata<-read.table(file="/mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/StG23/Metadata/TARA_MGs_Paired_With_Salazar_Data_Info+Metadata_Updated_with_Guidi.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE,check.names=FALSE)

all_div<-merge(all_div,metadata[,c("Sample_UID","Latitude_from_Table_W4","Longitude_from_Table_W4","Layer_from_Table_W1","Depth.nominal")],by.x="Sample_UID",all.x=TRUE)

colnames(all_div)[5]<-"Latitude"
colnames(all_div)[6]<-"Longitude"
colnames(all_div)[7]<-"Layer"
colnames(all_div)[8]<-"Depth"

summary(all_div)

all_div$Group<-as.factor("Predictions")
div_col_pal<-brewer.pal(9,"RdBu")[c(1:9)]
div_col_grad<-rev(colorRampPalette(div_col_pal)(n=99))

library(maps, lib.loc="/mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/Rlibs/")
world_map <- map_data("world")

all_div$Layer<-factor(all_div$Layer,levels=c("SRF","DCM","MIX","MES"))

figX<-ggplot(world_map, aes(x = long, y = lat, group = group))+
geom_polygon(fill="lightgray", colour = "lightgray")+
coord_cartesian(xlim=c(-180,180),ylim=c(-90,90))+
geom_jitter(all_div,mapping=aes(y=Latitude,x=Longitude,fill=Shannon_L2_Fold, group=Group),size=3,colour="black",shape=23,height=1,width=0)+
theme_bw()+xlab("Longitude")+ylab("Latitude")+
scale_fill_gradientn(colours = div_col_grad,limits = c(-0.3, 0.3))+
facet_wrap(. ~ Layer, ncol=2, nrow=2)

ggsave("ANN_Prok_Warming+3_Scenario_Sample_Shannon_Diversity_Fold_Maps.pdf",plot=figX,device=cairo_pdf,width=13,height=8,pointsize=8)


###Vir imp x host x metab
pred_info_df<-read.table(file=pred_info_file,sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE,check.names=FALSE)
colnames(pred_info_df)[1]<-"UID"
summary(pred_info_df)

annot_info_file<-"/mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/StG23/Viruses/Annotation/Genome_Info.tsv"
annot_df<-read.table(file=annot_info_file,sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE,check.names=FALSE)
colnames(annot_df)[1]<-"UID"
summary(annot_df)

imp_df<-read.table(file=importance_file,sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE,check.names=FALSE)
colnames(imp_df)[1]<-"UID"
summary(imp_df)
m_imp_df<-melt(imp_df,id=c("UID"),variable.name="Response_Variable",value.name="Relative_Importance")

m_imp_df<-m_imp_df[which(m_imp_df$Relative_Importance > 0),]

m_imp_df<-merge(m_imp_df,pred_info_df[,c("UID","Host_Phylum_PHIST")],by="UID",all.x=TRUE)

m_imp_df<-merge(m_imp_df,annot_df,by="UID",all.x=TRUE)

write.table(m_imp_df,file="VirRF_Predictor_Importance_With_Annotation.tsv",sep="\t",quote=FALSE,row.names=FALSE)

