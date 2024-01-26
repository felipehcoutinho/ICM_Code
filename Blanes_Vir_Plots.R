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

#zone_coloring<-c(brewer.pal(9,"YlOrRd")[c(4)],brewer.pal(9,"YlGnBu")[c(5,9)])
zone_coloring<-c(brewer.pal(9,"YlGnBu")[c(1,4,9)])
names(zone_coloring)<-c("Epipelagic","Mesopelagic","Bathypelagic")

taxon_order<-c("Acidobacteriota","Actinobacteriota","Bacteroidota","Bdellovibrionota","Chlamydiota","Chloroflexota","Cyanobacteria","Desulfobacterota","Desulfobacterota_D","Eremiobacterota","Gemmatimonadota","Margulisbacteria","Marinisomatota","Myxococcota","Nitrospinota","Patescibacteria","Planctomycetota","Poribacteria","Alphaproteobacteria","Gammaproteobacteria","SAR324","UBP7","Verrucomicrobiota","Halobacteriota","Hydrothermarchaeota","Nanoarchaeota","Thermoplasmatota","Thermoproteota")

phylum_coloring<-c(brewer.pal(9,"Set1"),brewer.pal(4,"Pastel1"),brewer.pal(8,"Dark2")[8],brewer.pal(5,"Set2")[5],brewer.pal(8,"Accent"),brewer.pal(6,"Dark2"))
names(phylum_coloring)<-c("Actinobacteriota","Alphaproteobacteria","Cyanobacteria","Crenarchaeota","Bacteroidota","Nanoarchaeota","Verrucomicrobiota","Marinisomatota","Myxococcota","Patescibacteria","Gammaproteobacteria","Acidobacteriota","UBP7","Poribacteria","Eremiobacterota","Margulisbacteria","Gemmatimonadota","Desulfobacterota","Nitrospinota","Desulfobacterota_D","Bdellovibrionota","Hydrothermarchaeota","Chlamydiota","SAR324","Halobacteriota","Chloroflexota","Thermoplasmatota","Planctomycetota","Thermoproteota")

phylum_custom_selection<-c("Actinobacteriota","Alphaproteobacteria","Cyanobacteria","Bacteroidota","Verrucomicrobiota","Marinisomatota","Gammaproteobacteria","Nitrospinota","SAR324","Chloroflexota","Thermoplasmatota","Planctomycetota","Thermoproteota")

###Raw phist preds
all_phist_preds<-read.table(file="/mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/Blanes_Vir_Nestor/PHIST_Predictions/PHIST_Output/positive_predictions.csv",sep=",",header=TRUE,quote="",comment="",stringsAsFactors=TRUE,check.names=FALSE)

summary(all_phist_preds)

Host_ID_data<-read.table(file="/mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/Blanes_Vir_Nestor/Blanes_MAG_Info.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE,check.names=FALSE)
colnames(Host_ID_data)[1]<-"Host_ID"

summary(Host_ID_data)

all_phist_preds<-read.table(file="/mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/Blanes_Vir_Nestor/PHIST_Predictions/PHIST_Output/positive_predictions.csv",sep=",",header=TRUE,quote="",comment="",stringsAsFactors=TRUE,check.names=FALSE)

filtered_phist_preds<-filter_phist_preds(phist_df=all_phist_preds,host_df=Host_ID_data,max_p_value=2.384e-14)

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

full_vir_scaff_data<-merge(vir_scaff_data,filtered_phist_preds,by="Virus_ID",all.x=TRUE)

###Seq info
#gVOD
metadata_file<-"/mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/Databases/gVOD/gVOD_Raw_Data.tsv"
scenarios_pred_abd_file<-"/mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/ICCOMM/gVOD/gVOD_V+P_ANN_Testing_Models_Inverse_Transformed_Predicted_Values.tsv"
pred_abd_file<-"/mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/ICCOMM/gVOD/gVOD_V+P_ANN_Training_Models_Inverse_Transformed_Predicted_Values.tsv"
prefix<-"ANN_gVOD"#"RF_Vir"#
###Load datasets
#scena_df<-read.table(file=scenarios_file,sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE,check.names=FALSE)

###Importance data
library(ggpubr)
library(reshape2)
library(ggplot2)
library(RColorBrewer)

#Prokaryotes
imp_df<-read.table(file="/mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/ICCOMM/ProkNoHptRF_Predictor_Importance.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE,check.names=FALSE)

m_imp<-melt(imp_df,id=c("Predictor_Variable"),variable.name="Response_Variable",value.name="Importance")

m_imp<-m_imp[which(m_imp$Importance != 0),]

pred_info_df<-read.table(file="/mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/Databases/OceanDNA/Sub_OceanDNA_Host_IDs_Info_Species_Rep_only.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE,check.names=FALSE)
colnames(pred_info_df)[1]<-"Taxon_UID"

m_imp<-merge(m_imp,subset(pred_info_df,select=c("Taxon_UID","Phylum")),by.x="Predictor_Variable",by.y="Taxon_UID",all.x=TRUE)

summary(m_imp)

imp_sums<-aggregate(Importance ~ Response_Variable + Phylum, data=m_imp, FUN=sum)


phylum_custom_selection<-c("Actinobacteriota","Alphaproteobacteria","Cyanobacteria","Bacteroidota","Verrucomicrobiota","Marinisomatota","Gammaproteobacteria","Nitrospinota","SAR324","Chloroflexota","Thermoplasmatota","Planctomycetota","Thermoproteota")

f_imp_sums<-imp_sums[which(imp_sums$Phylum %in% phylum_custom_selection),]

phylum_coloring<-c(brewer.pal(9,"Set1"),brewer.pal(4,"Pastel1"),brewer.pal(8,"Dark2")[8],brewer.pal(5,"Set2")[5],brewer.pal(8,"Accent"),brewer.pal(6,"Dark2"))
names(phylum_coloring)<-c("Actinobacteriota","Alphaproteobacteria","Cyanobacteria","Crenarchaeota","Bacteroidota","Nanoarchaeota","Verrucomicrobiota","Marinisomatota","Myxococcota","Patescibacteria","Gammaproteobacteria","Acidobacteriota","UBP7","Poribacteria","Eremiobacterota","Margulisbacteria","Gemmatimonadota","Desulfobacterota","Nitrospinota","Desulfobacterota_D","Bdellovibrionota","Hydrothermarchaeota","Chlamydiota","SAR324","Halobacteriota","Chloroflexota","Thermoplasmatota","Planctomycetota","Thermoproteota")

figY<-ggplot(data=f_imp_sums)+
        geom_bar(aes(x=Response_Variable, y= Importance , fill=Phylum, group=Phylum),alpha=0.9,position="dodge",stat="identity",colour="black")+
        scale_fill_manual(name="hylum",values=phylum_coloring)+
        theme_bw()+
        theme(legend.position="top", axis.text.x = element_text(angle = 90,hjust = 1,size=8))

ggsave("Prok_RF_Models_Imp_Bar.pdf",plot=figY,device=cairo_pdf,width=12,height=5, pointsize=8) 

###Viruses
imp_df<-read.table(file="/mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/ICCOMM/VirRF_Predictor_Gini_Importance.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE,check.names=FALSE)

m_imp<-melt(imp_df,id=c("Predictor_Variable"),variable.name="Response_Variable",value.name="Importance")

m_imp<-m_imp[which(m_imp$Importance != 0),]

pred_info_df<-read.table(file="/mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/StG23/Viruses/Info_Sequences/Full_Info_IMGVR_Scaffolds_Comp_75-100_Max_Conta_0_vOTU_Representatives.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE,check.names=FALSE)
colnames(pred_info_df)[1]<-"Taxon_UID"

m_imp<-merge(m_imp,subset(pred_info_df,select=c("Taxon_UID","Host_Phylum_PHIST")),by.x="Predictor_Variable",by.y="Taxon_UID",all.x=TRUE)

summary(m_imp)

imp_sums<-aggregate(Importance ~ Response_Variable + Host_Phylum_PHIST, data=m_imp, FUN=sum)


phylum_custom_selection<-c("Actinobacteriota","Alphaproteobacteria","Cyanobacteria","Bacteroidota","Verrucomicrobiota","Marinisomatota","Gammaproteobacteria","Nitrospinota","SAR324","Chloroflexota","Thermoplasmatota","Planctomycetota","Thermoproteota")


f_imp_sums<-imp_sums[which(imp_sums$Host_Phylum_PHIST %in% phylum_custom_selection),]


phylum_coloring<-c(brewer.pal(9,"Set1"),brewer.pal(4,"Pastel1"),brewer.pal(8,"Dark2")[8],brewer.pal(5,"Set2")[5],brewer.pal(8,"Accent"),brewer.pal(6,"Dark2"))
names(phylum_coloring)<-c("Actinobacteriota","Alphaproteobacteria","Cyanobacteria","Crenarchaeota","Bacteroidota","Nanoarchaeota","Verrucomicrobiota","Marinisomatota","Myxococcota","Patescibacteria","Gammaproteobacteria","Acidobacteriota","UBP7","Poribacteria","Eremiobacterota","Margulisbacteria","Gemmatimonadota","Desulfobacterota","Nitrospinota","Desulfobacterota_D","Bdellovibrionota","Hydrothermarchaeota","Chlamydiota","SAR324","Halobacteriota","Chloroflexota","Thermoplasmatota","Planctomycetota","Thermoproteota")

figY<-ggplot(data=f_imp_sums)+
        geom_bar(aes(x=Response_Variable, y= Importance , fill=Host_Phylum_PHIST, group=Host_Phylum_PHIST),alpha=0.9,position="dodge",stat="identity",colour="black")+
        scale_fill_manual(name="Host Phylum",values=phylum_coloring)+
        theme_bw()+
        theme(legend.position="top", axis.text.x = element_text(angle = 90,hjust = 1,size=8))

ggsave("Vir_RF_Models_Imp_Bar.pdf",plot=figY,device=cairo_pdf,width=12,height=5, pointsize=8) 

#Load and edit metadata
metadata_df<-read.table(file=metadata_file,sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE,check.names=FALSE)

colnames(metadata_df)[1]<-"Sample_UID"

colnames(metadata_df)[colnames(metadata_df) == "Ocean.region"]<-"Region"
colnames(metadata_df)[colnames(metadata_df) == "Layer_from_Table_W4"]<-"Layer"

if ((prefix == "ANN_Vir") | (prefix == "RF_Vir")) {
    layer_ids<-data.frame(do.call('rbind', strsplit(as.character(metadata_df$Sample_ID),'_',fixed=TRUE)))
    colnames(layer_ids)<-c("Station","Layer")
    metadata_df$Layer<-as.character(layer_ids$Layer)
    metadata_df$Layer[which(metadata_df$Layer == "dDCM")] <- "DCM"
    metadata_df$Layer<-as.factor(metadata_df$Layer)
}

metadata_df$Region<-as.factor(gsub("\\(.*","",metadata_df$Region,perl=TRUE))
region_custom_selection<-gsub("\\(.*","",region_custom_selection,perl=TRUE)

#Load and edit variable info (predictor or response depending on model type)
pred_info_df<-read.table(file=pred_info_file,sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE,check.names=FALSE)
colnames(pred_info_df)[1]<-"Taxon_UID"

###Parsing importance data and making SHAP importance violin plots
#Function to make SHAP importance violin plots
#Shap Importance violin plots
make_shap_imp_violinplots<-function(data_df=NA,valid_groups=c(),first_group_var=NA,second_group_var=NA,fig_width=45,fig_height=20) {
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

    if (debug == TRUE) {
        print("Summary of data_df")
        print(summary(data_df))
    }

    figX<-ggplot(data=data_df)+
    geom_violin(aes(x=Predictor, y=Importance, fill=Predictor),alpha=0.9)+
    #scale_y_continuous(limits=c(-100,100))+
    theme_bw()+
    theme(legend.position="top", axis.text.x = element_text(angle = 90,hjust = 1,size=8))+ #, strip.text.y.right = element_text(angle = 0, size =9),legend.title=element_text(size=9)
    scale_fill_manual(name="Predictor",values=env_var_coloring)+
    facet_wrap(. ~ Group_1, ncol=4)

    outname<-paste(prefix,rel_imp_prefix,"_Importance_by_",first_group_var,"_All_Samples_Violinlots.pdf",sep="")
    ggsave(outname,plot=figX,device=cairo_pdf,width=12,height=12,pointsize=8)

    if (!is.na(second_group_var)) {
        figY<-ggplot(data=data_df)+
        geom_violin(aes(x=Predictor, y=Importance, fill=Predictor),alpha=0.9)+
        #scale_y_continuous(limits=c(-100,100))+
        theme_bw()+
        theme(legend.position="top", axis.text.x = element_text(angle = 90,hjust = 1,size=8))+ #, strip.text.y.right = element_text(angle = 0, size =9),legend.title=element_text(size=9)
        scale_fill_manual(name="Predictor",values=env_var_coloring)+
        facet_grid(Group_2 ~ Group_1)

        outname<-paste(prefix,rel_imp_prefix,"_Importance_by_",first_group_var,"_and_",second_group_var,"_Violinlots.pdf",sep="")
        ggsave(outname,plot=figY,device=cairo_pdf,width=fig_width,height=fig_height, pointsize=8) #
    }
}

#Function to converte importance values to relative importance
calc_rel_imp<-function(raw_imp_df=NA) {
    print("Calculating relative importance")
    #Indentify numeric columns in abs_imp_df
    num_cols<-which(sapply(raw_imp_df,is.numeric))
    abs_imp_df<-abs(raw_imp_df[,num_cols])
    abs_col_max<-apply(abs_imp_df,1,max)
    rel_imp_df<-raw_imp_df
    rel_imp_df[,num_cols]<-raw_imp_df[,num_cols]/abs_col_max
    if (debug == TRUE) {
        print("Summary of relative importance")
        print(summary(rel_imp_df))
    }
    return(rel_imp_df)
}

#Read in the raw importance data
imp_df<-read.table(file=importance_file,sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE,check.names=FALSE)
colnames(imp_df)[1]<-"Sample_UID"

#ANN
colnames(imp_df)[colnames(imp_df) == "Response_Variable"]<-"Taxon_UID"

if (debug == TRUE) {
    print("Summary of raw importance")
    summary(imp_df)
}

#TRansform to relativ eimportances if specified
transf_to_rel_imp<-TRUE
rel_imp_prefix<-""

if (transf_to_rel_imp == TRUE) {
    imp_df<-calc_rel_imp(raw_imp_df=imp_df)
    rel_imp_prefix<-"_Relative"
}

#melt the importance, predictor to full_df
#ANN
#
full_df<-melt(imp_df,id=c("Sample_UID","Taxon_UID"),variable.name="Predictor",value.name="Importance")
#RF
#
full_df<-melt(imp_df,id=c("Sample_UID","Response_Variable"),variable.name="Taxon_UID",value.name="Importance")

#Build Taxon_Alias column in pred_info
# primary_group_var<-"Genus"
# secondary_group_var<-"Phylum"#
#primary_group_var<-"Host_Phylum_PHIST"#"Phylum"
primary_group_var<-"Phylum"
secondary_group_var<-NA

if (is.na(secondary_group_var)) {
    pred_info_df$Taxon_Alias<-pred_info_df[[primary_group_var]]
    #Now add the minimal subset of pred_info data to full_df
    full_df<-merge(full_df,subset(pred_info_df,select=c("Taxon_UID",primary_group_var,"Taxon_Alias")),by="Taxon_UID",all.x=TRUE)
} else {
    pred_info_df$Taxon_Alias<-as.factor(paste(pred_info_df[[secondary_group_var]],pred_info_df[[primary_group_var]],sep="|"))
    #Now add the minimal subset of pred_info data to full_df
    full_df<-merge(full_df,subset(pred_info_df,select=c("Taxon_UID",primary_group_var,secondary_group_var,"Taxon_Alias")),by="Taxon_UID",all.x=TRUE)
}

#Read in the sample metadata
if (debug == TRUE) {
    print("Summary of metadata")
    summary(metadata_df)
}

#Now add the minimal subset ofsample data to full_df
full_df<-merge(full_df,subset(metadata_df,select=c("Sample_UID", "Layer", "Region")),by="Sample_UID",all.x=TRUE,suffixes=c("_Importance","_Raw_Value"))

#call plot functions
#Prok
full_df$Layer<-factor(full_df$Layer,levels=c("SRF","DCM","MIX","MES"))
#Vir
full_df$Layer<-factor(full_df$Layer,levels=c("SRF","DCM","MES"))

###ANN
######Prok
make_shap_imp_violinplots(data_df=full_df[which(full_df$Taxon_Alias %in% phylum_custom_selection & full_df$Region %in% region_custom_selection),],first_group_var="Taxon_Alias",second_group_var="Region",fig_width=45,fig_height=15)

make_shap_imp_violinplots(data_df=full_df[which(full_df$Taxon_Alias %in% phylum_custom_selection & full_df$Region %in% region_custom_selection),],first_group_var="Taxon_Alias",second_group_var="Layer",fig_width=45,fig_height=15)

#make_shap_imp_violinplots(data_df=full_df[which(full_df$Genus %in% genus_custom_selection & full_df$Ocean.region %in% region_custom_selection),],first_group_var="Genus",second_group_var="Ocean.region",fig_width=45,fig_height=15)

#make_shap_imp_violinplots(data_df=full_df,valid_groups=phylum_custom_selection,first_group_var="Phylum",second_group_var="Ocean.region")
#make_shap_imp_violinplots(data_df=full_df[which(full_df$Genus %in% genus_custom_selection & full_df$Ocean.region %in% region_custom_selection),],first_group_var="Taxon_Alias",second_group_var="Ocean.region",fig_width=45,fig_height=15)

#make_shap_imp_violinplots(data_df=full_df[which(full_df$Genus %in% genus_custom_selection & full_df$Ocean.region %in% region_custom_selection),],first_group_var="Genus",second_group_var="Ocean.region",fig_width=45,fig_height=15)

#make_shap_imp_violinplots(data_df=full_df,valid_groups=phylum_custom_selection,first_group_var="Phylum",second_group_var="Ocean.region")

######Vir
make_shap_imp_violinplots(data_df=full_df[which(full_df$Host_Phylum_PHIST %in% phylum_custom_selection),],first_group_var="Taxon_Alias",second_group_var="Layer",fig_width=45,fig_height=10)

make_shap_imp_violinplots(data_df=full_df[which(full_df$Host_Phylum_PHIST %in% phylum_custom_selection & full_df$Region %in% region_custom_selection),],first_group_var="Taxon_Alias",second_group_var="Region",fig_width=45,fig_height=10)


###Perfo barplot
p_ann_perf_df<-read.table(file="/mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/ICCOMM/Prokaryotes/Niche_Models/ProkZ_ANN_Training_Models_Performance_Info.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE,check.names=FALSE)
p_ann_perf_df$Model_Type<-"ANN"
p_ann_perf_df$Response_Type<-"Prokaryotes"


v_ann_perf_df<-read.table(file="/mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/ICCOMM/Viruses/Niche_Models/VirZ_ANN_Training_Models_Performance_Info.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE,check.names=FALSE)
v_ann_perf_df$Model_Type<-"ANN"
v_ann_perf_df$Response_Type<-"Viruses"

full_perf_df<-as.data.frame(rbind(p_ann_perf_df,v_ann_perf_df))

figY<-ggplot(data=full_perf_df[which(full_perf_df$Pearson_R2 >= 0.5),])+
        geom_boxplot(aes(x=Response_Type, y= R2_Score , fill=Response_Type),alpha=0.9)+
        theme_bw()+
        theme(legend.position="top") #, , axis.text.x = element_text(angle = 90,hjust = 1,size=8), strip.text.y.right = element_text(angle = 0, size =9),legend.title=element_text(size=9)

ggsave("ANN_Models_Performance_Box.pdf",plot=figY,device=cairo_pdf,width=7,height=5, pointsize=8) #

#REF PERFO
p_rf_perf_df<-read.table(file="/mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/ICCOMM/Prokaryotes/RF/ProkZ_RF_Training_Models_Performance_Info.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE,check.names=FALSE)
p_rf_perf_df$Model_Type<-"RF"
p_rf_perf_df$Response_Type<-"Prokaryotes"

v_rf_perf_df<-read.table(file="/mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/ICCOMM/Viruses/RF/VirZ_RF_Training_Models_Performance_Info.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE,check.names=FALSE)
v_rf_perf_df$Model_Type<-"ANN"
v_rf_perf_df$Response_Type<-"Viruses"

full_perf_df<-as.data.frame(rbind(p_rf_perf_df,v_rf_perf_df))

full_perf_df$Response_Variable[which(full_perf_df$Response_Variable == "ChlorophyllA")]<-"Chlorophyll_A"
full_perf_df$Response_Variable[which(full_perf_df$Response_Variable == "Iron.5m")]<-"Iron_5m"
full_perf_df$Response_Variable[which(full_perf_df$Response_Variable == "Ammonium.5m")]<-"Ammonium_5m"

figY<-ggplot(data=full_perf_df)+
        geom_bar(aes(x=Response_Variable, y= R2_Score , fill=Response_Type),alpha=0.9,position="dodge",stat="identity",colour="black")+
        theme_bw()+
        theme(legend.position="top", axis.text.x = element_text(angle = 90,hjust = 1,size=8))

ggsave("RF_Models_Performance_Bar.pdf",plot=figY,device=cairo_pdf,width=7,height=5, pointsize=8) 

#RF
#data_df=full_df[which(full_df$Phylum %in% phylum_custom_selection & full_df$Ocean.region %in% region_custom_selection & full_df$Importance != 0),]
data_df=full_df[which(full_df$Importance != 0),]

 figY<-ggplot(data=data_df)+
        geom_boxplot(aes(fill=Phylum, y=Importance, x=Phylum))+
        #scale_y_continuous(limits=c(-100,100))+
        theme_bw()+
        theme(legend.position="top", axis.text.x = element_text(angle = 90,hjust = 1,size=8))+ #, strip.text.y.right = element_text(angle = 0, size =9),legend.title=element_text(size=9)
        scale_fill_manual(name="Phylum",values=phylum_coloring)+
         facet_grid(Response_Variable ~ .)
        #facet_grid(Response_Variable ~  Ocean.region)
        #facet_grid(Response_Variable Layer_from_Table_W4 ~  Ocean.region)

        #outname<-paste(prefix,rel_imp_prefix,"_Importance_by_",first_group_var,"_and_",second_group_var,"_Violinlots.pdf",sep="")
        ggsave("Top100_Prok_RF_SHAP_Rel_Imp.pdf",plot=figY,device=cairo_pdf,width=10,height=20, pointsize=8) #



###Generate Scenarios Data
metadata_df_s1<-metadata_df
metadata_df_s1$Scenario<-"Warming+3"
metadata_df_s1$Temperature<-metadata_df_s1$Temperature+3

metadata_df_s2<-metadata_df
metadata_df_s2$Scenario<-"Freshning+5"
metadata_df_s2$Salinity<-metadata_df_s2$Salinity+5

metadata_df_s3<-metadata_df
metadata_df_s3$Scenario<-"Warming+3 & Freshning+5"
metadata_df_s3$Salinity<-metadata_df_s3$Salinity+5
metadata_df_s3$Temperature<-metadata_df_s3$Temperature+3

metadata_df$Scenario<-"Current"

scenarios_metadata_df<-as.data.frame(rbind(metadata_df,metadata_df_s1,metadata_df_s2,metadata_df_s3))

#scenarios_metadata_df$Group<-paste(scenarios_metadata_df$Sample_UID,scenarios_metadata_df$Scenario,c(1:nrow(scenarios_metadata_df)),sep="|")
scenarios_metadata_df$No<-paste(scenarios_metadata_df$No,scenarios_metadata_df$Scenario,c(1:nrow(scenarios_metadata_df)),sep="|")

write.table(scenarios_metadata_df,file="Scenarios.tsv",sep="\t",quote=FALSE,row.names=FALSE)

###Trained x Test Predicted Values for Community diversity
#Testing set
test_set_preds_df<-read.table(file=scenarios_pred_abd_file,sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE,check.names=FALSE)

colnames(test_set_preds_df)[1]<-"Sample_Scenario_UID"

sample_uids<-data.frame(do.call('rbind', strsplit(as.character(test_set_preds_df$Sample_Scenario_UID),'|',fixed=TRUE)))

colnames(sample_uids)<-c("Sample_UID","Scenario","Sample_Posit")

sample_uids$Sample_UID<-as.factor(sample_uids$Sample_UID)
sample_uids$Scenario<-as.factor(sample_uids$Scenario)
sample_uids$Sample_Posit<-as.factor(sample_uids$Sample_Posit)

test_set_preds_df<-as.data.frame(cbind(sample_uids,test_set_preds_df))

m_ts_set<-melt(test_set_preds_df,id=c("Sample_Scenario_UID","Scenario","Sample_Posit","Sample_UID","Data_Type","Dataset"),variable.name="TaxonID",value.name="Abundance")

#Replace all negative values in testing set with 0
m_ts_set$Abundance[which(m_ts_set$Abundance < 0)]<-0

summary(m_ts_set)

#Training set
train_set_preds_df<-read.table(file=pred_abd_file,sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE,check.names=FALSE)

colnames(train_set_preds_df)[1]<-"Sample_UID"

m_tr_set<-melt(train_set_preds_df,id=c("Sample_UID","Data_Type","Dataset"),variable.name="TaxonID",value.name="Abundance")

#Replace all negative values in training set with 0
m_tr_set$Abundance[which(m_tr_set$Abundance < 0)]<-0

summary(m_tr_set)


###Trained x Test Predicted Values for response variables FC

#Merge trainign and testing set predictions
m_full_set<-merge(subset(m_tr_set,select=c("Sample_UID","TaxonID","Abundance")),subset(m_ts_set,select=c("Sample_Scenario_UID","Scenario","Sample_Posit","Sample_UID","TaxonID","Abundance")),by=c("Sample_UID","TaxonID"),all.x=TRUE,suffixes=c("_Train","_Test"))

dim(m_full_set)

m_full_set$Fold_Change<-((m_full_set$Abundance_Test+1)/(m_full_set$Abundance_Train+1))

#m_full_set<-m_full_set[1:5000,]

#m_full_set<-merge(m_full_set,metadata_df[,c("Sample_UID","Latitude_from_Table_W4","Longitude_from_Table_W4","Layer_from_Table_W1","Depth.nominal","Ocean.region")],by="Sample_UID",all.x=TRUE,suffixes=c("_Abd","_Metadata"))

m_full_set<-merge(m_full_set,metadata_df,by="Sample_UID",all.x=TRUE,suffixes=c("_Abd","_Metadata"))

dim(m_full_set)
#summary(m_full_set)

#gVOD
m_sub_set<-m_full_set[which(m_full_set$Scenario != "Current" & m_full_set[["Depth water [m]"]] <= 200),]

dim(m_sub_set)
summary(m_sub_set)

figX<-ggplot(world_map, aes(x = long, y = lat, group = group))+
geom_polygon(fill="lightgray", colour = "lightgray")+
coord_cartesian(xlim=c(-180,180),ylim=c(-90,90))+
geom_point(data=m_sub_set,mapping=aes(y=Latitude,x=Longitude,fill=log10(Fold_Change),group=NULL),size=3,colour="black",shape=23)+
theme_bw()+xlab("Longitude")+ylab("Latitude")+
scale_fill_gradientn(colours = div_col_grad,limits = c(-8, 8))+ #
facet_grid(TaxonID ~ Scenario)

ggsave("gVOD_Scenarios_FC.pdf",plot=figX,device=cairo_pdf,width=25,height=12,pointsize=8)

figX<-ggplot(world_map, aes(x = long, y = lat, group = group))+
geom_polygon(fill="lightgray", colour = "lightgray")+
coord_cartesian(xlim=c(-180,180),ylim=c(-90,90))+
geom_point(data=m_sub_set,mapping=aes(y=Latitude,x=Longitude,group=NULL),size=3,colour="black",shape=23)+
theme_bw()+xlab("Longitude")+ylab("Latitude")

ggsave("gVOD_Samples_Coords.pdf",plot=figX,device=cairo_pdf,width=25,height=12,pointsize=8)

write.table(m_sub_set[,c(1:10,75)],file="Scenario_Preds_processed.tsv",sep="\t",quote=FALSE,row.names=FALSE)

#m_sub_set$Group_3<-factor(m_sub_set$Group_3,levels=c("SRF","DCM","MIX","MES"))


primary_group_var<-"Phylum"
secondary_group_var<-"Scenario"
tertiary_group_var<-"Layer_from_Table_W1"

m_full_set<-merge(m_full_set,subset(pred_info_df,select=c("Taxon_UID",primary_group_var)),by.x="TaxonID",by.y="Taxon_UID",all.x=TRUE)

dim(m_full_set)

colnames(m_full_set)[colnames(m_full_set) == primary_group_var]<-"Group_1"
colnames(m_full_set)[colnames(m_full_set) == secondary_group_var]<-"Group_2"
colnames(m_full_set)[colnames(m_full_set) == tertiary_group_var]<-"Group_3"

m_full_set$Fold_Change<-(m_full_set$Abundance_Test+1)/(m_full_set$Abundance_Train+1)

#m_sub_set<-m_full_set[which(m_full_set$Group_1 %in% phylum_custom_selection & m_full_set$Group_3 %in% region_custom_selection),]

m_sub_set<-m_full_set[which(m_full_set$Group_1 %in% phylum_custom_selection & m_full_set$Group_2 != "Current"),]

#m_sub_set<-m_full_set[which(m_full_set$Group_3 %in% region_custom_selection),]

#Prok
figX<-ggplot(data=m_sub_set)+
geom_violin(aes(x=Group_2, y=log10(Fold_Change), fill=Group_2),alpha=0.9)+
#scale_y_continuous(limits=c(-100,100))+
theme_bw()+
theme(legend.position="top", axis.text.x = element_text(angle = 90,hjust = 1,size=8))+
facet_grid(Group_1 ~ Group_3)# #, strip.text.y.right = element_text(angle = 0, size =9),legend.title=element_text(size=9)
#scale_fill_manual(name=primary_group_var,values=phylum_coloring)
outname<-paste(prefix,"FC_by",primary_group_var,secondary_group_var,tertiary_group_var,"All_Samples_Violinlots.pdf",sep="_")

ggsave(outname,plot=figX,device=cairo_pdf,width=15,height=20,pointsize=8)


#Training set
library(vegan)

train_set_preds_df<-dcast(data=m_tr_set, formula = Sample_UID ~ TaxonID, fun.aggregate = NULL, fill=0, value.var = "Abundance") 
rownames(train_set_preds_df)<-train_set_preds_df$Sample_UID
train_set_preds_df$Sample_UID<-NULL

tr_div<-as.data.frame(diversity(train_set_preds_df, index = "shannon"))
tr_div$Sample_UID<-as.factor(rownames(tr_div))
colnames(tr_div)[1]<-"Shanon_Diversity_Index"

#Testing set
summary(m_ts_set)

test_set_preds_df<-dcast(data=m_ts_set, formula = Sample_Scenario_UID ~ TaxonID, fun.aggregate = NULL, fill=0, value.var = "Abundance") 
rownames(test_set_preds_df)<-test_set_preds_df$Sample_Scenario_UID
test_set_preds_df$Sample_Scenario_UID<-NULL

ts_div<-as.data.frame(diversity(test_set_preds_df, index = "shannon"))
ts_div$Sample_Scenario_UID<-as.factor(rownames(ts_div))
colnames(ts_div)[1]<-"Shanon_Diversity_Index"
summary(ts_div)

sample_uids<-data.frame(do.call('rbind', strsplit(as.character(ts_div$Sample_Scenario_UID),'|',fixed=TRUE)))
colnames(sample_uids)<-c("Sample_UID","Scenario","Sample_Posit")
sample_uids$Sample_UID<-as.factor(sample_uids$Sample_UID)
sample_uids$Scenario<-as.factor(sample_uids$Scenario)
sample_uids$Sample_Posit<-as.factor(sample_uids$Sample_Posit)

ts_div<-as.data.frame(cbind(sample_uids,ts_div))

all_div<-merge(tr_div,ts_div,by="Sample_UID",all.xy=TRUE,suffixes=c("_Train","_Test"))

all_div$Shannon_L2_Fold<-log2(all_div$Shanon_Diversity_Index_Test/all_div$Shanon_Diversity_Index_Train)


#Map of log2 shannon fold
#Prok
#all_div<-merge(all_div,metadata_df[,c("Sample_UID","Latitude_from_Table_W4","Longitude_from_Table_W4","Layer_from_Table_W1","Depth.nominal")],by.x="Sample_UID",all.x=TRUE)
colnames(all_div)[8]<-"Latitude"
colnames(all_div)[9]<-"Longitude"
colnames(all_div)[10]<-"Layer"
colnames(all_div)[11]<-"Depth"

#Vir
all_div<-merge(all_div,metadata_df,by="Sample_UID",all.x=TRUE)

summary(all_div)

all_div$Group<-as.factor("Predictions")
div_col_pal<-brewer.pal(9,"RdBu")[c(1:9)]
div_col_grad<-rev(colorRampPalette(div_col_pal)(n=99))

library(maps)
world_map <- map_data("world")

all_div<-all_div[which(all_div$Scenario != "Current" & all_div$Layer != "ZZZ"),]
all_div$Layer<-factor(all_div$Layer,levels=c("SRF","DCM","MIX","MES"))

figX<-ggplot(world_map, aes(x = long, y = lat, group = group))+
geom_polygon(fill="lightgray", colour = "lightgray")+
coord_cartesian(xlim=c(-180,180),ylim=c(-90,90))+
geom_jitter(all_div,mapping=aes(y=Latitude,x=Longitude,fill=Shannon_L2_Fold, group=Group),size=3,colour="black",shape=23,height=1,width=0)+
theme_bw()+xlab("Longitude")+ylab("Latitude")+
scale_fill_gradientn(colours = div_col_grad,limits = c(-0.6, 0.6))+
facet_grid(Layer ~  Scenario)

ggsave("Scenario_Sample_Shannon_Diversity_Fold_Maps.pdf",plot=figX,device=cairo_pdf,width=25,height=20,pointsize=8)


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



###Plot model performance x response var stats
#Read in the performance data
perf_df<-read.table(file=performance_file,sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE,check.names=FALSE)

#Read in the measured values stats
stats_df<-read.table(file=resp_stats_file,sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE,check.names=FALSE)

full_df<-merge(perf_df,stats_df,by.x="Response_Variable",by.y="UID",all.x=TRUE)

summary(full_df)

cor(full_df[,c(2,3,4,6,7)])
perfo_coloring<-rev(colorRampPalette(brewer.pal(11,"Spectral")[c(1:5,7:11)])(n=99))

figY<-ggplot(data=full_df)+
geom_point(aes(y=Mean_Abundance, x= Prevalence , colour=RMSE),shape=18,size=3)+
scale_y_log10()+
theme_bw()+
scale_colour_gradientn(name="RMSE",colours=perfo_coloring)

outname<-paste(prefix,"_PerfomancexMean_AbdxPrev_Scatterplot.pdf",sep="")

ggsave(outname,plot=figY,device=cairo_pdf,width=11,height=10, pointsize=8)

###Calc abundance stats summary
abd_df<-read.table(file=mes_abd_file,sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE,check.names=FALSE,row.names=1)

outname<-"/mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/ICCOMM/Prokaryotes/Niche_Models/Stats_Prok_Measured_Values.tsv"

calc_abund_stats<-function(abd_df=NA,outfile=NA) {
    abd_df<-abd_df[,which((colnames(abd_df) != "Dataset") & (colnames(abd_df) != "Data_Type") & (colnames(abd_df) != "Sample_UID"))]
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


###Calc abundance sums by group
calc_group_sums<-function(abd_df=NA,info_df=NA,first_group_var=NA) {
    m_abd_df<-melt(abd_df,id=c("Sample_UID"),variable.name="Taxon_UID",value.name="Abundance")

    m_abd_df<-merge(m_abd_df,info_df[,c("Taxon_UID",first_group_var)],by="Taxon_UID",all.x=TRUE)

    colnames(m_abd_df)[colnames(m_abd_df) == first_group_var]<-"Group"

    #Replace NA values in the Group column by "Unclassified"
    m_abd_df$Group<-as.character(m_abd_df$Group)
    m_abd_df$Group[which(m_abd_df$Group == "NA")] <- "Unclassified"
    m_abd_df$Group[which(m_abd_df$Group == "")] <- "Unclassified"
    m_abd_df$Group<-as.factor(m_abd_df$Group)

    if (debug == TRUE) {
        print("Summary of m_abd_df")
        print(summary(m_abd_df))
    }

    group_abd_df<-dcast(m_abd_df,Sample_UID ~ Group,value.var="Abundance",fun.aggregate=sum,fill=0)

    return(group_abd_df)
}

abd_df<-read.table(file=mes_prok_abd_file,sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE,check.names=FALSE)
colnames(abd_df)[1]<-"Sample_UID"

pred_info_df<-read.table(file=prok_pred_info_file,sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE,check.names=FALSE)
colnames(pred_info_df)[1]<-"Taxon_UID"

group_var<-"Genus"
group_abd_df<-calc_group_sums(abd_df=abd_df,info_df=pred_info_df,first_group_var=group_var)

group_stats_df<-calc_abund_stats(abd_df=group_abd_df,outfile="/mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/ICCOMM/Prokaryotes/Niche_Models/Stats_Genus_Abundance_Measured_Values.tsv")


###Current vs Projected abundance plots for single OTU
# resp_mes_file<-"/mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/StG23/Prokaryotes/Abundance/Transposed_RPKM_Abundance_OceanDNA_All_Species_Rep_Host_IDs_by_Host_ID_ID.tsv"

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

