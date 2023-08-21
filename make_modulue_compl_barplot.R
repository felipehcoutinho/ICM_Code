###Input options
args = commandArgs(trailingOnly=TRUE)
input_file<-args[1]
tax_level<-args[2]#
tax_name<-args[3]#
group_level<-args[4]#
min_mean_mond_comp<-args[5]
mw_test<-args[6]

use_presets<-FALSE
#use_presets<-TRUE
if (use_presets == TRUE) {
	input_file<-"/mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/Metabolic_Preds/Profiles_Malaspina_MAGs_METABOLIC_Module_Completeness_From_worksheet4.tsv"
	tax_level<-"Phylum"#"Domain"
	tax_name<-"Thermoplasmatota"#"Archaea"
	group_level<-"Zone"#"Phylum"
	min_mean_mond_comp<-0
	mw_test<-TRUE
}

###Libraries
library(ggplot2)
library(RColorBrewer)
library(reshape2)
library(stringr)
library(ggsignif)

#Fixed aesthetics
bar_width_val<-0.01

sub_taxon_coloring<-NA

if (group_level == "Zone") {
	subtaxon_coloring<-c(brewer.pal(9,"YlGnBu")[c(1,4,9)])
	names(subtaxon_coloring)<-c("Epipelagic","Mesopelagic","Bathypelagic")
}

###Taxon x Metabolic Module completenes plots
func_data<-read.table(file=input_file,sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=FALSE)

f_func_data<-func_data[which(func_data[[tax_level]] == tax_name),]

f_func_data[[group_level]][which(f_func_data[[group_level]] == "")]<-"Unclassified"

f_func_data$Group_Var<-f_func_data[[group_level]]

if (group_level == "Zone") {
	f_func_data$Group_Var<-factor(f_func_data$Group_Var,levels=c("Epipelagic","Mesopelagic","Bathypelagic"))
}

if (mw_test == TRUE) {
	f_func_data$Module<-paste(f_func_data$Module," (",f_func_data$Total_Steps,")",sep="")
	
	mag_count_data<-rbind()
	unique_taxa_list<-unique(f_func_data[[group_level]])
	for (tax in unique_taxa_list) {
		sub_data<-f_func_data[which(f_func_data[[group_level]] == tax),]
		magc<-length(unique(sub_data$MAG))
		mag_count_data<-rbind(mag_count_data,c(tax,magc))
	}
	mag_count_data<-as.data.frame(mag_count_data)
	colnames(mag_count_data)<-c("Group_Var","MAG_Count")
	
	f_func_data<-merge(f_func_data,mag_count_data,by="Group_Var",all.x=TRUE)
	
	f_func_data$Group_Var_w_Count<-paste(f_func_data$Group_Var," (",f_func_data$MAG_Count,")",sep="")
	
	mx_func_data<-as.data.frame(aggregate(f_func_data$Module_Completeness, by=list(Taxon=f_func_data[[group_level]],Module_UID=f_func_data$Module_UID,Module=f_func_data$Module,Module_Category=f_func_data$Module_Category), FUN=max))
	colnames(mx_func_data)[5]<-"Max_Module_Completeness"
	
	non_zero_mods<-unique(mx_func_data[which(mx_func_data$Max_Module_Completeness > 0),]$Module)
	
	sub_f_func_data<-f_func_data[which(f_func_data$Module %in% non_zero_mods),]
	
	sub_f_func_data$Module<-str_wrap(sub_f_func_data$Module,width=35)
	
	comps_list<-split(t(combn(as.vector(unique(sub_f_func_data$Group_Var)), 2)), seq(nrow(t(combn(as.vector(unique(sub_f_func_data$Group_Var)), 2)))))
	
	box_fig<-ggplot(data=sub_f_func_data,aes(x=Group_Var,y=Module_Completeness,fill=Group_Var))+geom_boxplot(position="dodge2",alpha=0.9,colour="Black",lwd=0.5)+geom_signif(comparisons = list(c("Epipelagic","Mesopelagic"),c("Mesopelagic","Bathypelagic"),c("Epipelagic","Bathypelagic")),map_signif_level = TRUE, textsize=2.5,test = "wilcox.test",step_increase=0.05,tip_length = 0.01,size=0.1,lwd=0.1)+scale_fill_manual(name=group_level,values=subtaxon_coloring)+theme_bw()+geom_hline(alpha=0.25,colour="black",size=0.3,linetype=2,yintercept=100)+ylab("Module Completeness")+xlab(group_level)+theme(axis.text.x = element_text(angle = 45,hjust = 1,size=6),legend.position="none")+facet_wrap(Module_Category ~ Module ,scales="free_x")+theme(strip.text.x = element_text(size=5,margin = margin(0,0,0,0, "cm")))#
	fig_name<-paste("MLP",tax_name,"by",group_level,"MW_Test","Metab_Comp_Boxplot.pdf",sep="_")
	ggsave(fig_name,plot=box_fig,device=cairo_pdf,width=30,height=45,pointsize=8)
	
	
	q("no")
}

f_func_data$Module<-paste(f_func_data$Module," (",f_func_data$Total_Steps,")",sep="")
#summary(f_func_data)

m_func_data<-as.data.frame(aggregate(f_func_data$Module_Completeness, by=list(Taxon=f_func_data[[group_level]],Module_UID=f_func_data$Module_UID,Module=f_func_data$Module,Module_Category=f_func_data$Module_Category), FUN=mean))

colnames(m_func_data)[5]<-"Mean_Module_Completeness"

mx_func_data<-as.data.frame(aggregate(f_func_data$Module_Completeness, by=list(Taxon=f_func_data[[group_level]],Module_UID=f_func_data$Module_UID,Module=f_func_data$Module,Module_Category=f_func_data$Module_Category), FUN=max))

colnames(mx_func_data)[5]<-"Max_Module_Completeness"

s_func_data<-as.data.frame(aggregate(f_func_data$Module_Completeness, by=list(Taxon=f_func_data[[group_level]],Module_UID=f_func_data$Module_UID,Module=f_func_data$Module,Module_Category=f_func_data$Module_Category), FUN=sd))

colnames(s_func_data)[5]<-"SD_Module_Completeness"
s_func_data$SD_Module_Completeness[is.na(s_func_data$SD_Module_Completeness)]<-0

mag_count_data<-rbind()
unique_taxa_list<-unique(f_func_data[[group_level]])
for (tax in unique_taxa_list) {
	sub_data<-f_func_data[which(f_func_data[[group_level]] == tax),]
	magc<-length(unique(sub_data$MAG))
	mag_count_data<-rbind(mag_count_data,c(tax,magc))
}
mag_count_data<-as.data.frame(mag_count_data)
colnames(mag_count_data)<-c("Taxon","MAG_Count")
#summary(mag_count_data)

mean_sd_func_data<-merge(m_func_data,s_func_data,by=c("Taxon","Module_UID","Module","Module_Category"),all.x=TRUE)
mean_sd_func_data<-merge(mean_sd_func_data,mx_func_data,by=c("Taxon","Module_UID","Module","Module_Category"),all.x=TRUE)
mean_sd_func_data<-merge(mean_sd_func_data,mag_count_data,by=c("Taxon"),all.x=TRUE)

#mean_sd_func_data$Taxon<-paste(mean_sd_func_data$Taxon," (",mean_sd_func_data$MAG_Count,")",sep="")
mean_sd_func_data$Taxon_MAG_Count<-paste(mean_sd_func_data$Taxon," (",mean_sd_func_data$MAG_Count,")",sep="")

mean_sd_func_data$Upper<-mean_sd_func_data$Mean_Module_Completeness + mean_sd_func_data$SD_Module_Completeness
mean_sd_func_data$Lower<-mean_sd_func_data$Mean_Module_Completeness - mean_sd_func_data$SD_Module_Completeness

#dim(mean_sd_func_data)
#summary(mean_sd_func_data)

f_mean_sd_func_data<-mean_sd_func_data[which(mean_sd_func_data$Mean_Module_Completeness > min_mean_mond_comp),]
#dim(f_mean_sd_func_data)
#summary(f_mean_sd_func_data)

f_mean_sd_func_data$Module_Category<-str_wrap(f_mean_sd_func_data$Module_Category,width=30)
f_mean_sd_func_data$Module<-str_wrap(f_mean_sd_func_data$Module,width=60)

subtaxon_coloring<-c(brewer.pal(8,"Accent"),brewer.pal(12,"Paired"),brewer.pal(8,"Set1"))[c(1:length(unique(f_mean_sd_func_data$Taxon)))]
names(subtaxon_coloring)<-unique(f_mean_sd_func_data$Taxon)


fig4<-ggplot(f_mean_sd_func_data,aes(x=Module,fill=Taxon,y=Mean_Module_Completeness,group=Taxon))+geom_bar(position=position_dodge2(width = bar_width_val, preserve = "single"),stat="identity",alpha=0.9,size=0.05,colour="Black")+geom_errorbar(position=position_dodge2(width = bar_width_val, preserve = "single"),stat = "identity",aes(ymin=Lower,ymax=Upper),na.rm=TRUE,size=0.03)+geom_point(aes(x=Module,fill=Taxon,y=Max_Module_Completeness,group=Taxon),position=position_dodge2(width=1,preserve="single"),shape=23,size=1)+geom_hline(alpha=0.8,colour="black",size=0.1,linetype=2,yintercept=100)+scale_fill_manual(name=group_level,values=subtaxon_coloring)+theme_bw()+theme(axis.text.x = element_text(angle = 75,hjust = 1,size=4),legend.position="top")+ylab(paste("% Mean ± SD Module Completeness ( Max value among all MAGs from",group_level,")",sep=""))+facet_grid(. ~ Module_Category,scales="free_x",space="free")+theme(strip.text.x.top = element_text(angle = 90))

fig_name<-paste("MLP",tax_name,"by",group_level,"Metab_Comp_Barplot.pdf",sep="_")

ideal_height<-min(40,max(c(length(unique(f_mean_sd_func_data$Module))/6.6,9)))
ideal_width<-min(45,max(c(length(unique(f_mean_sd_func_data$Taxon))*11,9)))

print(paste("Ideal width:",ideal_width,sep=" "))

ggsave(fig_name,plot=fig4,device=cairo_pdf,width=45,height=7,pointsize=8)

print(paste("Output file generated:",fig_name,sep=" "))

q("no")

#facet options: +facet_wrap(. ~ Module_Category,scales="free_x",ncol=6) ; +facet_grid(. ~ Module_Category,scales="free_x",space="free")
#geom_errorbar non aes options: ,size=0.1,width=0.1
#group=interaction(Taxon,Module),x=Module,
#geom_bar position options: position_dodge2(width = bar_width_val, preserve = "single") , position_dodge(width=bar_width_val) , position="dodge"
#geom_errorbar position options: position_dodge2(width = bar_width_val, padding = 0.5), position_dodge(width=bar_width_val)

if (group_level == "Phylum") {
	valid_taxa<-as.vector(unique(f_mean_sd_func_data$Taxon))
	print(paste("Setting colours for taxon: ",valid_taxa,sep=""))
	
	full_taxon_coloring<-c(brewer.pal(9,"Set1"),brewer.pal(4,"Pastel1"),brewer.pal(9,"Greys")[6],brewer.pal(5,"Set2")[5],brewer.pal(8,"Accent"),brewer.pal(6,"Dark2"))
	names(full_taxon_coloring)<-c("Actinobacteriota","Alphaproteobacteria","Cyanobacteria","Crenarchaeota","Bacteroidota","Nanoarchaeota","Verrucomicrobiota","Marinisomatota","Myxococcota","Patescibacteria","Gammaproteobacteria","Acidobacteriota","UBP7","Poribacteria","Eremiobacterota","Margulisbacteria","Gemmatimonadota","Desulfobacterota","Nitrospinota","Desulfobacterota_D","Bdellovibrionota","Hydrothermarchaeota","Chlamydiota","SAR324","Halobacteriota","Chloroflexota","Thermoplasmatota","Planctomycetota","Thermoproteota")
	
	full_taxon_order<-c("Acidobacteriota","Actinobacteriota","Bacteroidota","Bdellovibrionota","Chlamydiota","Chloroflexota","Cyanobacteria","Desulfobacterota","Desulfobacterota_D","Eremiobacterota","Gemmatimonadota","Margulisbacteria","Marinisomatota","Myxococcota","Nitrospinota","Patescibacteria","Planctomycetota","Poribacteria","Alphaproteobacteria","Gammaproteobacteria","SAR324","UBP7","Verrucomicrobiota","Halobacteriota","Hydrothermarchaeota","Nanoarchaeota","Thermoplasmatota","Thermoproteota")
	
	subtaxon_coloring<-full_taxon_coloring[which(full_taxon_coloring %in% valid_taxa)]
	
	print(paste("Using coloring: ",subtaxon_coloring,sep=""))
	
	sub_taxon_order<-full_taxon_order[which(full_taxon_order %in% unique(f_mean_sd_func_data$Taxon))]
	
	f_mean_sd_func_data$Taxon<-factor(f_mean_sd_func_data$Taxon,levels=sub_taxon_order)
	
	print(paste("Using taxon order: ",paste(sub_taxon_order,sep=""),sep=""))
	
	sub_taxon_order_df<-as.data.frame(cbind(sub_taxon_order))
	colnames(sub_taxon_order_df)[1]<-"Taxon"
	
	sub_taxon_mag_count_df<-subset(f_mean_sd_func_data[which(!duplicated(f_mean_sd_func_data$Taxon)),],select=c("Taxon","MAG_Count"))
	
	valid_taxa_mag_count_df<-merge(sub_taxon_order_df,sub_taxon_mag_count_df,by="Taxon",all.x=TRUE)
	valid_taxa_mag_count_unique_labels<-paste(valid_taxa_mag_count_df$Taxon," (",valid_taxa_mag_count_df$MAG_Count,")",sep="")
	
	fig4<-ggplot(f_mean_sd_func_data,aes(x=Module,fill=Taxon,y=Mean_Module_Completeness,group=Taxon))+geom_bar(position=position_dodge2(width = bar_width_val, preserve = "single"),stat="identity",alpha=0.9,size=0.05)+geom_errorbar(position=position_dodge2(width = bar_width_val, preserve = "single"),stat = "identity",aes(ymin=Lower,ymax=Upper),na.rm=TRUE,size=0.03)+geom_point(aes(x=Module,fill=Taxon,y=Max_Module_Completeness,group=Taxon),position=position_dodge2(width=1,preserve="single"),shape=23,size=1)+geom_hline(alpha=0.8,colour="black",size=0.1,linetype=2,yintercept=100)+scale_fill_manual(name=group_level,values=subtaxon_coloring,labels=valid_taxa_mag_count_unique_labels)+theme_bw()+theme(axis.text.x = element_text(angle = 75,hjust = 1,size=4),legend.position="top")+ylab(paste("% Mean ± SD Module Completeness (Dots = Max value among all MAGs from ",group_level,")",sep=""))+facet_grid(. ~ Module_Category,scales="free_x",space="free")+theme(strip.text.x.top = element_text(angle = 90))
}