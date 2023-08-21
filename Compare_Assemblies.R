##Assembly Info
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
