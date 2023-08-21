library(dplyr)
library(tidyr)
library(reshape2)

args = commandArgs(trailingOnly=TRUE)
scaff_file<-args[1]
cds_file<-args[2]
out_file<-args[3]

###Updated scaffold Info
scaff_data<-read.table(file=scaff_file,sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)
rownames(scaff_data)<-scaff_data$Sequence

#Add non redudant KO, Pathway and Metabolism info to scaffold info table
cds_info<-read.table(file=cds_file,sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)
colnames(cds_info)[1]<-"CDS"
cds_info$Sequence<-as.factor(gsub("_(\\d)+$","",cds_info$CDS,perl=TRUE))

categories<-c("KEGG_Best_Subject_Modules","KEGG_Best_Subject_Pathways","KEGG_Best_Subject")

for (categ in categories){
	cds_info$Category<-cds_info[[categ]]
	f_cds_info<-cds_info[!is.na(cds_info$Category),]
	f_cds_info<-f_cds_info[which(f_cds_info$Category != "set()"),]
	f_cds_info<-f_cds_info %>% separate_longer_delim(Category, delim = "', '")
	f_cds_info$Category<-as.factor(gsub("(\\{)|(\\})|'","",f_cds_info$Category,perl=TRUE))
	scaff_data[[categ]]<-NA
	for (sid in unique(scaff_data$Sequence)) {
		subdata<-f_cds_info[which(f_cds_info$Sequence == sid),]
		if(nrow(subdata) >0) {
            unique_vals_conc<-as.character(paste(as.vector(sort(unique(subdata$Category[!is.na(subdata$Category)]))),collapse="|"))
			scaff_data[sid,categ]<-unique_vals_conc
			#summary(subdata)
		}
	}
	scaff_data[[categ]]<-as.factor(scaff_data[[categ]])
}

#summary(scaff_data)
write.table(scaff_data,file=out_file,sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)
