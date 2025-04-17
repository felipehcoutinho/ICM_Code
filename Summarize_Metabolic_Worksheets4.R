library(reshape2)
###Metabolic module completeness x Sequence calculation
full_df<-rbind()

for (i in c(1:20)) {
    file_name<-paste("/mnt/smart/users/fcoutinho/Blanes_Vir_Nestor/Metabolic_Output/Metabolic_Outputs_Batch_",i,"/METABOLIC_result_each_spreadsheet/METABOLIC_result_worksheet4.tsv",sep="")
    
    func_data<-read.table(file=file_name,sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)
    dim(func_data)

    long_func_data<-melt(func_data,id=c("Module.step","Module","KO.id","Module.Category"))

    colnames(long_func_data)<-c("Module_Step","Module_Name","KO_IDs","Module_Category","Sequence","Module_Step_Presence")

    long_func_data$Sequence<-gsub(".Module.step.presence$","",long_func_data$Sequence,perl=TRUE)

    long_func_data$Sequence<-as.factor(long_func_data$Sequence)
    long_func_data$Module_Step_Presence<-as.factor(long_func_data$Module_Step_Presence)

    #summary(long_func_data)

    f_long_func_data<-long_func_data[which(long_func_data$Module_Step_Presence == "Present"),]

    #summary(f_long_func_data)
    full_df<-rbind(full_df,f_long_func_data)
    #summary(full_df)
}

summary(full_df)

write.table(full_df,file="/mnt/smart/users/fcoutinho/Blanes_Vir_Nestor/Metabolic_Output/Blanes_Viruses_Metabolic_Summary.tsv",sep="\t",append=FALSE,row.names=FALSE,col.names=TRUE,quote=FALSE)
