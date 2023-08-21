library(dplyr)
library(tidyr)
library(tibble)

args = commandArgs(trailingOnly=TRUE)
response_file<-args[1]
min_prev<-as.numeric(args[2])
output_file<-args[3]

raw_response_df<-read.table(file=response_file,sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE,row.names=1,check.names=FALSE)

#For each row in raw_response_df, calculate the number of samples with a value > 0
raw_response_df<-raw_response_df %>% mutate(non_zero_count=rowSums(.[2:ncol(raw_response_df)]!=0))

f_raw_response_df<-raw_response_df[which(raw_response_df$non_zero_count >= min_prev),]

f_raw_response_df$non_zero_count<-NULL

#summary(f_raw_response_df)

write.table(f_raw_response_df,file=output_file,sep="\t",append=FALSE,row.names=TRUE,col.names=NA,quote=FALSE)
