library(ranger)
library(scales)
library(Metrics)
#library(caret)
options(expressions = 5e5)

args = commandArgs(trailingOnly=TRUE)
response_file<-args[1]
response_variables<-strsplit(args[2],",")[[1]]
predictor_file<-args[3]
max_rmse<-as.numeric(args[4])
cor_cutoff<-as.numeric(args[5])
tree_num<-as.numeric(args[6]) #Number of trees to grow
var_tries<-as.numeric(args[7]) #Variables to try at each split when growing trees
node_size<-as.numeric(args[8]) #Node size for the probabilistic forest
threads<-as.numeric(args[9])

load_defaults<-FALSE
if(load_defaults == TRUE) {
	response_file<-"/mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/StG23/Metadata/Marine_Viral_Communities_Sample_Metadata_Updated_with_Guidi.tsv"#"/mnt/lustre/scratch/fcoutinho/StG_23/Metadata/TARA_Salazar_19_MGs_Metadata.tsv"
	response_variables<-"Mean.Flux.at.150m"#strsplit("Temperature,Oxygen,ChlorophyllA,Salinity,Iron.5m,Ammonium.5m",",")[[1]]
	predictor_file<-"/mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/StG23/Viruses/Abundance/StG23_Viruses_Min_Prev_50_Raw_Abundance.tsv"#"/mnt/lustre/scratch/fcoutinho/StG_23/Prokaryotes/Abundance/RPKM_Abundance_OceanDNA_All_Species_Rep_MAGs_by_Species_Cluster.tsv"
	max_rmse<-100000
	cor_cutoff<-0.5
	tree_num<-5000
	var_tries<-500
	node_size<-1
	threads<-47
}

#Load and filter predictor data (i.e. the sampletaxon abundance table)
raw_predictor_df<-read.table(file=predictor_file,sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE,check.names=FALSE,row.names=1)

z_predictor_df<-raw_predictor_df

t_z_predictor_df<-as.data.frame(t(z_predictor_df))

predictor_variables<-colnames(t_z_predictor_df)
#summary(z_predictor_df)

#Load and filter response data (i.e. the sample metdata table)
raw_response_df<-read.table(file=response_file,sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE,row.names=1,check.names=FALSE)

z_response_df<-subset(raw_response_df,select=response_variables)

fulldata<-merge(z_response_df,t_z_predictor_df,by="row.names",all.x=TRUE)

colnames(fulldata)[1]<-"MG_UID"
fulldata$MG_UID<-as.factor(fulldata$MG_UID)
rownames(fulldata)<-fulldata$MG_UID
#summary(fulldata[,1:15])

all_rf_stats<-c()
best_rf_stats<-c()
all_preds<-rbind()
all_best_stats<-c()
passed_importance_vals<-c()
best_rfs<-list()

###Build RF Models for each response variable
for (resp_var in response_variables) {
	#best_prod<-0
	best_rmse<-+Inf
	#best_pred<-NA
	best_RF<-NA
	#best_tmatrix<-NA
	#best_vmatrix<-NA
	all_cort<-c()
	#all_corv<-c()
	all_rmset<-c()
	#all_rmsev<-c()
	pass_cort<-c()
	#pass_corv<-c()
	pass_rmset<-c()
	#pass_rmsev<-c()
	best_stats<-c()
	
	subdata<-subset(fulldata,select=(c(resp_var,predictor_variables)))
	subdata<-subdata[complete.cases(subdata),]
	resp_df<-subdata[[resp_var]]
	pred_df<-subdata[,c(predictor_variables)]
	
	print(paste("Building RF model for Response variable:",resp_var,sep=" "))
	print("Subset DF dimensions:")
	print(dim(subdata))
	#print("Response DF dimensions:")
	#print(dim(resp_df))
	#print("Predictor DF dimensions:")
	#print(dim(pred_df))

	#trainIndex <- createDataPartition(z_sub_df$Z_Response_Variable, p = 0.8, list = FALSE)
	#z_train_df<-z_sub_df[trainIndex, ]
	#z_val_df<-z_sub_df[-trainIndex, ]

	set.seed(3435)	
	rf_model<-ranger(y=resp_df, x=pred_df, num.trees=tree_num, mtry=var_tries, write.forest=TRUE,probability=FALSE,importance="impurity",min.node.size=node_size,num.threads=threads)

	cor_t<-cor.test(rf_model$predictions,resp_df,method="pearson")

	#Calculate Root Mean Square Error
	rmse_t<-rmse(rf_model$predictions,resp_df)

	all_cort<-c(all_cort,cor_t$estimate)
	#all_corv<-c(all_corv,cor_v$estimate)

	all_rmset<-c(all_rmset,rmse_t)
		#all_rmsev<-c(all_rmsev,rmse_v)
		
	if ((cor_t$estimate < cor_cutoff) || (rmse_t > max_rmse)) {passed<-FALSE } else {passed<-TRUE}
		
	all_rf_stats<-rbind(all_rf_stats,c(resp_var,1,cor_t$estimate,rmse_t,passed))
		
	#Check if correlation and RMSE values are within allowed cutoff
	if (passed == FALSE) {
		pass_cort<-c(pass_cort,NA)
		#pass_corv<-c(pass_corv,NA)
		pass_rmset<-c(pass_rmset,NA)
		#pass_rmsev<-c(pass_rmsev,NA)
	} else {
		pass_cort<-c(pass_cort,cor_t$estimate)
		#pass_corv<-c(pass_corv,cor_v$estimate)
		pass_rmset<-c(pass_rmset,rmse_t)
		#pass_rmsev<-c(pass_rmsev,rmse_v)
	}

	#Keep track of the impotance values assigned to predictors for all the passed networks
	passed_importance<-as.data.frame(cbind(names(rf_model$variable.importance),rf_model$variable.importance))
	colnames(passed_importance)<-c("Predictor","importance")
	passed_importance$importance<-as.numeric(passed_importance$importance)
	passed_importance$Response<-resp_var
	#passed_importance$RMSEV<-rmse_t
	#passed_importance$CorV<-cor_t$estimate
	#passed_importance$Method<-"Impurity"
	#Convert absolute importance to percentages: Note that impurity method only yields positive importance values
	passed_importance$Relative_Importance<-(passed_importance$importance/sum(passed_importance$importance))*100
	passed_importance_vals<-rbind(passed_importance_vals,passed_importance)

	#Store best matrix and best ANN based ons RMSE
	if ((rmse_t <= best_rmse) & (cor_t$estimate >= cor_cutoff)) {
		best_rmse<-rmse_t
		best_RF<-rf_model
		best_stats<-passed_importance
	}
		
	all_best_stats<-rbind(all_best_stats,best_stats)
	if (is.na(best_RF)[1] == FALSE) {	
		best_rfs[[resp_var]]<-best_RF
		pred_data<-as.data.frame(cbind(rf_model$predictions,resp_df))
		colnames(pred_data)<-c("Predicted","Measured")
		pred_data$Response_Var<-resp_var
		pred_data$MG_UID<-as.factor(rownames(subdata))
		pred_data$Error<-sqrt((pred_data$Predicted-pred_data$Measured)**2)
		all_preds<-rbind(all_preds,subset(pred_data,select=c("Measured","Predicted","Response_Var","MG_UID")))
	}
}

colnames(all_rf_stats)<-c("Response","Iteration","CorV","RMSEv","Passed")
write.table(all_rf_stats,file="Reverse_Ecology_All_RF_Stats.tsv",sep="\t",append=FALSE,row.names=FALSE,col.names=TRUE,quote=FALSE)

write.table(passed_importance_vals,file="Reverse_Ecology_Best_RF_Importance.tsv",sep="\t",append=FALSE,row.names=FALSE,col.names=TRUE,quote=FALSE)

write.table(all_preds,file="Reverse_Ecology_Best_RF_Predictions.tsv",sep="\t",append=FALSE,row.names=FALSE,col.names=TRUE,quote=FALSE)

save.image("Reverse_Ecology_RF.RData")

quit(status=1)