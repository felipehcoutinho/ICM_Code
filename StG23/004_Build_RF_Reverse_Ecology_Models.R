library(ranger)
library(scales)
library(Metrics)
library(caret)
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
train_iters<-as.numeric(args[10])

load_defaults<-FALSE
if(load_defaults == TRUE) {
	response_file<-"/mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/StG23/Metadata/Marine_Viral_Communities_Sample_Metadata_Updated_with_Guidi.tsv"#"/mnt/lustre/scratch/fcoutinho/StG_23/Metadata/TARA_Salazar_19_MGs_Metadata.tsv"
	response_variables<-c("Mean.Flux.at.150m","NPP.8d.VGPM..mgC.m2.day.")#strsplit("Temperature,Oxygen,ChlorophyllA,Salinity,Iron.5m,Ammonium.5m",",")[[1]]
	predictor_file<-"/mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/StG23/Viruses/Abundance/StG23_Viruses_Min_Prev_50_Raw_Abundance.tsv"#"/mnt/lustre/scratch/fcoutinho/StG_23/Prokaryotes/Abundance/RPKM_Abundance_OceanDNA_All_Species_Rep_MAGs_by_Species_Cluster.tsv"
	max_rmse<-100000
	cor_cutoff<-0.5
	tree_num<-5000
	var_tries<-500
	node_size<-1
	threads<-47
	train_iters<-3
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
best_rfs<-list()
best_importance_vals<-rbind()

###Build RF Models for each response variable
for (resp_var in response_variables) {
	best_rmse<-+Inf
	best_rf_passed<-FALSE
	subdata<-subset(fulldata,select=(c(resp_var,predictor_variables)))
	subdata<-subdata[complete.cases(subdata),]

	full_resp_df<-subdata[[resp_var]]
	full_pred_df<-subdata[,c(predictor_variables)]
	
	print(paste("Building RF model for Response variable:",resp_var,sep=" "))
	print("Subset DF dimensions:")
	print(dim(subdata))

	#Create a data partition and separate the predictors and response df into training and validation subsets
	trainIndex <- createDataPartition(subdata[[resp_var]], p = 0.7, list = FALSE)
	train_resp_df<-subdata[trainIndex, resp_var]
	train_pred_df<-subdata[trainIndex, c(predictor_variables)]

	val_resp_df<-subdata[-trainIndex, resp_var]
	val_pred_df<-subdata[-trainIndex, c(predictor_variables) ]

	set.seed(3435)
	#Perform multiple training iterations
	for (ticount in c(1:train_iters)) {
		print(paste("Training iteration:",ticount,sep=" "))
		#Train a  model using the training data
		train_rf_model<-ranger(y=train_resp_df, x=train_pred_df, num.trees=tree_num, mtry=var_tries, write.forest=TRUE,probability=FALSE,importance="none",min.node.size=node_size,num.threads=threads)

		#Evaluate model performance by calculating perason correlation and RMSE for the taining set
		cor_t<-cor.test(train_rf_model$predictions,train_resp_df,method="pearson")

		rmse_t<-rmse(train_rf_model$predictions,train_resp_df)

		#Obtain predictions for the validation set with the predict function
		val_preds<-predict(train_rf_model, data = val_pred_df)

		#Evaluate model performance by calculating perason correlation and RMSE for the validation set
		cor_v<-cor.test(val_preds$predictions,val_resp_df,method="pearson")

		rmse_v<-rmse(val_preds$predictions,val_resp_df)

		#Train full model
		#set.seed(3435)	
		full_rf_model<-ranger(y=full_resp_df, x=full_pred_df, num.trees=tree_num, mtry=var_tries, write.forest=TRUE,probability=FALSE,importance="impurity",min.node.size=node_size,num.threads=threads)

		#Evaluate model performance by calculating perason correlation and RMSE for the full set
		cor_f<-cor.test(full_rf_model$predictions,full_resp_df,method="pearson")
		rmse_f<-rmse(full_rf_model$predictions,full_resp_df)
			
		if ((cor_f$estimate < cor_cutoff) || (rmse_f > max_rmse)) {passed<-FALSE } else {passed<-TRUE}
			
		all_rf_stats<-rbind(all_rf_stats,c(resp_var,ticount,cor_t$estimate,rmse_t,cor_v$estimate,rmse_v,cor_f$estimate,rmse_f,passed))

		#Keep track of the best model based on RMSE
		if (rmse_f <= best_rmse) {
			best_rmse<-rmse_f
			if ((cor_f$estimate >= cor_cutoff) & (rmse_f <= max_rmse)) {
				best_rf_passed<-TRUE
			}
			best_rfs[[resp_var]]<-full_rf_model
			pred_data<-as.data.frame(cbind(full_rf_model$predictions,full_resp_df))
			colnames(pred_data)<-c("Predicted","Measured")
			pred_data$Response_Var<-resp_var
			pred_data$MG_UID<-as.factor(rownames(subdata))
			pred_data$Error<-sqrt((pred_data$Predicted-pred_data$Measured)**2)
		}		
	}

	#Keep track of the impotance values assigned to predictors for the best full RF model only
	#if (is.null(best_rfs[[resp_var]])) {next} #Skip if no model was succesfuly built for this response variable
	best_full_RF_model<-best_rfs[[resp_var]]
	best_importance<-as.data.frame(cbind(names(best_full_RF_model$variable.importance),best_full_RF_model$variable.importance))
	colnames(best_importance)<-c("Predictor","importance")
	best_importance$importance<-as.numeric(best_importance$importance)
	best_importance$Response<-resp_var
	#Convert absolute importance to percentages: Note that impurity method only yields positive importance values
	best_importance$Relative_Importance<-(best_importance$importance/sum(best_importance$importance))*100
	best_importance_vals<-rbind(best_importance_vals,best_importance)
}

colnames(all_rf_stats)<-c("Response","Iteration","Pearson_Cor_Training","RMSE_Training","Pearson_Cor_Validation","RMSE_Validation","Pearson_Cor_Full","RMSE_Full","Passed")
write.table(all_rf_stats,file="Reverse_Ecology_All_RF_Stats.tsv",sep="\t",append=FALSE,row.names=FALSE,col.names=TRUE,quote=FALSE)

write.table(best_importance_vals,file="Reverse_Ecology_Best_RF_Importance.tsv",sep="\t",append=FALSE,row.names=FALSE,col.names=TRUE,quote=FALSE)

#write.table(all_preds,file="Reverse_Ecology_Best_RF_Predictions.tsv",sep="\t",append=FALSE,row.names=FALSE,col.names=TRUE,quote=FALSE)

save.image("Reverse_Ecology_RF.RData")

quit(status=1)