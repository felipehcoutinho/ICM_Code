library(scales)
library(foreach)
library(doParallel)
library(NeuralNetTools)
library(nnet)
library(Metrics)
library(dplyr)
library(caret)

load_defaults<-FALSE
if  (load_defaults == TRUE) {
	response_file<-"/mnt/lustre/scratch/fcoutinho/StG_23/Viruses/Abundance/RPKM_Abundance_NA.tsv"#"/mnt/lustre/scratch/fcoutinho/StG_23/Prokaryotes/Abundance/RPKM_Abundance_OceanDNA_All_Species_Rep_MAGs_by_Species_Cluster.tsv"#"/mnt/lustre/scratch/fcoutinho/StG_23/Viruses/Abundance/Raw_Abundance_Min_Prev_50_Samples_IMGVR_Scaffolds_Comp_75-100_Max_Conta_0_vOTU_Representatives.tsv"
	predictor_file<-"/mnt/lustre/scratch/fcoutinho/MaLME/Marine_Viral_Communities_Sample_Metadata.tsv"#"/mnt/lustre/scratch/fcoutinho/StG_23/Metadata/TARA_Salazar_19_MGs_Metadata.tsv"
	max_opt_tries<-100
	max_rmse<-0.5
	n_neurons<-5
	predictor_variables<-c("Temperature","Salinity","Oxygen","Chlorophyll_A","Iron_5m","Ammonium_5m")#c("Temperature","Salinity","Oxygen","ChlorophyllA","Iron.5m","Ammonium.5m")
	out_prefix<-"Viruses"
	min_resp_prev<-50
}

###Read in command line arguments
args = commandArgs(trailingOnly=TRUE)
response_file<-args[1]
predictor_file<-args[2]
n_neurons<-as.numeric(args[3])
max_opt_tries<-as.numeric(args[4])
max_rmse<-as.numeric(args[5])
min_pcc<-as.numeric(args[6])
predictor_variables<-strsplit(args[7],",")[[1]]
min_resp_prev<-as.numeric(args[8])
out_prefix<-args[9]


###Load and filter predictor data
raw_predictor_df<-read.table(file=predictor_file,sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE,check.names=FALSE)

colnames(raw_predictor_df)[1]<-"MG_UID"

raw_predictor_df<-subset(raw_predictor_df,select=c("MG_UID",predictor_variables))

raw_predictor_df<-raw_predictor_df[complete.cases(raw_predictor_df),]

sample_count<-nrow(raw_predictor_df)

#summary(raw_predictor_df)

#Build the scenarios of climate change for the raw predictors
raw_predictor_df$Scenario<-"Current"

warming_data<-raw_predictor_df
warming_data$Scenario<-"Warming"
warming_data$Temperature<-(warming_data$Temperature + 2)

freshening_data<-raw_predictor_df
freshening_data$Scenario<-"Freshening"
freshening_data$Salinity<-(freshening_data$Salinity - 5)

combined_data<-raw_predictor_df
combined_data$Scenario<-"Both"
combined_data$Temperature<-(combined_data$Temperature + 2)
combined_data$Salinity<-(combined_data$Salinity - 5)

raw_predictor_df<-rbind(raw_predictor_df,warming_data,freshening_data,combined_data)

scenario_names<-raw_predictor_df$Scenario
raw_predictor_df$Scenario<-NULL

MG_UIDs<-raw_predictor_df$MG_UID
raw_predictor_df$MG_UID<-NULL

z_predictor_df<-as.data.frame(scale(raw_predictor_df,center=TRUE,scale=TRUE))
z_predictor_df$Scenario<-as.factor(scenario_names)
z_predictor_df$MG_UID<-MG_UIDs

#print("Summary of Z trandformed predictors:")
#summary(z_predictor_df)

#Load and filter response data
raw_response_df<-read.table(file=response_file,sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE,row.names=1,check.names=FALSE)

matched_MG_UIDs<-colnames(raw_response_df)[colnames(raw_response_df) %in% z_predictor_df$MG_UID]

raw_response_df<-raw_response_df[,matched_MG_UIDs]

#summary(raw_response_df[,1:11])

t_raw_response_df<-as.data.frame(t(raw_response_df))

#Calculate the number on non-zero values for each column of t_raw_response_df
non_zero_counts<-apply(t_raw_response_df,2,function(x) {sum(x != 0)})
names(non_zero_counts)<-colnames(t_raw_response_df)
response_prev_perc<-(non_zero_counts/sample_count)*100
print("Summary of response variables prevalence percentage:")
summary(response_prev_perc)
passed_response_variables<-names(response_prev_perc)[which(response_prev_perc >= min_resp_prev)]

print(paste("Number of response variables passed minimum prevalence cutoff:",length(passed_response_variables),sep=" "))

response_means<-colMeans(t_raw_response_df)
names(response_means)<-colnames(t_raw_response_df)

prev_abd_df<-as.data.frame(cbind(names(response_prev_perc),response_prev_perc,response_means))
colnames(prev_abd_df)<-c("Response_Variable_Name","Response_Variable_Prevalence_Percentage","Response_Mean_Abundance")
prev_abd_df$Response_Variable_Prevalence_Percentage<-as.numeric(prev_abd_df$Response_Variable_Prevalence_Percentage)
prev_abd_df$Response_Mean_Abundance<-as.numeric(prev_abd_df$Response_Mean_Abundance)
prev_abd_df$Response_Variable_Name<-as.factor(prev_abd_df$Response_Variable_Name)
#summary(prev_abd_df)

#passed_response_variables<-passed_response_variables[1:10]

#Filter out columns in t_raw_response_df the percentage of non-zeo values is less than min_resp_prev
f_t_raw_response_df<-subset(t_raw_response_df,select=passed_response_variables)

f_t_z_response_df<-as.data.frame(scale(f_t_raw_response_df,center=TRUE,scale=TRUE))

f_t_z_response_df$MG_UID<-as.factor(rownames(f_t_z_response_df))

#dim(f_t_z_response_df)
#print("Summary of first 10 columns of the Z transformed, filtered, response variable data:")
#summary(f_t_z_response_df[,1:10])

shouter_val<-100
vars_done<-0
#my.cluster <- parallel::makeCluster(47, type = "PSOCK")
#doParallel::registerDoParallel(cl = my.cluster)
#all_results_list<-foreach(rvar = passed_response_variables) %dopar% {

z_full_df<-merge(f_t_z_response_df,z_predictor_df[which(z_predictor_df$Scenario == "Current"),], by="MG_UID")

print("Building ANNs")
all_results_list<-foreach(rvar = passed_response_variables) %do% {
	#Subset the data to only include the response variable and the predictor variables
	z_sub_df<-z_full_df[,c("MG_UID",rvar,predictor_variables)]
	z_sub_df$MG_UID<-NULL
	z_sub_df$Scenario<-NULL
	colnames(z_sub_df)[colnames(z_sub_df) == rvar]<-"Z_Response_Variable"
	
	#Split the data into a training set and a validation set
	trainIndex <- createDataPartition(z_sub_df$Z_Response_Variable, p = 0.8, list = FALSE)
	z_train_df<-z_sub_df[trainIndex, ]
	z_val_df<-z_sub_df[-trainIndex, ]

	#Keep track performance metrics for the best networks on the training and validation sets
	passed_t<-FALSE
	best_ANN_t<-NA
	best_rmse_t<-+Inf
	best_rmse_v<-+Inf
	best_pearson_cor_t<-NA
	best_pearson_cor_v<-NA

	#Keep track of the number of optimization tries performed for the current resonse variable
	optimization_num<-0
	set.seed(666)
	#Perform multiple optimizations of the ANN on the training set until a network is found that passes the performance evaluation criteria on the validation set, or reaches the maximum number of optimization iterations allwoed
	while ((passed_t == FALSE) & (optimization_num < max_opt_tries)) {
		#Keep count of how many optimization tries were performed for each response variable
		optimization_num<-optimization_num+1
		shouter<-paste("Response",rvar,"Training/Validation Set Iteration:",optimization_num,sep=" ")
		#print(shouter)

		#Train an ANN using the trainingtrainign set only
		t_net<-nnet(Z_Response_Variable ~ .,z_train_df,size=n_neurons,linout=TRUE,maxit=1000,decay=0.001,reltol=0.001,trace=FALSE)
		
		#If there is no variation in the fitted values or the response variable, the trainign failed, so we skip to the next iteration
		if ((sd(t_net$fitted.values) == 0) | (sd(z_train_df$Z_Response_Variable) == 0)) {
			next
		}

		#Calculate performance metrics (pearson correlation coefficient and root mean squared error) of the training network on the training set
		pearson_cor_t<-cor.test(z_train_df$Z_Response_Variable,t_net$fitted.values,method="pearson")
		rmse_t<-rmse(z_train_df$Z_Response_Variable,t_net$fitted.values)
		#OBtain the predictions of the validation set using the training network
		val_preds<-as.data.frame(predict(t_net,z_val_df,type=c("raw")))
		colnames(val_preds)[1]<-"Predicted_Z_Response_Variable"
		#Calculate performance metrics (pearson correlation coefficient and root mean squared error) of the training network on the validation set
		pearson_cor_v<-cor.test(z_val_df$Z_Response_Variable,val_preds$Predicted_Z_Response_Variable,method="pearson")
		rmse_v<-rmse(z_val_df$Z_Response_Variable,val_preds$Predicted_Z_Response_Variable)
		
		#Keep track of the best network on the training validation set. The best network is the one with the lowest RMSE on the validation set
		if (rmse_v < best_rmse_v) {
			best_ANN_t<-t_net
			best_rmse_t<-rmse_t
			best_pearson_cor_t<-pearson_cor_t$estimate
			best_rmse_v<-rmse_v
			best_pearson_cor_v<-pearson_cor_v$estimate
		}

		#A variable is considered to have passed training if the best network's RMSE on the validation set is less than the max_rmse and the best network's PCC on the validation set is greater than the min_pcc
		if ((best_rmse_v <= max_rmse) & (best_pearson_cor_v >= min_pcc)) {
			passed_t<-TRUE
		}

		#Print the network performance of the current network
		shouter<-paste("Response Variable:",rvar,"Training/Validation Set Iteration:",optimization_num,"RMSE:",round(rmse_t,digits=2),"PCC training:",round(pearson_cor_t$estimate,digits=2),"RMSE validation:",round(rmse_v,digits=2),"PCC validation:",round(pearson_cor_v$estimate,digits=2),sep=" ")
		print(shouter)
	}

	#Print the network performance of the best training network
	shouter<-paste(">>> Response Variable:",rvar,"best RMSE training:",round(best_rmse_t,digits=2),"best PCC training:",round(best_pearson_cor_t,digits=2),"best RMSE validation:",round(best_rmse_v,digits=2),"best PCC validation:",round(best_pearson_cor_v,digits=2),"Total Training/Validation Set Optimizations:",optimization_num,"Passed training:",passed_t,sep=" ")
	print(shouter)

	#Keep track performance metrics for the best networks on the full set
	passed_f<-FALSE
	best_ANN_f<-NA
	best_rmse_f<-+Inf
	best_pearson_cor_f<-NA
	#Restart the optimization counter
	optimization_num<-0
	set.seed(42069)
	#Perform multiple optimizations of the ANN on the full set until a network is found that passes the evaluation criteria on the training set (i.e. now the full set) or reaches the maximum number of optimization iterations allwoed
	while ((passed_f == FALSE) & (optimization_num < max_opt_tries)) {
		optimization_num<-optimization_num+1
		shouter<-paste("Response",rvar,"Full Set Iteration:",optimization_num,sep=" ")
		#print(shouter)
		full_net<-nnet(Z_Response_Variable ~ .,z_sub_df,size=n_neurons,linout=TRUE,maxit=1000,decay=0.001,reltol=0.001,trace=FALSE)
		if ((sd(full_net$fitted.values) == 0) | (sd(z_sub_df$Z_Response_Variable) == 0)) {next}

		#Calculate performance metrics (pearson correlation coefficient and root mean squared error) of the full network on the full set
		pearson_cor_full<-cor.test(z_sub_df$Z_Response_Variable,full_net$fitted.values,method="pearson")
		rmse_full<-rmse(z_sub_df$Z_Response_Variable,full_net$fitted.values)

		#Keep track of the best network on the fulln set. The best network is the one with the lowest RMSE on the full set
		if (rmse_full < best_rmse_f) {
			best_ANN_f<-full_net
			best_rmse_f<-rmse_full
			best_pearson_cor_full<-pearson_cor_full$estimate
		}

		#A variable is considered to have passed if the best full network's RMSE on the full set is less than the max_rmse and the best full network's PCC on the full set is equal or greater than min_pcc
		if ((best_rmse_f <= max_rmse) & (best_pearson_cor_full >= min_pcc)) {
			passed_f<-TRUE
		}

		#Print the network performance of the current network
		shouter<-paste("Response Variable:",rvar,"Full Set Iteration:",optimization_num,"RMSE full:",round(rmse_full,digits=2),"PCC Full:",round(pearson_cor_full$estimate,digits=2),"Passed Full:",passed_f,sep=" ")
		print(shouter)
	}

	#Print the network performance of the best full network
	shouter<-paste(">>> Response Variable:",rvar,"best RMSE full:",round(best_rmse_f,digits=2),"best PCC full:",round(best_pearson_cor_full,digits=2),"Total Full Set Optimizations:",optimization_num,"Passed Full:",passed_f,sep=" ")
	print(shouter)

	#Calculate the importance of the predictor variables in the best full network using the Olden method
	importance_olden<-NA
	importance_olden<-olden(best_ANN_f,bar_plot=FALSE)
	#Normalize the importance values to a range of -1 to 1
	importance_olden$Normalized_Importance<-importance_olden$importance/max((abs(importance_olden$importance)))
	importance_olden$Predictor_Variable<-rownames(importance_olden)

	#Calculate the predictions of the best full network on the scenarios of climate change
	scenario_predictions<-NA
	scenario_predictions<-as.data.frame(predict(best_ANN_f,z_predictor_df,type=c("raw")))
	colnames(scenario_predictions)[1]<-"Predicted_Z_Response_Variable"
	scenario_predictions$Scenario<-as.factor(scenario_names)
	scenario_predictions$Predicted_Raw_Response_Variable<-mean(t_raw_response_df[,rvar]) + (scenario_predictions$Predicted_Z_Response_Variable*sd(t_raw_response_df[,rvar]))
	scenario_predictions$Non_Negative_Predicted_Raw_Response_Variable<-scenario_predictions$Predicted_Raw_Response_Variable	
	scenario_predictions$Non_Negative_Predicted_Raw_Response_Variable[which(scenario_predictions$Non_Negative_Predicted_Raw_Response_Variable < 0)]<-0
	scenario_predictions$Response_Variable_Name<-rvar
	scenario_predictions$MG_UID<-MG_UIDs

	#Print the number of response variables that have been processed every 100 response variables
	vars_done<-vars_done+1
	if (vars_done == shouter_val) {
		shouter<-paste("ANNs for ",vars_done," response variables finished.",sep="")
		#print(shouter)
		shouter_val<-shouter_val+100
	}

	#Return a list of the results for the response variable
	list("Response_Variable_Name" = rvar, "Best_Full_ANN" = best_ANN_f, "Best_Full_ANN_RMSE" = best_rmse_f, "Best_Full_ANN_PCC" = best_pearson_cor_full, "Best_Full_ANN_Olden_Importance" =  importance_olden, "Best_Full_ANN_Scenario_Predictions" = scenario_predictions, "Best_Full_ANN_Passed" = passed_f, "At_least_One_TV_ANN_Passed" = passed_t)
}


#parallel::stopCluster(cl = my.cluster)

print("Finished building Networks")

###Print results to table files
get_preds<-function(x){
	return(x$Best_Full_ANN_Scenario_Predictions)
}
all_preds_df <-do.call(rbind,lapply(all_results_list,FUN=get_preds))
outfile<-paste(out_prefix,"_Parallel_Niche_Modelling_All_Scenario_Predictions.tsv",sep="")
write.table(all_preds_df,file=outfile,sep="\t",append=FALSE,row.names=FALSE,col.names=TRUE,quote=FALSE)

get_perf_metrics<-function(x){
	#print(paste("Fetching data for: ",x$Response_Variable_Name,sep=" "))
	ann_metrics<-c(x$Response_Variable_Name,x$Best_Full_ANN_PCC,x$Best_Full_ANN_RMSE,x$Best_Full_ANN_Passed,x$At_least_One_TV_ANN_Passed)
	names(ann_metrics)<-c("Response_Variable_Name","Best_Full_ANN_PCC","Best_Full_ANN_RMSE","Best_Full_ANN_Passed","At_least_One_TV_ANN_Passed")
	if (sum(is.na(c(x$Best_Full_ANN_Olden_Importance[predictor_variables,"importance"],x$Best_Full_ANN_Olden_Importance[predictor_variables,"Normalized_Importance"]))) == 0) {
		raw_imps<-x$Best_Full_ANN_Olden_Importance[predictor_variables,"importance"]#[["importance"]]
		norm_imps<-x$Best_Full_ANN_Olden_Importance[predictor_variables,"Normalized_Importance"]#[["Normalized_Importance"]]
	} else {
		raw_imps<-rep(NA,length(predictor_variables))
		norm_imps<-rep(NA,length(predictor_variables))
	}
	names(raw_imps)<-paste("Raw_Importance",predictor_variables,sep="_")
	names(norm_imps)<-paste("Normalized_Importance",predictor_variables,sep="_")
	rvar_prev<-prev_abd_df[x$Response_Variable_Name,"Response_Mean_Abundance"]
	names(rvar_prev)<-"Response_Variable_Mean_Abundance"
	rvar_mean_abd<-prev_abd_df[x$Response_Variable_Name,"Response_Variable_Prevalence_Percentage"]
	names(rvar_mean_abd)<-"Response_Variable_Prevalence_Percentage"
	ann_metrics<-c(ann_metrics,raw_imps,norm_imps,rvar_prev,rvar_mean_abd)
	#ann_metrics<-as.data.frame(merge(ann_metrics,prev_abd_df,by="Response_Variable_Name",all.x=TRUE))
	return(ann_metrics)
}

all_perf_df <-do.call(rbind,lapply(all_results_list,FUN=get_perf_metrics))
outfile<-paste(out_prefix,"_Parallel_Niche_Modelling_Performance_Metrics.tsv",sep="")
write.table(all_perf_df,file=outfile,sep="\t",append=FALSE,row.names=FALSE,col.names=TRUE,quote=FALSE)

outfile<-paste(out_prefix,"_Parallel_Niche_Modelling_Workspace.RData",sep="")
save.image(outfile)
quit(status=0)