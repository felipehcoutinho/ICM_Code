library(ggplot2)
library(gplots)
library(reshape2)
library(RColorBrewer)
library("GGally")
library(ggpubr, lib.loc="/mnt/lustre/bio/users/fcoutinho/Rlibs/")
require(maps)

library(scales)
library(foreach)
library(doParallel)
library(NeuralNetTools)
library(nnet)
library(Metrics)
library(dplyr)
library(caret)

args = commandArgs(trailingOnly=TRUE)
response_variables<-strsplit(args[1],",")[[1]]
predictor_variables<-strsplit(args[2],",")[[1]]
max_opt_tries<-as.numeric(args[3])
max_rmse<-as.numeric(args[4])
n_neurons<-as.numeric(args[5])
basic_plots<-args[6]

load_defaults<-FALSE
if  (load_defaults == TRUE) {
	response_variables<-c("Viral_Abundance")
	predictor_variables<-c("Depth","Prokaryote_Abundance","Temperature","Salinity")
	max_opt_tries<-1000
	max_rmse<-0.8
	n_neurons<-10
}


data<-read.table(file="/mnt/lustre/bio/users/fcoutinho/Databases/gVOD/gVOD_Raw_Data.tsv",header=TRUE,sep="\t",quote="",comment="",stringsAsFactors=TRUE,na.string=c("NA",""))
subdata<-subset(data,select=c("No","Depth.water..m.","Virus....ml...Particles.ml.","Prokaryotes....ml.","Temp...C.","Sal","pH"))
colnames(subdata)<-c("Sample_ID","Depth","Viral_Abundance","Prokaryote_Abundance","Temperature","Salinity","pH")
subdata$VMR<-subdata$Viral_Abundance/subdata$Prokaryote_Abundance
subdata<-subdata[which(subdata$Prokaryote_Abundance < 31082400 & subdata$Prokaryote_Abundance > 0),]
#summary(subdata)

###Step 2: Generate ANN models for the current and future scenarios
raw_predictor_df<-subdata
raw_predictor_df<-subset(raw_predictor_df,select=c("Sample_ID",response_variables,predictor_variables))
raw_predictor_df<-raw_predictor_df[complete.cases(raw_predictor_df),]
#summary(raw_predictor_df)

#Build the scenarios of climate change for the raw predictors
num_cols<-nums<-unlist(lapply(raw_predictor_df, is.numeric), use.names = FALSE)

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
#Replace negative salinity values in raw_predictor_df with 0
raw_predictor_df$Salinity[raw_predictor_df$Salinity < 0]<-0

scenario_names<-raw_predictor_df$Scenario
raw_predictor_df$Scenario<-NULL
Sample_IDs<-raw_predictor_df$Sample_ID
raw_predictor_df$Sample_ID<-NULL

z_predictor_df<-as.data.frame(scale(raw_predictor_df,center=TRUE,scale=TRUE))
z_predictor_df$Scenario<-as.factor(scenario_names)
z_predictor_df$Sample_ID<-Sample_IDs
summary(z_predictor_df)

#my.cluster <- parallel::makeCluster(47, type = "PSOCK")
#doParallel::registerDoParallel(cl = my.cluster)
#all_results_list<-foreach(rvar = response_variables) %dopar% {
all_results_list<-foreach(rvar = response_variables) %do% {
	z_full_df<-subset(z_predictor_df[which(z_predictor_df$Scenario == "Current"),],select=c(rvar,predictor_variables,"Sample_ID"))
	z_full_df$Sample_ID<-NULL
	z_full_df$Scenario<-NULL
	colnames(z_full_df)[colnames(z_full_df) == rvar]<-"Z_Response_Variable"
	
	trainIndex <- createDataPartition(z_full_df$Z_Response_Variable, p = 0.8, list = FALSE)
	z_train_df<-z_full_df[trainIndex, ]
	z_val_df<-z_full_df[-trainIndex, ]

	passed_tv<-FALSE
	best_ANN_t<-NA
	optimization_num<-0
	best_rmse_t<-+Inf
	best_rmse_v<-+Inf
	best_pearson_cor_t<-NA
	best_pearson_cor_v<-NA

	passed_f<-FALSE
	best_ANN_f<-NA
	optimization_num<-0
	best_rmse_f<-+Inf
	best_pearson_cor_f<-NA

	set.seed(666)
	while ((best_rmse_v > max_rmse) & (optimization_num < max_opt_tries)) {
		optimization_num<-optimization_num+1
		shouter<-paste("Response",rvar,"Iteration",optimization_num,sep=" ")
		print(shouter)
		t_net<-nnet(Z_Response_Variable ~ .,z_train_df,size=n_neurons,linout=TRUE,maxit=1000,decay=0.001,reltol=0.001,trace=FALSE)
		
		if ((sd(t_net$fitted.values) == 0) | (sd(z_train_df$Z_Response_Variable) == 0)) {next}

		pearson_cor_t<-cor.test(z_train_df$Z_Response_Variable,t_net$fitted.values,method="pearson")
		rmse_t<-rmse(z_train_df$Z_Response_Variable,t_net$fitted.values)
		val_preds<-as.data.frame(predict(t_net,z_val_df,type=c("raw")))
		colnames(val_preds)[1]<-"Predicted_Z_Response_Variable"
		pearson_cor_v<-cor.test(z_val_df$Z_Response_Variable,val_preds$Predicted_Z_Response_Variable,method="pearson")
		rmse_v<-rmse(z_val_df$Z_Response_Variable,val_preds$Predicted_Z_Response_Variable)
		
		if (rmse_v < best_rmse_v) {
			best_ANN_t<-t_net
			best_rmse_t<-rmse_t
			best_pearson_cor_t<-pearson_cor_t$estimate
			best_rmse_v<-rmse_v
			best_pearson_cor_v<-pearson_cor_v$estimate

			if (best_rmse_v <= max_rmse) {
				passed_tv<-TRUE
			}
		}
	}

	shouter<-paste("Response",rvar,"best RMSE training:",best_rmse_t,"best PCC training:",best_pearson_cor_t,"best RMSE validation:",best_rmse_v,"best PCC validation:",best_pearson_cor_v,"Opitimizations:",optimization_num,"Passed training:",passed_tv,sep=" ")
	print(shouter)

	optimization_num<-0
	set.seed(42069)

	while ((best_rmse_f > max_rmse) & (optimization_num < max_opt_tries)) {
		optimization_num<-optimization_num+1

		full_net<-nnet(Z_Response_Variable ~ .,z_full_df,size=n_neurons,linout=TRUE,maxit=1000,decay=0.001,reltol=0.001,trace=FALSE)

		if ((sd(full_net$fitted.values) == 0) | (sd(z_full_df$Z_Response_Variable) == 0)) {next}

		pearson_cor_full<-cor.test(z_full_df$Z_Response_Variable,full_net$fitted.values,method="pearson")
		rmse_full<-rmse(z_full_df$Z_Response_Variable,full_net$fitted.values)

		if (rmse_full < best_rmse_f) {
			best_ANN_f<-full_net
			best_rmse_f<-rmse_full
			best_pearson_cor_full<-pearson_cor_full$estimate
			if (best_rmse_f <= max_rmse) {
				passed_f<-TRUE
			}
		}
	}

	shouter<-paste("Response Variable:",rvar,"Best RMSE full:",round(best_rmse_f,digits=3),"Best PCC full:",round(best_pearson_cor_full,digits=3),"Opitimizations:",optimization_num,"Passed Full:",passed_f,sep=" ")
	print(shouter)

	importance_olden<-NA
	importance_olden<-olden(best_ANN_f,bar_plot=FALSE)
	importance_olden$Normalized_Importance<-importance_olden$importance/max((abs(importance_olden$importance)))
	importance_olden$Predictor_Variable<-rownames(importance_olden)

	scenario_predictions<-NA
	scenario_predictions<-as.data.frame(predict(best_ANN_f,z_predictor_df,type=c("raw")))
	colnames(scenario_predictions)[1]<-"Predicted_Z_Response_Variable"
	scenario_predictions$Scenario<-as.factor(scenario_names)	
	scenario_predictions$Predicted_Raw_Response_Variable<-(mean(as.numeric(raw_predictor_df[,rvar])) + (scenario_predictions$Predicted_Z_Response_Variable)*sd(as.numeric(raw_predictor_df[,rvar])))		
	scenario_predictions$Predicted_Raw_Response_Variable[which(scenario_predictions$Predicted_Raw_Response_Variable < 0)]<-0
	scenario_predictions$Response_Variable_Name<-rvar
	scenario_predictions$Sample_ID<-Sample_IDs

	list("Response_Variable_Name" = rvar, "Best_Full_ANN" = best_ANN_f, "Best_Full_ANN_RMSE" = best_rmse_t, "Best_Full_ANN_PCC" = best_pearson_cor_t, "Best_Full_ANN_Olden_Importance" =  importance_olden, "Best_Full_ANN_Scenario_Predictions" = scenario_predictions, "Best_Full_ANN_Passed" = passed_f, "At_least_One_TV_ANN_Passed" = passed_tv)
}


#parallel::stopCluster(cl = my.cluster)

print("Finished building Networks")

###Step 3: Print results to table files
get_preds<-function(x){
	return(x$Best_Full_ANN_Scenario_Predictions)
}
all_preds_df <-do.call(rbind,lapply(all_results_list,FUN=get_preds))
write.table(all_preds_df,file="gVOD_Niche_Modelling_All_Scenario_Predictions.tsv",sep="\t",append=FALSE,row.names=FALSE,col.names=TRUE,quote=FALSE)

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
	ann_metrics<-c(ann_metrics,raw_imps,norm_imps)
	return(ann_metrics)
}
all_perf_df <-do.call(rbind,lapply(all_results_list,FUN=get_perf_metrics))
write.table(all_perf_df,file="gVOD_Niche_Modelling_Performance_Metrics.tsv",sep="\t",append=FALSE,row.names=FALSE,col.names=TRUE,quote=FALSE)

save.image("gVOD_ANN_Niche_Modelling_Workspace.RData")

###Step 1: Generate plots describing the spatiotempotal nature of the data and plot correlogram
if (basic_plots == TRUE) {
	sta_count<-table(subdata$Sample_ID)
	passed_stas<-names(sta_count[which(sta_count >= 50)])
	sta_subdata<-subset(subdata,Sample_ID %in% passed_stas)
	fig_1D<-ggplot(sta_subdata,aes(y=log10(Viral_Abundance),x=log10(Prokaryote_Abundance)))+geom_point(alpha=0.5,shape=18)+geom_smooth(method = "lm",formula= y~x, se = FALSE)+stat_cor(method = "pearson", alternative = "two.sided",cor.coef.name = "R",size=3,label.y=4.2)+theme_bw()+facet_wrap(Sample_ID ~ ., scales="free")

	ggsave("gVOD_ProkxVir_Abd_Scatterplot_by_Sample_ID.pdf",plot=fig_1D,width=20,height=20,pointsize=8)

	
	subdata<-subdata[,c("Prokaryote_Abundance","Viral_Abundance","Depth","Temperature","Salinity","pH")]
	summary(subdata)

	fig_1A<-ggpairs(subdata,upper=list(continuous = wrap("cor", method = "spearman",size= 5)))
	ggsave("gVOD_Correlogram.pdf",plot=fig_1A,width=10,height=10,pointsize=8)

	mdata<-melt(subdata,id=c("Viral_Abundance"))
	colnames(mdata)<-c("Viral_Abundance","Variable","Value")
	fig_1B<-ggplot(mdata,aes(y=log10(Viral_Abundance),x=Value))+geom_point(alpha=0.5)+geom_smooth(method = "loess",formula= y~x, se = TRUE)+facet_wrap(Variable ~ ., scales="free_x",nrow=1)+theme_bw()
	ggsave("gVOD_VirusxEnv_Scatterplot.pdf",plot=fig_1B,width=15,height=5,pointsize=8)

	mdata<-melt(subdata,id=c("Prokaryote_Abundance"))
	colnames(mdata)<-c("Prokaryote_Abundance","Variable","Value")
	fig_1B<-ggplot(mdata,aes(y=log10(Prokaryote_Abundance),x=Value))+geom_point(alpha=0.5)+geom_smooth(method = "loess",formula= y~x, se = TRUE)+facet_wrap(Variable ~ ., scales="free_x",nrow=1)+theme_bw()
	ggsave("gVOD_ProkaryotexEnv_Scatterplot.pdf",plot=fig_1B,width=15,height=5,pointsize=8)
	
	subdata<-subset(data,select=c("No","Latitude","Longitude","Depth.water..m.","Virus....ml...Particles.ml.","Prokaryotes....ml.","Temp...C.","Sal","pH"))
	colnames(subdata)<-c("UID","Latitude","Longitude","Depth","Viral_Abundance","Prokaryote_Abundance","Temperature","Salinity","pH")
	subdata<-subdata[,c("UID","Latitude","Longitude","Prokaryote_Abundance","Viral_Abundance","Depth","Temperature","Salinity","pH")]

	col_grad<-rev(colorRampPalette(brewer.pal(9,"RdYlBu")[-6])(n=299))
	world_map <- map_data("world")
	
	subdata<-subdata[which(subdata$Depth >= 0 & subdata$Viral_Abundance >= 0),]
	summary(subdata)

	subdata$Zone<-"Surface"
	subdata$Zone[which(subdata$Depth > 200)]<-"Deep"
	subdata$Zone<-factor(subdata$Zone,levels=c("Surface","Deep"))

	fig_1C<-ggplot(world_map, aes(x = long, y = lat, group = group))+geom_polygon(fill="lightgray", colour = "lightgray")+geom_point(subdata,alpha=0.5,size=1,mapping=aes(y=Latitude,x=Longitude,group=UID,col=log10(Viral_Abundance)))+coord_cartesian(xlim=c(-180,180), ylim=c(-90, 90))+theme_bw()+xlab("Longitude")+ylab("Latitude")+scale_colour_gradientn(colours =col_grad,limits=c(0,9))+facet_wrap(Zone ~ .)+theme(legend.position="top")
	
	ggsave("gVOD_Current_Global_Viral_Abundance.pdf",plot=fig_1C,width=20,height=10,pointsize=8)
	
	fig_1D<-ggplot(world_map, aes(x = long, y = lat, group = group))+geom_polygon(fill="lightgray", colour = "lightgray")+geom_point(subdata,alpha=0.5,size=1,mapping=aes(y=Latitude,x=Longitude,group=UID,col=log10(Prokaryote_Abundance)))+coord_cartesian(xlim=c(-180,180), ylim=c(-90, 90))+theme_bw()+xlab("Longitude")+ylab("Latitude")+scale_colour_gradientn(colours =col_grad,limits=c(0,9))+facet_wrap(Zone ~ .)+theme(legend.position="top")
	
	ggsave("gVOD_Current_Global_Prokaryote_Abundance.pdf",plot=fig_1C,width=20,height=10,pointsize=8)
	#fig_1<-ggarrange(fig_1C, fig_1D,  nrow = 2, labels = "AUTO") 

	#ggsave("StG_23_gVOD_Exploratory_1.pdf",plot=fig_1,width=16,height=12,pointsize=8)

}


quit(status=1)