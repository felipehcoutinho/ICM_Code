#Read in the example dataset table as a Data Frame
example_df<-read.table(file="/mnt/lustre/scratch/fcoutinho/MaLME/Example_ANN_Input_Data_1.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)

#Print a summary of the dataframe
summary(example_df)

#Subset the table to only keep complete cases (i.e. Samples with non  NA values for all predictor and response variables)
sub_example_df<-example_df[complete.cases(example_df), ]

summary(sub_example_df)

#Load the ANN training library
library(nnet)

#Set a seed number to make your results reproducible
seed_num<-12223
set.seed(seed_num)

#Train the network
trained_net<-nnet(Shannon_Diversity ~ Temperature+Depth+Chlorophyll_A+Salinity+Ammonium_5m,sub_example_df,size=3,linout=TRUE,maxit=1000)

rmse(sub_example_df$Shannon_Diversity,trained_net$fitted.values)
cor.test(sub_example_df$Shannon_Diversity,trained_net$fitted.values, method="pearson")

#Check the contents of the output ANN object
ls(trained_net)

#Load the ANN interpretation library
library(NeuralNetTools)

#Calculate the importance of the predictors using the Olden method
importance_olden<-olden(trained_net,bar_plot=FALSE)

#Look at the outpt and identify the most important predictor
importance_olden

#Calculate RMSE
library(Metrics)
rmse(sub_example_df$Shannon_Diversity,trained_net$fitted.values)

#Calculate Pearson R between measured and predicted values
cor.test(sub_example_df$Shannon_Diversity,trained_net$fitted.values,method="pearson")

#Load the caret library
library(caret)

#Assign samples to the test and validation sets
set.seed(seed_num)

train_set_row_nums<-createDataPartition(sub_example_df$Shannon_Diversity,p=0.5,list=FALSE)

sub_example_df_train<-sub_example_df[train_set_row_nums,]

sub_example_df_valid<-sub_example_df[-train_set_row_nums,]

dim(sub_example_df_train)
dim(sub_example_df_valid)

summary(sub_example_df_train)
summary(sub_example_df_valid)

#Test if the differences in the values of the response variable are significantly different between train and test sets
wilcox.test(sub_example_df_train$Shannon_Diversity,sub_example_df_valid$Shannon_Diversity,paired=FALSE,exact=FALSE)

#Train a model using the training set only
set.seed(seed_num)

trained_net_ts<-nnet(Shannon_Diversity ~ Temperature+Depth+Chlorophyll_A+Salinity+Ammonium_5m,sub_example_df_train,size=3,linout=TRUE,maxit=1000)

#Evaluate the performance of the new model on the training set
rmse(sub_example_df_train$Shannon_Diversity,trained_net_ts$fitted.values)

#Calculate Pearson R between measured and predicted values of the model using the training set
cor.test(sub_example_df_train$Shannon_Diversity,trained_net_ts$fitted.values, method="pearson")

#Get predictions for the validation set
valid_preds<-predict(trained_net_ts,sub_example_df_valid,type=c("raw"))

summary(valid_preds)

#Evaluate the performance of the new model on the validation set
rmse(sub_example_df_valid$Shannon_Diversity,valid_preds)

cor.test(sub_example_df_valid$Shannon_Diversity,valid_preds)
