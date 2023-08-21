#Load the required packages
library(caret)
library(scales)
library(nnet)
library(NeuralNetTools)
library(Metrics)
library(ggplot2)
library(GGally)

#Read in the time-series data
full_ts_raw_data <- read.csv("/mnt/lustre/scratch/fcoutinho/xlopez/Blanes_TS/ANN_TS_input_data.csv", dec=",", header = TRUE, sep = ";",quote="",stringsAsFactors = TRUE)
#Print summary of full_ts_raw_data
summary(full_ts_raw_data)
dim(full_ts_raw_data)

#Apply the subset function to full_ts_raw_data to generate sub_full_ts_raw_data which only includes the columns with variables AbunVir, Salinity, Temperature, PO4, CHL, AbunBact, and Abun_HNF
sub_ts_raw_data <- subset(full_ts_raw_data, select = c(AbunVir, Salinity, AbunBact_HNA, PO4, CHL, Secchi_Disk))
#Remove rows with NA values usign complete.cases
comp_sub_ts_raw_data <- sub_ts_raw_data[complete.cases(sub_ts_raw_data),]
#Print summary of comp_sub_ts_raw_data
summary(comp_sub_ts_raw_data)
#z-transform the variables in comp_sub_ts_raw_data using the scale function and assign the output to z_comp_sub_ts_raw_data
comp_sub_ts_z_data <- as.data.frame(scale(comp_sub_ts_raw_data))
summary(comp_sub_ts_z_data)

set.seed(42069)
fitControl <- trainControl( method = "LGOCV", number = 10, repeats=10, savePredictions = "all", returnResamp="all", p=0.7,allowParallel=TRUE)

HyperPGrid<-expand.grid(size=c(1:10),decay=c(0,1E-1,1E-2,1E-3,1E-4,1E-5))

AnnTuning <- train(AbunVir ~ ., data = comp_sub_ts_z_data, method = "nnet",trControl = fitControl, tuneGrid = HyperPGrid,linout = TRUE, trace=FALSE, )

