#Load the required packages
library(caret)
library(scales)
library(nnet)
library(NeuralNetTools)
library(Metrics)
library(ggplot2)
library(reshape2)
library(GGally)
library(RColorBrewer)
#Read in the time-series data
full_ts_raw_data <- read.csv("/mnt/lustre/scratch/fcoutinho/xlopez/Blanes_TS/ANN_TS_input_data.csv", dec=",", header = TRUE, sep = ";",quote="",stringsAsFactors = TRUE)
full_ts_raw_data$Timepoint<-c(1:nrow(full_ts_raw_data))
full_ts_raw_data$Heteretrophic_Bacteria<-full_ts_raw_data$AbunBact
full_ts_raw_data$AbunBact<-full_ts_raw_data$AbunBact+full_ts_raw_data$Proclorococcus+full_ts_raw_data$Synechococcus
full_ts_raw_data$VMR<-full_ts_raw_data$AbunVir/full_ts_raw_data$AbunBact
summary(full_ts_raw_data)

cor_fn <- function(data, mapping, method="p", use="pairwise", ...){

              # grab data
              x <- eval_data_col(data, mapping$x)
              y <- eval_data_col(data, mapping$y)

              # calculate correlation
			  corr_val<-0
			  if ((length(x) > 5) & (length(y) > 5)) {
				corr_val <- cor(x, y, method=method, use="pairwise.complete.obs")
				}
              # calculate colour based on correlation value
              # Here I have set a correlation of minus one to blue, 
              # zero to white, and one to red 
              # Change this to suit: possibly extend to add as an argument of `my_fn`
              colFn <- colorRampPalette(rev(brewer.pal(9,"RdBu")), interpolate ='spline')
              fill <- colFn(100)[findInterval(corr_val, seq(-1, 1, length=100))]

              ggally_cor(data = data, size = 4, colour="black", mapping = mapping, ...) + 
                theme_void() +
                theme(panel.background = element_rect(fill=fill))
            }
			
scatter_fn <- function(data, mapping, method="p", use="pairwise", ...){
	
              # grab data
              x <- eval_data_col(data, mapping$x)
              y <- eval_data_col(data, mapping$y)

	ggally_smooth_loess(data = data, mapping=mapping, ...)+theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

scatter_fn_no_se <- function(data, mapping, method="p", use="pairwise", ...){
	
              # grab data
              x <- eval_data_col(data, mapping$x)
              y <- eval_data_col(data, mapping$y)

	ggally_smooth_loess(data = data, mapping=mapping, se=FALSE, ...)+theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

pdf("Blanes_Correlogram.pdf",width=35,height=35,pointsize=8)
plot<-ggpairs(full_ts_raw_data,columns=c("Temperature","Salinity", "Secchi_Disk","NO3","NH4" , "PO4","NO2","SiO4","CHL","AbunBact" , "AbunBact_HNA" , "Proclorococcus","Synechococcus" , "Abun_HNF" ,"Abun_PNF", "Abun_micromonas", "Abun_cryptomonas" , "AbunVir","Abun_V1","Abun_V2" ,"Abun_V3","Heteretrophic_Bacteria"),aes(alpha = 0.5),lower = list(continuous = scatter_fn),upper = list(continuous = cor_fn))
print(plot)
dev.off()

#Print summary of full_ts_raw_data
m_ts_data<-melt(full_ts_raw_data,id=c("SampleID","Year","Month","Timepoint"),variable.name="Variable",value.name="Value")

m_ts_data$Month<-factor(m_ts_data$Month,levels=c("January","February","March","April","May","June","July","August","September","October","November","December"))

ts_sct<-ggplot(m_ts_data,aes(x=Timepoint,y=Value))+
geom_smooth(method = "loess",se=FALSE, formula = y ~x)+
geom_point(size=3,aes(color=Month))+
theme_bw()+
facet_grid(Variable ~ ., scales="free_y")
ggsave("/mnt/lustre/scratch/fcoutinho/xlopez/Blanes_TS/Blanes_TS_Scatterplots.pdf",plot=ts_sct,width=15,height=25,pointsize=8)

ts_box<-ggplot(m_ts_data,aes(x=Month,y=log10(Value)))+
geom_boxplot(aes(fill=Month))+
theme_bw()+
theme(axis.text.x = element_text(angle = 45,hjust = 1))+
facet_grid(Variable ~ ., scales="free_y")
ggsave("/mnt/lustre/scratch/fcoutinho/xlopez/Blanes_TS/Blanes_TS_Boxplots_by_Month.pdf",plot=ts_box,width=7,height=35,pointsize=8)

ts_sct<-ggplot(m_ts_data,aes(x=Month,y=Value,group=Year))+
geom_smooth(method = "loess",se=FALSE, formula = y ~x,aes(colour=as.character(Year)),alpha=0.8)+
geom_point(aes(fill=Year),size=23,shape=23)+
theme_bw()+
facet_grid(Variable ~ ., scales="free_y")
ggsave("/mnt/lustre/scratch/fcoutinho/xlopez/Blanes_TS/Blanes_TS_Scatterplots_by_Month_and_Year.pdf",plot=ts_sct,width=9,height=35,pointsize=8)

#summary(full_ts_raw_data)
#dim(full_ts_raw_data)

#Apply the subset function to full_ts_raw_data to generate sub_full_ts_raw_data which only includes the columns with variables AbunVir, Salinity, Temperature, PO4, CHL, AbunBact, and Abun_HNF
sub_ts_raw_data <- subset(full_ts_raw_data, select = c(AbunVir, Salinity, AbunBact, PO4, CHL, Secchi_Disk, Proclorococcus, Synechococcus))
#Remove rows with NA values usign complete.cases
comp_sub_ts_raw_data <- sub_ts_raw_data[complete.cases(sub_ts_raw_data),]

comp_sub_ts_raw_data$AbunBact<-comp_sub_ts_raw_data$AbunBact+comp_sub_ts_raw_data$Proclorococcus+comp_sub_ts_raw_data$Synechococcus
dim(comp_sub_ts_raw_data)
summary(comp_sub_ts_raw_data)
#z-transform the variables in comp_sub_ts_raw_data using the scale function and assign the output to z_comp_sub_ts_raw_data
comp_sub_ts_z_data <- as.data.frame(scale(comp_sub_ts_raw_data))
#comp_sub_ts_z_data <- comp_sub_ts_raw_data
summary(comp_sub_ts_z_data)

#Split comp_sub_ts_raw_data into training and validation through createDataPartition based on the AbunVir variable.
set.seed(34+35)
trainIndex <- createDataPartition(comp_sub_ts_z_data$AbunVir, p = 0.7, list = FALSE)
train_sub_ts_z_data <- comp_sub_ts_z_data[trainIndex, ]
test_sub_ts_z_data <- comp_sub_ts_z_data[-trainIndex, ]
#Print summaries of train_sub_ts_z_data and test_sub_ts_z_data
summary(train_sub_ts_z_data)
summary(test_sub_ts_z_data)

#Create a neural network model using the nnet package using the training data, 5 neurons, maixmum of 1000 iterations, and a linear activation function
set.seed(42069)
tn_model <- nnet(AbunVir ~ Salinity+AbunBact+PO4+CHL+Secchi_Disk, data = train_sub_ts_z_data, size = 3, maxit = 1000, reltol=1E-5, linout = TRUE)#
#Evaluate model performance in the trainning data through pearson correlation
cor.test(tn_model$fitted.values, train_sub_ts_z_data$AbunVir,method="pearson")
rmse(tn_model$fitted.values, train_sub_ts_z_data$AbunVir)

#Evaluate model performance in the validation data through pearson correlation
test_pred <- predict(tn_model, test_sub_ts_z_data)
cor.test(test_pred, test_sub_ts_z_data$AbunVir,method="pearson")