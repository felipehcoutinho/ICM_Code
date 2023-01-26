library(ggplot2)
library(RColorBrewer)


com_data_file<-""
out_fig_name<-"Virome_Community.svg"
#com_df<-read.table(file=com_data_file,sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)

deg2rad <- function(deg) {
	rads<-((deg * pi ) / 180)
	return(rads)
}

get_rot_coord<-function(ent_df=NA,rot=0) {
	#As specified in (p replaced by phi, and q replaced by rho): https://math.stackexchange.com/questions/270194/how-to-find-the-vertices-angle-after-rotation
	rot_ent_df<-ent_df
	phi<-(min(rot_ent_df$x_coord)+max(rot_ent_df$x_coord)) / 2
	rho<-(min(rot_ent_df$y_coord)+max(rot_ent_df$y_coord)) / 2
	rads<-deg2rad(rot)
	sin_rot<-sin(rads)
	cos_rot<-cos(rads)
	print(paste("Phi:",phi,"Rho:",rho,"Rotation:",rot,"Radians:",rads,"sin:",sin_rot,"cos:",cos_rot,sep=" "))
	for (i in c(1:nrow(rot_ent_df))) {
		ori_x_coord<-rot_ent_df$x_coord[i]
		ori_y_coord<-rot_ent_df$y_coord[i]
		rot_ent_df$x_coord[i]<-((ori_x_coord - phi)*cos_rot) - ((ori_y_coord - rho)*sin_rot) + phi
		rot_ent_df$y_coord[i]<-((ori_x_coord - phi)*sin_rot) + ((ori_y_coord - rho)*cos_rot) + rho
		print(paste("Original (x,y):",ori_x_coord,",",ori_y_coord,sep=" "))
		print(paste("Rotated (x,y):",rot_ent_df$x_coord[i],",",rot_ent_df$y_coord[i],sep=" "))
	}
	print(rot_ent_df)
	return(rot_ent_df)
}

build_entity<-function(ent_name=NA,morphology="Unknown",col_var=NA,rot=0) {
	if (morphology == "Myoviridae") {
		capsid_df<-data.frame(x_coord=c(4,2,-2,-4,-2,2),y_coord=c(0,3.5,3.5,0,-3.5,-3.5),Subpart="Capsid",SUVI=c(1:6))
		tail_df<-data.frame(x_coord=c(-1.3,1.3,1.3,-1.3),y_coord=c(-3.5,-3.5,-12.5,-12.5),Subpart="Tail",SUVI=c(1:4))
		br_fiber_df<-data.frame(x_coord=c(1.05,3.75,6.45,6.20,3.75,1.30),y_coord=c(-12.5,-9,-12.5,-12.5,-9.3,-12.5),Subpart="BR_Fiber",SUVI=c(1:6))
		bl_fiber_df<-data.frame(x_coord=c(-6.45,-3.75,-1.05,-1.30,-3.75,-6.20),y_coord=c(-12.5,-9,-12.5,-12.5,-9.3,-12.5),Subpart="BL_Fiber",SUVI=c(1:6))
		ur_fiber_df<-data.frame(x_coord=c(1.05,3.75,6.95,6.70,3.75,1.30),y_coord=c(-12.5,-8,-12.5,-12.5,-8.3,-12.5),Subpart="UR_Fiber",SUVI=c(1:6))
		ul_fiber_df<-data.frame(x_coord=c(-6.95,-3.75,-1.05,-1.30,-3.75,-6.70),y_coord=c(-12.5,-8,-12.5,-12.5,-8.3,-12.5),Subpart="UL_Fiber",SUVI=c(1:6))
		baseplate_df<-data.frame(x_coord=c(-0.7,0.8,0.8,-0.7),y_coord=c(-12.5,-12.5,-13,-13),Subpart="Baseplate",SUVI=c(1:4))
		subparts_list<-list(capsid_df,ur_fiber_df,ul_fiber_df,br_fiber_df,bl_fiber_df,baseplate_df,tail_df)
	
	} else if (morphology=="Square") {
		capsid_df<-data.frame(x_coord=c(1,5,5,1),y_coord=c(5,5,1,1),Subpart="Capsid",SUVI=c("A","B","D","C"))
		subparts_list<-list(capsid_df)
	} else {
		print(paste(morphology," is not a valid morphology!",sep=""))
	}
	#subparts_list<-lapply(list(capsid_df,ur_fiber_df,ul_fiber_df,br_fiber_df,bl_fiber_df,baseplate_df,tail_df),FUN=get_rot_coord)#
	ent_df<-do.call(rbind,subparts_list)
	if (rot != 0) {
		ent_df<-get_rot_coord(ent_df=ent_df,rot=90)
	}
	ent_df$Morphology<-morphology
	ent_df$Entity<-ent_name
	ent_df$Col_Var<-col_var
	ent_df$Group<-paste(ent_df$Entity,ent_df$Subpart,sep="_")
	ent_df$UVI<-paste(ent_df$Subpart," ","v",ent_df$SUVI," ","(",round(ent_df$x_coord,digits=2),",",round(ent_df$y_coord,digits=2),")",sep="")
	return(ent_df)
}

build_sub_com<-function(x) {
	#return(rbind(build_entity(ent_name="Entity_1",col_var="Dummy_Colour_1"),build_entity(ent_name="Entity_2",col_var="Dummy_Colour_2")))
	sub_com_df<-build_entity(ent_name="Entity_2",col_var="Dummy_Colour_2",morphology="Myoviridae",rot=15)
	return(sub_com_df)
}

sub_com_df<-build_sub_com()

summary(sub_com_df)
head(sub_com_df)
table(sub_com_df$Subpart)

unique_subparts<-unique(sub_com_df$Subpart)
colours_num<-length(unique_subparts)
print(as.character(paste("Assigning ",colours_num," colours to: ",as.character(paste(as.vector(unique_subparts),sep=",",collapse=","))),sep="",collapse=","))
subpart_colours<-brewer.pal(9,"Set1")[c(1:colours_num)]
names(subpart_colours)<-unique_subparts

com_fig<-ggplot(sub_com_df, aes(x = x_coord, y = y_coord, group = Subpart))
com_fig<-com_fig+geom_polygon(alpha=0.7,colour="black",aes(fill=Subpart))
com_fig<-com_fig+theme_bw()
com_fig<-com_fig+geom_text(aes(label=UVI),position="identity",size=3, angle = 0)
com_fig<-com_fig+theme(legend.position="top")#"none"
#com_fig<-com_fig+scale_y_continuous(breaks=round(seq(min(sub_com_df$y_coord),max(sub_com_df$y_coord),.25),digits=2))+scale_x_continuous(breaks=round(seq(min(sub_com_df$x_coord),max(sub_com_df$x_coord),.5),digits=2))
com_fig<-com_fig+scale_fill_manual(name="Subpart",values=subpart_colours)
com_fig<-com_fig+geom_hline(yintercept = 0)+geom_vline(xintercept = 0)
#com_fig<-com_fig+facet_grid(Group ~ .,)
#svg(filename=out_fig_name,width=13.4,height=16,pointsize=8)
svg(filename=out_fig_name,width=12,height=12,pointsize=8)
#pdf("Virome_Community.pdf",width=10,height=20,pointsize=8)
print(com_fig)
dev.off()
