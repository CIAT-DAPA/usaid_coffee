############################
#############   Making the training table for the machine Learning Models
#############



rm(list=ls())

###
#            -----   Loading packages
###

x11()

libs=c("stringr","ggplot2","dplyr","reshape2","scales","sp", "raster","Hmisc","rgdal")
lapply(libs,require,character.only=T)

give.n <- function(x){
  return(c(y = max(x)*1.40, label = length(x))) 
}


###
#            -----   Setting Workspace
###
server="/dapadfs"

setwd(paste0("/",server,"/workspace_cluster_6/TRANSVERSAL_PROJECTS/MADR/COMPONENTE_2/USAID/"))

VI="EVI"
period="4_Months"
###
#     ---- Reading Slope Eventes

SlopeEvents=read.csv(paste0("ClimateComparison/VI_Signal/",VI,"/",period,"/IndividualValues",period,"_",VI,".csv"))


dfAllSlopeAv=SlopeEvents
valMean=data.frame(summarise(group_by(dfAllSlopeAv,DateRef),
                             meanVal=median(Pendiente,na.rm=T)))


dfAllSlopeAv$Semester2=factor(dfAllSlopeAv$DateRef,
                              levels=as.character(valMean[order(valMean$meanVal,decreasing = T),][,1]))

EventsWLessthan50=names(table(dfAllSlopeAv$Semester2))[table(dfAllSlopeAv$Semester2)<=100]

dfAllSlopeAvRed=dfAllSlopeAv[!dfAllSlopeAv$Semester2%in%EventsWLessthan50,]

dfAllSlopeAvRed$Period=dfAllSlopeAvRed$Semester2
ggplot(dfAllSlopeAvRed,aes(Semester2,Pendiente))+geom_boxplot()+
  labs(y="Slope",x="")+
  theme(text=element_text(size=8),
        axis.title.y=element_text(size = rel(1.2),colour = "#999999"),
        axis.title.x=element_text(size = rel(1.2),colour = "#888888"),
        axis.text.x  = element_text(angle=90, vjust=0.25,size = rel(1.3)),
        axis.text.y  = element_text(size = rel(1.5)),
        panel.background=element_rect(fill="white",colour = "black"),
        panel.grid.major = element_line(colour = "gray"))+
  stat_summary(fun.data = give.n, geom = "text", fun.y = median)

ggsave(paste0("ClimateComparison/VI_Signal/",VI,"/",period,"/VI_ThroughTime_Boxplot.png"),width = 26, height = 12, units = "cm")


#### Joining Climate month


load(file="ClimateComparison/CoffeeV1V2_With_Climate_Delta_2000_2016.RData")
head(dfAllSlopeAvRed)
head(DeltaClima[[1]])
x=str_split(dfAllSlopeAvRed$Period,pattern = " - ")[[1]]
y=str_split(x,pattern = "_")[[1]]
PeriodNum=lapply(str_split(dfAllSlopeAvRed$Period,pattern = " - "),function(x){
  Periods=str_split(x,pattern = "_")                                    
  Month1=as.numeric(which(str_sub(month.name,1,3)%in%Periods[[1]][length(Periods[[1]])]))
  Month2=as.numeric(which(str_sub(month.name,1,3)%in%Periods[[2]][length(Periods[[2]])]))
  Year1=as.numeric(Periods[[1]][1])
  Year2=as.numeric(Periods[[2]][1])
  return(cbind(paste0(Year1,"_",Month1),paste0(Year2,"_",Month2)))})

dfAllSlopeAvRed=cbind(dfAllSlopeAvRed,data.frame(do.call(rbind,PeriodNum)))

VarClima=DeltaClima[[1]]
i=1
j=18
Lack=5
Lacks=0:6
dfAllSlopeAvRed=dfAllSlopeAvRed[!dfAllSlopeAvRed$X1=="2000_6",]
dfAllSlopeAvRed=dfAllSlopeAvRed[!dfAllSlopeAvRed$X2=="2017_1",]
library(snowfall)
sfInit(parallel = T, cpus = 7)
sfExport("dfAllSlopeAvRed")
sfExport("DeltaClima")
sfExport("VI")
sfExport("period")

sfLibrary(plyr)
sfLibrary(stringr)
sfLapply(Lacks ,function(Lack){
  
  MonthClimate_Lote=lapply(1:nrow(dfAllSlopeAvRed),function(j){
    return(plyr::join_all(lapply(1:length(DeltaClima),function(i){
      cat(paste0(i,"_",j,"\n"))
      VarClima=DeltaClima[[i]]
      Start=which(names(VarClima)%in%as.character(dfAllSlopeAvRed$X1[j]))
      End=which(names(VarClima)%in%as.character(dfAllSlopeAvRed$X2[j]))
      IDVal=paste0(dfAllSlopeAvRed[j,"ID"],"t",
                   do.call(rbind,lapply(str_split(dfAllSlopeAvRed[j,]$Period,pattern = " - "),function(x){x[2]})))
      dfAux=cbind(IDVal,VarClima[VarClima$ID%in%dfAllSlopeAvRed$ID[j],names(VarClima)[(Start-Lack):(End-Lack)]])
      names(dfAux)=c("ID",paste0(names(DeltaClima)[i],"_",1:length(Start:End)))
      return(dfAux)
    }),by="ID"))
  })
  MonthClimate_Lote_df=as.data.frame(do.call(rbind,MonthClimate_Lote))
  
  ######
  
  dfAllSlopeAvRed$IDMatch=paste0(dfAllSlopeAvRed$ID,"t",
                                 do.call(rbind,lapply(str_split(dfAllSlopeAvRed$Period,pattern = " - "),function(x){x[2]})))
  
  names(MonthClimate_Lote_df)[1]="IDMatch"
  
  dfAllSlopeAvRed=plyr::join_all(list(dfAllSlopeAvRed,MonthClimate_Lote_df),by="IDMatch")
  write.csv(dfAllSlopeAvRed,file=paste0("ClimateComparison/Climatology_VISignal/",VI,"/DeltaClimeAssociation/DeltaValues",period,"_",VI,"_Coffee_L",Lack,".csv"))
  
})

sfStop()
####


############################
#############   
#############



rm(list=ls())

###
#            -----   Loading packages
###



libs=c("stringr","ggplot2","dplyr","reshape2","scales","sp", "raster","Hmisc","rgdal")
lapply(libs,require,character.only=T)

give.n <- function(x){
  return(c(y = max(x)*1.40, label = length(x))) 
}


###
#            -----   Setting Workspace
###


#setwd("D:/ToBackup/Documents/Projects/USAID/")
setwd("//dapadfs/workspace_cluster_6/TRANSVERSAL_PROJECTS/MADR/COMPONENTE_2/USAID/")

###
#     ---- Uploading the geo-coordenates reference

#CoffeePoints=read.csv("Analisis_V2/FusionTable/CSV/Coffee_Points_P1_P2_FT.csv")
CoffeePoints=read.csv("FusionTable/CSV/Coffee_Points_2015_2016Match.csv")

xy=cbind( CoffeePoints$Longitud, CoffeePoints$Latitud)
CoffeePoints_spRef = SpatialPointsDataFrame(coords = xy,data=CoffeePoints,
                                            proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))


###
#     ---- Reading Slope Eventes
Lacks=0:4
Lack=5
VI="NDVI"
MonthsN="4_Months"
ListValues2=lapply(Lacks, function(Lack){
  SlopeEvents=read.csv(paste0("ClimateComparison/Climatology_VISignal/",VI,"/DeltaClimeAssociation/DeltaValues",MonthsN,"_",VI,"_Coffee_L",Lack,".csv"))
  
  ID="IDMatch"
  Slope="Pendiente"
  
  SlopeEvents=SlopeEvents[!is.na(SlopeEvents$prec_1),]
  varNames1=c("prec")
  prec=apply(SlopeEvents[,grepl(names(SlopeEvents),pattern = varNames1)],1,mean,na.rm=T)
  varNames2=c("dtr")
  dtr=apply(SlopeEvents[,grepl(names(SlopeEvents),pattern = varNames2)],1,mean,na.rm=T)
  varNames3=c("tmean")
  tmean=apply(SlopeEvents[,grepl(names(SlopeEvents),pattern = varNames3)],1,mean,na.rm=T)
  varNames4=c("tmax")
  tmax=apply(SlopeEvents[,grepl(names(SlopeEvents),pattern = varNames4)],1,mean,na.rm=T)
  varNames5=c("tmin")
  tmin=apply(SlopeEvents[,grepl(names(SlopeEvents),pattern = varNames5)],1,mean,na.rm=T)
  
  SlopeDeltaDF=data.frame(ID=SlopeEvents[,ID],prec,dtr,tmean,tmin,tmax,Slope=SlopeEvents[,Slope])
  
  names(SlopeDeltaDF)[-c(1)]=paste0(names(SlopeDeltaDF)[-c(1)],"_L",Lack)
  return(SlopeDeltaDF)
})


dfEventsSlopeClim=plyr::join_all(ListValues2,by = "ID")

Percent95=quantile(dfEventsSlopeClim$Slope_L0,c(0.05,0.95),na.rm=T)
Percent975=quantile(dfEventsSlopeClim$Slope_L0,c(0.025,0.975),na.rm=T)

period="2014_Jul_2016_Oct"
setwd("//dapadfs/workspace_cluster_6/TRANSVERSAL_PROJECTS/MADR/COMPONENTE_2/USAID/")
dir.create(paste0("ClimateComparison/Climatology_VISignal/",VI,"/Exploratory_Graphs/",MonthsN,"/",period),showWarnings = F)
setwd(paste0("ClimateComparison/Climatology_VISignal/",VI,"/Exploratory_Graphs/",MonthsN,"/",period))

dateEnd=str_sub(dfEventsSlopeClim$ID,-8)
Year=str_sub(dateEnd,1,4)

dfEventsSlopeClim_Red_F1=dfEventsSlopeClim[as.numeric(as.character(Year))%in%c(2001,2004,2007,2010,2014,2015:2016),]


dateEnd_F2=str_sub(dfEventsSlopeClim_Red_F1$ID,-8)
Year_F2=str_sub(dateEnd_F2,1,4)
True2014=as.numeric(as.character(Year_F2))%in%2014 & dateEnd_F2%in% c("2014_Oct","2014_Nov","2014_Dec")
True2016=(as.numeric(as.character(Year_F2))%in%2016 & !dateEnd_F2%in% c("2016_Nov","2016_Dec"))
True2015=as.numeric(as.character(Year_F2))%in%2015
True2007=as.numeric(as.character(Year_F2))%in%2007
True2001=as.numeric(as.character(Year_F2))%in%2001
True2004=as.numeric(as.character(Year_F2))%in%2004
True2010=as.numeric(as.character(Year_F2))%in%2010


dfEventsSlopeClim_Red_F2=dfEventsSlopeClim_Red_F1[ True2014| True2015|  True2016 ,]
dfEventsSlopeClim=dfEventsSlopeClim_Red_F2

levels(factor(str_sub(dfEventsSlopeClim$ID,-8)))
varNames_All=c("prec","dtr","tmean","tmax","tmin")
varNames1="prec"
lapply(varNames_All,function(varNames1){
  prec=dfEventsSlopeClim[,grepl(names(dfEventsSlopeClim),pattern = varNames1)]
  prec=cbind(ID=dfEventsSlopeClim$ID,prec,Slope=dfEventsSlopeClim$Slope_L4)
  
  
  prec_graph=reshape(prec,idvar = "ID", varying = list(names(prec)[-c(1,ncol(prec))]),timevar = "Slack",
                     times = str_sub(names(prec)[-c(1,ncol(prec))],-2),direction = "long")
  
  prec_graph$Cat="Normal"
  
  
  
  prec_graph$Cat[prec_graph$Slope<=Percent95[1]]="Negative"
  prec_graph$Cat[prec_graph$Slope>=Percent95[2]]="Positive"
  
  ggplot(prec_graph,aes(prec_graph[,4]))+geom_density(alpha=0.6)+
    labs(x=varNames1,fill="Slope")+geom_vline(xintercept = 0,colour="coral3")+
    theme(text=element_text(size=14),
          axis.title.y=element_text(size = rel(1.3),colour = "#999999"),
          axis.title.x=element_text(size = rel(1.4),colour = "#888888"),
          axis.text.x  = element_text(angle=0, hjust=0.5),
          panel.background=element_rect(fill="white",colour = "black"),
          panel.grid.major = element_line(colour = "gray"))+
    facet_grid(Slack~.)
  ggsave(paste0("0_",VI,"_Density_",varNames1,"_5Months.png"),width = 18, height = 16, units = "cm")
  
  ggplot(prec_graph,aes(prec_graph[,4],Slope,colour=Cat))+geom_point(alpha=0.6)+
    labs(x=varNames1,fill="Slope")+geom_vline(xintercept = 0,colour="coral3")+
    geom_hline(yintercept = 0,colour="coral3")+
    theme(text=element_text(size=14),
          axis.title.y=element_text(size = rel(1.3),colour = "#999999"),
          axis.title.x=element_text(size = rel(1.4),colour = "#888888"),
          axis.text.x  = element_text(angle=0, hjust=0.5),
          panel.background=element_rect(fill="white",colour = "black"),
          panel.grid.major = element_line(colour = "gray"))+
    facet_grid(Slack~.)
  
  ggsave(paste0("0_",VI,"_Point_",varNames1,"_5Months.png"),width = 20, height = 16, units = "cm")
  
  
  ggplot(prec_graph,aes(prec_graph[,4]))+geom_density(aes(fill=Cat),alpha=0.6)+
    labs(x=varNames1,fill="Slope")+geom_vline(xintercept = 0,colour="coral3")+
    theme(text=element_text(size=14),
          axis.title.y=element_text(size = rel(1.3),colour = "#999999"),
          axis.title.x=element_text(size = rel(1.4),colour = "#888888"),
          axis.text.x  = element_text(angle=0, hjust=0.5),
          panel.background=element_rect(fill="white",colour = "black"),
          panel.grid.major = element_line(colour = "gray"))+
    facet_grid(Slack~.)
  
  ggsave(paste0("1_",VI,"_Density_",varNames1,"_5Months.png"),width = 20, height = 16, units = "cm")
  ggplot(prec_graph,aes(prec_graph[,4],fill=Cat))+geom_histogram(aes(y = ..count..),alpha=0.6,
                                                                 binwidth = 1,position="identity")+
    labs(x=varNames1,fill="Slope")+geom_vline(xintercept = 0,colour="coral3")+
    theme(text=element_text(size=14),
          axis.title.y=element_text(size = rel(1.3),colour = "#999999"),
          axis.title.x=element_text(size = rel(1.4),colour = "#888888"),
          axis.text.x  = element_text(angle=0, hjust=0.5),
          panel.background=element_rect(fill="white",colour = "black"),
          panel.grid.major = element_line(colour = "gray"))+
    facet_grid(Slack~.)
  
  ggsave(paste0("1_",VI,"_histogram_",varNames1,"_5Months.png"),width = 20, height = 16, units = "cm")
  
})

lapply(varNames_All,function(varNames1){
  prec=dfEventsSlopeClim[,grepl(names(dfEventsSlopeClim),pattern = varNames1)]
  prec=cbind(ID=dfEventsSlopeClim$ID,prec,Slope=dfEventsSlopeClim$Slope_L4)
  
  
  prec_graph=reshape(prec,idvar = "ID", varying = list(names(prec)[-c(1,ncol(prec))]),
                     times = str_sub(names(prec)[-c(1,ncol(prec))],-2),direction = "long")
  
  prec_graph$Cat="Normal"
  
  
  prec_graph$Cat[prec_graph$Slope<=Percent95[1]]="Negative"
  prec_graph$Cat[prec_graph$Slope>=Percent95[2]]="Positive"
  
  
  prec_graph$Cat[prec_graph$Slope<=Percent975[1]]="Extreme_Negative"
  prec_graph$Cat[prec_graph$Slope>=Percent975[2]]="Extreme_Positive"
  prec_graph=prec_graph[!prec_graph$Cat=="Normal",]
  ggplot(prec_graph,aes(prec_graph[,4]))+geom_density(aes(fill=Cat),alpha=0.6)+
    labs(x=varNames1)+geom_vline(xintercept = 0,colour="coral3")+
    theme(text=element_text(size=14),
          axis.title.y=element_text(size = rel(1.3),colour = "#999999"),
          axis.title.x=element_text(size = rel(1.4),colour = "#888888"),
          axis.text.x  = element_text(angle=0, hjust=0.5),
          panel.background=element_rect(fill="white",colour = "black"),
          panel.grid.major = element_line(colour = "gray"))+
    facet_grid(time~.)
  
  ggsave(paste0("2_",VI,"_Density_Extreme_",varNames1,"_5Months.png"),width = 20, height = 16, units = "cm")
  ggplot(prec_graph,aes(prec_graph[,4],fill=Cat))+geom_histogram(aes(y = ..count..),alpha=0.6,
                                                                 binwidth = 1,position="identity")+
    labs(x=varNames1,fill="Slope")+geom_vline(xintercept = 0,colour="coral3")+
    theme(text=element_text(size=14),
          axis.title.y=element_text(size = rel(1.3),colour = "#999999"),
          axis.title.x=element_text(size = rel(1.4),colour = "#888888"),
          axis.text.x  = element_text(angle=0, hjust=0.5),
          panel.background=element_rect(fill="white",colour = "black"),
          panel.grid.major = element_line(colour = "gray"))+
    facet_grid(time~.)
  
  ggsave(paste0("2_",VI,"_histogram_Extreme_",varNames1,"_5Months.png"),width = 20, height = 16, units = "cm")
  
})
varAux="tmean"
prec=dfEventsSlopeClim[,grepl(names(dfEventsSlopeClim),pattern = varAux)]
prec=cbind(ID=dfEventsSlopeClim$ID,prec,Slope=dfEventsSlopeClim$Slope_L4)


prec_graph=reshape(prec,idvar = "ID", varying = list(names(prec)[-c(1,ncol(prec))]),
                   times = str_sub(names(prec)[-c(1,ncol(prec))],-2),direction = "long")
prec_graph$Cat="Normal"

prec_graph$Cat[prec_graph$Slope<=Percent95[1]]="Negative"
prec_graph$Cat[prec_graph$Slope>=Percent95[2]]="Positive"
prec_graph=prec_graph[!prec_graph$Cat=="Normal",]
Cat2=array(dim = length(prec_graph$Slope))
Cat2[prec_graph[,4]>5]=paste0(" >5 ", varAux)
Cat2[prec_graph[,4]<=5]=paste0(" <=5 ", varAux)
varNames1="prec"


prec=dfEventsSlopeClim[,grepl(names(dfEventsSlopeClim),pattern = varNames1)]
prec=cbind(ID=dfEventsSlopeClim$ID,prec,Slope=dfEventsSlopeClim$Slope_L4)
prec_graph=reshape(prec,idvar = "ID", varying = list(names(prec)[-c(1,ncol(prec))]),
                   times = str_sub(names(prec)[-c(1,ncol(prec))],-2),direction = "long")
prec_graph$Cat="Normal"

prec_graph$Cat[prec_graph$Slope<=Percent95[1]]="Negative"
prec_graph$Cat[prec_graph$Slope>=Percent95[2]]="Positive"
prec_graph=prec_graph[!prec_graph$Cat=="Normal",]
prec_graph$Cat3=paste0(prec_graph$Cat,"_",Cat2)
ggplot(prec_graph,aes(prec_graph[,4]))+geom_density(aes(fill=Cat3),alpha=0.6)+
  labs(x=varNames1)+geom_vline(xintercept = 0,colour="coral3")+
  theme(text=element_text(size=14),
        axis.title.y=element_text(size = rel(1.3),colour = "#999999"),
        axis.title.x=element_text(size = rel(1.4),colour = "#888888"),
        axis.text.x  = element_text(angle=0, hjust=0.5),
        panel.background=element_rect(fill="white",colour = "black"),
        panel.grid.major = element_line(colour = "gray"))+
  facet_grid(time~.)

ggsave(paste0("3_",VI,"_Density_",varNames1,"_",varAux,"_5Months.png"),width = 20, height = 16, units = "cm")




prec_graph$Cat[prec_graph$Slope<=Percent975[1]]="Extreme_Negative"
prec_graph$Cat[prec_graph$Slope>=Percent975[2]]="Extreme_Positive"
prec_graph=prec_graph[!prec_graph$Cat=="Normal",]
prec_graph$Cat3=paste0(prec_graph$Cat,"_",Cat2)

ggplot(prec_graph,aes(prec_graph[,4]))+geom_density(aes(fill=Cat3),alpha=0.6)+
  labs(x=varNames1)+geom_vline(xintercept = 0,colour="coral3")+
  theme(text=element_text(size=14),
        axis.title.y=element_text(size = rel(1.3),colour = "#999999"),
        axis.title.x=element_text(size = rel(1.4),colour = "#888888"),
        axis.text.x  = element_text(angle=0, hjust=0.5),
        panel.background=element_rect(fill="white",colour = "black"),
        panel.grid.major = element_line(colour = "gray"))+
  facet_grid(time~.)

ggsave(paste0("3_",VI,"_Density_Extreme_",varNames1,"_",varAux,"_5Months.png"),width = 20, height = 16, units = "cm")


###########Tmax

varAux="tmax"
prec=dfEventsSlopeClim[,grepl(names(dfEventsSlopeClim),pattern = varAux)]
prec=cbind(ID=dfEventsSlopeClim$ID,prec,Slope=dfEventsSlopeClim$Slope_L4)


prec_graph=reshape(prec,idvar = "ID", varying = list(names(prec)[-c(1,ncol(prec))]),
                   times = str_sub(names(prec)[-c(1,ncol(prec))],-2),direction = "long")
prec_graph$Cat="Normal"

prec_graph$Cat[prec_graph$Slope<=Percent95[1]]="Negative"
prec_graph$Cat[prec_graph$Slope>=Percent95[2]]="Positive"
prec_graph=prec_graph[!prec_graph$Cat=="Normal",]
Cat2=array(dim = length(prec_graph$Slope))
Cat2[prec_graph[,4]>8]=paste0(" >8 ", varAux)
Cat2[prec_graph[,4]<=8]=paste0(" <=8 ", varAux)
varNames1="prec"


prec=dfEventsSlopeClim[,grepl(names(dfEventsSlopeClim),pattern = varNames1)]
prec=cbind(ID=dfEventsSlopeClim$ID,prec,Slope=dfEventsSlopeClim$Slope_L4)
prec_graph=reshape(prec,idvar = "ID", varying = list(names(prec)[-c(1,ncol(prec))]),
                   times = str_sub(names(prec)[-c(1,ncol(prec))],-2),direction = "long")
prec_graph$Cat="Normal"

prec_graph$Cat[prec_graph$Slope<=Percent95[1]]="Negative"
prec_graph$Cat[prec_graph$Slope>=Percent95[2]]="Positive"
prec_graph=prec_graph[!prec_graph$Cat=="Normal",]
prec_graph$Cat3=paste0(prec_graph$Cat,"_",Cat2)
ggplot(prec_graph,aes(prec_graph[,4]))+geom_density(aes(fill=Cat3),alpha=0.6)+
  labs(x=varNames1)+geom_vline(xintercept = 0,colour="coral3")+
  theme(text=element_text(size=14),
        axis.title.y=element_text(size = rel(1.3),colour = "#999999"),
        axis.title.x=element_text(size = rel(1.4),colour = "#888888"),
        axis.text.x  = element_text(angle=0, hjust=0.5),
        panel.background=element_rect(fill="white",colour = "black"),
        panel.grid.major = element_line(colour = "gray"))+
  facet_grid(time~.)

ggsave(paste0("3_",VI,"_Density_",varNames1,"_",varAux,"_5Months.png"),width = 20, height = 16, units = "cm")




prec_graph$Cat[prec_graph$Slope<=Percent975[1]]="Extreme_Negative"
prec_graph$Cat[prec_graph$Slope>=Percent975[2]]="Extreme_Positive"
prec_graph=prec_graph[!prec_graph$Cat=="Normal",]
prec_graph$Cat3=paste0(prec_graph$Cat,"_",Cat2)

ggplot(prec_graph,aes(prec_graph[,4]))+geom_density(aes(fill=Cat3),alpha=0.6)+
  labs(x=varNames1)+geom_vline(xintercept = 0,colour="coral3")+
  theme(text=element_text(size=14),
        axis.title.y=element_text(size = rel(1.3),colour = "#999999"),
        axis.title.x=element_text(size = rel(1.4),colour = "#888888"),
        axis.text.x  = element_text(angle=0, hjust=0.5),
        panel.background=element_rect(fill="white",colour = "black"),
        panel.grid.major = element_line(colour = "gray"))+
  facet_grid(time~.)

ggsave(paste0("3_",VI,"_Density_Extreme_",varNames1,"_",varAux,"_5Months.png"),width = 20, height = 16, units = "cm")


###########dtr

varAux="dtr"
prec=dfEventsSlopeClim[,grepl(names(dfEventsSlopeClim),pattern = varAux)]
prec=cbind(ID=dfEventsSlopeClim$ID,prec,Slope=dfEventsSlopeClim$Slope_L4)


prec_graph=reshape(prec,idvar = "ID", varying = list(names(prec)[-c(1,ncol(prec))]),
                   times = str_sub(names(prec)[-c(1,ncol(prec))],-2),direction = "long")
prec_graph$Cat="Normal"

prec_graph$Cat[prec_graph$Slope<=Percent95[1]]="Negative"
prec_graph$Cat[prec_graph$Slope>=Percent95[2]]="Positive"
prec_graph=prec_graph[!prec_graph$Cat=="Normal",]
Cat2=array(dim = length(prec_graph$Slope))
Cat2[prec_graph[,4]>4]=paste0(" >4 ", varAux)
Cat2[prec_graph[,4]<=4]=paste0(" <=4 ", varAux)
varNames1="prec"


prec=dfEventsSlopeClim[,grepl(names(dfEventsSlopeClim),pattern = varNames1)]
prec=cbind(ID=dfEventsSlopeClim$ID,prec,Slope=dfEventsSlopeClim$Slope_L4)
prec_graph=reshape(prec,idvar = "ID", varying = list(names(prec)[-c(1,ncol(prec))]),
                   times = str_sub(names(prec)[-c(1,ncol(prec))],-2),direction = "long")
prec_graph$Cat="Normal"

prec_graph$Cat[prec_graph$Slope<=Percent95[1]]="Negative"
prec_graph$Cat[prec_graph$Slope>=Percent95[2]]="Positive"
prec_graph=prec_graph[!prec_graph$Cat=="Normal",]
prec_graph$Cat3=paste0(prec_graph$Cat,"_",Cat2)
ggplot(prec_graph,aes(prec_graph[,4]))+geom_density(aes(fill=Cat3),alpha=0.6)+
  labs(x=varNames1)+geom_vline(xintercept = 0,colour="coral3")+
  theme(text=element_text(size=14),
        axis.title.y=element_text(size = rel(1.3),colour = "#999999"),
        axis.title.x=element_text(size = rel(1.4),colour = "#888888"),
        axis.text.x  = element_text(angle=0, hjust=0.5),
        panel.background=element_rect(fill="white",colour = "black"),
        panel.grid.major = element_line(colour = "gray"))+
  facet_grid(time~.)

ggsave(paste0("3_",VI,"_Density_",varNames1,"_",varAux,"_5Months.png"),width = 20, height = 16, units = "cm")




prec_graph$Cat[prec_graph$Slope<=Percent975[1]]="Extreme_Negative"
prec_graph$Cat[prec_graph$Slope>=Percent975[2]]="Extreme_Positive"
prec_graph=prec_graph[!prec_graph$Cat=="Normal",]
prec_graph$Cat3=paste0(prec_graph$Cat,"_",Cat2)

ggplot(prec_graph,aes(prec_graph[,4]))+geom_density(aes(fill=Cat3),alpha=0.6)+
  labs(x=varNames1)+geom_vline(xintercept = 0,colour="coral3")+
  theme(text=element_text(size=14),
        axis.title.y=element_text(size = rel(1.3),colour = "#999999"),
        axis.title.x=element_text(size = rel(1.4),colour = "#888888"),
        axis.text.x  = element_text(angle=0, hjust=0.5),
        panel.background=element_rect(fill="white",colour = "black"),
        panel.grid.major = element_line(colour = "gray"))+
  facet_grid(time~.)

ggsave(paste0("3_",VI,"_Density_Extreme_",varNames1,"_",varAux,"_5Months.png"),width = 20, height = 16, units = "cm")

######-----------------------------------------------
###########prec

varAux="prec"
prec=dfEventsSlopeClim[,grepl(names(dfEventsSlopeClim),pattern = varAux)]
prec=cbind(ID=dfEventsSlopeClim$ID,prec,Slope=dfEventsSlopeClim$Slope_L4)


prec_graph=reshape(prec,idvar = "ID", varying = list(names(prec)[-c(1,ncol(prec))]),
                   times = str_sub(names(prec)[-c(1,ncol(prec))],-2),direction = "long")
prec_graph$Cat="Normal"

prec_graph$Cat[prec_graph$Slope<=Percent95[1]]="Negative"
prec_graph$Cat[prec_graph$Slope>=Percent95[2]]="Positive"
#prec_graph=prec_graph[!prec_graph$Cat=="Normal",]
Cat2=array(dim = length(prec_graph$Slope))
Cat2[prec_graph[,4]>-20]=paste0(" >20 ", varAux)
Cat2[prec_graph[,4]<=-20]=paste0(" <=20 ", varAux)
varNames1="dtr"


prec=dfEventsSlopeClim[,grepl(names(dfEventsSlopeClim),pattern = varNames1)]
prec=cbind(ID=dfEventsSlopeClim$ID,prec,Slope=dfEventsSlopeClim$Slope_L4)
prec_graph=reshape(prec,idvar = "ID", varying = list(names(prec)[-c(1,ncol(prec))]),
                   times = str_sub(names(prec)[-c(1,ncol(prec))],-2),direction = "long")
prec_graph$Cat="Normal"

prec_graph$Cat[prec_graph$Slope<=Percent95[1]]="Negative"
prec_graph$Cat[prec_graph$Slope>=Percent95[2]]="Positive"
#prec_graph=prec_graph[!prec_graph$Cat=="Normal",]
prec_graph$Cat3=paste0(prec_graph$Cat,"_",Cat2)
ggplot(prec_graph,aes(prec_graph[,4]))+facet_grid(time~.)+
  geom_density(aes(fill=Cat3),alpha=0.6)+
  labs(x=varNames1)+geom_vline(xintercept = 0,colour="coral3")+
  geom_vline(xintercept = 3,colour="coral3")+
  theme(text=element_text(size=14),
        axis.title.y=element_text(size = rel(1.3),colour = "#999999"),
        axis.title.x=element_text(size = rel(1.4),colour = "#888888"),
        axis.text.x  = element_text(angle=0, hjust=0.5),
        panel.background=element_rect(fill="white",colour = "black"),
        panel.grid.major = element_line(colour = "gray"))


ggsave(paste0("3_",VI,"_Density_",varNames1,"_",varAux,"_5Months.png"),width = 20, height = 16, units = "cm")
