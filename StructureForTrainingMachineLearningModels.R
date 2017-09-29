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

setwd("/mnt/workspace_cluster_6/TRANSVERSAL_PROJECTS/MADR/COMPONENTE_2/USAID/")


###
#     ---- Reading Slope Eventes

VI="NDVI"
period="4_Months"

SlopeEvents=read.csv(paste0("ClimateComparison/VI_Signal/",VI,"/",period,"/IndividualValues",period,"_",VI,".csv"))

SlopeEvents=SlopeEvents[SlopeEvents$rsquared>=0.8,]
SlopeEvents=SlopeEvents[SlopeEvents$MeanBSI<=-0.28,]


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
  theme(text=element_text(size=10),
        axis.title.y=element_text(size = rel(1.2),colour = "#999999"),
        axis.title.x=element_text(size = rel(1.2),colour = "#888888"),
        axis.text.x  = element_text(angle=90, vjust=0.25,size = rel(1.3)),
        axis.text.y  = element_text(size = rel(1.5)),
        panel.background=element_rect(fill="white",colour = "black"),
        panel.grid.major = element_line(colour = "gray"))+
  stat_summary(fun.data = give.n, geom = "text", fun.y = median)


#### Joining Climate month

load(file="ClimateComparison/Climatology_VISignal/CoffeeV1V2_With_Climate_Delta_2000_2016_NewVar.RData")
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

DeltaClima[[1]]
VarClima=DeltaClima[[1]]
i=1
j=18
Lacks=0:4

dfAllSlopeAvRed=dfAllSlopeAvRed[!dfAllSlopeAvRed$X1=="2000_6",]
dfAllSlopeAvRed=dfAllSlopeAvRed[!dfAllSlopeAvRed$X2=="2017_1",]

lapply(Lacks ,function(Lack){
  library(snowfall)
  sfInit(parallel = T, cpus = 20)
  sfExport("dfAllSlopeAvRed")
  sfExport("DeltaClima")
  sfExport("VI")
  sfExport("period")
  sfLibrary(plyr)
  sfLibrary(stringr)
  MonthClimate_Lote=sfLapply(1:nrow(dfAllSlopeAvRed),function(j){
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
  sfStop()
  MonthClimate_Lote_df=as.data.frame(do.call(rbind,MonthClimate_Lote))
  
  ######
  
  dfAllSlopeAvRed$IDMatch=paste0(dfAllSlopeAvRed$ID,"t",
                                 do.call(rbind,lapply(str_split(dfAllSlopeAvRed$Period,pattern = " - "),function(x){x[2]})))
  
  names(MonthClimate_Lote_df)[1]="IDMatch"
  
  dfAllSlopeAvRed=plyr::join_all(list(dfAllSlopeAvRed,MonthClimate_Lote_df),by="IDMatch")
  
  write.csv(dfAllSlopeAvRed,file=paste0("ClimateComparison/Climatology_VISignal/",VI,"/DeltaClimeAssociation/DeltaValues",period,"_",VI,"_WithNewVar_Coffee_L",Lack,".csv"))
  
})


##################



SlopeEvents=read.csv("//dapadfs/workspace_cluster_6/TRANSVERSAL_PROJECTS/MADR/COMPONENTE_2/USAID/IndividualValues4Monthsper1Month2015_2016.csv",row.names = 1)
Lacks=0:6
Lack=1
ListValues2=lapply(Lacks, function(Lack){
  SlopeEvents=read.csv(paste0("ClimateComparison/DeltaValues4MonthsPer1Month_All_Coffee20152016_L",Lack,".csv"))
  
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

varNames1=c("prec")
prec=dfEventsSlopeClim[,grepl(names(dfEventsSlopeClim),pattern = varNames1)]
prec=cbind(ID=dfEventsSlopeClim$ID,prec,Slope=dfEventsSlopeClim$Slope_L4)


prec_graph=reshape(prec,idvar = "ID", varying = list(names(prec)[-c(1,ncol(prec))]),
                   times = str_sub(names(prec)[-c(1,ncol(prec))],-2),direction = "long")

prec_graph$Cat="Normal"

Percent95=quantile(prec_graph$Slope,c(0.05,0.95),na.rm=T)

prec_graph$Cat[prec_graph$Slope<=Percent95[1]]="Negative"
prec_graph$Cat[prec_graph$Slope>=Percent95[2]]="Positive"



ggplot(prec_graph,aes(prec_graph[,4]))+geom_density(alpha=0.6)+
  labs(x=varNames1)+
  theme(text=element_text(size=18),
        axis.title.y=element_text(size = rel(1.3),colour = "#999999"),
        axis.title.x=element_text(size = rel(1.4),colour = "#888888"),
        axis.text.x  = element_text(angle=0, hjust=0.5),
        panel.background=element_rect(fill="white",colour = "black"),
        panel.grid.major = element_line(colour = "gray"))+
  facet_wrap(~as.factor(time))
