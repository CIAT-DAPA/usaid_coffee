
# -------------> Counting the anomalies 



rm(list=ls())

###
#            -----   Loading packages
###

libs=c("maptools","rgdal","raster","sp","ggplot2","stringr","reshape2","ff")
sapply(libs,require,character.only=T)


CleaningTifNames= function(DirFol,Year=""){
  setwd(DirFol)
  filesIn=list.files(pattern="*.tif")
  filesIn=filesIn[grepl(filesIn,pattern = Year)]
  rx = gregexpr(".tif.*.",filesIn, perl = TRUE)
  valMTif=unique(str_sub(unlist(regmatches(filesIn, rx)),1))
  if(length(valMTif)>0){
    filesIn=filesIn[!grepl(valMTif,filesIn)]
  }
  return(filesIn)
}


ChangingPeriodAnotation=function(PeriodRef,divProg,PatternDiv="_"){
  Year1=str_sub(PeriodRef[1],1,4)
  MonthVal1=month.abb[((as.numeric(unlist(str_split(PeriodRef[1],PatternDiv))[2])*(12/divProg))-((12/divProg)-1))]
  MonthVal2=month.abb[((as.numeric(unlist(str_split(PeriodRef[length(PeriodRef)],PatternDiv))[2])*(12/divProg)))]
  Year2=str_sub(PeriodRef[length(PeriodRef)],1,4)
  return(paste0(Year1,"_",MonthVal1," - ",Year2,"_",MonthVal2))
}

dirFol="//dapadfs/workspace_cluster_6/TRANSVERSAL_PROJECTS/MADR/COMPONENTE_2/USAID/"
setwd(dirFol)

VI="NDVI"

CoffeePoints=read.csv("FusionTable/CSV/Coffee_Points_2015_2016Match.csv")
SlopeCat=read.csv(paste0("ClimateComparison/VI_Signal/",VI,"/4_Months/IndividualValues4_Months_",VI,".csv"))
head(SlopeCat)

SlopeCat$Lote=SlopeCat$ID
SlopeCat=plyr::join_all(list(SlopeCat,CoffeePoints[,1:3]),by="Lote")
SlopeCat$StartDate=as.Date(sapply(str_split(SlopeCat$DateRef,pattern = " - "),function(x){paste0(x[1],"_01")}),format="%Y_%b_%d")
SlopeCat$EndDate=as.Date(sapply(str_split(SlopeCat$DateRef,pattern = " - "),function(x){paste0(x[2],"_28")}),format="%Y_%b_%d")
RasterDirAnomalies="//dapadfs/workspace_cluster_6/TRANSVERSAL_PROJECTS/MADR/COMPONENTE_2/USAID/ClimateComparison/Climatology/Anomalies _Alt08_085/"
RasterDirAnomalies=paste(RasterDirAnomalies,CleaningTifNames(RasterDirAnomalies),sep="/")
SlopeCat$Pendiente


Percent80=quantile(SlopeCat$Pendiente,c(0.2,0.8),na.rm=T)
ggplot(SlopeCat, aes(Pendiente))+geom_density()+geom_vline(xintercept = Percent80)

SlopeCat$Extremo="Normal"
SlopeCat$Extremo[SlopeCat$Pendiente<=Percent80[1]]="Negative"
SlopeCat$Extremo[SlopeCat$Pendiente>=Percent80[2]]="Positive"


Proximity=c("Low","Up")

PosPattern=gregexpr("085/.*._20",RasterDirAnomalies,perl = T)
ClimateVariables=unique(str_sub(unlist(regmatches(RasterDirAnomalies,PosPattern)),6,-4))

ClimateVar=ClimateVariables[1]

x=Proximity[1]
ListCountEvents=lapply(Proximity,function(x){
  SymbolComp=switch(x,"Low"=">=","Up"="<")
  cat("Starting Multi values to points process: \n")
  NamesRasterLimits=RasterDirAnomalies[grepl(RasterDirAnomalies,pattern = x)]
  CountPerVar=lapply(ClimateVariables,function(ClimateVar){
    NamesRasterVar=NamesRasterLimits[grepl(NamesRasterLimits,pattern = paste0(ClimateVar,"_20"))]
    
    Sys.time()->start
    library(snowfall)
    ncores=20
    sfInit(parallel=T,cpus=ncores)
    sfLibrary(raster)
    sfLibrary(stringr)
    sfExport("SlopeCat")
    sfExport("x")
    RasterFiles=sfLapply(NamesRasterVar,function(y){
      m=raster(y)
      crs(m) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
      rx = gregexpr(paste0("20.*.",str_sub(x,1,2)),y, perl = TRUE)
      yearMonth=unique(str_sub(unlist(regmatches(y, rx)),1))
      yearMonth=str_sub(yearMonth,1,-4)
      year=str_sub(yearMonth,1,4)
      month=str_sub(yearMonth,6)
      DateImage=as.Date(paste0(year,"-",month,"-10"),format= "%Y-%m-%d")
      SlopeDFRed=SlopeCat[as.numeric(SlopeCat$StartDate) < as.numeric(DateImage) & as.numeric(DateImage) < as.numeric(SlopeCat$EndDate),]
      SlopeDFRed2=SlopeDFRed[SlopeDFRed$Extremo=="Normal",]
      xy=cbind( SlopeDFRed2$Longitud, SlopeDFRed2$Latitud)
      plot(m)
      m[!1:length(m)%in%cellFromXY(m,xy)]=NA
      plot(m)
      return(m)})
    sfStop()
    
    PosPattern=gregexpr(paste0("",str_sub(ClimateVar,-2),".*_",str_sub(x,1,1)),NamesRasterVar,perl = T)
    NamesDates=str_sub(unlist(regmatches(NamesRasterVar,PosPattern)),4,-3)
    names(RasterFiles)=NamesDates
    NameList=names(RasterFiles)[3]
    PixelsWAnom=sapply(names(RasterFiles), function(NameList){
      length(which(eval(parse(text = paste0("!RasterFiles[[NameList]][]" ,SymbolComp, "0")))))
    })
    print(Sys.time()-start)
    YearMonth=sapply(strsplit(names(PixelsWAnom),"_"),function(x){t(x)})
    Year=YearMonth[1,]
    Month=YearMonth[2,]
    YearMonth=data.frame(Year,Month=as.numeric(Month))
    YearMonth=YearMonth[order(YearMonth$Month),]
    YearMonth=YearMonth[order(as.numeric(YearMonth$Year)),]
    NamesDateRef=paste0(YearMonth[,1], "_",YearMonth[,2] )
    PixelsWAnom=PixelsWAnom[NamesDateRef]
    return(PixelsWAnom)
    
  })
  names(CountPerVar)=ClimateVariables
  dfCounting=as.data.frame(t(do.call(rbind,CountPerVar)))
  dfCounting$Period=row.names(dfCounting)
  return(dfCounting)
})

names(ListCountEvents)=Proximity

dirFol="//dapadfs/workspace_cluster_6/TRANSVERSAL_PROJECTS/MADR/COMPONENTE_2/USAID/"
setwd(dirFol)

save(ListCountEvents,file=paste0("ClimateComparison/Climatology_VISignal/",VI,"/CountingAnomalies/",VI,"_CountingAlternative_Normal0208_080_085.RData"))
load(file=paste0("ClimateComparison/Climatology_VISignal/",VI,"/CountingAnomalies/",VI,"_CountingAlternative_Positive_080_085.RData"))

divPer=12

index=Proximity[1]
CountingPeriod=lapply(Proximity,function(index){
  DFIndex=ListCountEvents[[index]]
  month=as.numeric(str_sub(DFIndex$Period,6))
  year=as.numeric(str_sub(DFIndex$Period,1,4))
  Period=ceiling(month/(12/divPer))
  PeriodAccum=((year-min(year))*divPer)+Period
  LoopPeriods=unique(PeriodAccum[order(PeriodAccum)])
  dfExport=data.frame()
  dfDateRef=array()
  for(i in LoopPeriods[-c((length(LoopPeriods)-2):length(LoopPeriods))]){
    SelectRows=PeriodAccum%in% (i: (i+(3)))
    NumberPixels=sapply(DFIndex[SelectRows,-ncol(DFIndex)],mean,na.rm=T)
    DateRef=ChangingPeriodAnotation(DFIndex[SelectRows,ncol(DFIndex)],divPer,"_")
    dfExport=rbind(dfExport,NumberPixels)
    dfDateRef=if(i==1){dfDateRef=DateRef}else{dfDateRef=c(dfDateRef,DateRef)}
  }
  names(dfExport)=names(DFIndex)[-ncol(DFIndex)]
  dfExport$DateRef=dfDateRef
  return(dfExport)
})

names(CountingPeriod)=Proximity



GraphicList=lapply(1:2, function(i){
  graphLow=melt(CountingPeriod[[i]],id="DateRef")
  graphLow$value=as.numeric(as.character(graphLow$value))
  as.factor(graphLow$DateRef)
  graphLow
  graphLow$RefStart=sapply(strsplit(unique(graphLow$DateRef)," - "),function(x){x[[1]]})
  graphLow$RefEnd=sapply(strsplit(unique(graphLow$DateRef)," - "),function(x){x[[2]]})
  YearMonth=data.frame(do.call(rbind,lapply(strsplit(graphLow$RefEnd,"_"),function(x){x})))
  graphLow$YearEnd=YearMonth[,1]
  graphLow$MonthEnd=sapply(YearMonth[,2],function(x){which(month.abb%in%x)})
  graphLow=graphLow[order(graphLow$MonthEnd),]
  graphLow=graphLow[order(graphLow$YearEnd),]
  
  graphLow$DateRef=factor(graphLow$DateRef,levels=unique(graphLow$DateRef))
  graphLow$Limit=switch(i,"Low","Up")
  graphLow$ID=paste0(graphLow$DateRef,"V_",graphLow$variable)
  return(graphLow)
})

graphDF=do.call(rbind,GraphicList)

#write.csv(graphDF,paste0("ClimateComparison/Climatology_VISignal/",VI,"/CountingAnomalies/DatosGrafica.csv"))


graphLow=GraphicList[[2]]
graphLow=graphDF
graphLow$value=as.numeric(as.character(graphLow$value))

m=ggplot(graphLow,aes(DateRef,value,fill=Limit))+geom_bar(stat = "identity")+
  labs(x="",y="Promedio del Conteo")+
  theme(text=element_text(size=16),
        axis.title.y=element_text(size = rel(1),colour = "#999999"),
        axis.title.x=element_text(size = rel(1),colour = "#888888"),
        axis.text.x  = element_text(angle=90, hjust=1,vjust=0.5),
        panel.background=element_rect(fill="white",colour = "black"),
        panel.grid.major = element_line(colour = "gray"))+
  facet_wrap(~ variable,ncol = 1,scales = "free_y")
m

ggsave(paste0("ClimateComparison/Climatology_VISignal/",VI,"/CountingAnomalies/PNG/",VI,"_countingAnomaliesNegativeCat.png"),m,width =24 ,height = 14)
