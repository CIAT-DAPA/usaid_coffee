
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


CoffeePoints=read.csv("FusionTable/CSV/Coffee_Points_2015_2016Match.csv")
xy=cbind( CoffeePoints$Longitud, CoffeePoints$Latitud)
CoffeePoints_spRef = SpatialPointsDataFrame(coords = xy,data=CoffeePoints,
                                            proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))




RasterDirAnomalies="//dapadfs/workspace_cluster_6/TRANSVERSAL_PROJECTS/MADR/COMPONENTE_2/USAID/ClimateComparison/Climatology/Anomalies"
RasterDirAnomalies=paste(RasterDirAnomalies,CleaningTifNames(RasterDirAnomalies),sep="/")


Proximity=c("Low","Up")

PosPattern=gregexpr("ies/.*._20",RasterDirAnomalies,perl = T)
ClimateVariables=unique(str_sub(unlist(regmatches(RasterDirAnomalies,PosPattern)),5,-4))

ClimateVar=ClimateVariables[1]

x=Proximity[1]
ListCountEvents=lapply(Proximity,function(x){
  SymbolComp=switch(x,"Low"=">=","Up"="<")
  cat("Starting Multi values to points process: \n")
  NamesRasterLimits=RasterDirAnomalies[grepl(RasterDirAnomalies,pattern = x)]
  ClimateVar=ClimateVariables[1]
  CountPerVar=lapply(ClimateVariables,function(ClimateVar){
    NamesRasterVar=NamesRasterLimits[grepl(NamesRasterLimits,pattern = paste0(ClimateVar,"_20"))]
    
    Sys.time()->start
    library(snowfall)
    ncores=2
    sfInit(parallel=T,cpus=ncores)
    sfLibrary(raster)
    RasterFiles=sfLapply(NamesRasterVar,function(x){
      m=raster(x)
      crs(m) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
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

save(ListCountEvents,file="ClimateComparison/Climatology/Counting1RIQ.RData")
load(file="ClimateComparison/Climatology/Counting1_5RIQ.RData")

divPer=12


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

graphLow=reshape::melt(CountingPeriod[[1]],id="DateRef")
graphLow$value=as.numeric(as.character(graphLow$value))

m=ggplot(graphLow,aes(DateRef,value))+geom_bar(stat = "identity")+
  labs(x="",y="Promedio del Conteo")+
  theme(text=element_text(size=10),
        axis.title.y=element_text(size = rel(1),colour = "#999999"),
        axis.title.x=element_text(size = rel(1),colour = "#888888"),
        axis.text.x  = element_text(angle=90, hjust=1,vjust=0.5),
        panel.background=element_rect(fill="white",colour = "black"),
        panel.grid.major = element_line(colour = "gray"))+
  facet_wrap(~ variable,ncol = 1,scales = "free_y")
m

