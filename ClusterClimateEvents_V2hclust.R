###### Functions

rm(list=ls())

######

ReadingPackages=function(package){
  lapply(package,function(x){
    if(!require(x, character.only = T)){
      install.packages(x);require(x, character.only = T)
    }else{
      require(x, character.only = T)
    }
  })
}


libs= c("reshape2","dtwclust","dtw","gtools","agricolae","stringr","snowfall","TSclust","TSdist","cluster")

ReadingPackages(libs)

##### Setting the workspace

server="/dapadfs"
server="mnt"
dirfol=paste0("/",server,"/workspace_cluster_6/TRANSVERSAL_PROJECTS/MADR/COMPONENTE_2/USAID/")
source(paste0("/",server,"/workspace_cluster_6/TRANSVERSAL_PROJECTS/MADR/COMPONENTE_2/USAID/ClimateComparison/Cluster/Script/funciones_cluster_temporal.R"))

setwd(dirfol)


##### Reading the file

varEstudio=c("prec")
Lacks=0:4
Lack=1
VI="NDVI"
NameID="IDMatch"
namePendiente="Pendiente"
period="4_Months"
sfInit(parallel = T, cpus = 5)
sfExport("NameID")
sfExport("namePendiente")
sfExport("dirfol")
sfExport("varEstudio")
sfExport("distDtwMV")
sfExport("period")
sfExport("VI")
sfExport("multiplotClust")
sfLibrary(reshape2)
sfLibrary(ggplot2)
sfLibrary(dtw)
sfLibrary(dtwclust)
sfLibrary(TSclust)
sfLibrary(TSdist)
sfLibrary(cluster)

sfLapply(Lacks,function(Lack){
  setwd(dirfol)
  EventsFile= paste0("ClimateComparison/Climatology_VISignal/",VI,"/DeltaClimeAssociation/DeltaValues",period,"_",VI,"_WithNewVar_Coffee_V2_L",Lack,".csv")
  DataBase=read.csv(EventsFile,row.names = 1)
  DataBase=DataBase[!is.na(DataBase$prec_1),]
  NewFolder1=paste0("Lag_",Lack)
  dir.create(paste0("ClimateComparison/Cluster/",VI),showWarnings = F)
  dir.create(paste0("ClimateComparison/Cluster/",VI,"/",period,"/"),showWarnings = F)
  dir.create(paste0("ClimateComparison/Cluster/",VI,"/",period,"/",NewFolder1),showWarnings = F)
  setwd(paste0("ClimateComparison/Cluster/",VI,"/",period,"/",NewFolder1))
  
  counti=1
  ListPosNeg=list()
  lapply(1:2, function(counti){
    NewFolder2=switch(counti,"NegativeSlope","PositiveSlope")
    setwd(paste0(dirfol,"ClimateComparison/Cluster/",VI,"/",period,"/",NewFolder1,"/"))
    dir.create(NewFolder2,showWarnings = F)
    setwd(NewFolder2)
    RowSelection=switch(counti, DataBase$CatSlope=="Negative",DataBase$CatSlope=="Positive")
    
    DataBase_Aux=DataBase[RowSelection,]
    DataBase_Aux
    head(DataBase_Aux)
    
    
    NamesColumns=names(DataBase_Aux)[grepl(paste0(varEstudio,collapse = "|"),names(DataBase_Aux))]
    DTW_DF=DataBase_Aux[,c(NameID,NamesColumns)]
    
    ListVarLong=lapply(varEstudio,function(varName){
      LongDF=melt(DTW_DF,id.vars = NameID,measure.vars = NamesColumns[grepl(NamesColumns,pattern = varName)])[,c(1,3)]
    })
    DTW_DF_Join=do.call(cbind,ListVarLong)[,c(1,seq(2,length(varEstudio)*2,by = 2))]
    names(DTW_DF_Join)[c(-1)]=varEstudio
    
    tail(DTW_DF_Join[order(DTW_DF_Join[,NameID]),])
    tail(DTW_DF[order(DTW_DF[,NameID]),])
    
    evenF= split(DTW_DF_Join,as.character(DTW_DF_Join[,NameID]))
    evenF[[1]]
    length(evenF)
    x=evenF[[1]]
    x[,1]
    
    data3=lapply(evenF,function(x)(x[,2]))
    data <- zscore(data3)
    
    fc = tsclust(data, type = "h", k = 15L,
                 distance = "dtw_basic",
                 control = hierarchical_control(method = diana),
                 args = tsclust_args(dist = list(window.size = 18L)))
    
    
    #CREATE DTW DISTANCE
    NewFolder3=paste0(varEstudio,collapse = "_")
    dir.create(NewFolder3,showWarnings = F)
    cat(getwd())
    setwd(NewFolder3)
    
    IdEvent <- names(fc@cluster)
    EventsClass=data.frame(IdEvent,fc@cluster)
    
    listEventsDate=lapply(evenF,function(x){return(as.data.frame(cbind(Dates=rep(1:dim(x)[1]),x)))})
    membAll=fc@cluster
    nlevents=listEventsDate
    head(nlevents)
    minVal  <- apply(do.call(rbind,lapply(nlevents,function(x){apply(x[,-c(2)],2,min,na.rm = T)})),2,min,na.rm = T)
    maxVal  <- apply(do.call(rbind,lapply(nlevents,function(x){apply(x[,-c(2)],2,max,na.rm = T)})),2,max,na.rm = T)
    
    limts <- as.data.frame(rbind(minVal,maxVal))
    
    NomVarxAxix="Dates"
    nomVarClust=varEstudio
    listMeans=list()
    
    for(i in 1:length(unique(membAll))){
      png(paste0(NewFolder1,"_",NewFolder2,"_DTWcluster",i,".png"),width =16, height = 10,res=200,units = 'in')
      layout(cbind(1:length(nomVarClust)))
      subBas <- nlevents[which(membAll==i)]
      
      listMeansVar=list()
      for(j in nomVarClust)
      {
        limtsV=limts[[j]]
        listMeansVar[[j]]=multiplotClust(listG = subBas,NomVarx=NomVarxAxix,NomVary = j,limts =limtsV,textNum = T,mainTitle = paste0("Cluster_",i),SaveMean = T) #### Activar para guardar 
      }
      dev.off()
      listMeans[[i]]=listMeansVar
    }
    ListPosNeg[[NewFolder2]]=listMeans
    save(listMeans, file="ClimateSeries_Means.RData")
    save(fc, file="dtwClust_data.RData")
    
    write.csv(EventsClass,"ClassifiedEvents.csv")
    
    
  })
})

sfStop()
#load("distMatrixCluster.RData")
Lack=0
Lacks=1:4
Clustering=T
listLack=list()
for(Lack in Lacks){
  setwd(dirfol)
  NewFolder1=paste0("Lag_",Lack)
  setwd(paste0("ClimateComparison/Cluster/",VI,"/",period,"/",NewFolder1))
  ListPosNeg=list()
  for(counti in 1:2){
    
    NewFolder2=switch(counti,"NegativeSlope","PositiveSlope")
    setwd(paste0(dirfol,"ClimateComparison/Cluster/",VI,"/",period,"/",NewFolder1,"/",NewFolder2,"/",paste0(varEstudio,collapse = "_")))
    cat(paste0(getwd(),"\n"))
    load(file = "distMatrixCluster.RData")
    load(file = "Normal_Events.RData")
    load(file = "List_Events.RData")
    if(Clustering){
      hClustEvents <- hirarCluster(distMatrix = distAllMatrix)
      IdEvent <- names(evenN)
      EventsClass=data.frame(IdEvent,hClustEvents)
      write.csv(EventsClass,"eventsClasificated.csv")
    }else{
      hClustEvents=read.csv("eventsClasificated.csv")[,3]
    }
    x=evenF[[1]]
    listEventsDate=lapply(evenF,function(x){return(as.data.frame(cbind(Dates=rep(1:dim(x)[1]),x)))})
    membAll=hClustEvents
    nlevents=listEventsDate
    head(nlevents)
    minVal  <- apply(do.call(rbind,lapply(nlevents,function(x){apply(x[,-c(2)],2,min,na.rm = T)})),2,min,na.rm = T)
    maxVal  <- apply(do.call(rbind,lapply(nlevents,function(x){apply(x[,-c(2)],2,max,na.rm = T)})),2,max,na.rm = T)
    minVal[2]=-45
    maxVal[2]=45
    limts <- as.data.frame(rbind(minVal,maxVal))
    
    NomVarxAxix="Dates"
    nomVarClust=varEstudio
    listMeans=list()
    
    for(i in 1:length(unique(membAll))){
      png(paste0(NewFolder1,"_",NewFolder2,"_DTWcluster",i,".png"),width =16, height = 10,res=200,units = 'in')
      layout(cbind(1:length(nomVarClust)))
      subBas <- nlevents[which(membAll==i)]
      
      listMeansVar=list()
      for(j in nomVarClust)
      {
        limtsV=limts[[j]]
        listMeansVar[[j]]=multiplotClust(listG = subBas,NomVarx=NomVarxAxix,NomVary = j,limts =limtsV,textNum = T,mainTitle = paste0("Cluster_",i),SaveMean = T) #### Activar para guardar 
      }
      dev.off()
      listMeans[[i]]=listMeansVar
    }
    ListPosNeg[[NewFolder2]]=listMeans
  }
  listLack[[Lack]]=ListPosNeg
}

ggplot(EventsClass,aes(Date,hClustEvents))+geom_point()+
  theme(text=element_text(size=15),
        axis.title.y=element_text(size = rel(1.3),colour = "#999999"),
        axis.title.x=element_text(size = rel(1.2),colour = "#888888"),
        axis.text.x  = element_text(angle=90, hjust=1,vjust = 0.5, size =rel(0.8)),
        panel.background=element_rect(fill="white",colour = "black"),
        panel.grid.major = element_line(colour = "gray"))+
  scale_y_continuous(limits=c(1,(max(EventsClass$hClustEvents)+1)),breaks=seq(1,max(EventsClass$hClustEvents),by = 1))+
  stat_summary(fun.data = give.n, geom = "text", fun.y = median)
