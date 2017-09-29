###### Functions

rm(list=ls())
#source("D:/ToBackup/Documents/Scripts/R/Clustering/funciones_cluster_temporal.R")
libs= c("reshape2","dtwclust","dtw","gtools","agricolae","stringr","snowfall")

lapply(libs, require, character.only=T)

##############################################################

library("dtwclust")

require("TSclust")
require("TSdist")

require("cluster")
x=evenF[[1]]
data2=lapply(evenF,function(x)ts(x[,2]))
data3=lapply(evenF,function(x)(x[,2]))

data <- zscore(data2)

pc_k <- tsclust(data2, type = "p", k = 4L:20L,
                distance = "dtw_basic", centroid = "pam",
                args = tsclust_args(dist = list(window.size = 20L)),
                seed = 93)

names(pc_k) <- paste0("k_", 4L:20L)
sapply(pc_k, cvi, type = "internal")


fc = tsclust(data3, type = "h", k = 28L,
        distance = "dtw_basic",
        control = hierarchical_control(method = diana),
        args = tsclust_args(dist = list(window.size = 18L)))



plot(fc, series = data3, type = "series")
plot(fc, type = "centroids")
table(fc@cluster)

pc_ks <- tsclust(data3, k = 10L,
                 distance = "sbd", centroid = "shape",
                 seed = 8, trace = TRUE)

plot(pc_ks, series = data3, type = "series")
table(pc_ks@cluster)

pc_dtw <- tsclust(data3, k = 10L,
                  distance = "dtw_basic", centroid = "dba",
                  trace = TRUE, seed = 8,
                  norm = "L2", window.size = 20L,
                  args = tsclust_args(cent = list(trace = TRUE)))

plot(pc_dtw, series = data3, type = "series")
table(pc_dtw@cluster)


pc_tp <-tsclust(data3, k = 10L, type = "t",
                seed = 8, trace = TRUE,
                control = tadpole_control(dc = 1.5,
                                          window.size = 3L))
plot(pc_tp, series = data3, type = "series")
table(pc_tp@cluster)

pc_dtwlb <- tsclust(data3, k = 10L,
                    distance = "dtw_lb", centroid = "dba",
                    trace = TRUE, seed = 8,
                    norm = "L2", window.size = 3L,
                    args = tsclust_args(cent = list(trace = TRUE)))

acf_fun <- function(dat, ...) {lapply(dat, function(x) {
  as.numeric(acf(x, lag.max = 50, plot = FALSE)$acf)
})
}


fc <- tsclust(data3, type = "f", k = 10L,
              preproc = acf_fun, distance = "L2",
              seed = 42)


pc_k_Eval <- tsclust(data, type = "p", k = 4L:15L,
                distance = "dtw_basic", centroid = "pam",
                args = tsclust_args(dist = list(window.size = 20L)),
                seed = 93)


names(pc_k_Eval) <- paste0("k_",  4L:15L)
sapply(pc_k_Eval, cvi, type = "internal")

plot(fc, series = data3, type = "series")
plot(fc, type = "centroids")

plot(pc_dtwlb, series = data3, type = "series")
table(pc_dtwlb@cluster)

sapply(list(DTW = pc_dtw, DTW_LB = pc_dtwlb, kShape = pc_ks, TADPole = pc_tp, DTW_Basic=fc),
       cvi, b = names(data3), type = "VI")

names(CharTraj)[60:100]

CharTrajLabels[60L:100L]

IdEvent <- names(fc@cluster)
EventsClass=data.frame(IdEvent,fc@cluster)


listEventsDate=lapply(evenF,function(x){return(as.data.frame(cbind(Dates=rep(1:dim(x)[1]),x)))})
membAll=fc@cluster
nlevents=listEventsDate
head(nlevents)
minVal  <- apply(do.call(rbind,lapply(nlevents,function(x){apply(x[,-c(2)],2,min,na.rm = T)})),2,min,na.rm = T)
maxVal  <- apply(do.call(rbind,lapply(nlevents,function(x){apply(x[,-c(2)],2,max,na.rm = T)})),2,max,na.rm = T)
minVal[2]=-90
maxVal[2]=130
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





pc_k <- tsclust(listObje[1:300], type = "p", k = 2L:6L,
                distance = "dtw_basic", centroid = "pam",
                args = tsclust_args(dist = list(window.size = 20L)),
                seed = 93)


data <- reinterpolate(CharTraj, new.length = max(lengths(CharTraj)))

pc_dtw <- tsclust(data[1:50], k = 4L,
                  distance = "dtw_basic", centroid = "dba",
                  trace = TRUE, seed = 8,
                  norm = "L2", window.size = 20L,
                  args = tsclust_args(cent = list(trace = TRUE)))


distDtwMV=function(listObje)
{
  require(dtw)
  len      <-  length(listObje)
  listDist <- matrix(0,nrow=len,ncol=len)
  
  disDtw <- array(0,len)
  # pb <- winProgressBar(title="Progress bar", label="0% done", min=0, max=100, initial=0)
  pb <- txtProgressBar(min = 0, max = len, style = 3)
  for(i in 1:len){ 
    
    for(j in i:len){
      listDist[j,i]<- dtw(listObje[[i]],listObje[[j]])$distance
    }
    setTxtProgressBar(pb, i)
    #  info <- sprintf("%d%% done", round((i/(len)*100)))
    # setWinProgressBar(pb, i/(len)*100, label=info)
  }
  close(pb)
  #close(pb)
  rownames <- names(listObje)
  colnames <- names(listObje)
  return(as.dist(listDist))
}


##### Setting the workspace

server="/dapadfs"
server="mnt"
dirfol=paste0("/",server,"/workspace_cluster_6/TRANSVERSAL_PROJECTS/MADR/COMPONENTE_2/USAID/")
source(paste0("/",server,"/workspace_cluster_6/TRANSVERSAL_PROJECTS/MADR/COMPONENTE_2/USAID/ClimateComparison/Cluster/Script/funciones_cluster_temporal.R"))

setwd(dirfol)


##### Reading the file

varEstudio=c("prec")
Lacks=0:4
Lack=0
VI="NDVI"
NameID="IDMatch"
namePendiente="Pendiente"
period="4_Months"
sfInit(parallel = T, cpus = 7)
sfExport("NameID")
sfExport("namePendiente")
sfExport("dirfol")
sfExport("varEstudio")
sfExport("distDtwMV")
sfExport("period")
sfExport("VI")
sfLibrary(reshape2)
sfLibrary(ggplot2)
sfLibrary(dtw)


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
    
    evenF[["L_660010008612t2016_Nov"]]
    minVal  <- apply(do.call(rbind,lapply(evenF,function(x){apply(x,2,min,na.rm = T)})),2,min,na.rm = T)
    maxVal  <- apply(do.call(rbind,lapply(evenF,function(x){apply(x,2,max,na.rm = T)})),2,max,na.rm = T)
    
    x=evenF[[1]]
    evenN=lapply(evenF,function(x){
      x=data.frame(sapply(2:(length(varEstudio)+1),function(y){
        (x[,y]-(as.numeric(minVal[y])))/(as.numeric(maxVal[y])-(as.numeric(minVal[y])))
      }))
      names(x)=varEstudio
      return(x)
    })
    length(evenN)
    tsnleventsN <- lapply(evenN,ts)
    
    #CREATE DTW DISTANCE
    NewFolder3=paste0(varEstudio,collapse = "_")
    dir.create(NewFolder3,showWarnings = F)
    cat(getwd())
    distAllMatrix <- distDtwMV(listObje = tsnleventsN)
    setwd(NewFolder3)
    write.csv(DataBase_Aux[,c(NameID,NamesColumns,namePendiente)],"SlopeValues3.csv")
    save(evenN,file = "Normal_Events.RData")
    save(evenF,file = "List_Events.RData")
    save(distAllMatrix,file = "distMatrixCluster.RData")
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
