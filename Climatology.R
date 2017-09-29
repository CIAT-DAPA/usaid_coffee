######### Calculating the Climatologic Normal


rm(list=ls())

###
#            -----   Loading packages
###

libs=c("maptools","rgdal","raster","sp","stringr","ff","dplyr","ggplot2","reshape2")
lapply(libs,require,character.only=T)

DataShape=readShapeSpatial("//dapadfs/workspace_cluster_6/TRANSVERSAL_PROJECTS/MADR/COMPONENTE_2/USAID/SHP/Risaralda.shp")
proj4string(DataShape) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")


###
#            -----   Reading the raster Files
###

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


RasterDir2010="//dapadfs/Workspace_cluster_9/USAID_Project/Product_6_resilient_coffee/02-monthly-interpolations/outputs_yearly_v2/average"
FileNames2010=paste(RasterDir2010,CleaningTifNames(RasterDir2010),sep="/")

RasterDir2016="//dapadfs/Workspace_cluster_9/USAID_Project/Product_6_resilient_coffee/02-monthly-interpolations/outputs_yearly_v2_2011_2016/average"
FileNames2016=paste(RasterDir2016,CleaningTifNames(RasterDir2016),sep="/")

CVariables=unique(unlist(lapply(str_split(CleaningTifNames(RasterDir2016), pattern="_"),function(x){x[1]})))
VarI=CVariables[1]

lapply(CVariables,function(VarI){
  NameVarfiles=c(FileNames2010,FileNames2016)[grepl(x = c(FileNames2010,FileNames2016),pattern = VarI)]
  library(snowfall)
  ncores=20
  sfInit(parallel=T,cpus=ncores)
  sfLibrary(raster)
  RasterFilesDTR=sfLapply(NameVarfiles,function(x){
    m=raster(x)
    crs(m) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
    return(m)})
  crop(RasterFilesDTR[[1]], DataShape)
  sfStop()
  RasterExtentCor=lapply(RasterFilesDTR,function(x){
    if(!x@extent[1,]==RasterFilesDTR[[1]]@extent[1,]){
      return(crop(x, DataShape))
    }else{return(x)}
  })
  x=RasterExtentCor[[1]]
  NamesRaster=unlist(lapply(RasterExtentCor,function(x){
    names(x)
  }))
  names(RasterExtentCor)=NamesRaster
  lapply(unique(str_sub(NamesRaster,-2)), function(Month_){
    RasterMonth=NamesRaster[grepl(str_sub(NamesRaster,-2),pattern =Month_)]
    STRasterFiles=stack(RasterExtentCor[names(RasterExtentCor)%in%RasterMonth])
    MeanRaster=mean(STRasterFiles)
    if(str_detect(Month_,pattern = "_")){
      rx = gregexpr("_.",Month_, perl = TRUE)
      valMTif=unique(str_sub(unlist(regmatches(Month_, rx)),2))
    }else{
      valMTif= Month_
    }
    nameFile=paste0("//dapadfs/workspace_cluster_6/TRANSVERSAL_PROJECTS/MADR/COMPONENTE_2/USAID/ClimateComparison/Climatology/",VarI,"_",valMTif,".tif")
    writeRaster(MeanRaster, filename=nameFile, format="GTiff", overwrite=T)
  })
})



