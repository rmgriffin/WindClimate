## Use this to start every program.  This clears out previous information from memory
rm(list=ls())

## Initalize renv for library lockfile
library(renv)
#renv::init()

## Packages
#Sys.setenv(RENV_PATHS_RTOOLS = "C:/rtools40/") # https://github.com/rstudio/renv/issues/225

PKG <- c("raster","sf","tidyverse","rgdal","ncdf4","RColorBrewer","lattice","googledrive","tmap","chron","furrr","nngeo")

for (p in PKG) {
  if(!require(p,character.only = TRUE)) {  
    install.packages(p)
    require(p,character.only = TRUE)}
}
rm(p,PKG)

## Snapshot of libraries used
renv::snapshot()

## Downloading supporting data
# Download from google drive to directory "Data"
setwd("~/Github/WindClimate")
dir.create(file.path('Data'), recursive = TRUE)
folder_url<-"https://drive.google.com/open?id=1Yc9RiqXXj0qja2DaQaXBwyi267mMZ4Yd"
folder<-drive_get(as_id(folder_url))
files<-drive_ls(folder)
dl<-function(files){
  walk(files, ~ drive_download(as_id(.x), overwrite = TRUE))
}
setwd("./Data")
system.time(map(files$id,dl))
setwd("..")
rm(files, folder, folder_url, dl)

## Download file https://na-cordex.org/data.html
options(timeout=100000) # Keeps download.file from returning a timeout for long-downloading files
download.file("ftp://anonymous:anonymous@ftp2.psl.noaa.gov/pub/Public/dswales/wrfout_d01_2007.nc", destfile = "./Data/wrfout_d01_2007.nc",)

## Actions with a NetCDF
# Reading data
df<-nc_open("./Data/wrfout_d01_2007.nc")
df1<-nc_open("./Data/sfcWind.hist.GFDL-ESM2M.RegCM4.day.NAM-22i.raw.nc")
dfb<-brick("./Data/wrfout_d01_2007.nc")
print(df)

p1<-tm_shape(dfb[[1]]) +
  tm_raster(palette = "Greens", colorNA = NULL, title = "Wind Speed (m/s)") + 
  tm_layout(legend.outside = TRUE)
p1

# Exploring variables
y<-ncvar_get(df,"lat")
ny<-dim(y)
head(y)

x<-ncvar_get(df,"lon")
nx<-dim(x)
head(x)

print(c(nx,ny))

t<-ncvar_get(df,"year")
head(t)

tunits<-ncatt_get(df,"time","units")
nt<-dim(t)
nt
tunits

# Getting data
ws<-ncvar_get(df,"sfcWind")
dim(ws)

# Creating a dataframe from one slice
xy<-as.matrix(expand.grid(x,y))
dim(xy)
slice1<-ws[,,1]
vec1<-as.vector(slice1)
dfs1<-data.frame(cbind(xy,vec1))
names(dfs1)<-c("lon","lat","ws1") # names(tmp_df01) <- c("lon","lat",paste(dname,as.character(m), sep="_"))