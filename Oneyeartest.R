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
# dfb<-brick("./Data/wrfout_d01_2007.nc", crs = 4326)
# dfb<-projectRaster(dfb, crs = 3857)
cst<-st_read("./Data/global_polygon.gpkg")
cst<-st_transform(cst,crs = st_crs(3857))
eeza<-st_read("./Data/eez_atlantic.gpkg") # US EEZ vector https://www.marineregions.org/gazetteer.php?p=details&id=8456
eeza<-st_transform(eeza,st_crs(3857))
eeza<-eeza[1:2]
print(df)

p1<-tm_shape(dfb[[1003]]) +
  tm_raster(palette = "Greens", colorNA = NULL, title = "Wind Speed (m/s)") + 
  # tm_shape(cst) +
  # tm_borders(col = "black", lty = "dashed") +
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

year<-ncvar_get(df,"year")
head(year)

day<-ncvar_get(df,"day")
head(day)

hour<-ncvar_get(df,"hour") # One observation every 3 hours
head(hour)

z<-ncvar_get(df,"level")
head(z)

tunits<-ncatt_get(df,"time","units")
nt<-dim(t)
nt
tunits

# Getting data
wsu<-as.vector(ncvar_get(df,"u"))
dim(wsu)
head(wsu)
wsv<-as.vector(ncvar_get(df,"v"))
dim(wsv)
head(wsv)


# Testing domain
x2<-as.matrix(expand.grid(x))
y2<-as.matrix(expand.grid(y))
xy<-as.data.frame(cbind(x2,y2))
names(xy)<-c("lon","lat")
# xy<-st_as_sf(xy, coords = c("lon", "lat"),crs=4269, remove = FALSE)
# xy<-st_transform(xy,st_crs(3857))

ggplot() + # Quick plot of data points
  geom_sf(data = cst) +
  geom_sf(data = xy)

# Creating a dataframe from whole array

expand.grid.df<-function(...) Reduce(function(...) merge(..., by=NULL), list(...)) # https://stackoverflow.com/questions/11693599/alternative-to-expand-grid-for-data-frames
system.time(xyz<-expand.grid.df(xy,z,day))
df<-cbind(xyz,wsu,wsv)