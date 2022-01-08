## Use this to start every program.  This clears out previous information from memory
rm(list=ls())

## Initalize renv for library lockfile
library(renv)
#renv::init()

## Packages
#Sys.setenv(RENV_PATHS_RTOOLS = "C:/rtools40/") # https://github.com/rstudio/renv/issues/225

PKG <- c("raster","sf","tidyverse","rgdal","ncdf4","RColorBrewer","lattice","googledrive","tmap")

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

#testing following https://pjbartlein.github.io/REarthSysSci/netCDF.html
## Download file https://na-cordex.org/data.html
# Useful reference for NA CORDEX https://doi.org/10.1016/j.cliser.2021.100233
options(timeout=100000) # Keeps download.file from returning a timeout for long-downloading files
download.file("https://tds.ucar.edu/thredds/fileServer/datazone/cordex/data/raw/NAM-22i/day/RegCM4/GFDL-ESM2M/hist/sfcWind/sfcWind.hist.GFDL-ESM2M.RegCM4.day.NAM-22i.raw.nc", destfile = "./Data/sfcWind.hist.GFDL-ESM2M.RegCM4.day.NAM-22i.raw.nc",)

## Actions with a NetCDF
# Reading data
df<-nc_open("./Data/sfcWind.hist.GFDL-ESM2M.RegCM4.day.NAM-22i.raw.nc")
dfb<-brick("./Data/sfcWind.hist.GFDL-ESM2M.RegCM4.day.NAM-22i.raw.nc")
print(df)

# Exploring variables
y<-ncvar_get(df,"lat")
ny<-dim(y)
head(y)

x<-ncvar_get(df,"lon")
nx<-dim(x)
head(x)

print(c(nx,ny))

t<-ncvar_get(df,"time")
head(t)

tunits<-ncatt_get(df,"time","units")
nt<-dim(t)
nt
tunits

# Getting data
ws<-ncvar_get(df,"sfcWind")
dim(ws)

# Visualizing data (one slice)
slice1<-ws[,,1]
image(x,y,slice1, col=rev(brewer.pal(10,"RdBu"))) # Plot 1
grid<-expand.grid(lon=x, lat=y)
cutpts<-seq(0,30,3)
levelplot(slice1 ~ lon * lat, data=grid, at=cutpts, cuts=9, pretty=T, 
          col.regions=(rev(brewer.pal(11,"RdBu")))) # Plot 2

# Creating a dataframe from one slice
xy<-as.matrix(expand.grid(x,y))
dim(xy)
vec1<-as.vector(slice1)
dfs1<-data.frame(cbind(xy,vec1))
names(dfs1)<-c("lon","lat","ws1") # names(tmp_df01) <- c("lon","lat",paste(dname,as.character(m), sep="_"))

# Creating a dataframe from whole array
dfws<-as.vector(ws)
dfws<-matrix(dfws,nrow=nx*ny,ncol=nt)
dfws<-data.frame(cbind(xy,dfws)) # x by y by t dataframe
names(dfws)<-c("lon","lat",paste("t",as.character(t), sep=""))

## Geospatial analysis
# US EEZ - Atlantic Coast
eeza<-st_read("./Data/eez_atlantic.gpkg") # US EEZ vector https://www.marineregions.org/gazetteer.php?p=details&id=8456
eeza<-st_transform(eeza,st_crs(dfb))

# Mask raster brick by EEZ
system.time(dfbeez<-mask(dfb,eeza))

# Plot raster and eez
p1<-tm_shape(eeza) +
  tm_borders(col = "black", lty = "dashed")+
  tm_shape(dfbeez[[3]]) +
  tm_raster(palette = "Greens", colorNA = NULL, title = "Wind Speed (m/s)") + 
  tm_layout(legend.outside = TRUE)
p1

# # EEZA in meters projected and buffered
# eezaUTM18N<-st_transform(eeza,crs = 32618)
# eezaUTM18N<-st_buffer(eezaUTM18N,dist = 35000)
# st_write(eezaUTM18N,"./eezaUTM18N.gpkg")

