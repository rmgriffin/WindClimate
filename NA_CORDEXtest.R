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

# Creating a dataframe from one slice
xy<-as.matrix(expand.grid(x,y))
dim(xy)
slice1<-ws[,,1]
vec1<-as.vector(slice1)
dfs1<-data.frame(cbind(xy,vec1))
names(dfs1)<-c("lon","lat","ws1") # names(tmp_df01) <- c("lon","lat",paste(dname,as.character(m), sep="_"))

## Geospatial analysis
# Creating a dataframe from whole array
dfws<-as.vector(ws)
rm(ws)
dfws<-matrix(dfws,nrow=nx*ny,ncol=nt)
dfws<-data.frame(cbind(xy,dfws)) # x by y by t dataframe
names(dfws)<-c("lon","lat",paste("t",as.character(t), sep=""))
# Lots of NAs due to missing data, remove all rows that only have NAs
dfws<-dfws[complete.cases(dfws),]
# Prep for dealing with time units
tustr<-strsplit(tunits$value, " ")
tdstr<-strsplit(unlist(tustr)[3], "-")
tmonth<-as.integer(unlist(tdstr)[2])
tday<-as.integer(unlist(tdstr)[3])
tyear<-as.integer(unlist(tdstr)[1])

# Power generation calculation
W<-18 # Number of turbines
TT<-3.6 # Turbine rated power  
RD<-120 # Rotor diameter
A<-0.97 # Wind Availability (% as fraction)
EL<-0.98 # Wind Energy Losses (% as fraction)

dates<-as.data.frame(chron(t,origin=c(tmonth, tday, tyear))) # Cannot transform whole dataset wide to long, takes up too much memory, breaking into annual time chunks, running analysis, and reassembling
dates$year<-years(dates$`chron(t, origin = c(tmonth, tday, tyear))`)
dates$year<-as.numeric(dates$year)
xy<-dfws[,1:2]
dfws<-dfws[,3:ncol(dfws)]
switches<-c(1,1+which(diff(dates$year)!=0)) # https://stackoverflow.com/questions/20896242/finding-the-index-of-first-changes-in-the-elements-of-a-vector
switches2<-c(switches[2:length(switches)],length(dates$year))
powerbreak<-function(arg_1, arg_2) {
  dfwss<-dfws[,arg_1:(arg_2-1)]
  dfwss<-cbind(xy,dfwss)
  system.time(dfwsl<-dfwss %>% 
                pivot_longer(cols = 3:ncol(dfwss), names_to = "t", names_prefix = "t", values_to = "ws"))
  dfwsl$t<-as.numeric(dfwsl$t)
  system.time(dfwsl$t<-chron(dfwsl$t,origin=c(tmonth, tday, tyear))) # Date
  system.time(dfwsl$year<-years(dfwsl$t)) # Extracting year factor
  
  dfwsl$kWh<-24*(TT*1000)*(.087*dfwsl$ws-((TT*1000)/RD^2))*A*EL*W # Daily kWh
  dfwsannual<-dfwsl %>% group_by(lon, lat, year) %>% summarise(kWh = sum(kWh)) # Annual aggregation
  dfwsannual<-as.data.frame(dfwsannual)
  dfwsannual$year<-as.character(dfwsannual$year)
  return(dfwsannual)
}


plan(multisession, workers = 2)
options(future.globals.maxSize= 891289600000)
system.time(dfy<-future_map2_dfr(switches,switches2,powerbreak,.progress = TRUE)) # Energy generation at annual time step for each location
rm(dfws)

# Costs
dfy<-st_as_sf(dfy, coords = c("lon", "lat"), crs=st_crs(dfb), remove = FALSE) # Convert annual energy generation df to sf object, don't remove coordinate columns
dfy<-st_transform(dfy,crs = st_crs(3857)) # Reprojecting all layers to EPSG:3857 WGS/Pseudo Mercator - seems consistent with TNC marine mapper
cst<-st_read("./Data/global_polygon.gpkg")
cst<-st_transform(cst,crs = st_crs(3857))
eeza<-st_read("./Data/eez_atlantic.gpkg") # US EEZ vector https://www.marineregions.org/gazetteer.php?p=details&id=8456
eeza<-st_transform(eeza,st_crs(3857))
eeza<-eeza[1:2]

system.time(dfy<-st_join(dfy,eeza)) # Keeping only those observations that are within the US Atlantic EEZ
dfy<-dfy %>% drop_na()

eeza10kmbuff<-st_buffer(eeza,dist = 10000)
csteeza10kmbuff<-st_intersection(cst,eeza10kmbuff) # Creating a limited nearshore polygon to search over for nearest landing point for wind points

system.time(dfy$nearest<-st_nearest_feature(dfy,csteeza10kmbuff))

#csteeza10kmbuff<-csteeza10kmbuff[1,] # Only mainland connection points

system.time(dfy$dist<-as.vector(st_distance(dfy, csteeza10kmbuff[dfy$nearest,], by_element=TRUE)))



ggplot() +
  geom_sf(data = csteeza10kmbuff) +
  geom_sf(data = xy)
  
system.time(xy2<-nngeo::st_nn(xy,csteeza10kmbuff,returnDist = TRUE))

xy3<-cbind(xy2$nn,xy2$dist)
# LCOE
# Revenue
# NPV

## Plotting

# Mask raster brick by EEZ
system.time(dfbeez<-mask(dfb,eeza))

# Visualizing data (one slice)
image(x,y,slice1, col=rev(brewer.pal(10,"RdBu"))) # Plot 1
grid<-expand.grid(lon=x, lat=y)
cutpts<-seq(0,30,3)
levelplot(slice1 ~ lon * lat, data=grid, at=cutpts, cuts=9, pretty=T, 
          col.regions=(rev(brewer.pal(11,"RdBu")))) # Plot 2

# Plot raster and eez
p1<-tm_shape(eeza) +
  tm_borders(col = "black", lty = "dashed")+
  tm_shape(dfbeez[[3]]) +
  tm_raster(palette = "Greens", colorNA = NULL, title = "Wind Speed (m/s)") + 
  tm_layout(legend.outside = TRUE)
p1

writeRaster(dfbeez[[3]],"windspeed.tif")
# # EEZA in meters projected and buffered
# eezaUTM18N<-st_transform(eeza,crs = 32618)
# eezaUTM18N<-st_buffer(eezaUTM18N,dist = 35000)
# st_write(eezaUTM18N,"./eezaUTM18N.gpkg")

