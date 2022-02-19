## Use this to start every program.  This clears out previous information from memory
rm(list=ls())

## Initalize renv for library lockfile
library(renv)
#renv::init()
#renv::activate() # Run if starting from cloned repository to load package dependencies

## Packages
#Sys.setenv(RENV_PATHS_RTOOLS = "C:/rtools40/") # https://github.com/rstudio/renv/issues/225

PKG <- c("collapse","raster","sf","tidyverse","rgdal","ncdf4","RColorBrewer","lattice","googledrive","tmap","chron","furrr","nngeo")

for (p in PKG) {
  if(!require(p,character.only = TRUE)) {  
    install.packages(p)
    require(p,character.only = TRUE)}
}
rm(p,PKG)

options(collapse_mask = "manip") # https://twitter.com/grant_mcdermott/status/1493400952878952448?s=20&t=H9w3Azc04_LlJnqqk-56gQ

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
#dfb<-brick("./Data/sfcWind.hist.GFDL-ESM2M.RegCM4.day.NAM-22i.raw.nc")
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

## Geospatial analysis
# Creating a dataframe from whole array
xy<-as.matrix(expand.grid(x,y))
dim(xy)
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
  dfwsannual<-dfwsl %>% group_by(lon, lat, year) %>% summarise(kWh = sum(kWh)) # Annual aggregation of kWh
  ws<-dfwsl %>% group_by(lon, lat, year) %>% summarise(ws = mean(ws)) # Annual mean wind speed
  dfwsannual<-as.data.frame(dfwsannual)
  ws<-as.data.frame(ws)
  dfwsannual$year<-as.character(dfwsannual$year)
  dfwsannual$ws<-ws$ws
  return(dfwsannual)
}

plan(multisession, workers = 2)
options(future.globals.maxSize= 891289600000)
system.time(dfy<-future_map2_dfr(switches,switches2,powerbreak)) # Energy generation at annual time step for each location
rm(dfws)

# Spatializing data
dfy<-st_as_sf(dfy, coords = c("lon", "lat"), crs=4326, remove = FALSE) # Convert annual energy generation df to sf object, don't remove coordinate columns
dfy<-st_transform(dfy,crs = st_crs(3857)) # Reprojecting all layers to EPSG:3857 WGS/Pseudo Mercator - seems consistent with TNC marine mapper
cst<-st_read("./Data/global_polygon.gpkg")
cst<-st_transform(cst,crs = st_crs(3857))
eeza<-st_read("./Data/eez_atlantic.gpkg") # US EEZ vector https://www.marineregions.org/gazetteer.php?p=details&id=8456
eeza<-st_transform(eeza,st_crs(3857))
eeza<-eeza[1:2]

system.time(dfy<-st_join(dfy,eeza)) # Keeping only those observations that are within the US Atlantic EEZ
dfy<-dfy %>% drop_na()
dfy<-dfy %>% 
  select (-c(geoname, mrgid))

eeza10kmbuff<-st_buffer(eeza,dist = 10000)
csteeza10kmbuff<-st_intersection(cst,eeza10kmbuff) # Creating a limited nearshore polygon to search over for nearest landing point for wind points

# Figures for power generation and wind speed
dfy %>% st_drop_geometry(.) %>% distinct(lon,lat) %>% nrow(.) # distinct combinations of lat/lon 
dfy$index<-interaction(dfy$lon,dfy$lat, drop = TRUE)

ggplot(dfy %>% filter(index %in% sample(levels(dfy$index), size = 12)), aes(x = year, y=kWh, group=index)) + # filter draws a random sample of factors to plot 
  geom_line(aes(color=index), size=1) +
  #geom_point(aes(color=index), size=1) +
  #scale_color_gradient() + 
  theme_classic() +
  labs(y = "kWh", x = "Year")

ggplot(dfy %>% filter(index %in% sample(levels(dfy$index), size = 12)), aes(x = year, y=ws, group=index)) + 
  geom_line(aes(color=index), size=1) +
  #geom_point(aes(color=index), size=1) +
  #scale_color_gradient() + 
  theme_classic() +
  labs(y = "m/s", x = "Year")

# Costs
xy<-as.matrix(expand.grid(x,y)) # Unique wind data points 
xy<-as.data.frame(xy)
names(xy)<-c("lon","lat")
xy<-st_as_sf(xy, coords = c("lon", "lat"), crs=4326, remove = FALSE) 
xy<-st_transform(xy,crs = st_crs(3857))
system.time(xy<-st_join(xy,eeza)) # Keeping only those observations that are within the US Atlantic EEZ
xy<-xy %>% drop_na()
xy<-xy %>% 
  select (-c(geoname, mrgid))

system.time(xy$nearest<-st_nearest_feature(xy,csteeza10kmbuff)) # Finding nearest feature to observation points
system.time(xy$dist<-as.vector(st_distance(xy, csteeza10kmbuff[xy$nearest,], by_element=TRUE))) # Distance to nearest feature
#csteeza10kmbuff<-csteeza10kmbuff[1,] # Only mainland connection points

# dfy<-dfy %>% 
#   left_join(.,st_drop_geometry(xy), by = c("lon","lat")) # Merge using unique wind data points is much quicker than finding distance of nearest for all points for all years

#system.time(dfy$nearest<-st_nearest_feature(dfy,csteeza10kmbuff))
#system.time(dfy$dist<-as.vector(st_distance(dfy, csteeza10kmbuff[dfy$nearest,], by_element=TRUE)))

xy$trnsc<-ifelse(xy$dist>60000,810000*W*TT+1360000*xy$dist/1000,1090000*W*TT+890000*xy$dist/1000) # Transmission capital cost, uses different functions less/more than 60km

TS<-6410000 # Wind 3.6MW Turbine Unit Cost 
TL<-10600000 # Wind 5.0MW Turbine Unit Cost
IC<-305000 # Wind Infield Cable Cost per km
MF<-1860000 # Wind Monopile Foundation Unit Cost
JF<-2060000 # Wind Jacketed Foundation Unit Cost
TI<-.20 # Wind Installation Cost as a Percentage of CAPEX
TM<-.08 # Wind Miscellaneous Costs as a Percentage of CAPEX
TO<-.035 # Wind Operations and Management Costs as a Percentage of CAPEX per Year
TD<-.133 # Wind Weighted Average Cost of Capital (High Discount Rate (Levitt, 2011))
D<-.070 # Wind Decomissioning (occurs at time WDT)
WDT<-30 # Wind decomissioning year

xy$wce<-W*(if (TT==3.6) TS else TL)+W*(if (TT==3.6) MF else JF)+W*.91*IC + xy$trnsc # Capex of wind farm equipment, including transmission. Foundations are 3.6MW = Monopile, 5.0MW = Jacketed
xy$wcapex<-xy$wce/(1-TI-TM) # Total capex with installation and misc costs

xy$womc<-xy$wcapex*TO # Annual O&M Costs 
xy$pv_womc<-(xy$womc/TD)*(1-(1/((1+TD)^WDT))) # Present value of wind O&M costs (assuming O&M is annually constant)
xy$pv_d<-(D*xy$wcapex)/((1+TD)^WDT)
xy$pv_costs<-xy$wcapex+xy$pv_womc+xy$pv_d

## Start here

# LCOE
# Revenue
# NPV
# Add cut-in and cut-out speeds here and below in the powerbreak function
# Normalize by area

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

