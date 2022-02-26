rm(list=ls())

## Initalize renv for library lockfile
library(renv)
#renv::init()

## Packages
#Sys.setenv(RENV_PATHS_RTOOLS = "C:/rtools40/") # https://github.com/rstudio/renv/issues/225

PKG<-c("raster","sf","tidyverse","rgdal","ncdf4","RColorBrewer","lattice","googledrive","tmap","chron","furrr","nngeo","EnvStats","gganimate","transformr","ggridges", "patchwork", "lmtest")

for (p in PKG) {
  if(!require(p,character.only = TRUE)) {  
    install.packages(p)
    require(p,character.only = TRUE)}
}
rm(p,PKG)

## Snapshot of libraries used
renv::snapshot()

options(collapse_mask = "manip") # https://twitter.com/grant_mcdermott/status/1493400952878952448?s=20&t=H9w3Azc04_LlJnqqk-56gQ

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

# Download files
a<-paste0("ftp://anonymous:anonymous@ftp2.psl.noaa.gov/pub/Public/dswales/wrfout_d01_",seq(2020,2070,1),".nc") # List of files on server. Only downloading files for 2020 - 2070, runs from 2007 - 2080
b<-seq(2020,2070,1) # years to download (ranges from 2007 - 2080)
c<-paste0("./Data/ncdf/wrfout_d01_",b,".nc") # List of names to save locally as

dl<-function(param1,param2){
  options(timeout=100000) # Keeps download.file from returning a timeout for long-downloading files
  download.file(url = param1, destfile = param2)}

plan(multisession(workers = nbrOfFreeWorkers())) # Parallel downloading
system.time(future_map2(a,c,dl)) # In RStudio, this may hang even though all files are downloaded

# Reading in non-wind data
cst<-st_read("./Data/global_polygon.gpkg")
cst<-st_transform(cst,crs = st_crs(3857))
eeza<-st_read("./Data/eez_atlantic.gpkg") # US EEZ vector https://www.marineregions.org/gazetteer.php?p=details&id=8456
eeza<-st_transform(eeza,st_crs(3857))
eeza<-eeza[1:2]

# Preparing data frames to append data to 
df<-nc_open("./Data/ncdf/wrfout_d01_2020.nc") # Uses a random layer, could be any year
#print(df)

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
days<-colMeans(matrix(day, nrow=8))

hour<-ncvar_get(df,"hour") # One observation every 3 hours
head(hour)

z<-ncvar_get(df,"level")
head(z)

x2<-as.matrix(expand.grid(x))
y2<-as.matrix(expand.grid(y))
xy<-as.data.frame(cbind(x2,y2))
rm(x2,y2)
names(xy)<-c("lon","lat")
xy<-st_as_sf(xy, coords = c("lon", "lat"),crs=4269, remove = FALSE) 
xy<-st_transform(xy,st_crs(3857)) # Grid of all point observations from the modeling

expand.grid.df<-function(...) Reduce(function(...) merge(..., by=NULL), list(...)) # https://stackoverflow.com/questions/11693599/alternative-to-expand-grid-for-data-frames
system.time(xyz<-expand.grid.df(xy,z,seq(1,365,1))) # Creating a long format dataframe that is # of grid points * # of hub heights * # of days
names(xyz)<-c("lon","lat","hub","day","geometry")

system.time(xy<-st_join(xy,eeza)) # Keeping only those observations that are within the US Atlantic EEZ
xy<-xy %>% drop_na()
xy<-xy %>% 
  select (-c(geoname, mrgid))

# ggplot() + # Plot of observation locations within domain
#   geom_sf(data = cst) +
#   geom_sf(data = xy)

# Function that measures shape and scale parameters of a Weibull distribution, and mean wind speeds, for a specified set of netcdfs, and aggregates values into a dataframe
annualwind<-function(param1,param2,param3){ 
  # param1: downloaded netcdfs file location list 
  # param2: vector of corresponding years represented in downloaded netcdfs, equal in length to param1
  # param3: desired hub height vector for netcdfs (options are 109, 135, 161, 187, 213, 239 meters), equal in length to param1
  df<-nc_open(param1)

  # Getting wind speed data
  wsu<-as.vector(ncvar_get(df,"u"))
  wsv<-as.vector(ncvar_get(df,"v"))

  # Creating a dataframe from whole array
  ws<-sqrt(wsu^2+wsv^2)
  #hist(ws, xlim = c(0, 100), breaks = seq(0,max(ws)+2,2))
  #max(ws)
  #length(ws[ws > 32.7])/length(ws) # Ratio of observations greater than hurricane speed
  wsday<-colMeans(matrix(ws, nrow=8)) # A vector of daily mean wind speeds, summarized from a vector of 3 hr wind speeds

  xyz2<-xyz
  xyz2$wsday<-wsday
  system.time(xyz2<-xyz2 %>% 
    filter(.,hub==param3)) # Filtering by chosen hub height

  system.time(xyz2<-st_join(xyz2,eeza)) # Keeping only those observations that are within the US Atlantic EEZ
  xyz2<-xyz2 %>% drop_na()
  xyz2<-xyz2 %>% 
    select (-c(geoname, mrgid))

  xyz2$llindex<-interaction(xyz2$lon,xyz2$lat,drop = TRUE) # Unique index value for each point
  
  llilevels<-levels(xyz2$llindex) # List of index identifiers

  multiweibull<-function(param4){
    xyz2<-xyz2 %>% filter(llindex == param4) # Subset for each point
    ll<-xyz2[1,1:2] %>% st_drop_geometry() # Grabbing lat lon
    mwresult<-eweibull(xyz2$wsday) # Estimate weibull parameters
    mwss<-as.data.frame(t(mwresult$parameters))
    mws<-as.data.frame(mean(xyz2$wsday)) # Mean annual wind speed
    names(mws)<-"mean_ws_yr"
    mwss<-cbind(ll,mwss,mws) # Attaching lat lon to Weibull parameters
    mwss$year<-param2
    return(mwss) # Return parameters
  }

  system.time(mw<-future_map_dfr(llilevels,multiweibull)) # Data frame of shape and scale parameters for all points
  return(mw)
}

ncs<-list.files("./Data/ncdf/", full.names = TRUE)
hub<-rep(109, times = length(ncs))

list3<-list(ncs,b,hub) # Prepping lists for pmap
system.time(wparam<-pmap_dfr(list3,annualwind)) # Resistant to parallel processing due to dataframe size. Takes ~7hrs to run

#write.csv2(wparam,"wparam.csv") # Save to have on hand in case to save time

# Animation of change in wind speed through time
wparam<-st_as_sf(wparam, coords = c("lon", "lat"),crs=4269, remove = FALSE) 
wparam<-st_transform(wparam,st_crs(3857))

getPalette<-colorRampPalette(brewer.pal(9, "Blues")) # Function, change if different palette from Colorbrewer is desired
cols<-getPalette(11) 
scales::show_col(cols)

wparam$mean_ws_yr_f<-cut(wparam$mean_ws_yr,
                          breaks = c(seq(2,12,1))) # Manually specifying breaks for symbology

windtime<-ggplot() +
  geom_sf(data=wparam, aes(color=mean_ws_yr_f, group=year), size=2) + # Don't use geom_point for sf data! Group here is key to just have changes appear in place, versus moving around # [wparam$year<=2022,]
  scale_color_manual(values = cols, labels = c("2-3","3-4","4-5","5-6","6-7","7-8","8-9","9-10","10-11","11-12"), name = "Mean Annual \nWind Speed (m/s) \n{closest_state}", drop = FALSE) +
  geom_sf(data = eeza, fill = "transparent", color = "black", inherit.aes = FALSE) +
  ggthemes::theme_map() +
  theme(legend.position = c(0.9,0.3), legend.key.size = unit(.5, 'cm'), text=element_text(size=16), panel.border = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "null")) + #legend.position = c(0.9,0.5)
  transition_states(year)
  #labs(title = "Year: {closest_state}")

animate(windtime, nframes = 112, end_pause = 10, height=900, width=1200, renderer = gifski_renderer("./windtime_2020_2070.gif"))

wrange<-wparam %>% 
  group_by(llindex) %>% 
  summarise(mean = mean(mean_ws_yr), sd = sd(mean_ws_yr))

wrange_plot<-ggplot(data=wparam, aes(x = mean_ws_yr, y = year2)) +
  geom_boxplot()

winddisttime<-ggplot(data=wparam, aes(x = mean_ws_yr, y = year2)) +
  geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01)

wrangeyear<-wparam %>% 
  group_by(year) %>% 
  summarise(mean = mean(mean_ws_yr), sd = sd(mean_ws_yr))

wrangeyear_plot_mean<-ggplot() +
  geom_line(data=wrangeyear, aes(x = year, y = mean))

wrangeyear_plot_sd<-ggplot() +
  geom_line(data=wrangeyear, aes(x = year, y = sd))

wrangeyear_plot_mean + wrangeyear_plot_sd # Patchwork

OLS_wind_year<-lm(mean_ws_yr ~ year, data = wparam)
summary(OLS_wind_year) # Slight evidence of an increase in wind speed through time (sig at 10% level)
bptest(OLS_wind_year) # Unlikely heteroskedastic through time (insig at 5% level using a Breusch Pagan test)



# # Valuation options
# 1. Change in (LCOE, NPV) locations through time, assuming constant technology, using USD 2020
# 2. NPV of all locations starting now, with a rebuild at decommissioning date
# 3. Change in (LCOE, NPV) locations through time, assuming changing technology, using USD 2020
# 4. Best locations to build out
# 

# Steps
# For each year file:
# Calculate Weibull parameters for each location
# Calculate annual mean wind speed
# Calculate annual wind energy production
# Row bind to XY dataframe  
# Return dataframe






