rm(list=ls())

## Initalize renv for library lockfile
library(renv)
#renv::init()

## Packages
#Sys.setenv(RENV_PATHS_RTOOLS = "C:/rtools40/") # https://github.com/rstudio/renv/issues/225

PKG<-c("raster","sf","tidyverse","rgdal","ncdf4","RColorBrewer","lattice","googledrive","tmap","chron","furrr","nngeo","EnvStats","gganimate","transformr","ggridges", "patchwork", "lmtest", "reshape2")

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

# Download netcdf files
a<-paste0("ftp://anonymous:anonymous@ftp2.psl.noaa.gov/pub/Public/dswales/wrfout_d01_",seq(2020,2099,1),".nc") # List of files on server. Only downloading files for 2020 - 2070, runs from 2007 - 2080
b<-seq(2020,2099,1) # years to download (ranges from 2007 - 2080)
c<-paste0("./Data/ncdf/wrfout_d01_",b,".nc") # List of names to save locally as

dl2<-function(param1,param2){
  options(timeout=100000) # Keeps download.file from returning a timeout for long-downloading files
  download.file(url = param1, destfile = param2)}

plan(multisession, workers = availableCores()) # Parallel downloading
#system.time(future_map2(a,c,dl2)) # In RStudio, this may hang even though all files are downloaded

# Reading in non-netcdf data
cst<-st_read("./Data/global_polygon.gpkg")
cst<-st_transform(cst,crs = st_crs(3857))
eeza<-st_read("./Data/eez_atlantic.gpkg") # US EEZ vector https://www.marineregions.org/gazetteer.php?p=details&id=8456
eeza<-st_transform(eeza,st_crs(3857))
eeza<-eeza[1:2]

# Preparing data frames to append data to 
df<-nc_open("./Data/ncdf/wrfout_d01_2020.nc") # Uses a random layer, could be any year
#dfu<-brick("./Data/ncdf/wrfout_d01_2020.nc", varname = "u", level = 1) # Don't want to use raster as I have irregularly spaced points

dfu<-ncvar_get(df, "u", start=c(1,1,1,1))
system.time(dfu<-melt(dfu))
system.time(dfu<-dfu %>% 
              filter(.,Var3==1))
dfu<-rename(dfu,"u"="value")

dfv<-ncvar_get(df, "v", start=c(1,1,1,1))
system.time(dfv<-melt(dfv))
system.time(dfv<-dfv %>% 
              filter(.,Var3==1))
dfv<-rename(dfv,"v"="value")

df<-cbind(dfu,dfv$v)
df<-rename(df,"v"="dfv$v")

rm(dfu,dfv)

y<-ncvar_get(df,"lat")
y<-as.data.frame(y)
x<-ncvar_get(df,"lon")
x<-as.data.frame(x)

y2<-melt(y)
y2$variable<-gsub("[a-zA-Z ]", "", y2$variable) # Remove prefix
y2$Var1<-rep(seq(1,203,1),194)
y2<-rename(y2,"Var2"="variable","y"="value")
y2$Var2<-as.numeric(y2$Var2)

x2<-melt(x)
x2$variable<-gsub("[a-zA-Z ]", "", x2$variable) # Remove prefix
x2$Var1<-rep(seq(1,203,1),194)
x2<-rename(x2,"Var2"="variable","x"="value")
x2$Var2<-as.numeric(x2$Var2)

df<-left_join(df,y2)
df<-left_join(df,x2)

df$ws<-sqrt(df$u^2+df$v^2)
df<-df %>% 
  select (Var4,y,x,ws)

# Look at one time slice
pts<-203*194
pts2<-2*pts
df2<-df[1:pts,]
df2<-st_as_sf(df2, coords = c("x", "y"),crs=4269, remove = FALSE) #4269
df2<-st_transform(df2,st_crs(3857))

system.time(df2<-df2 %>% 
              filter(.,ws!=0))

ggplot() +
  geom_sf(data=df2, aes(color=ws, geometry = geometry), size=2) + 
  #geom_sf(data = eeza, fill = "transparent", color = "black", inherit.aes = FALSE) +
  ggthemes::theme_map()










x2<-as.matrix(expand.grid(x))
y2<-as.matrix(expand.grid(y))
xy<-as.data.frame(cbind(x2,y2))
rm(x2,y2)
names(xy)<-c("lon","lat")
# xy<-st_as_sf(xy, coords = c("lon", "lat"),crs=4269, remove = FALSE) 
# xy<-st_transform(xy,st_crs(3857)) # Grid of all point observations from the modeling

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
  # param1: netcdf file location
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
  yrdays<-ncvar_get(df,"day")
  yrdays<-colMeans(matrix(yrdays, nrow=8))
  yrdays<-length(yrdays) # Number of days per year
  
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
    # param4: index that identifies unique groups to estimate weibull parameters for
    xyz2<-xyz2 %>% filter(llindex == param4) # Subset for each point
    ll<-xyz2[1,1:2] %>% st_drop_geometry() # Grabbing lat lon
    mwresult<-eweibull(xyz2$wsday) # Estimate weibull parameters
    mwss<-as.data.frame(t(mwresult$parameters))
    mws<-as.data.frame(mean(xyz2$wsday)) # Mean annual wind speed
    names(mws)<-"mean_ws_yr"
    mwss<-cbind(ll,mwss,mws) # Attaching lat lon to Weibull parameters
    mwss$year<-param2
    mwss$days<-yrdays # Accounting for leap years (this data doesn't appear to include them though)
    return(mwss) # Return parameters
  }

  system.time(mw<-future_map_dfr(llilevels,multiweibull)) # Data frame of shape and scale parameters for all points
  return(mw)
}

ncs<-list.files("./Data/ncdf/", full.names = TRUE)
hub<-rep(109, times = length(ncs))

list3<-list(ncs[1],b[1],hub[1]) # Prepping lists for pmap
system.time(wp<-pmap_dfr(list3,annualwind)) # Weibull parameter estimation for all netcdf files. Resistant to parallel processing due to dataframe size. Takes ~8hrs to run
rm(xyz)
#write.csv2(wp,"wp.csv") # Save to have on hand in case to save time
#wp<-read.csv2("wp.csv")

# Animation of change in wind speed through time
wp<-st_as_sf(wp, coords = c("lon", "lat"),crs=4269, remove = FALSE) 
wp<-st_transform(wp,st_crs(3857))

quantile(wp$mean_ws_yr, probs = seq(0, 1, 1/10))
getPalette<-colorRampPalette(brewer.pal(9, "Blues")) # Function, change if different palette from Colorbrewer is desired
cols<-getPalette(11) 
scales::show_col(cols)

wp$mean_ws_yr_f<-cut(wp$mean_ws_yr, breaks = c(seq(2,12,1))) # Manually specifying breaks for symbology

windtime<-ggplot() +
  geom_sf(data=wp[wp$year<=2070,], aes(color=mean_ws_yr_f, group=year, geometry = geometry), size=2) + # Don't use geom_point for sf data! Group here is key to just have changes appear in place, versus moving around # [wp$year<=2022,]
  scale_color_manual(values = cols, labels = c("2-3","3-4","4-5","5-6","6-7","7-8","8-9","9-10","10-11","11-12"), name = "Mean Annual \nWind Speed (m/s) \n{closest_state}", drop = FALSE) +
  geom_sf(data = eeza, fill = "transparent", color = "black", inherit.aes = FALSE) +
  ggthemes::theme_map() +
  theme(legend.position = c(0.9,0.3), legend.key.size = unit(.5, 'cm'), text=element_text(size=16), panel.border = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "null")) + #legend.position = c(0.9,0.5)
  transition_states(year)
  #labs(title = "Year: {closest_state}")

animate(windtime, nframes = 112, end_pause = 10, height=900, width=1200, renderer = gifski_renderer("./windtime_2020_2070.gif")) # Had to fiddle with nframes to make sure all years were displayed

# Other wind speed plots
wp$llindex<-interaction(wp$lon,wp$lat,drop = TRUE) # Unique index value for each point

wstats<-wp %>% 
  group_by(llindex) %>% 
  summarise(mean = mean(mean_ws_yr), sd = sd(mean_ws_yr))

wstats_plot<-ggplot(data=wp, aes(y = mean_ws_yr, x = as.factor(year))) +
  geom_boxplot()

winddisttime<-ggplot(data=wp[wp$year<=2023,], aes(x = mean_ws_yr, y = as.factor(year))) +
  geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01)

wstatsyear<-wp %>% 
  group_by(year) %>% 
  summarise(mean = mean(mean_ws_yr), sd = sd(mean_ws_yr))

wstatsyear_plot_mean<-ggplot() +
  geom_line(data=wstatsyear, aes(x = year, y = mean))

wstatsyear_plot_sd<-ggplot() +
  geom_line(data=wstatsyear, aes(x = year, y = sd))

wstatsyear_plot_mean + wstatsyear_plot_sd # Patchwork

# Exploratory regression to look at trends in wind speed through time in the region
lm_wind_year<-lm(mean_ws_yr ~ year, data = wp) # Linear regression model
summary(lm_wind_year) # Slight evidence of an increase in wind speed through time (sig at 10% level)
bptest(lm_wind_year) # Unlikely heteroskedastic through time (insig at 5% level using a Breusch Pagan test)

# Costs
eeza10kmbuff<-st_buffer(eeza,dist = 10000)
csteeza10kmbuff<-st_intersection(cst,eeza10kmbuff) # Creating a limited nearshore polygon to search over for nearest landing point for wind points

system.time(xy$nearest<-st_nearest_feature(xy,csteeza10kmbuff)) # Finding nearest feature to observation points. Merge using unique wind data points is much quicker than finding distance of nearest for all points for all years.
system.time(xy$dist<-as.vector(st_distance(xy, csteeza10kmbuff[xy$nearest,], by_element=TRUE))) # Distance to nearest feature
#csteeza10kmbuff<-csteeza10kmbuff[1,] # Only mainland connection points

W<-80 # Number of turbines
TT<-3.6 # Turbine rated power  
RD<-120 # Rotor diameter
A<-0.97 # Wind availability (% as fraction)
EL<-0.98 # Wind energy losses (% as fraction)
TS<-6410000 # Wind 3.6MW turbine unit cost 
TL<-10600000 # Wind 5.0MW turbine unit cost
IC<-305000 # Wind infield cable cost per km
MF<-1860000 # Wind monopile foundation unit cost
JF<-2060000 # Wind jacketed foundation unit cost
TI<-.20 # Wind installation cost as a percentage of CAPEX
TM<-.08 # Wind miscellaneous costs as a percentage of CAPEX
TO<-.035 # Wind operations and management costs as a percentage of CAPEX per year
TD<-.133 # Wind weighted average cost of capital (high discount rate (Levitt, 2011))
D<-.070 # Wind decommissioning (occurs at time WDT)
WDT<-30 # Wind decommissioning year
AD<-1.225 # Air density in kg/m3 at sea level
WI<-4 # Cut-in wind speed
WR<-12.5 # Rated wind speed
WO<-25 # Cut-out wind speed

xy$trnsc<-ifelse(xy$dist>60000,810000*W*TT+1360000*xy$dist/1000,1090000*W*TT+890000*xy$dist/1000) # Transmission capital cost, uses different functions less/more than 60km. Values are hardcoded in 2012 USD from https://invest-userguide.readthedocs.io/en/latest/wind_energy.html
xy$wce<-W*(if (TT==3.6) TS else TL)+W*(if (TT==3.6) MF else JF)+W*.91*IC + xy$trnsc # Capex of wind farm equipment, including transmission. Foundations are 3.6MW = Monopile, 5.0MW = Jacketed
xy$wcapex<-xy$wce/(1-TI-TM) # Total capex with installation and misc costs

xy$womc<-xy$wcapex*TO # Annual O&M Costs 
xy$pv_womc<-(xy$womc/TD)*(1-(1/((1+TD)^WDT))) # Present value of wind O&M costs (assuming O&M is annually constant)
xy$pv_d<-(D*xy$wcapex)/((1+TD)^WDT) # Present value of decommissioning cost
xy$pv_costs<-xy$wcapex+xy$pv_womc+xy$pv_d # Total present value of costs for a wind farm built at this time (not adjusted for future inflation)

# Energy generation
weib<-function(x, shape, scale) {(shape/scale)*((x/scale)^(shape-1))*(exp(-1*((x/scale)^shape)))} # Weibull distribution
# integrate(weib, shape = 2, scale = 6, lower = 0, upper = Inf) # Test function is coded and working correctly, should equal 1
weib2<-function(x, shape, scale, cut_in, wrated) {((x-cut_in)/(wrated-cut_in))*(shape/scale)*((x/scale)^(shape-1))*(exp(-1*((x/scale)^shape)))} # Weibull distribution with polynomial 

result<-vector("list",nrow(wp)) # https://stackoverflow.com/questions/50705751/writing-a-loop-and-the-integrate-function
for (i in 1:(nrow(wp))){
  result[[i]] <- integrate(weib, shape = wp$shape[i], scale = wp$scale[i], lower = WR, upper = WO)$value # $value saves only the estimated value from the integral
}
result<-as.vector(unlist(result))

result2<-vector("list", nrow(wp))
for (i in 1:(nrow(wp))){
  result2[[i]] <- integrate(weib2, shape = wp$shape[i], scale = wp$scale[i], lower = WI, upper = WR, cut_in = WI, wrated = WR)$value
}
result2<-as.vector(unlist(result2))

wp$result<-result
wp$result2<-result2
wp$kWh<-wp$days*24*((AD-(1.194*10^(-4))*hub[1])/AD)*TT*(wp$result+wp$result2)*EL*1000*W # Total kWh of farm per year 
#wp$kWh2<-wp$days*24*(TT*1000)*(.087*wp$mean_ws_yr-((TT*1000)/RD^2))*A*EL*W # Validation of above calculation, using alternate method

windowsum<-function(param1,param2){ # Function to calculate lifetime energy production for a wind farm started in year "param1" for site "param2." Returns total kWh, discounted kWh for LCOE, lat/lon, and year
  # param1: year
  # param2: index
  
  wp2<-wp[wp$llindex==param2&wp$year>=param1&wp$year<param1+WDT,] # Extracting relevant years and location from dataframe
  #wp2<-wp[wp$llindex=="-82.2669067382812.23.8599700927734"&wp$year>=2020&wp$year<2020+WDT,]
  
  tdf<-wp2[1,] %>% select(lon,lat) # Extracting lat lon
  tdf$kWhD<-sum(wp2$kWh) # Sum of total energy generation
  wp2$df<-(1/((1+TD)^(wp2$year-param1))) # Discount factor through time
  tdf$kWhDLCOE<-sum(wp2$kWh*wp2$df) # Denominator of LCOE equation
  tdf$year<-param1
  return(tdf)
}

yearslist<-sort(rep(seq(b[1],b[length(b)],1),length(levels(wp$llindex))))
llindexlist<-rep(levels(wp$llindex),length(seq(b[1],b[length(b)],1)))

system.time(kWhD<-future_map2_dfr(yearslist,llindexlist,windowsum))
kWhD<-kWhD %>% st_drop_geometry()
wp<-left_join(wp,kWhD,by = c("lon","lat","year"))
rm(kWhD)

# LCOE
xydf<-xy %>% st_drop_geometry()
wp<-merge(wp,xydf,by = c("lon","lat"))  # Merge costs (left_join doesn't work, presumably because of this https://stackoverflow.com/questions/61170525/simple-left-join-isnt-working-on-numeric-vector-to-join-by)
wp$LCOE<-wp$pv_costs/wp$kWhDLCOE # Calc LCOE
rm(xydf)

# Animation of LCOE through time
getPalette<-colorRampPalette(brewer.pal(9, "Blues")) # Function, change if different palette from Colorbrewer is desired
cols<-getPalette(11) 
#scales::show_col(cols)

breaks10LCOE<-quantile(wp$LCOE, probs = seq(0, 1, 1/10))
wp$LCOE_f<-cut(wp$LCOE, breaks = breaks10LCOE) # Manually specifying breaks for symbology

LCOEtime<-ggplot() +
  geom_sf(data=wp, aes(color=LCOE_f, group=year, geometry = geometry), size=2) + # Don't use geom_point for sf data! Group here is key to just have changes appear in place, versus moving around # [wp$year<=2022,]
  scale_color_manual(values = cols, labels = levels(wp$LCOE_f), name = "Levelized Cost \nof Energy ($/kWh) \n{closest_state}", drop = FALSE) +
  geom_sf(data = eeza, fill = "transparent", color = "black", inherit.aes = FALSE) +
  ggthemes::theme_map() +
  theme(legend.position = c(0.9,0.3), legend.key.size = unit(.5, 'cm'), text=element_text(size=16), panel.border = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "null")) + #legend.position = c(0.9,0.5)
  transition_states(year)
#labs(title = "Year: {closest_state}")

animate(LCOEtime, nframes = 112, end_pause = 10, height=900, width=1200, renderer = gifski_renderer("./lcoetime_2020_2070.gif")) # Had to fiddle with nframes to make sure all years were displayed


# # Valuation options
# 1. Change in (LCOE) locations through time, assuming constant technology, using USD 2020 (normalize by area to capture economies of scale)
# 2. Net present value calculation of build out, inclusive of externalities?
# #  Below require information about changing technology and market conditions
# 3. NPV of all locations starting now, with a rebuild at decommissioning date
# 4. Change in (LCOE, NPV) locations through time, assuming changing technology, using USD 2020
# 5. Best locations to build out
# 6. 







