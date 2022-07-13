rm(list=ls())

## Initalize renv for library lockfile
library(renv)
#renv::init()

## Packages
#Sys.setenv(RENV_PATHS_RTOOLS = "C:/rtools40/") # https://github.com/rstudio/renv/issues/225

PKG<-c("raster","terra","sf","tidyverse","rgdal","ncdf4","RNetCDF","RColorBrewer","lattice","googledrive","tmap","chron","furrr","nngeo","EnvStats","gganimate","transformr","ggridges", "patchwork", "lmtest", "reshape2", "exactextractr")

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

# a<-paste0("ftp://anonymous:anonymous@ftp2.psl.noaa.gov/pub/Public/dswales/wrfout_d01_",seq(2020,2099,1),".nc") # List of files on server. Only downloading files for 2020 - 2070, runs from 2007 - 2080
# b<-seq(2020,2099,1) # years to download (ranges from 2007 - 2080)
# c<-paste0("./Data/ncdf/wrfout_d01_",b,".nc") # List of names to save locally as
# 
# dl2<-function(param1,param2){
#   options(timeout=100000) # Keeps download.file from returning a timeout for long-downloading files
#   download.file(url = param1, destfile = param2)}
# 
# plan(multisession, workers = availableCores()) # Parallel downloading
# #system.time(future_map2(a,c,dl2)) # In RStudio, this may hang even though all files are downloaded

options(timeout=100000)
download.file(url = "ftp://anonymous@ftp.cdc.noaa.gov/pub/Public/jscott/RGriffin/wrf.winds.hgt.NEUS.zip", destfile = "./Data/ncdf/wrf.winds.hgt.NEUS.zip")

# Unzip code function -----------------------------------------------------
decompress_file <- function(file, directory, .file_cache = FALSE) {
  
  if (.file_cache == TRUE) {
    print("decompression skipped")
  } else {
    
    # Set working directory for decompression
    # simplifies unzip directory location behavior
    wd <- getwd()
    setwd(directory)
    
    # Run decompression
    decompression <-
      system2("unzip",
              args = c("-o", # include override flag
                       file),
              stdout = TRUE)
    
    # uncomment to delete archive once decompressed
    file.remove(file) 
    
    # Reset working directory
    setwd(wd); rm(wd)
    
    # Test for success criteria
    # change the search depending on 
    # your implementation
    if (grepl("Warning message", tail(decompression, 1))) {
      print(decompression)
    }
  }
}

system.time(decompress_file("./Data/ncdf/wrf.winds.hgt.NEUS.zip","./Data/ncdf/"))

# Reading in non-netcdf data
cst<-st_read("./Data/global_polygon.gpkg")
cst<-st_transform(cst,crs = st_crs(3857))
eeza<-st_read("./Data/eez_atlantic.gpkg") # US EEZ vector https://www.marineregions.org/gazetteer.php?p=details&id=8456
eeza<-st_transform(eeza,st_crs(3857))
eeza<-eeza[1:2]

# Reading in netcdf data and assembling a dataframe of relevant years
annualwind<-function(param1){ 
  # param1: netcdf file locations
  
  nc<-nc_open(param1) # Uses a random layer, could be any year
  
  df<-ncvar_get(nc, "speed", start=c(1,1,1,1))
  system.time(df<-melt(df))
  df<-rename(df,"ws"="value")
  
  y<-ncvar_get(nc,"lat")
  y<-as.data.frame(y)
  x<-ncvar_get(nc,"lon")
  x<-as.data.frame(x)
  
  y2<-melt(y)
  y2$variable<-gsub("[a-zA-Z ]", "", y2$variable) # Remove prefix
  y2$Var1<-rep(seq(1,37,1),58)
  y2<-rename(y2,"Var2"="variable","y"="value")
  y2$Var2<-as.numeric(y2$Var2)
  
  x2<-melt(x)
  x2$variable<-gsub("[a-zA-Z ]", "", x2$variable) # Remove prefix
  x2$Var1<-rep(seq(1,37,1),58)
  x2<-rename(x2,"Var2"="variable","x"="value")
  x2$Var2<-as.numeric(x2$Var2)
  
  df<-left_join(df,y2)
  df<-left_join(df,x2)
  
  lev<-cbind(as.data.frame(nc$dim$level$vals),seq(1,3,1))
  colnames(lev)
  lev<-rename(lev,"Var3"="seq(1, 3, 1)","hubhgt"="nc$dim$level$vals")
  df<-left_join(df,lev, by = "Var3")
  
  tm<-cbind(as.data.frame(utcal.nc(nc$dim$time$units, value = nc$dim$time$vals)),seq(1,length(nc$dim$time$vals),1)) # https://stackoverflow.com/questions/66376301/convert-from-minutes-since-origin-to-date-time-for-normal-people-yyyy-mm-dd-hh
  colnames(tm)[ncol(tm)]
  tm<-rename(tm,"Var4"=colnames(tm)[ncol(tm)])
  df<-left_join(df,tm, by = "Var4")
  
  df<-df %>% 
    select (-c(Var1, Var2, Var3, Var4, minute, second, month, day, hour)) %>% 
    filter(year==median(year)) # drops carryover days due to leap year from prior year
  
  # # Plot to check import
  # df2<-st_as_sf(df, coords = c("x", "y"),crs=4269, remove = FALSE) #4269
  # df2<-st_transform(df2,st_crs(3857))
  # 
  # getPalette<-colorRampPalette(brewer.pal(9, "RdBu")) # Function, change if different palette from Colorbrewer is desired
  # cols<-getPalette(50)
  # cols<-rev(cols) # Reverse list elements
  # #scales::show_col(cols)
  # breaks4<-seq(0,15,1)
  # df2$wsbreak<-cut(df2$ws, breaks=breaks4)
  # 
  # ggplot() +
  #   geom_sf(data=df2[df2$month==2&df2$day==1&df2$hour==3,], aes(color=wsbreak, geometry = geometry), size=2) +
  #   scale_color_manual(values = cols, labels = levels(df2$wsbreak), name = "Wind Speed (m/s)", drop = FALSE, guide = guide_legend(nrow = 2)) +
  #   geom_sf(data = eeza, fill = "transparent", color = "black", inherit.aes = FALSE) +
  #   ggthemes::theme_map() +
  #   theme(legend.position="bottom", legend.key.size = unit(0, 'cm'), legend.direction = "horizontal") +
  #   guides(fill = guide_legend(label.position = "bottom", title.hjust = 0.5)) +
  #   facet_wrap( ~ hubhgt, ncol=3)
  # 
  # df_1<-df2 %>% 
  #   filter(month==1) %>% # Only data from January
  #   group_by(y,x,hubhgt) %>% 
  #   summarise(meanws=mean(ws))
  
  # Weibull parameters for each lat/long/month/hubhgt combo
  #system.time(df$llindex<-interaction(df$x,df$y,df$hubhgt,df$year,drop = TRUE))
  #system.time(df$llindex<-as.factor(paste(df$x,df$y,df$hubhgt,df$year))) # Faster than "interaction" to create a grouping variable
  system.time(df$llindex<-as.factor(df$x*df$y*df$hubhgt)) # Faster than pasting to create a grouping variable (watch out for non-unique though)
  llilevels<-levels(df$llindex) # List of index identifiers
  
  multiweibull<-function(param2){
    # param2: index that identifies unique groups to estimate weibull parameters for
    df2<-df %>% filter(llindex == param2) # Subset for each point
    ll<-df2[1,2:3] # Grabbing lat lon
    mwresult<-eweibull(df2$ws) # Estimate weibull parameters
    mwss<-as.data.frame(t(mwresult$parameters))
    mws<-as.data.frame(mean(df2$ws)) # Mean wind speed
    names(mws)<-"mean_ws"
    mwss<-cbind(ll,mwss,mws) # Attaching lat lon to Weibull parameters
    mwss$year<-df2$year[1]
    mwss$days<-nrow(df2)/8 # Number of days
    mwss$hubhgt<-df2$hubhgt[1]
    
    return(mwss) # Return parameters
  }
  
  plan(multisession, workers = 3)
  options(future.globals.maxSize= 891289600000)
  system.time(mw<-future_map_dfr(llilevels,multiweibull)) # Data frame of shape and scale parameters for all points
  
  return(mw)
}

ncs<-list.files("./Data/ncdf/", full.names = TRUE)

plan(multisession, workers = 3)
options(future.globals.maxSize= 891289600000)
system.time(df<-future_map_dfr(ncs,annualwind))

write.csv2(df,"wp.csv") # Save to have on hand in case to save time
#df<-read.csv2("wp.csv")

# Checking years and time
yrtime<-function(param3){
  nc<-nc_open(param3)
  d<-data.frame(SVal=numeric(),
                EVal=numeric(),
                File=character(),
                Obs=numeric(),
                stringsAsFactors=FALSE)
  newrow<-data.frame(SVal=nc$dim$time$vals[1],EVal=nc$dim$time$vals[length(nc$dim$time$vals)],File=param3,Obs=length(nc$dim$time$vals))
  d<-rbind(d,newrow)
  tm<-as.data.frame(utcal.nc("minutes since 1947-07-01 00:00:00", value = nc$dim$time$vals))
  d<-cbind(d,tm[1,])

  return(d)
}

system.time(yt<-map_dfr(ncs,yrtime))

# Comparison of historical to future

df5022<-df %>% filter(year<=2022)
df2200<-df %>% filter(year>2022)

df5022a<-df5022 %>% 
  group_by(y,x,hubhgt) %>% 
  summarise(meanws=mean(mean_ws))

df2200a<-df2200 %>% 
  group_by(y,x,hubhgt) %>% 
  summarise(meanws=mean(mean_ws))

df5022a<-rename(df5022a,"meanws19502022"="meanws")
df2200a<-rename(df2200a,"meanws20222100"="meanws")

dfdelta<-left_join(df5022a,df2200a, by = c("y","x","hubhgt"))
dfdelta$wspctchg<-((dfdelta$meanws20222100-dfdelta$meanws19502022)/dfdelta$meanws19502022)*100

# Plot of change in mean wind speeds between historical and future periods
dfdelta<-st_as_sf(dfdelta, coords = c("x", "y"),crs=4269, remove = FALSE) 
dfdelta<-st_transform(dfdelta,st_crs(3857))

quantile(dfdelta$wspctchg, probs = seq(0, 1, 1/10))
getPalette<-colorRampPalette(brewer.pal(9, "RdBu")) # Function, change if different palette from Colorbrewer is desired
cols<-getPalette(16) 
cols<-rev(cols) # Reverse list elements
cols<-cols[1:11]
scales::show_col(cols)

dfdelta$wspctchgbreak<-cut(dfdelta$wspctchg, breaks = c(seq(-4,1.5,.5))) # Manually specifying breaks for symbology

ggplot() +
  geom_sf(data=dfdelta[dfdelta$hubhgt==109,], aes(color=wspctchgbreak, geometry = geometry), size=2) + # Don't use geom_point for sf data! Group here is key to just have changes appear in place, versus moving around # [wp$year<=2022,]
  scale_color_manual(values = cols, name = "Mean Annual \nWind Speed (m/s)", drop = FALSE) +
  geom_sf(data = eeza, fill = "transparent", color = "black", inherit.aes = FALSE) +
  ggthemes::theme_map() +
  theme(legend.position = c(0.9,0.3), legend.key.size = unit(.5, 'cm'), text=element_text(size=16), panel.border = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "null")) 

# st_write(dfdelta,"wp.gpkg")

rm(df2200,df2200a,df5022,df5022a,dfdelta)

# Animation of change in wind speed through time

wp<-st_as_sf(df, coords = c("x", "y"),crs=4269, remove = FALSE) 
wp<-st_transform(wp,st_crs(3857))

quantile(wp$mean_ws, probs = seq(0, 1, 1/10))
getPalette<-colorRampPalette(brewer.pal(9, "Blues")) # Function, change if different palette from Colorbrewer is desired
cols<-getPalette(15) 
scales::show_col(cols)

wp$mean_ws_f<-cut(wp$mean_ws, breaks = c(seq(5,12,0.5))) # Manually specifying breaks for symbology
labs<-c("5.0-5.5","5.5-6.0","6.0-6.5","6.5-7.0","7.0-7.5","7.5-8.0","8.0-8.5","8.5-9.0","9.0-9.5","9.5-10.0","10.0-10.5","10.5-11.0","11.0-11.5","11.5-12.0")

windtime<-ggplot() +
  geom_sf(data=wp[wp$year<=1952,], aes(color=mean_ws_f, group=year, geometry = geometry), size=2) + # Don't use geom_point for sf data! Group here is key to just have changes appear in place, versus moving around # [wp$year<=2022,]
  scale_color_manual(values = cols, labels = c("2-3","3-4","4-5","5-6","6-7","7-8","8-9","9-10","10-11","11-12"), name = "Mean Annual \nWind Speed (m/s) \n{closest_state}", drop = FALSE) +
  geom_sf(data = eeza, fill = "transparent", color = "black", inherit.aes = FALSE) +
  ggthemes::theme_map() +
  theme(legend.position = c(0.9,0.3), legend.key.size = unit(.5, 'cm'), text=element_text(size=16), panel.border = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "null")) + #legend.position = c(0.9,0.5)
  transition_states(year) +
  labs(title = "Year: {closest_state}")

# ggplot() + # Test for above animation
#   geom_sf(data=wp[wp$year==1980,], aes(color=mean_ws_f, geometry = geometry), size=2) + # Don't use geom_point for sf data! Group here is key to just have changes appear in place, versus moving around # [wp$year<=2022,]
#   scale_color_manual(values = cols, labels = labs, name = "Mean Annual \nWind Speed (m/s)", drop = FALSE) +
#   geom_sf(data = eeza, fill = "transparent", color = "black", inherit.aes = FALSE) +
#   ggthemes::theme_map() +
#   theme(legend.position = c(0.9,0.3), legend.key.size = unit(.5, 'cm'), text=element_text(size=16), panel.border = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "null")) #legend.position = c(0.9,0.5)

animate(windtime, nframes = 112, end_pause = 10, height=900, width=1200, renderer = gifski_renderer("./windtime_2020_2070.gif")) # Had to fiddle with nframes to make sure all years were displayed

# Other wind speed plots
wp$llindex<-interaction(wp$x,wp$y,drop = TRUE) # Unique index value for each point

wstats<-wp %>% 
  group_by(llindex) %>% 
  summarise(mean = mean(mean_ws), sd = sd(mean_ws))

wstats_plot<-ggplot(data=wp, aes(y = mean_ws, x = as.factor(year))) +
  geom_boxplot() 

winddisttime<-ggplot(data=wp[wp$year<=1983,], aes(x = mean_ws, y = as.factor(year))) +
  geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01) # Plot of distribution of mean less useful than raw data (or sampled data from Weibull)

wstatsyear<-wp %>% 
  group_by(year) %>% 
  summarise(mean = mean(mean_ws), sd = sd(mean_ws))

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

xy<-distinct_at(df,vars(x,y)) # Creating a vector layer with one observation per point
xy<-st_as_sf(xy, coords = c("x", "y"),crs=4269, remove = FALSE) 
xy<-st_transform(xy,st_crs(3857))

system.time(xy$nearest<-st_nearest_feature(xy,csteeza10kmbuff)) # Finding nearest feature to observation points. Merge using unique wind data points is much quicker than finding distance of nearest for all points for all years.
system.time(xy$dist<-as.vector(st_distance(xy, csteeza10kmbuff[xy$nearest,], by_element=TRUE))) # Distance to nearest feature
#csteeza10kmbuff<-csteeza10kmbuff[1,] # Only mainland connection points

globaldem<-rast("./Data/global_dem.tif") # Using terra for raster following https://geocompr.robinlovelace.net/reproj-geo-data.html
globaldem<-project(globaldem, "EPSG:3857", method = "bilinear")

xy$depth<-extract(globaldem,vect(xy))$global_dem # https://tmieno2.github.io/R-as-GIS-for-Economists/extracting-values-from-raster-layers-for-vector-data.html#extracting-to-points

W<-80 # Number of turbines
TT<-5.0 # Turbine rated power  
RD<-126 # Rotor diameter
A<-0.97 # Wind availability (% as fraction)
EL<-0.98 # Wind energy losses (% as fraction)
TS<-6410000 # Wind 3.6MW turbine unit cost 
TL<-10600000 # Wind 5.0MW turbine unit cost
IC<-305000 # Wind infield cable cost per km
MF<-1860000 # Wind monopile foundation unit cost
JF<-2060000 # Wind jacketed foundation unit cost
FF<-8173000 # Floating foundation unit cost (for ~5MW turbine - https://doi.org/10.3390/jmse8110835)
TI<-.20 # Wind installation cost as a percentage of CAPEX
TM<-.08 # Wind miscellaneous costs as a percentage of CAPEX
TO<-.035 # Wind operations and management costs as a percentage of CAPEX per year
TD<-.133 # Wind weighted average cost of capital (high discount rate (Levitt, 2011))
D<-.070 # Wind decommissioning (occurs at time WDT)
WDT<-30 # Wind decommissioning year
AD<-1.225 # Air density in kg/m3 at sea level
WI<-3 # Cut-in wind speed
WR<-12.5 # Rated wind speed
WO<-30 # Cut-out wind speed
HH<-109 # {109, 187, 239}

xy$trnsc<-ifelse(xy$dist>60000,810000*W*TT+1360000*xy$dist/1000,1090000*W*TT+890000*xy$dist/1000) # Transmission capital cost, uses different functions less/more than 60km. Values are hardcoded in 2012 USD from https://invest-userguide.readthedocs.io/en/latest/wind_energy.html
xy$wce<-W*(if (TT==3.6) TS else TL)+W*.91*IC + xy$trnsc + W*ifelse(xy$depth<=-60,FF,MF) # Capex of wind farm equipment, including transmission
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
wp$kWh<-365*24*((AD-(1.194*10^(-4))*wp$hubhgt)/AD)*TT*(wp$result+wp$result2)*EL*1000*W # Total kWh of farm per year (doesn't use wp$days instead of 365, but should if I can fix data importing for each year to get all days) 
#wp$kWh2<-wp$days*24*(TT*1000)*(.087*wp$mean_ws_yr-((TT*1000)/RD^2))*A*EL*W # Validation of above calculation, using alternate method

wp<-wp[wp$hubhgt==HH,]

windowsum<-function(param1,param2){ # Function to calculate lifetime energy production for a wind farm started in year "param1" for site "param2." Returns total kWh, discounted kWh for LCOE, lat/lon, and year
  # param1: year
  # param2: index
  
  wp2<-wp[wp$llindex==param2&wp$year>=param1&wp$year<param1+WDT,] # Extracting relevant years and location from dataframe
  #wp2<-wp[wp$llindex=="-82.2669067382812.23.8599700927734"&wp$year>=2020&wp$year<2020+WDT,]
  
  tdf<-wp2[1,] %>% select(x,y) # Extracting lat lon
  tdf$kWhD<-sum(wp2$kWh) # Sum of total energy generation
  wp2$df<-(1/((1+TD)^(wp2$year-param1))) # Discount factor through time
  tdf$kWhDLCOE<-sum(wp2$kWh*wp2$df) # Denominator of LCOE equation
  tdf$year<-param1
  return(tdf)
}

b<-unique(wp$year)

yearslist<-sort(rep(seq(b[1],b[length(b)],1),length(levels(wp$llindex))))

llindexlist<-rep(levels(wp$llindex),length(seq(b[1],b[length(b)],1)))

# plan(multisession, workers = 3)
# options(future.globals.maxSize= 891289600000)
system.time(kWhD<-map2_dfr(yearslist,llindexlist,windowsum))
kWhD<-kWhD %>% st_drop_geometry()
wp<-left_join(wp,kWhD,by = c("x","y","year"))
rm(kWhD)

# LCOE
xydf<-xy %>% st_drop_geometry()
wp<-merge(wp,xydf,by = c("x","y"))  # Merge costs (left_join doesn't work, presumably because of this https://stackoverflow.com/questions/61170525/simple-left-join-isnt-working-on-numeric-vector-to-join-by)
wp$LCOE<-wp$pv_costs/wp$kWhDLCOE # Calc LCOE
rm(xydf)

# Only within EEZ
wpeez<-st_join(wp,eeza)
wpeez<-wpeez %>% drop_na()

# Animation of LCOE through time
getPalette<-colorRampPalette(brewer.pal(9, "Blues")) # Function, change if different palette from Colorbrewer is desired
cols<-getPalette(11) 
#scales::show_col(cols)

breaks10LCOE<-quantile(wpeez$LCOE, probs = seq(0, 1, 1/10))
wpeez$LCOE_f<-cut(wpeez$LCOE, breaks = breaks10LCOE) # Manually specifying breaks for symbology

ggplot() +
  geom_sf(data=wpeez[wpeez$year==2022,], aes(color=LCOE_f, geometry = geometry), size=2) + # Don't use geom_point for sf data! Group here is key to just have changes appear in place, versus moving around # [wp$year<=2022,]
  scale_color_manual(values = cols, labels = levels(wpeez$LCOE_f), name = "Levelized Cost \nof Energy ($/kWh)", drop = FALSE) +
  geom_sf(data = eeza, fill = "transparent", color = "black", inherit.aes = FALSE) +
  ggthemes::theme_map() +
  theme(legend.position = c(0.9,0.3), legend.key.size = unit(.5, 'cm'), text=element_text(size=16), panel.border = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "null"))  #legend.position = c(0.9,0.5)

LCOEtime<-ggplot() +
  geom_sf(data=wpeez[wpeez$year<=2070,], aes(color=LCOE_f, group=year, geometry = geometry), size=2) + # Don't use geom_point for sf data! Group here is key to just have changes appear in place, versus moving around # [wp$year<=2022,]
  scale_color_manual(values = cols, labels = levels(wpeez$LCOE_f), name = "Levelized Cost \nof Energy ($/kWh) \n{closest_state}", drop = FALSE) +
  geom_sf(data = eeza, fill = "transparent", color = "black", inherit.aes = FALSE) +
  ggthemes::theme_map() +
  theme(legend.position = c(0.9,0.3), legend.key.size = unit(.5, 'cm'), text=element_text(size=16), panel.border = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "null")) + #legend.position = c(0.9,0.5)
  transition_states(year)
#labs(title = "Year: {closest_state}")

animate(LCOEtime, nframes = 121, end_pause = 10, height=900, width=1200, renderer = gifski_renderer("./lcoetime_1950_2070.gif")) # Had to fiddle with nframes to make sure all years were displayed

wp5MW<-wpeez %>% 
  select(x,y,shape,scale,mean_ws,year,hubhgt,kWh,depth,LCOE)

write.csv(wp5MW,"wp5MW.csv")
# # Valuation options
# 1. Change in (LCOE) locations through time, assuming constant technology, using USD 2020 (normalize by area to capture economies of scale)
# 2. Net present value calculation of build out, inclusive of externalities?
# #  Below require information about changing technology and market conditions
# 3. NPV of all locations starting now, with a rebuild at decommissioning date
# 4. Change in (LCOE, NPV) locations through time, assuming changing technology, using USD 2020
# 5. Best locations to build out
# 6. 







