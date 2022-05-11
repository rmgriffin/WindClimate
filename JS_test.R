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

options(collapse_mask = "manip")

nc<-nc_open("wrf.winds.mon.clim.1979-2021.nc")

df<-ncvar_get(nc, "speed", start=c(1,1,1,1))
system.time(df<-melt(df))
df<-rename(df,"ws"="value")

y<-ncvar_get(nc,"lat")
y<-as.data.frame(y)
x<-ncvar_get(nc,"lon")
x<-as.data.frame(x)

y2<-melt(y)
y2$variable<-gsub("[a-zA-Z ]", "", y2$variable) # Remove prefix
y2$Var1<-rep(seq(1,121,1),121)
y2<-rename(y2,"Var2"="variable","y"="value")
y2$Var2<-as.numeric(y2$Var2)

x2<-melt(x)
x2$variable<-gsub("[a-zA-Z ]", "", x2$variable) # Remove prefix
x2$Var1<-rep(seq(1,121,1),121)
x2<-rename(x2,"Var2"="variable","x"="value")
x2$Var2<-as.numeric(x2$Var2)

df<-left_join(df,y2)
df<-left_join(df,x2)

lev<-cbind(as.data.frame(nc$dim$level$vals),seq(1,7,1))
colnames(lev)
lev<-rename(lev,"Var3"="seq(1, 7, 1)","hubhgt"="nc$dim$level$vals")

df<-left_join(df,lev, by = "Var3")

# df2<-df %>% 
#   group_by(y,x,Var3) %>% 
#   summarise(meanws=mean(ws))

df2<-df %>% 
  filter(Var4==1) %>% # Only data from January
  group_by(y,x,hubhgt) %>% 
  summarise(meanws=mean(ws))

# Reading in non-netcdf data
cst<-st_read("./Data/global_polygon.gpkg")
cst<-st_transform(cst,crs = st_crs(3857))
eeza<-st_read("./Data/eez_atlantic.gpkg") # US EEZ vector https://www.marineregions.org/gazetteer.php?p=details&id=8456
eeza<-st_transform(eeza,st_crs(3857))
eeza<-eeza[1:2]

# Assigning a projection to data
df2<-st_as_sf(df2, coords = c("x", "y"),crs=4269, remove = FALSE) #4269
df2<-st_transform(df2,st_crs(3857))

getPalette<-colorRampPalette(brewer.pal(9, "RdBu")) # Function, change if different palette from Colorbrewer is desired
cols<-getPalette(50)
cols<-rev(cols) # Reverse list elements
#scales::show_col(cols)
breaks4<-seq(.5,12.75,.25)
df2$meanwsf<-cut(df2$meanws, breaks=breaks4)


ggplot() +
  geom_sf(data=df2, aes(color=meanwsf, geometry = geometry), size=2) +
  scale_color_manual(values = cols, labels = levels(df2$meanwsf), name = "Mean Wind Speed (m/s)", drop = FALSE, guide = guide_legend(nrow = 2)) +
  geom_sf(data = eeza, fill = "transparent", color = "black", inherit.aes = FALSE) +
  ggthemes::theme_map() +
  theme(legend.position="bottom", legend.key.size = unit(0, 'cm'), legend.direction = "horizontal") +
  guides(fill = guide_legend(label.position = "bottom", title.hjust = 0.5)) +
  facet_wrap( ~ hubhgt, ncol=4)