## Use this to start every program.  This clears out previous information from memory
rm(list=ls())

## Initalize renv for library lockfile
library(renv)
#renv::init()

## Packages
#Sys.setenv(RENV_PATHS_RTOOLS = "C:/rtools40/") # https://github.com/rstudio/renv/issues/225

# sudo apt-get install cargo # Needed for gifski package in ubuntu, run from terminal

PKG <- c("tidyverse","googledrive","gganimate","RColorBrewer","scales", "ggthemes","gifski")

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

## Read data
df<-read.csv("./Data/NE Lobster Commercial - Sheet1.csv")
df$Metric.tons<-as.numeric(gsub(",", "", df$Metric.tons)) # https://stackoverflow.com/questions/1523126/how-to-read-data-when-some-numbers-contain-commas-as-thousand-separator
df$State<-as.factor(df$State)

## Plot data (following https://www.r-bloggers.com/2021/05/animated-graph-gif-with-gganimate-ggplot/)
getPalette<-colorRampPalette(brewer.pal(6, "RdBu")) # function, change if different palette is desired
cols<-getPalette(6)
show_col(cols)
cols<-rev(cols)
cols<-c("MAINE" = cols[1] , "NEW HAMPSHIRE" = cols[2], "MASSACHUSETTS" = cols[3], "RHODE ISLAND" = cols[4], "CONNECTICUT" = cols[5], "NEW YORK" = cols[6])
# Remove two values (0 and NA) from dataset
df<-df %>% 
  na.exclude() %>% 
  filter(.,Metric.tons!=0)
# Percent change https://stackoverflow.com/questions/48196552/calculate-percentage-change-in-r-using-dplyr
df<-df %>%
  group_by(State) %>% 
  arrange(Year, .by_group = TRUE) %>%
  mutate(pct_change = ((Metric.tons/lag(Metric.tons) - 1)*100))
# Lagged percent change
df<-df %>%
  group_by(State) %>% 
  arrange(Year, .by_group = TRUE) %>%
  mutate(lag5_pct_change = ((lag(pct_change,n=5)+lag(pct_change,n=4)+lag(pct_change,n=3)+lag(pct_change,n=2)+lag(pct_change,n=1))/5))
df<-df %>%
  group_by(State) %>% 
  arrange(Year, .by_group = TRUE) %>%
  mutate(lag10_pct_change = ((lag(pct_change,n=5)+lag(pct_change,n=4)+lag(pct_change,n=3)+lag(pct_change,n=2)+lag(pct_change,n=1)+lag(pct_change,n=10)+lag(pct_change,n=9)+lag(pct_change,n=8)+lag(pct_change,n=7)+lag(pct_change,n=6))/10))

# Plot
lobplot<-ggplot(data = df[df$State!="NEW JERSEY"& df$Year>=1990,], aes(x = Year, y=lag10_pct_change, group=State)) + #df[df$State!="MAINE" & df$Year>=2000,]
  geom_line(aes(color=State), size=3) +
  geom_point(aes(color=State), size=1) +
  scale_color_manual(values = cols) + theme_tufte() + ylim(-20,30) + theme(text = element_text(size=20)) +
  geom_hline(yintercept = 0) +
  labs(y = "Percent change in lobster harvest (mt, 10 year average)", x = "Year") +
  transition_reveal(Year)

animate(lobplot, duration = 8, height=1080, width=1080, end_pause = 30, renderer = gifski_renderer("./lobster_change_1990.gif"))
