## Program Name: Windaqua
## Description: Simulates Monte Carlo outcomes of NPV for Mussel and Wind Energy firms in the German North Sea
## Date Started: 1/19/2013
## Date Modified: 2/27/2014

## Use this to start every program.  This clears out previous information from memory
rm(list=ls())
## Restricts numeric output 
options(digits=15)


### Packages
## Triangular distributions
#install.packages("triangle")
library("triangle")
## Alpha transparency for plotting
library("scales")
## Plotting
#install.packages("ggplot2")
library("ggplot2")


##### Scenarios
## Number of Turbines
W<-18
## Type of Turbine (3.6MW or 5.0MW) - foundation type is automatically configured
TT<-3.6
## Total Megawatt Capacity of Farm
MW<-W*TT
## Rotor Diameter
RD<-ifelse(TT==3.6,107,ifelse(TT==5.0,120,0))
## Time Horizon (years) - fixed for this analysis due to mussel equipment replacement
T<-20
## Number of 4 Plot groups. Also use for number of vessels (4 Plots = 1 Vessel)
F<-round(W/20)
## Energy Price, EUR per kWh
WP<- rep(.15,T) # .15 euro/kWh for t= 1, ...,12
WP[13]<-.059 # partial year extension of .15 euro/kWh rate based on distance from shore
WP[14:T]<-.035 # .035 euro/kWh for t= 12, ...,T

##### Random Variable Creation
# Number of runs for Monte Carlo
n<-1000
# Set Seed
set.seed(1000)


##### Parameterizing the Model
# Wind Speed in meters/sec
WS<-rtriangle(n,a=9.5,b=10.0,c=9.8)
# Wind Availability (%)
A<-rtriangle(n,a=.90,b=.97,c=.95)
# Wind Energy Losses (%)
EL<-rtriangle(n,a=.77,b=.98,c=.83)
# Wind 3.6MW Turbine Unit Cost 
TS<-rtriangle(n,a=5540000,b=6410000,c=6030000)
# Wind 5.0MW Turbine Unit Cost
TL<-rtriangle(n,a=10400000,b=10600000,c=10500000)
# Wind Infield Cable Cost per km
IC<-rtriangle(n,a=151000,b=305000,c=185000)
# Wind Monopile Foundation Unit Cost
MF<-rtriangle(n,a=1500000,b=1860000,c=1720000)
# Wind Jacketed Foundation Unit Cost
JF<-rtriangle(n,a=1820000,b=2060000,c=1940000)
# Wind Installaion Cost as a Percentage of CAPEX
TI<-rtriangle(n,a=.07,b=.20,c=.17)
# Wind Miscellaneous Costs as a Percentage of CAPEX
TM<-rtriangle(n,a=.02,b=.08,c=.05)
# Wind Operations and Management Costs as a Percentage of CAPEX per Year
TO<-rtriangle(n,a=.03,b=.035,c=.0325)
# Wind Weighted Average Cost of Capital
TD<-rtriangle(n,a=.096,b=.133,c=.116) #High Discount Rate (Levitt, 2011)
#TD<-rtriangle(n,a=.06,b=.09,c=.08)
# Wind Decomissioning (occurs at time T)
D<-rtriangle(n,a=.026,b=.070,c=.037)
# Mussel Biomass Harvested per 4 plots per year
B<-rtriangle(n,a=4757000,b=7135500,c=5946250)
# Mussel Discount Rate
MD<-rtriangle(n,a=.096,b=.133,c=.116) #High Discount Rate (Levitt, 2011)
#MD<-rtriangle(n,a=.06,b=.09,c=.08)
# Mussel Fuel Cost per Vessel per Year
FC<-rtriangle(n,a=92067,b=111128,c=101867)
# Mussel Infrastructure Cost
I<-rtriangle(n,a=3684392,b=3868612,c=3776502)
# Mussel Wages
MWAGE<-rtriangle(n,a=130904,b=134831,c=132867)
## Mussel Price, EUR per kWh (assumed fixed throughout horizon)
MP<-matrix(rnorm(n*T,0,.442205),T,n,TRUE)+(1.5819892+seq(0,T-1,1)*.048384)


## Other parameters
# Mussel License Costs (per 4 Plots)
L<-1102
# Mussel Land Facility Costs (per 4 Plots)
LF<-1653691
# Mussel Vessel Costs (one vessel per 4 Plots)
V<-4409843
# Mussel Motor Overhaul (occurs once per vessel in T=10)
M<-424447
# Mussel O&M and Misc Costs per year (15%/4years from Buck et al 2010)
MOM<-.0375
# Wind Transmission Distance in km
WD<-27

##### Wind Energy Costs (Euros)
### Transmissions Costs (Euros)
## AC Transmission Cost from Estimated Regression
AC<-MW*607997.1+WD*1022494 + rnorm(n,0,61995000)

## CAPEX Equipment (add AC if including transmission costs)
## Foundations are: 3.6MW = Monopile, 5.0MW = Jacketed
WCE<-W*(if (TT==3.6) TS else TL)+W*(if (TT==3.6) MF else JF)+W*.91*IC +AC
## Total CAPEX with Installation and Misc Costs
WCAPEX<-WCE/(1-TI-TM)
### PV of O&M Costs
# O&M Costs Vector
WOMC<-WCAPEX*TO
# O&M Costs Matrix for T years
WOMCT<-matrix(WOMC,T,n,TRUE)
## Discount Factors Matrix for T years
# Discount Rate Draws Matrix
WDF1<-matrix(TD,T,n,TRUE)
# Time Indicator Vector for T years
WDF2<-seq(1,T,1)
# Creating Discount Factor Matrix [Txn]
WDF<-1/((1+WDF1)^WDF2)
# Discounted O&M Values Matrix [Txn]
WDOM<-WOMCT*WDF
# Present Value of Discounted O&M for all n Draws
WPVOM<-colSums(WDOM)
## Total Present Value of Costs 
# Decomissioning Cost PV
WPVD<-D*WCAPEX/((1+TD)^T)
# Total Wind Cost PV
WCPV<-WCAPEX+WPVOM+WPVD

##### Wind Energy Generation (in kWh) and Revenue (in Euros)
# Energy Per Year Per Farm
kWh<-365.25*24*(TT*1000)*(.087*WS-((TT*1000)/RD^2))*A*EL*W
### Present Value of Revenue Stream to time T
# Energy Matrix for T years and n draws
WR<-matrix(kWh,T,n,TRUE)
# Revenue Matrix (Tx1)(Txn)[Txn]
WR1<-WP*WR
# Discounted Revenue Matrix [Txn]
WR2<-WR1*WDF
# Present Value of Revenue for all n Draws
WRPV<-colSums(WR2)

##### Net Present Value of Wind Farm
WNPV<-WRPV-WCPV
# Histogram of n Runs for Wind Net Present Value
hist(WNPV)

##### Mussel Farm Costs
### CAPEX Expenditures
# Year Zero (Only half of the plots are installed for the first year)
MCE1<-(V/2+I/2+L/2+LF/2)*F
# Present Value of Year One CAPEX
MPVCE2<-((V/2+I/2+L/2+LF/2)*F)/((1+MD)^1)
# Total PV of Initial CAPEX
MCAPEX0<-MCE1+MPVCE2
# Present Value of Vessel Motor Overhaul
MPVCE10<-(M*F)/((1+MD)^10)
# Present Value of Land Facility Maintenance Cost years 15-20
MPVLF<-(F*LF)/((1+MD)^16)+(F*LF)/((1+MD)^17)+(F*LF)/((1+MD)^18)+(F*LF)/((1+MD)^19)+(F*LF)/((1+MD)^20)
## Present Value of Replacement Plots (once every 4 years - half are done one year, the other half the next)
# Replacement Cost Matrix for Plots (Half are replaced at a time)
MPRCM<-matrix((F*I)/2,T,n,TRUE)
# Isolating Cost Periods for Plots
MPRCM[1:4,]=MPRCM[7:8,]=MPRCM[11:12,]=MPRCM[15:16,]=MPRCM[19:20,]=0
# Discount Rate Draws Matrix
MDF1<-matrix(MD,T,n,TRUE)
# Time Indicator Vector for T years
MDF2<-seq(1,T,1)
# Creating Discount Factor Matrix [Txn]
MDF<-1/((1+MDF1)^MDF2)
# Discounted Plot Replacement Matrix [Txn]
MDPRCM<-MPRCM*MDF
# Present Value of Plot Replacement for all n Draws
MPVPRC<-colSums(MDPRCM)
## Total Present Value of Fixed Costs
MPVFC<-MPVPRC+MCAPEX0+MPVCE10+MPVLF

### Variable Costs
## Fuel Costs
# Fuel Cost Matrix (first and last periods require half the vessel effort)
MFCM<-matrix((F*FC),T,n,TRUE)
MFCM[1,]=MFCM[T,]=(F*FC)/2
# Present Value of Fuel Costs
MPVF<-colSums(MFCM*MDF)
## Present Value of O&M and Misc Costs
MPVOM<-matrix(MOM*MCAPEX0,T,n,TRUE)
# Dividing first period by 2
Id<-matrix(1,T,n,TRUE)
Id[,1]<-.5
MPVOM<-MPVOM*Id
# NPV of O&M
MPVOM2<-colSums(MPVOM*MDF)
## Wages
MWM<-matrix((F*MWAGE),T,n,TRUE)
# Present Value of Wage Cost
MPVW<-colSums(MWM*MDF)
## Present Value of all Variable Costs
MPVVC<-MPVF+MPVOM2+MPVW
## Net Present Value of Mussel Costs
MPVC<-MPVVC+MPVFC

##### Mussel Production and Revenues
### Mussel Production Matrix (in kg), [Txn]
MBIO<-matrix((F*B),T,n,TRUE)
# Mussel production is zero in first period and last period (at least for T=20)
MBIO[1,]=0
# Mussel production is half in the last period (at least for T=20)  
MBIO[T,]=F*B/2
### Mussel Present Value of Revenue
MPVR<-colSums(MP*MDF*MBIO)

##### Mussel Net Present Value
MNPV<-MPVR-MPVC

##### Graphical Interpretation
# Setting up dataframe
analysis<-data.frame(WNPV,MNPV) #data frame
names(analysis)<-c("WNPV","MNPV") #variable names
### Graphs 
# Histograms of NPV
qplot(WNPV,data=analysis, geom="histogram")
ggplot(analysis,aes(MNPV)) + 
    geom_density(fill=alpha("darkblue",.5),color="black",binwidth=1000000) + 
    geom_density(data=analysis,aes(WNPV), fill=alpha("red",.5),color="black",binwidth=1000000)
## Box plots of NPV
# Need new dataframe! 
analysis2<-data.frame(c(MNPV,WNPV))
analysis2[1:n,2]<-"Mussels"
analysis2[n:(2*n),2]<-"Wind"
names(analysis2)<-c("NPV","Source")

ggplot(analysis2,aes(Source,NPV)) + geom_boxplot() 
qplot(Source, NPV, data=analysis2, geom="boxplot")

##### Data Analysis
## Mean NPV
mean(WNPV)
mean(MNPV)
## Standard Deviation of Monte Carlo Net Present Value Simulations
sd(WNPV)
sd(MNPV)

##### Cost sharing
### Wind
# Summary stats for Total Costs, Total Revenue, O&M
mean(WCAPEX)/1000000
sd(WCAPEX)/1000000
mean(WRPV)/1000000
sd(WRPV)/1000000
mean(WPVOM)/1000000
sd(WPVOM)/1000000
mean(WPVD)/1000000
sd(WPVD)/1000000
# O&M as a percent of total costs in present value
WPVOMCS<-WPVOM/WCPV
mean(WPVOMCS)
### Mussel 
mean(MPVR)/1000000
sd(MPVR)/1000000
#Vessel expenditures and motor overhaul
mean((V/2)*F+((V/2)*F)/(1+MD)+MPVCE10)/1000000
sd((V/2)*F+((V/2)*F)/(1+MD)+MPVCE10)/1000000
#Fuel Costs
mean(MPVF)/1000000
sd(MPVF)/1000000
#O&M Costs
mean(MPVOM2)/1000000
sd(MPVOM2)/1000000
#Wages
mean(MPVW)/1000000
sd(MPVW)/1000000
#All other costs
mean((MPVC-((V/2)*F+((V/2)*F)/(1+MD)+MPVCE10)-MPVF-MPVOM2-MPVW))/1000000
sd((MPVC-((V/2)*F+((V/2)*F)/(1+MD)+MPVCE10)-MPVF-MPVOM2-MPVW))/1000000

##### Mean cost of downtime in lost revenue
(mean(WR)/365)/80

##### Ttest that the pop mean of WNPV is sig different from zero
t.test(WNPV,mu=0) #Ho: mu=0