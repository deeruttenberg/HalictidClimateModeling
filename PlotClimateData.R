library(RMAWGEN)
library(RGENERATEPREC)
library(lubridate)
data(trentino)
NJData <- read.csv("https://www.ncei.noaa.gov/orders/cdo/2592648.csv")
NJData <- NJData[732:18993,]
PRECIPITATION$NJ <- NJData$PRCP
TEMPERATURE_MAX$NJ <- NJData$TMAX
TEMPERATURE_MIN$NJ <- NJData$TMIN

ColA <- read.csv("https://www.ncei.noaa.gov/orders/cdo/2592607.csv")
ColA <- ColA[ColA$STATION == "USR0000CTAY",]
ColA$Month <- format(as.Date(ColA$DATE), "%m")
ColA$Day <- format(as.Date(ColA$DATE), "%d")
ColA$Year <- format(as.Date(ColA$DATE), "%y")
year_min <- 1969
year_max <- 2000
for(i in 1:nrow(ColA)){
  ColA$Year[i] <- as.numeric(ColA$Year[i]) + 1900
  if(ColA$Year[i] < 1950){
    ColA$Year[i] <- as.numeric(ColA$Year[i]) + 100
  }
}
for(i in 1:18262){
  y <- as.numeric(TEMPERATURE_MAX$year[i])
  m <- as.numeric(TEMPERATURE_MAX$month[i])
  d <- as.numeric(TEMPERATURE_MAX$day[i])
  Row <- ColA[as.numeric(ColA$Year)==y & as.numeric(ColA$Day)==d & as.numeric(ColA$Month)==m,]
  if(nrow(Row) != 0){
    TEMPERATURE_MAX$COLA[i] <- Row$TMAX[1]
    TEMPERATURE_MIN$COLA[i] <- Row$TMIN[1]
  }
  else{
    TEMPERATURE_MAX$COLA[i] <- NA
    TEMPERATURE_MIN$COLA[i] <- NA
  }
}

ColB <- read.csv("https://www.ncei.noaa.gov/orders/cdo/2592607.csv")
ColB$Month <- format(as.Date(ColB$DATE), "%m")
ColB$Day <- format(as.Date(ColB$DATE), "%d")
ColB$Year <- format(as.Date(ColB$DATE), "%y")
ColB <- ColB[ColB$STATION == "USC00051959",]

year_min <- 1961
year_max <- 2001
for(i in 1:11872){
  ColB$Year[i] <- as.numeric(ColB$Year[i]) + 1900
  if(ColB$Year[i] < 1950){
    ColB$Year[i] <- as.numeric(ColB$Year[i]) + 100
  }
}
for(i in 1:18262){
  y <- as.numeric(TEMPERATURE_MAX$year[i])
  m <- as.numeric(TEMPERATURE_MAX$month[i])
  d <- as.numeric(TEMPERATURE_MAX$day[i])
  Row <- ColB[as.numeric(ColB$Year)==y & as.numeric(ColB$Day)==d & as.numeric(ColB$Month)==m,]
  if(nrow(Row) != 0){
    TEMPERATURE_MAX$COLB[i] <- Row$TMAX[1]
    TEMPERATURE_MIN$COLB[i] <- Row$TMIN[1]
  }
  else{
    TEMPERATURE_MAX$COLB[i] <- NA
    TEMPERATURE_MIN$COLB[i] <- NA
  }
}

origin <- paste(year_min,1,1,sep="-")
period <- TEMPERATURE_MAX$year>=year_min & TEMPERATURE_MAX$year<=year_max
Tx_mes <- TEMPERATURE_MAX[period,]
Tn_mes <- TEMPERATURE_MIN[period,]
station <- c("COLA")


col <- rainbow(n=12,start=0.1,end=0.9)
col[6:1] <- col[1:6]
col[7:12] <- col[12:7]
plot_sample(x=TEMPERATURE_MAX$COLA,
            sample="monthly",
            origin=origin,axes=FALSE,xlab="Tn
            [degC]",ylab="x",abline=NULL,col=col)
year_max_sim <- 2000
year_min_sim <- 2000
n_GPCA_iter <- 5
n_GPCA_iteration_residuals <- 5
p <- 1
nscenario=1