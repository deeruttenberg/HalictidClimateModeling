library(RMAWGEN)
library(RGENERATEPREC)
library(lubridate)
library(plyr)
library(readr)
data(trentino)
setwd("~/Downloads/")
files <- list.files(pattern = "station_data-19650204-19751214.csv", full.names = T,recursive=T) #List of files to concatenate
BritWeath <- ldply(files, read_csv)
Suss <- BritWeath[BritWeath$src_id == 803,]
Suss$Month <- format(as.Date(Suss$ob_end_time), "%m")
Suss$Day <- format(as.Date(Suss$ob_end_time), "%d")
Suss$Year <- format(as.Date(Suss$ob_end_time), "%y")
year_min <- 1966
year_max <- 1974
for(i in 1:nrow(Suss)){
  Suss$Year[i] <- as.numeric(Suss$Year[i]) + 1900
  if(Suss$Year[i] < 1950){
    Suss$Year[i] <- as.numeric(Suss$Year[i]) + 100
  }
}
for(i in 1:18262){
  y <- as.numeric(TEMPERATURE_MAX$year[i])
  m <- as.numeric(TEMPERATURE_MAX$month[i])
  d <- as.numeric(TEMPERATURE_MAX$day[i])
  Row <- Suss[as.numeric(Suss$Year)==y & as.numeric(Suss$Day)==d & as.numeric(Suss$Month)==m,]
  if(nrow(Row) != 0){
    TEMPERATURE_MAX$Suss[i] <- Row$max_air_temp[1]
    TEMPERATURE_MIN$Suss[i] <- Row$min_air_temp[1]
  }
  else{
    TEMPERATURE_MAX$Suss[i] <- NA
    TEMPERATURE_MIN$Suss[i] <- NA
  }
}
Hexh <- BritWeath[BritWeath$src_id == 297,]
Hexh$Month <- format(as.Date(Hexh$ob_end_time), "%m")
Hexh$Day <- format(as.Date(Hexh$ob_end_time), "%d")
Hexh$Year <- format(as.Date(Hexh$ob_end_time), "%y")
year_min <- 1965
year_max <- 1975
for(i in 1:nrow(Hexh)){
  Hexh$Year[i] <- as.numeric(Hexh$Year[i]) + 1900
  if(Hexh$Year[i] < 1950){
    Hexh$Year[i] <- as.numeric(Hexh$Year[i]) + 100
  }
}
for(i in 1:18262){
  y <- as.numeric(TEMPERATURE_MAX$year[i])
  m <- as.numeric(TEMPERATURE_MAX$month[i])
  d <- as.numeric(TEMPERATURE_MAX$day[i])
  Row <- Hexh[as.numeric(Hexh$Year)==y & as.numeric(Hexh$Day)==d & as.numeric(Hexh$Month)==m,]
  if(nrow(Row) != 0){
    TEMPERATURE_MAX$Hexh[i] <- Row$max_air_temp[1]
    TEMPERATURE_MIN$Hexh[i] <- Row$min_air_temp[1]
  }
  else{
    TEMPERATURE_MAX$Hexh[i] <- NA
    TEMPERATURE_MIN$Hexh[i] <- NA
  }
}
Inve <- BritWeath[BritWeath$src_id == 115,]
Inve$Month <- format(as.Date(Inve$ob_end_time), "%m")
Inve$Day <- format(as.Date(Inve$ob_end_time), "%d")
Inve$Year <- format(as.Date(Inve$ob_end_time), "%y")
year_min <- 1965
year_max <- 1975
for(i in 1:nrow(Inve)){
  Inve$Year[i] <- as.numeric(Inve$Year[i]) + 1900
  if(Inve$Year[i] < 1950){
    Inve$Year[i] <- as.numeric(Inve$Year[i]) + 100
  }
}
for(i in 1:18262){
  y <- as.numeric(TEMPERATURE_MAX$year[i])
  m <- as.numeric(TEMPERATURE_MAX$month[i])
  d <- as.numeric(TEMPERATURE_MAX$day[i])
  Row <- Inve[as.numeric(Inve$Year)==y & as.numeric(Inve$Day)==d & as.numeric(Inve$Month)==m,]
  if(nrow(Row) != 0){
    TEMPERATURE_MAX$Inve[i] <- Row$max_air_temp[1]
    TEMPERATURE_MIN$Inve[i] <- Row$min_air_temp[1]
  }
  else{
    TEMPERATURE_MAX$Inve[i] <- NA
    TEMPERATURE_MIN$Inve[i] <- NA
  }
}
  Dart <- BritWeath[BritWeath$src_id == 1350,]
Dart$Month <- format(as.Date(Dart$ob_end_time), "%m")
Dart$Day <- format(as.Date(Dart$ob_end_time), "%d")
Dart$Year <- format(as.Date(Dart$ob_end_time), "%y")
year_min <- 1965
year_max <- 1975
for(i in 1:nrow(Dart)){
  Dart$Year[i] <- as.numeric(Dart$Year[i]) + 1900
  if(Dart$Year[i] < 1950){
    Dart$Year[i] <- as.numeric(Dart$Year[i]) + 100
  }
}
for(i in 1:18262){
  y <- as.numeric(TEMPERATURE_MAX$year[i])
  m <- as.numeric(TEMPERATURE_MAX$month[i])
  d <- as.numeric(TEMPERATURE_MAX$day[i])
  Row <- Dart[as.numeric(Dart$Year)==y & as.numeric(Dart$Day)==d & as.numeric(Dart$Month)==m,]
  if(nrow(Row) != 0){
    TEMPERATURE_MAX$Dart[i] <- Row$max_air_temp[1]
    TEMPERATURE_MIN$Dart[i] <- Row$min_air_temp[1]
  }
  else{
    TEMPERATURE_MAX$Dart[i] <- NA
    TEMPERATURE_MIN$Dart[i] <- NA
  }
}

#SimulateData
origin <- paste(year_min,1,1,sep="-")
period <- TEMPERATURE_MAX$year>=year_min & TEMPERATURE_MAX$year<=year_max
Tx_mes <- TEMPERATURE_MAX[period,]
Tn_mes <- TEMPERATURE_MIN[period,]
station <- c("Hexh")


col <- rainbow(n=12,start=0.1,end=0.9)
col[6:1] <- col[1:6]
col[7:12] <- col[12:7]
plot_sample(x=TEMPERATURE_MAX$Dart,
            sample="monthly",
            origin=origin,axes=FALSE,xlab="Tn
            [degC]",ylab="x",abline=NULL,col=col)
year_max_sim <- 1970
year_min_sim <- 1970
n_GPCA_iter <- 5
n_GPCA_iteration_residuals <- 5
p <- 1
nscenario=1

for(i in 1:100){
  set.seed(i)
  generation00 <-ComprehensiveTemperatureGenerator(station=station, Tx_all=Tx_mes,Tn_all=Tn_mes,year_min=year_min,year_max=year_max,
                                                   p=p,n_GPCA_iteration=n_GPCA_iter,n_GPCA_iteration_residuals=n_GPCA_iteration_residuals, sample="monthly",year_min_sim=year_min_sim,year_max_sim=year_max_sim)
  year <- generation00$output$Tx_gen
  yearmin <- generation00$output$Tn_gen
  for(i in 5:365){
    if(year[i,] > 14 & year[i-1,] > 14 & year[i-2,] > 14 & year[i-3,] > 14 & year[i-4,] > 14){
      Start <- i
      break
    }
  }
  for(i in Start+50:365){
    if(year[i,] < 14 & year[i-1,] < 14 & year[i-2,] < 14 & year[i-3,] < 14 & year[i-4,] < 14){
      Finish <- i
      break
    }
  }
  print(Finish-Start)
}





