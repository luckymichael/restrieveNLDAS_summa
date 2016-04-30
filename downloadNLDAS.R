#!/sw/r-3.2.4/bin/Rscript
##!/usr/bin/Rscript
###!/opt/R/bin/Rscript --vanilla
#setwd("/Volumes/UUI/SUMMA")
#setwd("/Users/michaelou/OneDrive/!@postdoc/SUMMA/columbia")
#setwd("D:/!Cloud/OneDrive/!@postdoc/SUMMA/columbia")
#setwd("/home/mgou/uwhydro/summaProj/summaData/summa_columbia/annual_forcing")
#/home/mgou/summaProj/summaData/summa_columbia/annual_forcing
setwd("/usr/lusers/mgou/uwhydro/summaProj/summaData/summa_columbia/annual_forcing")
library("parallel")



#args <- c("2013-05-01","2013-05-01","32")
args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("Usage: columbia_forcing startDate endDate nProc", call.=FALSE)
} else if (length(args)==1) {
  # default output file and single thread
  args[2] = args[1]
  args[3] = "1"
}
# define the simulation perid
start.date <- as.Date(args[1])
end.date <-   as.Date(args[2])
nproc = as.integer(args[3])


# calculate the Julian day of the year
julianDay <- function(this.day) {
  # Julian day
  daysInMonth <- c(31,28,31,30,31,30,31,31,30,31,30,31)
  yy <- as.integer(format(this.day,"%Y"))
  mm <- as.integer(format(this.day,"%m"))
  dd <- as.integer(format(this.day,"%d"))
  if (((yy %% 4 == 0) & (yy %% 100 != 0)) | (yy %% 400 == 0)) daysInMonth[2] = 29
  if (mm>1) dd <- dd + sum(daysInMonth[1:(mm-1)])
  return(dd)
}


# define the url of NLDAS data (prefix) 
# sample:     ftp://hydro1.sci.gsfc.nasa.gov/data/s4pa/NLDAS/NLDAS_FORA0125_H.002/1994/005/NLDAS_FORA0125_H.A19940105.0000.002.grb
#             ftp://hydro1.sci.gsfc.nasa.gov/data/s4pa/NLDAS/NLDAS_FORA0125_H.002/2005/001/NLDAS_FORA0125_H.A200501010000.002.grb
#url.part1 <- "ftp://hydro1.sci.gsfc.nasa.gov/data/s4pa/NLDAS/NLDAS_FORA0125_H.002/"


# function to parse the date for a hourly loop to download NLDAS data
getThisDayGrib <- function(this.day){
  # Julian day
  dd <- julianDay(this.day) 
  
  jday <- sprintf("%03i",dd)
  
  #find file not exist
  hours <- sprintf("%02i",0:23)
  NLfname <- paste0("NLDAS_FORA0125_H.A",format(this.day,"%Y%m%d."),    hours,"00.002.grb")
  to.down <- NLfname[ ! file.exists(paste0("gribs/",NLfname)) ]
  for (down in to.down) {
    NLurl <- paste0("ftp://hydro1.sci.gsfc.nasa.gov/data/s4pa/NLDAS/NLDAS_FORA0125_H.002/",format(this.day,"%Y"),"/",jday,"/",down)
    if (capabilities("libcurl")) download.file(NLurl,down,quiet=FALSE,method="libcurl")  else download.file(NLurl,down,quiet=FALSE)
  }
  
  
  print(paste("Finish", format(this.day,"%Y%m%d")))
}


# loop to create all folders
for (mm in format(seq(start.date,end.date,by="month"),"%Y%m")){
  if (! dir.exists(mm)) dir.create(mm,recursive = TRUE)
}



# request the available cores
cl <- makeCluster(nproc, type="FORK")

# export shared variables and functions
clusterExport(cl, c("julianDay"))

### excute ###
all.date <- seq(start.date,end.date,by="day")
clusterMap(cl, getThisDayGrib, this.day = all.date, .scheduling = 'dynamic')

# close cluster run
stopCluster(cl)
stop("finished download grib")
