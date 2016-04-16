#!/opt/R/bin/Rscript --vanilla
#setwd("/Volumes/UUI/SUMMA")
#setwd("/Users/michaelou/OneDrive/!@postdoc/SUMMA/columbia")
#setwd("D:/!Cloud/OneDrive/!@postdoc/SUMMA/columbia")
setwd("/home/mgou/summaData/summa_columbia/annual_forcing")
library("foreign")
library("foreach");library("doParallel")
library("rgdal")
library("ncdf4")
library("parallel")
library("zoo")

args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[2] = args[1]
}
# define the simulation perid
start.date <- as.Date(args[1])
end.date <-   as.Date(args[2])
min.frac.area <- 5000 # filter those fractions smaller than this value square meters


NLDAS.nx <- 464; NLDAS.ny <- 224

# dir()

# read hru-nldas fraction dbf
frac <- read.dbf("HRU_NLDAS.dbf")
# clean and filter small fractions
frac <- frac[!(is.na(frac$hru_id2) | is.na(frac$fracArea)),]
frac <- subset(frac, frac$fracArea > min.frac.area)
# calculate the areal weight
hru.area <- aggregate(frac$fracArea, list(hruID <- frac$hru_id2), FUN = sum)
nhru <- nrow(hru.area)

frac$frac <- frac$fracArea/hru.area$x[match(frac$hru_id2,hru.area[,1])]

# calculate the 1D index of the cell ( from west to east then north to south !!! )
frac$idx <- NLDAS.nx*(NLDAS.ny-frac$NLDAS_Y) + frac$NLDAS_X



# define dimensions of output netcdf
dimList <- list(
  dimHRU  <-ncdim_def(name="hruId", units="", vals=hru.area[,1]),
  dimTime <- ncdim_def(name="time",  units="days since 1900-01-01", vals=0, unlim=TRUE, calendar="standard", longname="Observation time")
)


# define variables of output netcdf
varList <-list(
  nc.TMP <- ncvar_def( name="airtemp",  units="K",       dim=dimList, missval=-999., longname="Air temperature at the measurement height", prec="double"),
  nc.SPFH <- ncvar_def( name="spechum", units="kg/kg",   dim=dimList, missval=-999., longname="Specific humidity at the measurement height", prec="double"),
  nc.PRES <- ncvar_def( name="airpres", units="Pa",      dim=dimList, missval=-999., longname="Pressure at the measurement height", prec="double"),
  nc.WIND <- ncvar_def( name="windspd", units="m/s",     dim=dimList, missval=-999., longname="Wind speed at the measurement height", prec="double"),
  nc.DLWRF <- ncvar_def(name="LWRadAtm",units="W/m^2",   dim=dimList, missval=-999., longname="Longwave radiation at the upper boundary", prec="double"),
  nc.PCP <- ncvar_def( name="pptrate",  units="kg/m^2/s",dim=dimList, missval=-999., longname="Precipitation rate", prec="double"),
  nc.DSWRF <- ncvar_def(name="SWRadAtm",units="W/m^2",   dim=dimList, missval=-999., longname="Short radiation at the upper boundary", prec="double")
)





convertNLDAS2NC <- function(NLurl,NLfname,NCfname,t){
  # this function download the NLDAS file and convert it to a NC file for HRUs
  # NLurl   -- the link of the NLDAS grib file
  # NLfname -- a string of the file name to save the NLDAS grib file on the local disk
  # NCfname -- a string of the file name of the output netCDF
  # t       -- the time since of 1900-01-01 00:00 for the variables represented in this NLDAS file

  # download the NLDAS file
  if (capabilities("libcurl")) download.file(NLurl,NLfname,quiet=FALSE,method="libcurl")  else download.file(NLurl,NLfname,quiet=FALSE)

  # read the grib file, alternative method  using rNOMADS package: grib <- ReadGrib(file.name = NLfname, levels = 1,variables = c("TMP2m"))
  grib <- readGDAL(NLfname)
  # if want to plot use: image(grib, attr=1)
  
  
  # NLDAS grib file strcture 
  # be careful of the dimensions of the file: from east to west (changing first) then from north to south !!!
  # Col 1. Record number
  # Col 2. Position in bytes
  # Col 3. Date (YYYYMMDDHH)
  # Col 4. Parameter name
  # Col 5. Type of level/layer (grib PDS octet 10)
  # Col 6. KPDS5, KPDS6, KPDS7 (grib PDS octets 9, 10, 11-12)
  # Col 7. Forecasts, analysis, etc.
  # Col 8. Description of parameter type
  # 1:0:D=2001010118:TMP:2 m above gnd:kpds=11,105,2:anl:"Temperature [K]
  # 2:143796:D=2001010118:SPFH:2 m above gnd:kpds=51,105,2:anl:"Specific humidity [kg/kg]
  # 3:317756:D=2001010118:PRES:sfc:kpds=1,1,0:anl:"Pressure [Pa]
  # 4:491716:D=2001010118:UGRD:10 m above gnd:kpds=33,105,10:anl:"Zonal wind speed [m/s]
  # 5:615402:D=2001010118:VGRD:10 m above gnd:kpds=34,105,10:anl:"Meridional wind speed [m/s]
  # 6:739088:D=2001010118:DLWRF:sfc:kpds=205,1,0:anl:"LW radiation flux downwards (surface) [W/m^2]
  # 7:902994:D=2001010117:CONVfrac:sfc:kpds=153,1,0:0-1hr acc:"Fraction of total precipitation that is convective [unitless]
  # 8:1006570:D=2001010118:CAPE:180-0 mb above gnd:kpds=157,116,46080:anl:"Convective Available Potential Energy [J/kg]
  # 9:1180530:D=2001010117:PEVAP:sfc:kpds=228,1,0:0-1hr acc:"Potential evaporation [kg/m^2]
  # 10:1324326:D=2001010117:APCP:sfc:kpds=61,1,0:0-1hr acc:"Precipitation hourly total [kg/m^2]
  # 11:1498286:D=2001010118:DSWRF:sfc:kpds=204,1,0:anl:"SW radiation flux downwards (surface) [W/m^2]
  
  # fill missing values
  if (anyNA(grib$band1)) grib$band1 <- na.approx(grib$band1, rule = 2)
  if (anyNA(grib$band2)) grib$band2 <- na.approx(grib$band2, rule = 2)
  if (anyNA(grib$band3)) grib$band3 <- na.approx(grib$band3, rule = 2)
  if (anyNA(grib$band4)) grib$band4 <- na.approx(grib$band4, rule = 2)
  if (anyNA(grib$band5)) grib$band5 <- na.approx(grib$band5, rule = 2)
  if (anyNA(grib$band6)) grib$band6 <- na.approx(grib$band6, rule = 2)
  if (anyNA(grib$band10)) grib$band10 <- na.approx(grib$band10, rule = 2)
  if (anyNA(grib$band11)) grib$band11 <- na.approx(grib$band11, rule = 2)
  
  # extract the variables needed for SUMMA HRUs
  nl <- data.frame(airtemp = grib$band1[frac$idx],
                   spechum = grib$band2[frac$idx],
                   airpres = grib$band3[frac$idx],
                   windspd = sqrt(grib$band4[frac$idx]^2+grib$band5[frac$idx]^2),
                   LWRadAtm = grib$band6[frac$idx],
                   pptrate = grib$band10[frac$idx]/3600.0,
                   SWRadAtm = grib$band11[frac$idx])
  
  # multiply the fraction rate of each fraction  
  nl <- nl*frac$frac
  
  # sum the fraction to hru level
  hru.var <- aggregate(nl, list(hruID <- frac$hru_id2), FUN = sum)

  if (nrow(hru.var) != nhru) {
    print(paste0("convertNLDAS2NC fails at ", NLfname))
    stop
  }
  
  # output varaible names
  varnames <- names(hru.var)

  # creat and write the netCDF file
  nc <- nc_create(NCfname, vars=varList)
  ncvar_put(nc,"time",t)
  for (i in 2:8){
    ncvar_put(nc,varnames[i],hru.var[,i])
  }
  nc_close(nc)
  # delete the grib file
  file.remove(NLfname)
  return(0)
  
}
  
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

# day time
ref.date <- as.Date("1900-01-01")
hours <- sprintf("%02i",0:23)
all.date <- seq(start.date,end.date,by="day")

# define the url of NLDAS data (prefix) 
# sample:     ftp://hydro1.sci.gsfc.nasa.gov/data/s4pa/NLDAS/NLDAS_FORA0125_H.002/1994/005/NLDAS_FORA0125_H.A19940105.0000.002.grb
#             ftp://hydro1.sci.gsfc.nasa.gov/data/s4pa/NLDAS/NLDAS_FORA0125_H.002/2005/001/NLDAS_FORA0125_H.A200501010000.002.grb
#url.part1 <- "ftp://hydro1.sci.gsfc.nasa.gov/data/s4pa/NLDAS/NLDAS_FORA0125_H.002/"


# function to parse the date for a hourly loop to download NLDAS data
getThisDay <- function(this.day){
  ddiff <- as.numeric(this.day-as.Date("1900-01-01"))
  # Julian day
  dd <- julianDay(this.day) 
  
  jday <- sprintf("%03i",dd)
  for (hh in 0:23){
    # specify the NLDAS file names NLDAS_FORA0125_H.A<YYYYMMDD>.<HH>00.002.grb
    hrs <- sprintf("%02i",hh)
    NLfname <- paste0("NLDAS_FORA0125_H.A",format(this.day,"%Y%m%d."),    hrs,"00.002.grb")
    NCfname <- paste0(format(this.day,"%Y%m/"),format(this.day,"%Y%m%d."),hrs,".nc")
    # obtain the link
    NLurl <- paste0("ftp://hydro1.sci.gsfc.nasa.gov/data/s4pa/NLDAS/NLDAS_FORA0125_H.002/",format(this.day,"%Y"),"/",jday,"/",NLfname)
    t <- ddiff+hh/24
    convertNLDAS2NC(NLurl,NLfname,NCfname,t)
  }
  
  print(paste("Finish", format(this.day,"%Y%m%d")))
}


# loop to create all folders
for (mm in format(seq(start.date,end.date,by="month"),"%Y%m")){
  if (! dir.exists(mm)) dir.create(mm,recursive = TRUE)
}



# request the available cores
cl <- makeCluster(16, type="FORK")

# export shared variables and functions
clusterExport(cl, c("frac", "dimList", "varList", "nhru", "convertNLDAS2NC"))
clusterEvalQ(cl, library(rgdal))
clusterEvalQ(cl, library(ncdf4))
clusterEvalQ(cl, library(zoo))


### excute ###
all.date <- seq(start.date,end.date,by="day")
clusterMap(cl, getThisDay, this.day = all.date, .scheduling = 'dynamic')

# close cluster run
stopCluster(cl)

# final check if any netcdf file is missing (single thread)
foreach(this.day=all.date) %do% {
  ddiff <- as.numeric(this.day-as.Date("1900-01-01"))
  dd <- julianDay(this.day)
  foreach(hh = 0:23) %do% {
    hrs <- sprintf("%02i",hh)
    NCfname <- paste0(format(this.day,"%Y%m/"),format(this.day,"%Y%m%d."),hrs,".nc")
    if (! file.exists(NCfname)){
      NLfname <- paste0("NLDAS_FORA0125_H.A",format(this.day,"%Y%m%d."),    hrs,"00.002.grb")
      jday <- sprintf("%03i",dd)
      NLurl <- paste0("ftp://hydro1.sci.gsfc.nasa.gov/data/s4pa/NLDAS/NLDAS_FORA0125_H.002/",format(this.day,"%Y"),"/",jday,"/",NLfname)
      convertNLDAS2NC(NLurl,NLfname,NCfname, ddiff+hh/24.0)
    }
  }
  
  # check if every file is there
  if (length(dir(format(this.day,"%Y%m"),pattern=paste0(format(this.day,"%Y%m%d."),"*.nc"))) != 24) print(paste("Missing data on ",format(this.day,"%Y%m%d.")))
}  
