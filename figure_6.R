##############################################################################################
# Map projections from GAMs predicting RMR/MMR from metabolic experimental data (Lonthair et al. in prep)
# Contact bmuhling@ucsc.edu
##############################################################################################

library(ggplot2)
library(ncdf4)
library(mgcv)
library(reshape2)
library(terra)
library(lubridate)
library(lubridate)

# Load the GAMs
load("./caClimate/models/gamRMR.rda")
load("./caClimate/models/gamMMR.rda")

# Loop through the 3km SST from 3 ESMs
fs <- list.files("E:/projections/roms3km", pattern = ".nc") # 333, 111 per ESM

# # ncdf4 method to explore file if you want
# nc <- nc_open(paste0("E:/projections/roms3km/", fs[i]))
# print(nc)
# lon <- ncvar_get(nc, "xrho")
# lat <- ncvar_get(nc, "yrho")
# nc_close(nc)

for(i in 2:length(fs)) {
  r <-  terra::rast(paste0("E:/projections/roms3km/", fs[i]), nlyr <- 3) 
  # dim(r) # should be 346 lon, 331 lat, 365 time for 3 layers (sst, chl, zoo)
  suppressWarnings(rm(outStack, outList, percRMR))
  # Loop through 365 days
  for (j in 1:dim(r)[3]) { 
    Temp <- r[[j]]
    names(Temp) <- "Temp"
    wgtRMR <- Temp
    values(wgtRMR) <- 43 # is a 15cm fish
    names(wgtRMR) <- "wgtRMR"
    # Combine
    r_list <- list(Temp, wgtRMR)
    r_c <- rast(r_list) 
    predRMR <- terra::predict(r_c, gamRMR, type = "response")
    # plot(predRMR)
    predMMR <- terra::predict(r_c, gamMMR, type = "response")
    # plot(predMMR)
    percRMR <- predRMR / predMMR
    # plot(percRMR) # I think is ok? Low, but higher in north
    names(percRMR) <- "percRMR"
    # Stack daily predictions
    if(exists("outStack")) {
      outList <- list(outStack, percRMR)
      outStack <- rast(outList)
    } else {
      outStack <- percRMR
    }
    # Output some sort of progress counter
    print(paste0("Day ", j, " is complete"))
  }
  # Save outputs
  year <- year(as.Date("1900-01-01") + (terra::depth(r)[20] / 86400)) # Using 20 just to make sure we're not right at year start/end
  esm <- ifelse(grepl("gfdl", fs[i]), "gfdl",
                ifelse(grepl("ipsl", fs[i]), "ipsl", "hadl"))
  writeRaster(outStack, paste0("E:/sdms/romsDownscaled3km/sardMO2/sardPredRMRasPercMMR_", year, "_", esm, ".tif")) # 83MB: quite large
} 

# Re-load and plot a esm/year/date if desired:
had1990 <- rast("/vsizip/E:/sdms/romsDownscaled3km/sardMO2/sardPredRMRasPercMMR_hadl.zip/sardPredRMRasPercMMR_1990_hadl.tif")
had2100 <- rast("/vsizip/E:/sdms/romsDownscaled3km/sardMO2/sardPredRMRasPercMMR_hadl.zip/sardPredRMRasPercMMR_2100_hadl.tif")

# Get the yearday for the date you want
yday(as.Date("1990-09-07"))
# Create rasters and plot
r1990 <- had1990[[250]]
r2100 <- had2100[[250]]
plot(r1990)
plot(r2100)
# Save as single rasters if desired
writeRaster(r1990, "./caClimate/outputs/sardPredRMRasPercMMR_1990_hadl.tif")
writeRaster(r2100, "./caClimate/outputs/sardPredRMRasPercMMR_2100_hadl.tif")
