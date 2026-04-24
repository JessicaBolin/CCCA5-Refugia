################################################################################################
# Build adult sardine species distribution model using environmental predictors from the 3km WC15 ROMS
# Predictors are SST, upper 50m zooplankton biomass, and spawning stock biomass 
# Chlorophyll didn't improve AUC, so dropped from models
# Contact bmuhling@ucsc.edu
################################################################################################

# Load libraries
library(sdmTMB)
library(sp)
library(ggforce) # enables plot_anisotropy
library(pROC)
# remotes::install_github("pbs-assess/sdmTMBextra", dependencies = TRUE)
library(sdmTMBextra) # for adding barrier mesh
library(dplyr)
library(sf)
library(geodata) # To get coastlines
library(ggplot2)
library(visreg)
library(ggeffects)
library(patchwork)

# Load data (not public)
dat <- readRDS("./caClimate/data/sardObsEnvExtracted.rds")
# Log transform chlorophyll
dat$logChl <- log(dat$chl)
# Change 4 early October samples to September
dat$mo <- ifelse(dat$mo == 10, 9, dat$mo)

# Project data
dat2 <- dat
coordinates(dat2) <- c("lon", "lat")
proj4string(dat2) <- CRS("+proj=longlat +datum=WGS84")
datP <- spTransform(dat2, CRS("+proj=utm +zone=10 ellps=WGS84")) 
# Convert back to dataframe
datP <- as.data.frame(datP)
datP$X <- datP$coords.x1 / 1000 # Changes from UTM m to UTM km
datP$Y <- datP$coords.x2 / 1000

# Test/train split
yrs <- data.frame("yr" = unique(datP$yr))
set.seed(2)
index <- sample(1:nrow(yrs), round(0.75 * nrow(yrs))) 
trainYrs <- yrs[index,]
testYrs <- yrs[-index,]
train <- subset(datP, yr %in% trainYrs) 
test <- subset(datP, yr %in% testYrs) 
# table(train$yr) 
# table(test$yr)

# Next construct the mesh: see https://pbs-assess.github.io/sdmTMB/articles/basic-intro.html:
# Based on testing various mesh sizes and testing against withheld data, cutoff = 75 is reasonable
meshTrain <- make_mesh(train, xy_cols = c("X", "Y"), cutoff = 75)  
# plot(meshTrain) 
# Add barrier to mesh for coast
# First use geodata::gadm to get a coastline 
us <- gadm(country = c('United States'), level = 0, path = "./caClimate/data") 
mx <- gadm(country = c('Mexico'), level = 0, path = "./caClimate/data")
# Join and trim
usmx <- rbind(us, mx)
bbox <- st_bbox(c(xmin = -130, xmax = -114, ymin = 30, ymax = 49))
usmxCrop <- sf::st_as_sf(usmx) %>% st_crop(bbox)
# Project to be same projection as mesh
usmxProj <- st_transform(usmxCrop, crs = "+proj=utm +zone=10 ellps=WGS84") 
# Add barrier
mesh2 <- add_barrier_mesh(meshTrain, usmxProj, proj_scaling = 1000)
# plot(mesh2)

# For exporting projections: save a version of mesh2 that can be used to subset environmental netcdfs
# (Although nearly all of 3km ROMS projections are within the mesh)
meshLocns <- data.frame(mesh2$mesh_sf)
# plot(meshLocns$V1, meshLocns$V2)
# Hull around outermost points
lonlat <- meshLocns[, 1:2]
ch <- chull(lonlat[c("V1", "V2")])
coords <- lonlat[c(ch, ch[1]), ]
colnames(coords) <- c("lon", "lat")
# Project to lon/lat (WGS84)
coords$lon <- coords$lon * 1000
coords$lat <- coords$lat * 1000
coordinates(coords) <- c("lon", "lat")
proj4string(coords) <- CRS("+proj=utm +zone=10 ellps=WGS84")
coords2 <- spTransform(coords, CRS("+proj=longlat +datum=WGS84"))
coords2 <- as.data.frame(coords2)
saveRDS(coords2, "./caClimate/data/meshOutlineWGS84.rds")

#################################################################################################
# Now build the models: one with spatiotemporal effects, one without
set.seed(1)
fit0 <- sdmTMB(sardPA ~ s(sst) + s(zoo50) + s(ssb, k = 3), 
               data = train, spatial = "off", 
               family = binomial(link = "logit"), anisotropy = FALSE)
set.seed(1)
fit1 <- sdmTMB(sardPA ~ s(sst) + s(zoo50) + s(ssb, k = 3), 
               data = train, mesh = mesh2, time = "mo", spatiotemporal = "ar1", 
               family = binomial(link = "logit"), anisotropy = FALSE)

# Check model skill, fit, and convergence
sanity(fit0) 
sanity(fit1) 

# Predict (outputs a new df)
pTrain0 <- predict(fit0, newdata = train, type = "response")
pTest0 <- predict(fit0, newdata = test, type = "response")
pTrain1 <- predict(fit1, newdata = train, type = "response")
pTest1 <- predict(fit1, newdata = test, type = "response")
# Calculate overall AUCs 
auc(pTrain0$sardPA, pTrain0$est, direction = "<", quiet = TRUE) 
auc(pTest0$sardPA, pTest0$est, direction = "<", quiet = TRUE) 
auc(pTrain1$sardPA, pTrain1$est, direction = "<", quiet = TRUE) 
auc(pTest1$sardPA, pTest1$est, direction = "<", quiet = TRUE) 

# Save models
out0 <- list("fit0" = fit0, "train" = pTrain0, "test" = pTest0)
save(out0, file = "./caClimate/models/sdmTMB_adults_noST_wholeYr_split_FinalModel.rda")
out1 <- list("fit1" = fit1, "train" = pTrain1, "test" = pTest1)
save(out1, file = "./caClimate/models/sdmTMB_adults_wST_wholeYr_split_FinalModel_cutoff75.rda") 

# Simple partial plots
# SST
sstPartial0 <- ggpredict(fit0, terms = "sst [all]") |> plot() + xlab("SST") + ylab("Probability of Sardine Occurrence (%)") + 
  ggtitle("No Spatiotemporal Effects")
sstPartial0 
sstPartial1 <- ggpredict(fit1, terms = "sst [all]", typical = "median") |> plot() + xlab("SST") + ylab("Probability of Sardine Occurrence (%)") + 
  ggtitle("With Spatiotemporal Effects")
sstPartial1 

# Zoo50
zooPartial0 <- ggpredict(fit0, terms = "zoo50 [all]") |> plot() + xlab("Upper 50m Zooplankton") + 
  ylab("Probability of Sardine Occurrence (%)") + ggtitle("")
zooPartial0 
zooPartial1 <- ggpredict(fit1, terms = "zoo50 [all]", typical = "median") |> plot() + xlab("Upper 50m Zooplankton") + 
  ylab("Probability of Sardine Occurrence (%)") + ggtitle("")
zooPartial1 

# SSB
ssbPartial0 <- ggpredict(fit0, terms = "ssb [all]") |> plot() + xlab("Spawning Stock Biomass (mt)") + 
  ylab("Probability of Sardine Occurrence (%)") + ggtitle("")
ssbPartial0 
ssbPartial1 <- ggpredict(fit1, terms = "ssb [all]", typical = "median") |> plot() + xlab("Spawning Stock Biomass (mt)") + 
  ylab("Probability of Sardine Occurrence (%)") + ggtitle("")
ssbPartial1 

# Combine figures and export at high resolution
p <- (sstPartial0 + sstPartial1) / (zooPartial0 + zooPartial1)  / ( ssbPartial0 + ssbPartial1)
p
ggsave(filename = "./caClimate/plots/fig1_sdmtmb_partials.tiff", plot = p, dpi = 600, compression = "lzw",
       width = 8, height = 10)
